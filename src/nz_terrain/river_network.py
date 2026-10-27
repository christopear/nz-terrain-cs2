"""Experimental mapped-network selection with per-reach headwater trimming.

Width is an explicitly heuristic proxy: coefficient * sqrt(upstream mapped km).
It is not a survey and is not a hydrological discharge model.
"""
import json
from collections import Counter
import networkx as nx
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from pyproj import Transformer
from shapely.geometry import shape,LineString,Point,Polygon
from shapely.ops import transform,unary_union,substring
from . import terrain


def parts(geometry):
    if geometry.is_empty:return []
    if geometry.geom_type=='LineString':return [geometry]
    return [part for g in geometry.geoms for part in parts(g)] if hasattr(geometry,'geoms') else []


def prepare_graph(config,location,z,baseline):
    from .rivers import sample
    affine=terrain.grid(location,14336)
    project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
    boundary=Polygon([affine*p for p in [(0,0),(4096,0),(4096,4096),(0,4096)]])
    featurefile=terrain.output_file(baseline,'lake-polygons.geojson')
    water=unary_union([transform(project,shape(f['geometry'])) for f in json.loads(featurefile.read_text())['features']])
    features=json.loads((terrain.BASE/'river-configs'/config['routes']).read_text())['features']
    records=[]
    for f in features:
        geometry=transform(project,shape(f['geometry'])).intersection(boundary).difference(water)
        for line in parts(geometry):
            if line.length>1:records.append({'line':line,'name':f['properties'].get('name') or 'Unnamed stream',
                                            'trusted_flow':f['properties'].get('flow_direction')=='upstream_to_downstream'})
    endpoints=np.array([r['line'].coords[i] for r in records for i in [0,-1]])
    clusters=nx.Graph();clusters.add_nodes_from(range(len(endpoints)))
    clusters.add_edges_from(cKDTree(endpoints).query_pairs(config['endpoint_snap_m']))
    degree={}
    for nodes in nx.connected_components(clusters):
        indices=list(nodes);mean=endpoints[indices].mean(axis=0)
        for i in indices:endpoints[i]=mean;degree[i]=len(indices)
    for i,record in enumerate(records):
        coords=list(record['line'].coords);coords[0]=endpoints[2*i];coords[-1]=endpoints[2*i+1];record['line']=LineString(coords);record['splits']=[]
    # Join a dangling tributary to the nearest mapped channel, within a bounded
    # tolerance. Retain and report offsets so alignment is reviewable.
    snaps=[]
    for i,xy in enumerate(endpoints):
        if degree[i]>1:continue
        point=Point(xy);own=i//2
        candidates=[(r['line'].distance(point),j) for j,r in enumerate(records) if j!=own]
        if not candidates:continue
        distance,j=min(candidates)
        if distance>config['tributary_snap_m']:continue
        target=records[j]['line'];s=target.project(point)
        if s<config['endpoint_snap_m'] or s>target.length-config['endpoint_snap_m']:continue
        match=target.interpolate(s)
        levels=sample(z,np.array([xy,match.coords[0]]),affine)
        if abs(levels[0]-levels[1])>5:continue
        coords=list(records[own]['line'].coords);coords[0 if i%2==0 else -1]=match.coords[0]
        records[own]['line']=LineString(coords);records[j]['splits'].append(match)
        snaps.append({'distance_m':distance,'river':records[own]['name']})
    graph=nx.Graph()
    for record in records:
        line=record['line'];ss=sorted(set([0.,line.length]+[line.project(p) for p in record['splits']]))
        for start,end in zip(ss,ss[1:]):
            segment=substring(line,start,end)
            if segment.length<.1:continue
            u=tuple(np.round(segment.coords[0],3));v=tuple(np.round(segment.coords[-1],3))
            if u==v:continue
            if graph.has_edge(u,v) and graph[u][v]['length']<=segment.length:continue
            graph.add_edge(u,v,line=segment,length=segment.length,name=record['name'],
                           downstream=v if record['trusted_flow'] else None)
    return graph,water,boundary,snaps


def analyse(config,location,z,baseline):
    from .rivers import sample,solve_profile
    graph,water,boundary,snaps=prepare_graph(config,location,z,baseline)
    affine=terrain.grid(location,14336)
    water_levels=[]
    if baseline is not None:
        levels={x['id']:x['water_level_m'] for x in json.loads(terrain.output_file(baseline,'water-report.json').read_text())['lakes']}
        project=Transformer.from_crs(4326,terrain.projection(location['longitude'],location['latitude']),always_xy=True).transform
        for f in json.loads(terrain.output_file(baseline,'lake-polygons.geojson').read_text())['features']:
            identifier=str(f['properties']['t50_fid'])
            if identifier in levels:water_levels.append((transform(project,shape(f['geometry'])),levels[identifier]))
    def water_level_at(xy):
        point=Point(xy)
        matches=[(g.distance(point),level) for g,level in water_levels if g.distance(point)<5]
        return min(matches)[1] if matches else None
    profiles={};lines={};node_beds={};decisions=[];counter=0;ignored_cycles=0
    for component in sorted(nx.connected_components(graph),key=len,reverse=True):
        sub=graph.subgraph(component)
        tips=[n for n in sub if sub.degree(n)==1] or list(sub)
        directed_tips=[n for n in tips if any(e.get('downstream')==n for e in sub[n].values())]
        if directed_tips:tips=directed_tips
        heights=sample(z,np.array(tips),affine)
        root=tips[int(np.argmin(heights))]
        distances,paths=nx.single_source_dijkstra(sub,root,weight='length')
        parent={n:path[-2] for n,path in paths.items() if len(path)>1}
        ignored_cycles+=sub.number_of_edges()-len(parent)
        children={n:[] for n in sub}
        for n,p in parent.items():children[p].append(n)
        upstream={n:0. for n in sub}
        for n in sub:
            if children[n]:continue
            name=next(iter(sub[n].values()))['name']
            if Point(n).distance(water)<20:
                assumed=35 if 'Hollyford' in name else 12 if 'Awarua' in name else 8
                upstream[n]=(assumed/config['width_coefficient'])**2*1000
            elif Point(n).distance(boundary.boundary)<20:
                upstream[n]=(config['minimum_river_width_m']/config['width_coefficient'])**2*1000
        for n in sorted(sub,key=lambda n:distances[n],reverse=True):
            if n in parent:upstream[parent[n]]+=upstream[n]+sub[n][parent[n]]['length']
        leaves=sorted([n for n in sub if not children[n] and n!=root],key=lambda n:distances[n],reverse=True)
        for leaf in leaves:
            nodes=[leaf]
            while nodes[-1] in parent:
                n=parent[nodes[-1]];nodes.append(n)
                if n in node_beds:break
            coords=[];positions=[0.];names=Counter();total=0.
            for u,v in zip(nodes,nodes[1:]):
                e=sub[u][v];points=list(e['line'].coords)
                if np.linalg.norm(np.array(points[0])-u)>np.linalg.norm(np.array(points[-1])-u):points.reverse()
                points[0]=u;points[-1]=v
                coords.extend(points if not coords else points[1:]);length=LineString(points).length
                total+=length;positions.append(total);names[e['name']]+=length
            if not coords:continue
            line=LineString(coords);label=names.most_common(1)[0][0]
            override=config.get('named_overrides',{}).get(label,{})
            cut_limit=override.get('max_cut_m',config['max_cut_m'])
            s=np.unique(np.r_[np.linspace(0,line.length,int(np.ceil(line.length/config['sample_spacing_m']))+1),positions])
            accumulated=np.interp(s,positions,[upstream[n] for n in nodes])
            widths=config['width_coefficient']*np.sqrt(accumulated/1000)
            if 'width_m' in override:widths=np.full_like(widths,override['width_m'])
            accepted=np.flatnonzero(widths>=config['minimum_river_width_m'])
            if not len(accepted):
                decisions.append({'name':label,'length_m':line.length,'status':'below estimated width cutoff'});continue
            first=int(accepted[0]);width_trim=float(s[first])
            xy=np.array([line.interpolate(d).coords[0] for d in s]);centre=sample(z,xy,affine)
            # Filter on uniform metre spacing, not the uneven vertex samples.
            uniform=np.arange(0,line.length+5,5.)
            reference=np.interp(s,uniform,gaussian_filter1d(np.interp(uniform,s,centre),config['surface_smoothing_m']/5))
            reference=np.maximum(reference,0)
            start_level=water_level_at(xy[0]);end_level=water_level_at(xy[-1])
            for endpoint,level,distance in [(0,start_level,s),(-1,end_level,line.length-s)]:
                if level is not None:
                    weight=np.clip(1-distance/200,0,1);weight=weight*weight*(3-2*weight)
                    reference+=weight*(level-reference[endpoint])
            depth=np.clip(config['default_depth_m']+.05*widths,config['default_depth_m'],3)
            if 'water_level_m' in override:reference[:]=override['water_level_m']
            if 'depth_m' in override:depth[:]=override['depth_m']
            outlet=node_beds.get(nodes[-1],float(reference[-1]-depth[-1]))
            upper=np.minimum(reference-depth,centre-.25)
            lower=np.maximum(reference-cut_limit,centre-cut_limit+.1)
            required=outlet;profile_trim=0.
            # Find the longest feasible downstream suffix; never reject a whole
            # river merely because its upper mapped route crosses a DEM ridge.
            for i in range(len(s)-2,first-1,-1):
                required=max(required,lower[i])
                if required>upper[i]+1e-6:
                    first=i+1;profile_trim=float(s[first])-width_trim;break
            if first>=len(s)-1 or line.length-s[first]<config['minimum_retained_length_m']:
                decisions.append({'name':label,'length_m':line.length,'status':'insufficient feasible downstream length','width_trim_m':width_trim,'profile_trim_m':profile_trim});continue
            ss=s[first:];ref=reference[first:];dep=depth[first:]
            bed=solve_profile(ref,ss,dep,cut_limit,config['preferred_bed_slope'],outlet,source=centre[first:])
            kept=substring(line,float(ss[0]),line.length)
            # A profile cut can start inside a graph edge. A later tributary
            # may then trim back to the same edge, despite shared node anchors.
            # Do not re-carve/count that already retained downstream suffix.
            if lines and kept.difference(unary_union(list(lines.values())).buffer(.1)).length<config['minimum_retained_length_m']:
                decisions.append({'name':label,'status':'downstream suffix already covered',
                                  'original_length_m':line.length,'width_trim_m':width_trim,
                                  'profile_trim_m':profile_trim})
                continue
            counter+=1;key=f'reach-{counter:03d}'
            # One raster cell is insufficient for a reliable cross-section;
            # selection cutoff is 5 m, while modelled bottoms use at least 7 m.
            bottom=np.clip(widths[first:]*.6,7,config['maximum_bottom_width_m'])
            banks=np.maximum(bottom+8,widths[first:]+4)
            banks=np.minimum(banks,60)
            profiles[key]={'s':ss-ss[0],'xy':xy[first:],'centre':centre[first:],'reference':ref,'bed':bed,
                           'bottom_width':bottom,'bank_width':banks,'label':label,'estimated_width':widths[first:]}
            lines[key]=kept
            config['rivers'][key]={'bottom_width_m':float(np.median(bottom)),'bank_width_m':float(np.median(banks)),'depth_m':float(np.median(dep)),'name':label,'max_cut_m':cut_limit}
            for n,d in zip(nodes,positions):
                if d>=ss[0]:node_beds[n]=float(np.interp(d,ss,bed))
            decisions.append({'id':key,'name':label,'status':'carved','original_length_m':line.length,'retained_length_m':kept.length,
                              'inlet_water_level_anchor_m':start_level,'outlet_water_level_anchor_m':end_level,
                              'width_trim_m':width_trim,'profile_trim_m':profile_trim,'estimated_width_range_m':[float(widths[first:].min()),float(widths[first:].max())]})
    config['network_report']={'width_method':'5 * sqrt(upstream mapped channel length in km); configurable coefficient. Assumed width at lake/border inflows. Not measured widths.',
                              'selection_cutoff_m':config['minimum_river_width_m'],'minimum_modelled_bottom_m':7,
                              'decisions':decisions,'bounded_endpoint_snaps':snaps,'non_tree_edges_omitted':ignored_cycles,
                              'notes':['Branches below estimated width are trimmed upstream; downstream sections remain.',
                                       'Upper sections needing more than maximum cut are also trimmed and reported.',
                                       'Flow direction inferred from lowest component endpoint and shortest connected route, not surveyed flow direction.',
                                       'Mapped lake interiors excluded from river routing; lagoon/lake target handles their basin.']}
    if not profiles:raise ValueError('No feasible river reaches selected')
    print(f'Network: {len(profiles)} carved reaches, {sum(l.length for l in lines.values())/1000:.1f} km',flush=True)
    return profiles,lines,{}
