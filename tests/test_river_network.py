import unittest
from unittest.mock import patch
import networkx as nx
import numpy as np
from shapely.geometry import LineString,Polygon,GeometryCollection
from nz_terrain import river_network

class NetworkSelectionTests(unittest.TestCase):
    def run_network(self,trusted=False):
        graph=nx.Graph();nodes=[(-1000.,0.),(0.,0.),(1000.,0.)]
        for a,b in zip(nodes,nodes[1:]):
            graph.add_edge(a,b,line=LineString([a,b]),length=1000.,name='Test stream',downstream=b if trusted else None)
        boundary=Polygon([(-7000,-7000),(7000,-7000),(7000,7000),(-7000,7000)])
        config={'sample_spacing_m':10,'preferred_bed_slope':.001,'max_cut_m':10,'surface_smoothing_m':20,
                'minimum_river_width_m':5,'width_coefficient':5,'default_depth_m':1.5,
                'maximum_bottom_width_m':24,'minimum_retained_length_m':30,'rivers':{}}
        if trusted:config['named_overrides']={'Test stream':{'width_m':10,'water_level_m':1,'depth_m':2,'max_cut_m':20}}
        location={'longitude':168.,'latitude':-44.,'rotation_deg':0}
        # Synthetic source heights are known analytically, isolating routing and
        # selection from raster interpolation. The second case deliberately has
        # bad surface elevations that increase downstream.
        def sample(z,xy,affine):return (8+xy[:,0]*.002) if trusted else (30-xy[:,0]*.005)
        with patch('nz_terrain.river_network.prepare_graph',return_value=(graph,GeometryCollection(),boundary,[])),patch('nz_terrain.rivers.sample',side_effect=sample):
            profiles,lines,junctions=river_network.analyse(config,location,np.zeros((1,1)),None)
        return config,profiles,lines

    def test_width_cutoff_keeps_downstream_half_of_same_river(self):
        config,profiles,lines=self.run_network()
        self.assertEqual(len(lines),1)
        line=next(iter(lines.values()))
        self.assertAlmostEqual(line.length,1000,delta=11)
        self.assertEqual(tuple(line.coords[-1]),(1000.,0.))
        self.assertGreater(next(iter(profiles.values()))['estimated_width'].min(),4.99)
        self.assertGreater(config['network_report']['decisions'][0]['width_trim_m'],990)

    def test_trusted_outlet_direction_overrides_inverted_dem(self):
        config,profiles,lines=self.run_network(trusted=True)
        line=next(iter(lines.values()))
        self.assertEqual(tuple(line.coords[0]),(-1000.,0.))
        self.assertEqual(tuple(line.coords[-1]),(1000.,0.))
        bed=next(iter(profiles.values()))['bed']
        self.assertTrue(np.allclose(bed,-1))

    def test_trimmed_tributaries_do_not_duplicate_shared_downstream_suffix(self):
        graph=nx.Graph();join=(0.,0.);root=(2000.,0.)
        for a,b in [((-2000.,1000.),join),((-2000.,-1000.),join),(join,root)]:
            line=LineString([a,b])
            graph.add_edge(a,b,line=line,length=line.length,name='Shared river',downstream=None)
        boundary=Polygon([(-7000,-7000),(7000,-7000),(7000,7000),(-7000,7000)])
        config={'sample_spacing_m':10,'preferred_bed_slope':.001,'max_cut_m':10,'surface_smoothing_m':20,
                'minimum_river_width_m':5,'width_coefficient':5,'default_depth_m':1.5,
                'maximum_bottom_width_m':24,'minimum_retained_length_m':30,'rivers':{}}
        location={'longitude':168.,'latitude':-44.,'rotation_deg':0}
        def sample(z,xy,affine):
            x=xy[:,0]
            return 30-.005*x+50*np.exp(-((x-900)/50)**2)
        with patch('nz_terrain.river_network.prepare_graph',return_value=(graph,GeometryCollection(),boundary,[])),patch('nz_terrain.rivers.sample',side_effect=sample):
            profiles,lines,_=river_network.analyse(config,location,np.zeros((1,1)),None)
        self.assertEqual(len(profiles),1)
        self.assertGreater(next(iter(lines.values())).coords[0][0],0)
        self.assertIn('downstream suffix already covered',[d['status'] for d in config['network_report']['decisions']])

if __name__=='__main__':unittest.main()
