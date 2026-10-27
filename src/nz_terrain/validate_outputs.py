"""Validate completed packages and write an overview of all configured maps."""
import hashlib,json,zipfile
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from . import terrain
from .paths import preview_font_path


def validate():
    locations=json.loads((terrain.BASE/'locations.json').read_text())
    summary={};thumbs=[]
    for slug,location in locations.items():
        base=terrain.BASE/'output'/slug;candidate=terrain.BASE/'output'/f'{slug}-rivers'
        manifest=json.loads(terrain.output_file(base,'manifest.json').read_text())
        report=json.loads(terrain.output_file(candidate,'river-report.json').read_text())
        assert manifest['location']==report['location']==location,slug
        baseline_manifest=json.loads(terrain.output_file(candidate,'baseline-manifest.json').read_text())
        assert baseline_manifest==manifest,(slug,'stale baseline metadata')
        original={};carved={};checks={}
        for kind in terrain.WIDTHS:
            p=terrain.output_file(base,f'{slug}-{kind}.png');q=terrain.output_file(candidate,f'{slug}-rivers-{kind}.png')
            original[kind]=np.array(Image.open(p));carved[kind]=np.array(Image.open(q))
            for a in [original[kind],carved[kind]]:assert a.dtype==np.uint16 and a.shape==(4096,4096)
            checks[kind]={'baseline_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rivers_sha256':hashlib.sha256(q.read_bytes()).hexdigest()}
            assert checks[kind]['baseline_sha256']==manifest['grids'][kind]['sha256']
        assert report['baseline_sha256']==checks['playable']['baseline_sha256']
        assert np.all(carved['playable']<=original['playable']),slug
        changed=carved['playable']!=original['playable']
        assert int(changed.sum())==report['changed_pixels']
        affected=changed.reshape(1024,4,1024,4).any(axis=(1,3))
        means=carved['playable'].astype('float64').reshape(1024,4,1024,4).mean(axis=(1,3))
        centre=carved['world'][1536:2560,1536:2560]
        error=float(np.abs(centre[affected]-means[affected]).max()) if affected.any() else 0
        assert error<=.5,(slug,error)
        expected=original['world'].copy();expected[1536:2560,1536:2560][affected]=np.rint(means[affected]).astype('uint16')
        assert np.array_equal(expected,carved['world']),(slug,'unexpected scenery changes')
        for directory,label in [(base,slug),(candidate,slug+'-rivers')]:
            with zipfile.ZipFile(directory/f'{label}-heightmaps.zip') as archive:
                assert archive.testzip() is None
                members=archive.namelist()
                assert all(not n.startswith('/') and '..' not in Path(n).parts for n in members)
                for kind in terrain.WIDTHS:
                    name=f'heightmaps/{label}-{kind}.png'
                    assert archive.read(name)==(directory/name).read_bytes()
        water=json.loads(terrain.output_file(base,'water-report.json').read_text())
        summary[slug]={'location':location,'river_reaches':len(report['rivers']),
                       'river_length_km':sum(r['length_m'] for r in report['rivers'].values())/1000,
                       'lakes_in_playable':sum('playable' in x['grids'] for x in water['lakes']),
                       'lakes_in_world':sum('world' in x['grids'] for x in water['lakes']),
                       'uphill_steps_over_0_2m':report['raster_checks']['uphill_steps_over_0_2m'],
                       'maximum_centre_rounding_units':error,'pngs':checks,'zip_integrity':'passed'}
        image=Image.open(candidate/'previews/playable-preview.jpg').convert('RGB');image.thumbnail((500,500))
        tile=Image.new('RGB',(520,555),'#172534');tile.paste(image,((520-image.width)//2,45))
        draw=ImageDraw.Draw(tile);font=ImageFont.truetype(preview_font_path(),19)
        draw.text((12,12),f"{location['name']} · {location['rotation_deg']}°",font=font,fill='white');thumbs.append(tile)
        print(slug,summary[slug]['river_reaches'],'reaches;',summary[slug]['lakes_in_playable'],'playable lakes; validated',flush=True)
    canvas=Image.new('RGB',(1560,555*((len(thumbs)+2)//3)),'#172534')
    for i,tile in enumerate(thumbs):canvas.paste(tile,((i%3)*520,(i//3)*555))
    canvas.save(terrain.BASE/'output/locations-preview.jpg',quality=92)
    (terrain.BASE/'output/validation-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

def main():
    validate()


if __name__ == "__main__":
    main()
