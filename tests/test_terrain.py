import json
import unittest
from unittest.mock import patch
from tempfile import TemporaryDirectory
from pathlib import Path
import numpy as np
from pyproj import Geod
from nz_terrain import terrain

class GeographyTests(unittest.TestCase):
    def setUp(self):
        self.config=json.loads((terrain.BASE/'locations.json').read_text())

    def test_rotated_grid_preserves_centre_scale_and_bearing(self):
        for slug,c in self.config.items():
            with self.subTest(location=slug):
                a=terrain.grid(c,14336)
                self.assertTrue(np.allclose(a*(2048,2048),(0,0),atol=1e-8))
                self.assertAlmostEqual(abs(a.determinant),3.5**2)
                up=np.array(a*(2048,0));bearing=np.degrees(np.arctan2(up[0],up[1]))%360
                self.assertAlmostEqual(bearing,c['rotation_deg']%360)
                playable=terrain.grid(c,14336);world=terrain.grid(c,57344)
                for x,y in [(0,0),(4096,0),(4096,4096),(0,4096)]:
                    self.assertTrue(np.allclose(playable*(x,y),world*(1536+x/4,1536+y/4)))

    def test_sheet_prefilter_against_cached_stac_extents(self):
        transform=terrain.Transformer.from_crs(4326,2193,always_xy=True)
        checked=0
        for d in json.loads((Path(__file__).parent/'linz-sheet-fixtures.json').read_text()):
            if 'geometry' not in d:continue
            predicted=terrain.sheet_bounds(d.get('id','')+'.json')
            if predicted is None:continue
            # STAC geometries may be clipped to capture areas, so they must
            # fit inside their named sheet/tile; they need not fill it.
            coords=d['geometry']['coordinates']
            if d['geometry']['type']!='Polygon':continue
            for lon,lat in coords[0]:
                x,y=transform.transform(lon,lat)
                self.assertGreaterEqual(x,predicted[0]-.1)
                self.assertLessEqual(x,predicted[2]+.1)
                self.assertGreaterEqual(y,predicted[1]-.1)
                self.assertLessEqual(y,predicted[3]+.1)
            checked+=1
        self.assertEqual(checked,4)

    def test_depth_profile_and_encoding(self):
        self.assertAlmostEqual(30*np.tanh(300/600),13.8635,places=3)
        self.assertAlmostEqual(30*np.tanh(600/600),22.8478,places=3)
        a=terrain.encode(np.array([-30,0,1000],dtype='float32'))
        self.assertEqual(a.dtype,np.uint16)
        self.assertLess(a[0],terrain.SEA)
        self.assertEqual(a[1],terrain.SEA)
        with self.assertRaises(ValueError):terrain.encode(np.array([np.nan]))
        with self.assertRaises(ValueError):terrain.encode(np.array([4000]))

    def test_ocean_fill_requires_absence_from_both_tile_catalogues(self):
        stem='Copernicus_DSM_COG_10_S43_00_E169_00_DEM'
        item={'url':f'https://example.org/{stem}/{stem}.tif','bbox':[169,-43,170,-42]}
        with TemporaryDirectory() as directory,patch('nz_terrain.terrain.CACHE',Path(directory)):
            thirty=Path(directory)/'copernicus-30m-tileList.txt'
            ninety=Path(directory)/'copernicus-90m-tileList.txt'
            thirty.write_text(stem+'\n');ninety.write_text('')
            self.assertFalse(terrain.verified_ocean_tile(item))
            thirty.write_text('');ninety.write_text(stem.replace('_COG_10_','_COG_30_')+'\n')
            self.assertFalse(terrain.verified_ocean_tile(item))
            ninety.write_text('')
            self.assertTrue(terrain.verified_ocean_tile(item))
            self.assertEqual(len(item['ocean_tile_verified_against']),2)

    def test_verified_ocean_fill_is_confined_to_its_geographic_tile(self):
        config={'longitude':168.,'latitude':-44.,'rotation_deg':0}
        result=terrain.ocean_tile_grid({'bbox':[168,-45,169,-43]},config,320,size=32)
        self.assertTrue(np.isnan(result[:,:16]).all())
        self.assertTrue((result[:,16:]==0).all())

if __name__=='__main__':unittest.main()
