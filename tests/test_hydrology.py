import unittest
import numpy as np
from rasterio.transform import from_bounds, Affine
from shapely.geometry import Polygon
from nz_terrain import hydrology

class WaterBedTests(unittest.TestCase):
    def test_absolute_targets_overlap_without_adding_and_are_idempotent(self):
        original=np.full((5,5),100.,dtype='float32')
        lake=np.full_like(original,-99.) # surface +1 m, depth 100 m
        river=np.full_like(original,-2.)
        lm=np.zeros((5,5),bool);lm[1:4,1:4]=True
        rm=np.zeros_like(lm);rm[2,:]=True
        targets=[(lake,lm),(river,rm)]
        result=hydrology.combine_targets(original,targets)
        self.assertEqual(result[2,2],-99)
        self.assertEqual(result[2,0],-2)
        self.assertTrue(np.array_equal(result[~(lm|rm)],original[~(lm|rm)]))
        self.assertTrue(np.array_equal(result,hydrology.combine_targets(result,targets)))
        self.assertTrue(np.array_equal(result,hydrology.combine_targets(original,list(reversed(targets)))))

    def test_islands_remain_land_and_rotated_polygon_is_continuous(self):
        lake=Polygon([(-100,-100),(100,-100),(100,100),(-100,100)],
                     holes=[[(-20,-20),(20,-20),(20,20),(-20,20)]])
        affine=Affine.rotation(30)*from_bounds(-150,-150,150,150,100,100)
        window=hydrology.lake_window(lake,affine,100)
        dst,src=hydrology.interior_window(window,100)
        mask=np.zeros((100,100),bool);mask[dst]=window[2][src]
        self.assertFalse(mask[50,50]);self.assertTrue(mask[50,70]);self.assertFalse(mask[0,0])
        depth=hydrology.basin_depth(mask,3,30,600,14)
        self.assertEqual(depth[50,50],0)
        self.assertGreater(depth[50,70],0)
        self.assertLess(depth.max(),30)

    def test_depth_uses_metres_and_parameter_is_ceiling(self):
        mask=np.zeros((201,201),bool);mask[1:-1,1:-1]=True
        d=hydrology.basin_depth(mask,10,109,600,0)
        self.assertAlmostEqual(d[100,100],109*np.tanh(995/600),places=4)
        self.assertTrue(np.all(d>=0));self.assertTrue(np.all(d<=109))
        with self.assertRaises(ValueError):hydrology.basin_depth(mask,10,-30,600,0)

if __name__=='__main__':unittest.main()
