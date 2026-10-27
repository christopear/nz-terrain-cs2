import unittest
import numpy as np
from nz_terrain import rivers

class RiverTests(unittest.TestCase):
    def test_insufficient_drop_relaxes_slope_without_deep_trench(self):
        s=np.arange(0,10001,100.)
        reference=10-s/1000
        bed=rivers.solve_profile(reference,s,3,8,.005,-3)
        self.assertTrue(np.all(np.diff(bed)<=1e-7))
        self.assertAlmostEqual(bed[-1],-3)
        self.assertTrue(np.all(reference-bed<=8+1e-6))
        self.assertLess((bed[0]-bed[-1])/10000,.005)

    def test_tributary_outlet_anchors_and_bad_constraints_fail(self):
        s=np.arange(0,1001,10.);reference=20-s*.01
        bed=rivers.solve_profile(reference,s,2,8,.0025,5)
        self.assertAlmostEqual(bed[-1],5)
        self.assertTrue(np.all(np.diff(bed)<=1e-7))
        with self.assertRaises(ValueError):
            rivers.solve_profile(np.array([10.,0.,20.,0.]),np.arange(4.)*100,2,5,.001,-2)

    def test_carve_preserves_outside_and_caps_cut(self):
        location={'longitude':172.8,'latitude':-40.85,'rotation_deg':30}
        affine=rivers.terrain.grid(location,14336)
        xy=np.array([affine*(12,32),affine*(50,32)])
        profiles={'test':{'xy':xy,'bed':np.array([5.,4.])}}
        config={'max_cut_m':2,'rivers':{'test':{'bank_width_m':28,'bottom_width_m':10}}}
        z=np.full((64,64),10.,dtype='float32')
        result,mask=rivers.carve(z,profiles,config,location)
        self.assertTrue(np.array_equal(result[~mask],z[~mask]))
        self.assertLessEqual(float((z-result).max()),2)
        self.assertGreater(np.count_nonzero(result<z),0)

if __name__=='__main__':unittest.main()
