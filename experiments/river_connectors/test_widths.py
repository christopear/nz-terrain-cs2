import unittest
import numpy as np
from shapely.geometry import LineString
from widths import tapered_corridor


class WidthTests(unittest.TestCase):
    def test_measured_sections_and_reversal(self):
        line=LineString([(0,0),(1000,0)])
        poly,distances,widths,run=tapered_corridor(line,(28,18),200,mode="end-transitions")
        for station,expected in [(0,28),(100,25.5),(200,23),(500,23),(800,23),(900,20.5),(1000,18)]:
            section=LineString([(station,-50),(station,50)])
            self.assertAlmostEqual(poly.intersection(section).length,expected,places=6)
        reverse,*_=tapered_corridor(LineString([(1000,0),(0,0)]),(18,28),200,mode="end-transitions")
        self.assertLess(poly.symmetric_difference(reverse).area,1e-7)

    def test_linear_widths_over_full_reach(self):
        poly,s,w,run=tapered_corridor(LineString([(0,0),(1000,0)]),(28,18),mode='linear')
        for station,expected in [(0,28),(100,27),(500,23),(900,19),(1000,18)]:
            section=LineString([(station,-50),(station,50)])
            self.assertAlmostEqual(poly.intersection(section).length,expected,places=6)
        self.assertTrue(np.allclose(np.diff(w),-.01))

    def test_short_reach_ramps_meet_without_overshoot(self):
        poly,s,w,run=tapered_corridor(LineString([(0,0),(100,0)]),(28,18),200,mode="end-transitions")
        self.assertEqual(run,50)
        self.assertAlmostEqual(w[0],28)
        self.assertAlmostEqual(w[-1],18)
        self.assertAlmostEqual(w[len(w)//2],23)
        self.assertTrue(np.all(np.diff(w)<=1e-9))
        self.assertTrue(poly.is_valid)

if __name__=='__main__':unittest.main()
