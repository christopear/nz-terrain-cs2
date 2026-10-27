import unittest
import numpy as np
from affine import Affine
from shapely.geometry import LineString,box
from network_trial import one_anchor,cut_by_station,sample_dem
from widths import tapered_corridor
from build_heightmaps import nearest_stations
import shapely

class NetworkTrialTests(unittest.TestCase):
    def test_side_entry_measures_river_width_not_length(self):
        receiving=box(-500,0,500,100)
        curve,width,angle,samples=one_anchor(LineString([(0,-1000),(0,0)]),receiving)
        self.assertAlmostEqual(width,100,places=6)
        self.assertTrue(curve.is_simple)
        self.assertTrue(receiving.covers(shapely.Point(curve.coords[-1])))
        self.assertLess(abs(curve.coords[-1][0]),1e-8)

    def test_width_and_grade_exclusions_make_real_gaps(self):
        line=LineString([(0,0),(100,0)])
        _,s,w,_=tapered_corridor(line,(.5,10),spacing=1)
        grade=np.where((s>=70)&(s<=80),20.,0.)
        keep=(w>=5)&(grade<15)
        polygon=cut_by_station(line,w,s,keep)
        self.assertFalse(polygon.covers(shapely.Point(20,0)))
        self.assertTrue(polygon.covers(shapely.Point(60,0)))
        self.assertFalse(polygon.covers(shapely.Point(75,0)))
        self.assertTrue(polygon.covers(shapely.Point(90,0)))

    def test_nearest_station_uses_segment_projection(self):
        curve=LineString([(0,0),(10,0),(10,10)])
        points=shapely.points([[3,2],[12,7]])
        np.testing.assert_allclose(nearest_stations(points,curve,np.array([0,10,20])),[3,17])

    def test_dem_sampling_at_pixel_centres(self):
        image=np.arange(9.).reshape(3,3)
        np.testing.assert_allclose(sample_dem(image,Affine.identity(),[[.5,.5],[1.5,1.5]]),[0,4])

if __name__=='__main__':unittest.main()
