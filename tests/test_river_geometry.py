import unittest
import numpy as np
from affine import Affine
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union
from nz_terrain.river.geometry import polygon_width, sample_surface, transect


class RiverGeometryTests(unittest.TestCase):
    def test_width_follows_tapered_banks_without_a_narrow_cap(self):
        river=Polygon([(0,-20),(1000,-150),(1000,150),(0,20)])
        line=LineString([(0,0),(1000,0)])
        for distance,expected in [(250,105),(750,235)]:
            section=transect(line,distance,500)
            width,status=polygon_width(section,line.interpolate(distance),river)
            self.assertEqual(status,'mapped_polygon')
            self.assertAlmostEqual(width,expected)

    def test_islands_and_separate_braids_are_not_counted_as_water(self):
        river=Polygon([(-100,-100),(100,-100),(100,100),(-100,100)],
                      holes=[[(-10,-20),(10,-20),(10,20),(-10,20)]])
        section=LineString([(0,-200),(0,200)])
        self.assertEqual(polygon_width(section,Point(0,-50),river),(80.,'mapped_polygon'))
        self.assertEqual(polygon_width(section,Point(0,0),river),(None,'centre_outside_polygon'))
        braids=unary_union([box(-100,-50,100,-30),box(-100,10,100,50)])
        self.assertEqual(polygon_width(section,Point(0,20),braids),(40.,'mapped_polygon'))

    def test_unknown_and_censored_widths_are_not_zero_or_fallback_widths(self):
        river=box(-100,-100,100,100)
        self.assertEqual(polygon_width(LineString([(0,-50),(0,50)]),Point(0,0),river),
                         (None,'transect_truncated'))
        self.assertEqual(polygon_width(LineString([(200,-200),(200,200)]),Point(200,0),river),
                         (None,'centre_outside_polygon'))

    def test_diagonal_route_has_a_perpendicular_metre_transect(self):
        line=LineString([(0,0),(100,100)])
        section=transect(line,line.length/2,500)
        self.assertAlmostEqual(section.length,1000)
        delta=np.diff(np.array(section.coords),axis=0)[0]
        self.assertAlmostEqual(float(delta@np.array([100,100])),0)
        self.assertAlmostEqual(section.distance(line.interpolate(line.length/2)),0)

    def test_off_grid_terrain_is_unknown_instead_of_clamped(self):
        surface=np.arange(16,dtype='float32').reshape(4,4)
        xy=np.array([[.5,.5],[1.,1.],[-1.,.5],[5.,.5]])
        values=sample_surface(surface,xy,Affine.identity())
        self.assertAlmostEqual(values[0],0)
        self.assertAlmostEqual(values[1],2.5)
        self.assertTrue(np.isnan(values[2:]).all())


if __name__=='__main__':unittest.main()
