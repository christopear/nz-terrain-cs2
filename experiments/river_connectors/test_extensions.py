import unittest
from shapely.geometry import LineString,box
from trial_mckenzie import extend,extend_into_polygons

class ExtensionTests(unittest.TestCase):
    def test_straight_extension_preserves_observed_line(self):
        line=LineString([(0,0),(100,0)])
        route=extend(line)
        self.assertEqual(list(route.coords),[(-75.,0.),(0.,0.),(100.,0.),(175.,0.)])
        self.assertEqual(route.intersection(line).length,line.length)

    def test_inward_fallback_only_changes_invented_sections(self):
        line=LineString([(0,0),(100,0)])
        # Starting point lies on a corner; the leftward extension must turn inward.
        polys=[box(-100,-100,0,-.01),box(100,-20,200,20)]
        route,turns=extend_into_polygons(line,polys)
        self.assertNotEqual(turns[0],0)
        self.assertEqual(turns[1],0)
        self.assertEqual(route.intersection(line).length,line.length)
        self.assertGreater(route.intersection(polys[0]).length,35)

if __name__=='__main__':unittest.main()
