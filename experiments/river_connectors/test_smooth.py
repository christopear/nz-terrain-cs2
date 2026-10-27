"""Run: uv run python -m unittest discover -s experiments/river_connectors -p 'test_*.py'."""
import unittest
import numpy as np
from shapely.geometry import LineString,box
from connector import connect
from smooth import fit_curve


class SmoothCurveTests(unittest.TestCase):
    def test_offset_line_centres_on_banks_and_reverses(self):
        polygons=[box(-100,-10,0,10),box(200,-20,300,20)]
        curves=[]
        for points in [[(-100,5),(300,5)],[(300,5),(-100,5)]]:
            route=LineString(points)
            result=connect(route,*polygons)
            curve,corridor,ends,_=fit_curve(route,result)
            curves.append(curve)
            self.assertTrue(curve.is_simple and corridor.is_valid)
            self.assertTrue(all(corridor.intersects(p) for p in polygons))
            for end in ends:
                self.assertAlmostEqual(end[1][1],0,places=8)
                self.assertAlmostEqual(abs(end[2][0]),1,places=8)
            xy=np.array(curve.coords)
            for delta in (xy[1]-xy[0],xy[-1]-xy[-2]):
                self.assertLess(abs(delta[1]/delta[0]),.02)
            self.assertLess(max(abs(xy[:,1])),6)
        self.assertLess(curves[0].hausdorff_distance(curves[1]),1e-8)

if __name__=='__main__': unittest.main()
