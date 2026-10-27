import unittest

from affine import Affine
import numpy as np
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon, box
from water_carver import BedParameters, Grid, generate_bed


class BedTests(unittest.TestCase):
    def grid(self, shape=(10, 10), transform=Affine.identity()):
        return Grid(shape, transform, 'EPSG:2193')

    def test_rectangle_has_analytic_depth_and_absolute_elevation(self):
        grid = self.grid()
        result = generate_bed(box(0, 0, 10, 10), grid, water_level_m=100,
                              parameters=BedParameters(3, 1, 'linear'))
        row, col = np.indices(grid.shape)
        expected = np.minimum.reduce([row+.5, col+.5, 9.5-row, 9.5-col, np.full(grid.shape, 3)])
        np.testing.assert_array_equal(result.depth_m, expected)
        np.testing.assert_array_equal(result.elevation_m, 100-expected)
        self.assertTrue(result.mask.all())

    def test_islands_and_disconnected_water_preserve_land(self):
        lake = Polygon([(0,0),(10,0),(10,10),(0,10)], holes=[[(4,4),(6,4),(6,6),(4,6)]])
        water = MultiPolygon([lake, box(12, 0, 15, 3)])
        result = generate_bed(water, self.grid((10, 16)), water_level_m=0,
                              parameters=BedParameters(10, 1, 'linear'))
        self.assertFalse(result.mask[4, 4])
        self.assertTrue(np.isnan(result.elevation_m[4, 4]))
        self.assertEqual(result.depth_m[4, 3], .5)
        self.assertFalse(result.mask[0, 11])
        self.assertTrue(result.mask[0, 12])

    def test_open_ends_do_not_create_sills(self):
        water = box(0, 0, 20, 10)
        banks = MultiLineString([[(0,0),(20,0)], [(0,10),(20,10)]])
        result = generate_bed(water, self.grid((10,20)), water_level_m=0, banks=banks,
                              parameters=BedParameters(10, 1, 'linear'))
        np.testing.assert_array_equal(result.depth_m[:,0], result.depth_m[:,10])
        self.assertEqual(result.depth_m[4, 0], 4.5)

    def test_map_edges_are_not_shores_and_chunking_is_invariant(self):
        water = box(-100,-100,100,100)
        p = BedParameters(200, 1, 'linear')
        a = generate_bed(water, self.grid(), water_level_m=0, parameters=p, chunk_rows=1)
        b = generate_bed(water, self.grid(), water_level_m=0, parameters=p, chunk_rows=128)
        self.assertEqual(a.depth_m[0,0], 99.5)
        np.testing.assert_array_equal(a.elevation_m,b.elevation_m)

    def test_rotated_grid_uses_world_distances_and_levels(self):
        affine = Affine.translation(30,40) * Affine.rotation(27) * Affine.scale(2,3)
        water = Polygon([affine*p for p in [(0,0),(10,0),(10,10),(0,10)]])
        levels = np.arange(100, dtype=float).reshape(10,10)
        before = levels.copy()
        result = generate_bed(water, self.grid(transform=affine), water_level_m=levels,
                              parameters=BedParameters(100,1,'linear'))
        row,col = np.indices((10,10))
        expected = np.minimum.reduce([(row+.5)*3,(9.5-row)*3,(col+.5)*2,(9.5-col)*2])
        np.testing.assert_allclose(result.depth_m, expected, atol=1e-12)
        np.testing.assert_allclose(result.elevation_m, levels-expected, atol=1e-12)
        np.testing.assert_array_equal(levels,before)

    def test_profiles_are_monotone_bounded_and_have_requested_initial_slope(self):
        for name in ('linear','tanh','exponential'):
            with self.subTest(profile=name):
                p=BedParameters(30,.05,name)
                values=p.depths(np.linspace(0,10000,1000))
                self.assertEqual(values[0],0)
                self.assertTrue((np.diff(values)>=0).all())
                self.assertTrue((values<=30).all())
                self.assertAlmostEqual(float(p.depths(1e-5)/1e-5),.05,places=8)

    def test_rejects_invalid_units_geometry_banks_and_parameters(self):
        for crs in ('EPSG:4326','EPSG:2263'):
            with self.assertRaises(ValueError):Grid((10,10),Affine.identity(),crs)
        for p in (dict(max_depth_m=0),dict(bank_gradient=float('nan')),dict(profile='unknown')):
            with self.assertRaises(ValueError):BedParameters(**p)
        for geometry in (Polygon(),LineString([(0,0),(1,1)]),Polygon([(0,0),(2,2),(0,2),(2,0)])):
            with self.assertRaises(ValueError):generate_bed(geometry,self.grid(),water_level_m=0)
        with self.assertRaises(ValueError):
            generate_bed(box(0,0,10,10),self.grid(),water_level_m=0,banks=LineString([(1,1),(2,2)]))
        with self.assertRaises(ValueError):
            generate_bed(box(0,0,10,10),self.grid(),water_level_m=np.zeros((2,2)))
        with self.assertRaises(ValueError):
            generate_bed(box(0,0,10,10),self.grid(),water_level_m=np.nan)

    def test_levels_may_be_missing_only_outside_water(self):
        level=np.full((10,10),np.nan);level[:2,:2]=50
        result=generate_bed(box(0,0,2,2),self.grid(),water_level_m=level)
        self.assertEqual(result.mask.sum(),4)
        self.assertTrue(np.isfinite(result.elevation_m[result.mask]).all())
        level[0,0]=np.nan
        with self.assertRaises(ValueError):generate_bed(box(0,0,2,2),self.grid(),water_level_m=level)

    def test_outside_polygon_returns_empty_mask(self):
        result=generate_bed(box(100,100,110,110),self.grid(),water_level_m=0)
        self.assertFalse(result.mask.any())
        self.assertTrue(np.isnan(result.depth_m).all())


if __name__ == '__main__':
    unittest.main()
