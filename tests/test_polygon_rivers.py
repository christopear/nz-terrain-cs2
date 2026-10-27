"""Checks for polygon-only river integration and the selected profile."""
import unittest
import numpy as np
from shapely.geometry import Polygon, LineString
from nz_terrain.river.polygon_build import retained_banks
from water_carver.river import profile, smooth_creases


class PolygonRiverTests(unittest.TestCase):
    def test_openings_preserve_island_and_side_banks(self):
        water=Polygon([(0,0),(100,0),(100,20),(0,20)],holes=[[(40,8),(60,8),(60,12),(40,12)]])
        banks=retained_banks(water,[LineString([(0,0),(0,20)]),LineString([(100,0),(100,20)])])
        self.assertAlmostEqual(banks.length,water.boundary.length-40)
        self.assertAlmostEqual(banks.intersection(LineString(water.interiors[0].coords)).length,48)
        with self.assertRaises(ValueError):
            retained_banks(water,[LineString([(10,0),(10,20)])])

    def test_preferred_depth_and_no_central_hump(self):
        for width,lower,upper in [(5,.5,1.5),(10,2,4),(40,8,10)]:
            x=np.linspace(0,width,round(width/.1)+1);d=np.minimum(x,width-x)
            raw=profile(d,2,10,.05,1)
            z,_=smooth_creases(raw,d,.1,1)
            self.assertGreater(z.max(),lower);self.assertLessEqual(z.max(),upper)
            np.testing.assert_allclose(z,z[::-1],atol=1e-12)
            self.assertTrue(np.all(np.diff(z[:len(z)//2+1])>=-1e-12))
            np.testing.assert_array_equal(z[d<=.5],raw[d<=.5])

    def test_land_does_not_transmit_smoothing(self):
        d=np.array([0.,1.,2.,np.nan,0.,1.,2.])
        a,_=smooth_creases(profile(d,2,10,.05,1),d,1,1)
        depth=profile(d,2,10,.05,1);depth[4:]=9
        b,_=smooth_creases(depth,d,1,1)
        np.testing.assert_array_equal(a[:3],b[:3])
        self.assertTrue(np.isnan(b[3]))


class IndexedDistanceTests(unittest.TestCase):
    def test_index_matches_direct_bank_distances_with_island_and_opening(self):
        from affine import Affine
        from water_carver import BedParameters, Grid, generate_bed
        from nz_terrain.river.polygon_build import depth_field
        water=Polygon([(0,0),(16,0),(16,16),(0,16)],holes=[[(6,6),(10,6),(10,10),(6,10)]])
        banks=retained_banks(water,[LineString([(0,0),(0,16)])])
        settings=dict(shore_transition_m=2,max_depth_m=10,shore_gradient=.05,interior_gradient=1,
                      smoothing=dict(scale_m=1,shore_lock_m=.5,shore_blend_m=1.5))
        actual=depth_field(water,banks,Affine.identity(),(16,16),'EPSG:2193',settings)
        fine=Affine.scale(.25)*Affine.translation(-8,-8)
        exact=generate_bed(water,Grid((80,80),fine,'EPSG:2193'),water_level_m=0,
                           parameters=BedParameters(100000,1,'linear'),banks=banks)
        z,_=smooth_creases(profile(exact.depth_m,2,10,.05,1),exact.depth_m,.25,1)
        z=z[8:-8,8:-8].reshape(16,4,16,4)
        count=np.isfinite(z).sum(axis=(1,3))
        expected=np.divide(np.nansum(z,axis=(1,3)),count,out=np.full((16,16),np.nan),where=count>0)
        np.testing.assert_allclose(actual,expected,atol=1e-12,equal_nan=True)


if __name__=='__main__':unittest.main()
