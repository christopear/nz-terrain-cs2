import unittest
from shapely.geometry import box
from joins import smooth_joins


class JoinTests(unittest.TestCase):
    def test_addition_only_and_reconstructable(self):
        original=box(0,0,100,20)
        connector=box(-100,-2,10,22)
        result=smooth_joins([original],connector)
        self.assertGreater(result.additions.area,0)
        self.assertTrue(result.final.is_valid)
        self.assertLess(result.before.difference(result.final).area,1e-8)
        self.assertLess(result.additions.difference(result.seam_zone).area,1e-8)
        self.assertLess(original.union(result.synthetic_patch).symmetric_difference(result.final).area,1e-8)
        # The far, unrelated end is untouched.
        self.assertLess(result.before.symmetric_difference(result.final).intersection(box(50,-20,120,40)).area,1e-8)

    def test_disconnected_polygons_have_no_seam_changes(self):
        result=smooth_joins([box(0,0,20,20)],box(50,0,70,20))
        self.assertTrue(result.additions.is_empty)

if __name__=='__main__':unittest.main()
