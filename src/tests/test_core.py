import math
import unittest
from fractions import Fraction as F

import sphere

def close(x, y):
    return abs(x - y) < 1e-15

class TestCore(unittest.TestCase):
    def setUp(self):
        self.p1 = sphere.Point(F(1), F(0), F(0))
        self.p2 = sphere.Point(F(0), F(2), F(0))
        self.p3 = sphere.Point(F(0), F(0), F(3))
        self.p4 = sphere.Point(F(4), F(3), F(0))
        self.p5 = sphere.Point(F(-4), F(3), F(0))
        self.p6 = sphere.Point(F(1), F(2), F(2))
        self.p7 = sphere.Point(F(-1), F(2), F(2))
        self.p8 = sphere.Point(F(1), F(-1, 10), F(2))
        self.p9 = sphere.Point(F(0), F(-1), F(1))
        self.points = [self.p1, self.p2, self.p3, self.p4, self.p5, self.p6,
                       self.p6, self.p8, self.p9]
        self.seg1 = sphere.Segment(self.p1, self.p2)
        self.seg2 = sphere.Segment(self.p3, self.p4)
        self.seg3 = sphere.Segment(self.p4, self.p5)
        self.seg4 = sphere.Segment(self.p2, self.p3)
        self.seg5 = sphere.Segment(self.p6, self.p7)
        self.seg6 = sphere.Segment(sphere.Point(F(0), F(1), F(2)),
                                   sphere.Point(F(0), F(2), F(1)))
        self.segs = [self.seg1, self.seg2, self.seg3, self.seg4, self.seg5,
                                self.seg6]
        self.c1 = sphere.GreatCircle.through(self.p1, self.p2)
        self.c2 = sphere.GreatCircle.through(self.p1, self.p3)
        self.c3 = sphere.GreatCircle.through(self.p1, self.p4)
        self.poly1 = sphere.Polygon([self.p3, self.p1, self.p4, self.p2],
                                        self.p5)
        self.poly2 = sphere.Polygon([self.p3, self.p4, self.p5],
                                        -self.p5)
        self.polys = [self.poly1, self.poly2]
    def test_point_distance(self):
        self.assertEqual(self.p1.distance_to(self.p2), math.sqrt(2))
        self.assertEqual(self.p1.distance_to(self.p1), 0)
    def test_equality(self):
        self.assertFalse(self.p1 == self.p2)
        self.assertTrue(self.p1 == self.p1)
        self.assertTrue(self.p1 != self.p2)
        self.assertFalse(self.p1 != self.p1)
    def test_segment(self):
        self.assertTrue(close(self.seg1.get_arc(), math.pi / 2))
    def test_cross(self):
        self.assertEqual(sphere.cross_product(self.p1, self.p2), self.p3)
        self.assertEqual(sphere.cross_product(self.p2, self.p1), -self.p3)
        self.assertEqual(sphere.cross_product(self.p1, self.p4), self.p3)
    def test_orthogonal(self):
        self.assertTrue(self.p1.orthogonal_to(self.p2))
        self.assertTrue(self.p1.orthogonal_to(self.p3))
        self.assertFalse(self.p1.orthogonal_to(self.p4))
    def test_greatcircle(self):
        self.assertTrue(self.c1.contains(self.p1))
        self.assertTrue(self.c1.contains(self.p2))
        self.assertTrue(self.c1.contains(self.p4))
        self.assertTrue(self.c1.contains(self.p5))
        self.assertFalse(self.c1.contains(self.p3))
        self.assertTrue(self.c1 == self.c1)
        self.assertFalse(self.c1 != self.c1)
        self.assertTrue(self.c1 != self.c2)
        self.assertFalse(self.c1 == self.c2)
        self.assertFalse(self.c1 != self.c3)
        self.assertTrue(self.c1 == self.c3)
        self.assertTrue(self.c1.intersection(self.c2).collinear(self.p1))
    def test_segment(self):
        self.assertTrue(self.seg1.contains(self.p1))
        self.assertTrue(self.seg1.contains(self.p2))
        self.assertTrue(self.seg1.contains(self.p4))
        self.assertFalse(self.seg1.contains(self.p5))
        self.assertFalse(self.seg1.contains(self.p5, True))
        self.assertFalse(self.seg1.contains(self.p5))
    def test_segment_intersection(self):
        self.assertTrue(self.seg5.intersection(self.seg4))
        self.assertTrue(self.seg1.intersection(self.seg2))
        self.assertTrue(self.seg1.intersection(self.seg5) is None)
        self.assertTrue(self.seg1.intersection(self.seg2))
    def test_polygon_contains(self):
        # TODO add tests for all (or at least more) cornercases
        self.assertTrue(self.poly1.contains(self.p6))
        self.assertFalse(self.poly1.contains(self.p5))
        self.assertFalse(self.poly1.contains(self.p7))
        self.assertFalse(self.poly1.contains(-self.p5))
        self.assertFalse(self.poly1.contains(self.p8))
    def test_polygon_intersection(self):
        # TODO add tests for all (or at least more) cornercases
        ret = self.poly1.intersection(self.poly2)
        self.assertEqual(len(ret), 1)
        self.assertEqual(ret[0],
                        sphere.Polygon([self.p2, self.p3, self.p4], self.p1))

    def test_internal_point(self):
        internal_point1 = self.poly1.find_internal_point()
        internal_point2 = self.poly2.find_internal_point()
        self.assertTrue(self.poly1.contains(internal_point1))
        self.assertTrue(self.poly2.contains(internal_point2))
        self.assertFalse(self.poly1.border_contains(internal_point1))
        self.assertFalse(self.poly2.border_contains(internal_point2))
    def test_union(self):
        ret = self.poly1.intersection(self.poly2)
        # TODO add test
    def test_collinear(self):
        self.assertTrue(sphere.e_z.collinear(sphere.e_z))
        self.assertTrue(sphere.e_z.collinear(-sphere.e_z))
        self.assertFalse(sphere.e_z.collinear(sphere.e_x))
        self.assertFalse(sphere.e_z.collinear(-sphere.e_x))
    def test_orthonormal_basis_for(self):
        for point in self.points:
            basis = sphere.orthonormal_basis_for(point)
            self.assertEqual(basis[0], point)
            self.assertTrue(basis[1].orthogonal_to(basis[0]))
            self.assertTrue(basis[1].orthogonal_to(basis[2]))
            self.assertTrue(basis[0].orthogonal_to(basis[2]))
    def test_angle_of(self):
        self.assertTrue(close(0,
            self.c1.point_at(self.c1.angle_of(self.p1)).distance_to(self.p1)))
        self.assertTrue(close(0,
            self.c1.point_at(self.c1.angle_of(self.p2)).distance_to(self.p2)))
        self.assertTrue(close(0,
            self.c2.point_at(self.c2.angle_of(self.p1)).distance_to(self.p1)))
        self.assertTrue(close(0,
            self.c2.point_at(self.c2.angle_of(self.p3)).distance_to(self.p3)))
        self.assertTrue(close(0,
            self.c3.point_at(self.c3.angle_of(self.p1)).distance_to(self.p1)))
        self.assertTrue(close(0,
            self.c3.point_at(self.c3.angle_of(self.p4)).distance_to(self.p4)))
    def test_open_split(self):
        polys = self.poly1.open_split()
        # TODO add test
    def test_northern_orthocomplement(self):
        import sphere.viewer
        sphere.viewer.view(self.poly2.northern_orthocomplement())

        # TODO add test

if __name__ == '__main__':
    unittest.main()
