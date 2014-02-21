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
        self.seg1 = sphere.Segment(self.p1, self.p2)
        self.seg2 = sphere.Segment(self.p3, self.p4)
        self.seg3 = sphere.Segment(self.p4, self.p5)
        self.seg4 = sphere.Segment(self.p2, self.p3)
        self.seg5 = sphere.Segment(self.p6, self.p7)
        self.c1 = sphere.GreatCircle.through(self.p1, self.p2)
        self.c2 = sphere.GreatCircle.through(self.p1, self.p3)
        self.c3 = sphere.GreatCircle.through(self.p1, self.p4)
        self.poly1 = sphere.Polygon([self.p3, self.p1, self.p4, self.p2],
                                        self.p5)
        self.poly2 = sphere.Polygon([self.p3, self.p4, self.p5],
                                        -self.p5)
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
        print self.poly1.intersection(self.poly2)

if __name__ == '__main__':
    unittest.main()
