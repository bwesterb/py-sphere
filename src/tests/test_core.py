import math
import unittest
from fractions import Fraction as F

import sphere

def close(x, y):
    return abs(x - y) < 1e-15

class TestCore(unittest.TestCase):
    def setUp(self):
        self.p1 = sphere.Point(F(1), F(0), F(0))
        self.p2 = sphere.Point(F(0), F(1), F(0))
        self.p3 = sphere.Point(F(0), F(0), F(1))
        self.p4 = sphere.Point(F(4, 5), F(3, 5), F(0))
        self.p5 = sphere.Point(F(-4, 5), F(3, 5), F(0))
        self.seg = sphere.Segment(self.p1, self.p2)
        self.c1 = sphere.GreatCircle.through(self.p1, self.p2)
        self.c2 = sphere.GreatCircle.through(self.p1, self.p3)
        self.c3 = sphere.GreatCircle.through(self.p1, self.p4)
    def test_point_distance(self):
        self.assertEqual(self.p1.distance_to(self.p2), math.sqrt(2))
        self.assertEqual(self.p1.distance_to(self.p1), 0)
    def test_equality(self):
        self.assertFalse(self.p1 == self.p2)
        self.assertTrue(self.p1 == self.p1)
    def test_segment(self):
        self.assertTrue(close(self.seg.get_arc(), math.pi / 2))
    def test_cross(self):
        self.assertEqual(sphere.cross(self.p1, self.p2), self.p3)
        self.assertEqual(sphere.cross(self.p2, self.p1), -self.p3)
        self.assertEqual(sphere.cross(self.p1, self.p4), self.p3)
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
        self.assertTrue(self.c1.intersect(self.c2).collinear(self.p1))
    def test_segment(self):
        self.assertTrue(self.seg.contains(self.p1))
        self.assertTrue(self.seg.contains(self.p2))
        self.assertTrue(self.seg.contains(self.p4))
        self.assertFalse(self.seg.contains(self.p5))
        self.assertFalse(self.seg.contains(self.p5, True))
        self.assertFalse(self.seg.contains(self.p5))

if __name__ == '__main__':
    unittest.main()
