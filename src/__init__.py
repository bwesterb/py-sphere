""" The unit sphere. """
import math

class Polygon:
    """ A polygon on the sphere. """
    def __init__(self, vertices, external_point):
        """ Given vertices v1, ..., vn and an external point, creates
            the polygon v1 - v2 - ... - vn - v1 such that the given
            external point is outside and the segment between vi and vj
            is the short segment between vi and vj. """
        assert all([v != external_point for v in vertices])
        assert len(frozenset(vertices)) == len(vertices)
        self.vertices = vertices
        self.external_point = external_point
    def __repr__(self):
        return "<sphere.Polygon %s without %s>" % (
                    self.vertices, self.external_point)
    def __eq__(self, other):
        raise NotImplementedError # TODO
    def __ne__(self, other):
        raise NotImplementedError # TODO

class Segment:
    """ A segment of a Great Circle with arc less than pi radians.  """
    def __init__(self, point1, point2):
        """ Creates a segment by its endpoints. """
        self.point1 = point1
        self.point2 = point2

    def get_arc(self):
        """ Returns the arc in radians covered by the segment. """
        return 2 * math.asin(.5 * self.point1.distance_to(self.point2))
    def get_greatcircle(self):
        """ Returns the great circle of which this arc is a part. """
        return GreatCircle(self.point1, self.point2)
    def contains(self, point, assume_on_greatcircle=False):
        """ Returns whether the point is contained in the segment.

            If it is already known that the point is on the same great
            circle, (but not whether it is on the arc), set
            the flag assume_on_greatcircle to improve performance. """
        if (not assume_on_greatcircle and
                not self.get_greatcircle().contains(point)):
            return False
        if point == self.point1 or point == self.point2:
            return True
        # TODO add reference for this lemma.
        #      Thanks to Bram Westerbaan for the suggestion
        return (self.point1.distance_to(point) ** 2 +
                self.point2.distance_to(point) ** 2
                    <= self.point1.distance_to(self.point2) ** 2)

class GreatCircle:
    """ A circle that splits the sphere in two.  Equivalently: the
        intersection of a plane with the sphere. """
    def __init__(self, point1, point2):
        """ Creates a great circle given by two distinct point. """
        assert point1 != point2
        self.point1 = point1
        self.point2 = point2
    def contains(self, point):
        """ Checks whether point is contained in this great circle. """
        return cross(self.point1, self.point2).orthogonal_to(point)

    def __repr__(self):
        return "<sphere.GreatCircle %s %s>" % (self.point1, self.point2)

    def __eq__(self, other):
        return (self.point1, self.point2) == (other.point1, other.point2)
    def __ne__(self, other):
        return (self.point1, self.point2) != (other.point1, other.point2)

class Point:
    """ A point on the unit sphere. """
    def __init__(self, x, y, z, normalize=False):
        """ Creates a point by its 3D coordinates """
        if normalize:
            n = math.sqrt(x**2 + y**2 + z**2)
            x, y, z = x/n, y/n, z/n
        else:
            assert x**2 + y**2 + z**2 == 1
        self.x, self.y, self.z = x, y, z

    def __str__(self):
        return "(%s, %s, %s)" % (self.x, self.y, self.z)

    def __repr__(self):
        return "<sphere.Point (%s %s %s)>" % (self.x, self.y, self.z)

    def distance_to(self, other):
        """ Returns the distance to another point. """
        return math.sqrt((self.x - other.x)**2 +
                         (self.y - other.y)**2 +
                         (self.z - other.z)**2)
    def orthogonal_to(self, other):
        """ Returns whether this point is orthogonal to other. """
        return self.x * other.x + self.y * other.y + self.z * other.z == 0

    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)
    def __ne__(self, other):
        return (self.x, self.y, self.z) != (other.x, other.y, other.z)
    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)

def cross(point1, point2):
    """ Returns the normalized cross product of point1 and point2. """
    return Point(point1.y * point2.z - point1.z * point2.y,
                 point1.z * point2.x - point1.x * point2.z,
                 point1.x * point2.y - point1.y * point2.x,
                 normalize=True)
