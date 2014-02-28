""" The unit sphere. """
import math
import random
import fractions

import numpy
import numpy.linalg

def sign(x):
    """ Returns 0 if x is zero.  1 if x is strictly positive.  -1, otherwise,
        when x is strictly negative. """
    if x == 0:
        return 0
    if x >= 0:
        return 1
    return -1

class Error(Exception):
    pass
class SameGreatCircle(Error):
    pass
class DegeneratePolygon(Error):
    pass
class CouldNotSplitThere(Error):
    pass

class Region:
    """ A union of (open) disjoint polygons. """
    def __init__(self, polygons):
        # TODO check whether the polygons are disjoint
        self.polygons = polygons

    def intersection(self, other):
        """ Returns the intersection of this region with the other. """
        ret = []
        for poly1 in self.polygons:
            for poly2 in other.polygons:
                ret.extend(poly1.intersection(poly2))
        return Region(ret)

    def union(self, other):
        """ Returns the union of this region with another. """
        return self.complement().intersection(other.complement()).complement()

    def complement(self):
        """ Returns the complement of this region. """
        if not self.polygons:
            return full_sphere
        ret = None
        for poly in self.polygons:
            if ret is None:
                ret = Region([poly.complement()])
                continue
            ret = ret.intersection(Region([poly.complement()]))
        return ret

    def open_split(self):
        """ Splits this region into two regions such that the union of
            their interior is the union of the interior of this region. """
        if not self.polygons:
            return []
        if len(self.polygons) > 1:
            return (Region([self.polygons[0]]), Region(self.polygons[1:]))
        return Region(self.polygons[0].open_split())

    def __repr__(self):
        return '<Region %s>' % self.polygons

    def __eq__(self, other):
        raise NotImplementedError
    def __ne__(self, other):
        raise NotImplementedError

class Polygon:
    """ A polygon on the sphere. """
    def __init__(self, vertices, external_point):
        """ Given vertices v1, ..., vn and an external point, creates
            the polygon v1 - v2 - ... - vn - v1 such that the given
            external point is outside and the segment between vi and vj
            is the short segment between vi and vj. """
        assert all([v != external_point for v in vertices])
        assert len(frozenset(vertices)) == len(vertices)
        assert len(vertices) >= 3
        self.vertices = vertices
        self.external_point = external_point
        self.segments = []
        for i in xrange(len(self.vertices)):
            if i == 0:
                first_point = self.vertices[-1]
            else:
                first_point = self.vertices[i - 1]
            self.segments.append(Segment(first_point, self.vertices[i]))

    def union(self, other):
        """ Returns the union of two polygons. """
        return Region([self]).union(Region([other])).polygons

    def intersection(self, other):
        """ Finds the intersection of this polygon with another.  Will
            return a list of polygons. """
        # (I) First, check for the trivial cases:
        # (I.1) this polygon is contained in the other.
        self_in_other = True
        for vertex in self.vertices:
            if not other.contains(vertex):
                self_in_other = False
                break
        if self_in_other:
            return [self]
        # (I.2) the other polygon is contained in this.
        other_in_self = True
        for vertex in other.vertices:
            if not self.contains(vertex):
                other_in_self = False
                break
        if other_in_self:
            return [other]
        # (II) Now create copies of this and the other polygon, adding the
        # intersection points.
        # (II.1) First, find the intersection points.
        intersection_points = set()
        self_new_points = {}
        other_new_points = {}
        for i in xrange(len(self.segments)):
            self_new_points[i] = set()
        for i in xrange(len(other.segments)):
            other_new_points[i] = set()
        for self_seg_i, self_seg in enumerate(self.segments):
            # TODO improve performance.  Bounding box?
            for other_seg_i, other_seg in enumerate(other.segments):
                try:
                    intersection_point = self_seg.intersection(other_seg)
                except SameGreatCircle:
                    for point in (other_seg.point1, other_seg.point2):
                        if not self_seg.contains(point, True):
                            continue
                        if self_seg.is_endpoint(point):
                            continue
                        self_new_points[self_seg_i].add(point)
                    for point in (self_seg.point1, self_seg.point2):
                        if not other_seg.contains(point, True):
                            continue
                        if other_seg.is_endpoint(point):
                            continue
                        other_new_points[other_seg_i].add(point)
                    continue
                if intersection_point is None:
                    continue
                intersection_points.add(intersection_point)
                if not self_seg.is_endpoint(intersection_point):
                    self_new_points[self_seg_i].add(intersection_point)
                if not other_seg.is_endpoint(intersection_point):
                    other_new_points[other_seg_i].add(intersection_point)
        # (II.2) If there are no intersection points, return empty list.
        if not intersection_points:
            return []
        # (II.3) Sort intersection points by order on segments and
        #        create the vertex lists
        self2_vertices = []
        other2_vertices = []
        for i in xrange(len(self.segments)):
            self_new_points[i] = list(self_new_points[i])
            self_new_points[i].sort(lambda x,y:
                    -1 if Segment(self.segments[i].point1, y).contains(x)
                                else 1)
            self2_vertices.append(self.segments[i].point1)
            self2_vertices.extend(self_new_points[i])
        for i in xrange(len(other.segments)):
            other_new_points[i] = list(other_new_points[i])
            other_new_points[i].sort(lambda x,y:
                    -1 if Segment(other.segments[i].point1, y).contains(x)
                                else 1)
            other2_vertices.append(other.segments[i].point1)
            other2_vertices.extend(other_new_points[i])
        # (II.4) Create the point to segment look-up-tables
        point_to_self2 = {}
        point_to_other2 = {}
        for i, point in enumerate(self2_vertices):
            point_to_self2[point] = i
        for i, point in enumerate(other2_vertices):
            point_to_other2[point] = i
        # (III) Now, walk over the edges of this polygon and check at which
        #    intersection points the other polygon actually crosses this one.
        crossing_points = set()
        previously_outside = False
        previously_inside = False
        for i in xrange(len(self2_vertices) + 1):
            i = i % len(self2_vertices)
            if self2_vertices[i] not in intersection_points:
                continue
            next_segment = Segment(self2_vertices[i],
                    self2_vertices[(i + 1) % len(self2_vertices)])
            point = next_segment.get_center()
            # TODO faster lookup with point_to_other2
            if other.border_contains(point):
                continue
            # TODO unify these two cases
            if other.contains(point, assume_not_on_border=True):
                # First case: we're going inside the other polygon
                if previously_inside:
                    continue
                if not previously_outside:
                    previously_inside = True
                    continue
                crossing_points.add(self2_vertices[i])
                previously_inside = True
                previously_outside = False
                continue
            # The other case: we're going outside the other polygon
            if previously_outside:
                continue
            if not previously_inside:
                previously_outside = True
                continue
            crossing_points.add(self2_vertices[i])
            previously_inside = False
            previously_outside = True
        # (IV) Now, extract each polygon of the intersection
        ret = []
        while crossing_points:
            vertices = []
            # (IV.1) Pick a new crossing point we did not consider yet.
            first_point = crossing_points.pop()
            # (IV.2) Find the polygon and direction on that polygon to walk,
            # such that we enter the intersection.
            # TODO unify code for the two cases.
            ok = False
            self_direction = 0
            other_direction = 0
            for direction in (1, -1):
                probe = Segment(self2_vertices[point_to_self2[first_point]],
                        self2_vertices[(point_to_self2[first_point] + direction)
                                % len(self2_vertices)]).get_center()
                if (not other.border_contains(probe)
                        and other.contains(probe, assume_not_on_border=True)):
                    on_self = True
                    ok = True
                    self_direction = direction
                    break
            if not ok:
                for direction in (1, -1):
                    probe = Segment(other2_vertices[point_to_other2[
                                    first_point]],
                            other2_vertices[(point_to_other2[first_point]
                                        + direction)
                                    % len(other2_vertices)]).get_center()
                    if (not self.border_contains(probe)
                            and self.contains(probe, assume_not_on_border=True)):
                        on_self = False
                        other_direction = direction
                        ok = True
                        break
            assert ok
            previous_point = first_point
            # (IV.3) Now, walk.
            # TODO beautify code a bit
            while True:
                vertices.append(previous_point)
                # TODO improve performance
                if on_self:
                    if self_direction == 0:
                        # TODO improve performance
                        ok = False
                        for self_direction in (1, -1):
                            probe = Segment(self2_vertices[point_to_self2[
                                            previous_point]], self2_vertices[(
                                    point_to_self2[previous_point]
                                                + self_direction)
                                            % len(self2_vertices)]).get_center()
                            if (not other.border_contains(probe)
                                    and other.contains(probe,
                                        assume_not_on_border=True)):
                                ok = True
                                break
                        assert ok
                    next_point = self2_vertices[(point_to_self2[previous_point]
                                    + self_direction) % len(self2_vertices)]
                else:
                    if other_direction == 0:
                        ok = False
                        for other_direction in (1, -1):
                            probe = Segment(other2_vertices[point_to_other2[
                                            previous_point]], other2_vertices[(
                                    point_to_other2[previous_point]
                                                + other_direction)
                                            % len(other2_vertices)]).get_center()
                            if (not self.border_contains(probe)
                                    and self.contains(probe,
                                        assume_not_on_border=True)):
                                ok = True
                                break
                        assert ok
                    next_point = other2_vertices[(point_to_other2[
                                            previous_point] + other_direction)
                                            % len(other2_vertices)]
                if next_point == first_point:
                    # Done!
                    break
                if next_point in crossing_points:
                    crossing_points.remove(next_point)
                    on_self = not on_self
                previous_point = next_point
            ret.append(Polygon(vertices, self.external_point))
        return ret

    def border_contains(self, point):
        """ Checks whether the border of the polygon contains the given
            point. """
        for segment in self.segments:
            if segment.contains(point):
                return True
        return False

    def contains(self, point, assume_not_on_border=False):
        """ Checks whether the polygon contains the given point. """
        # TODO prevent iteration over every segment.  Bounding box?
        # First check whether the point is on one of the segments.
        if not assume_not_on_border and self.border_contains(point):
            return True
        # Or is the fixed external point.
        if point == self.external_point:
            return False
        # Now we consider the segment from the point to the fixed external
        # point.  Then we count the proper intersections with the polygon.  If
        # this count is odd, then the point is inside the polygon.  Otherwise
        # it is outside.
        if point != -self.external_point:
            ray = Segment(point, self.external_point)
        else:
            # If the point and the external point are antipodal, then we
            # pick any other point, check whether it is an external point
            # and use it to determine the containment.
            N = 2
            while True:
                N *= 2
                other_point = Point(random.randint(1, N),
                                    random.randint(0, N),
                                    random.randint(0, N))
                if other_point.collinear(self.external_point):
                    continue
                ok = True
                for segment in self.segments:
                    if segment.contains(other_point):
                        ok = False
                        break
                if not ok:
                    continue
                break
            other_point_is_external = not self.contains(other_point)
            other_polygon = Polygon(self.vertices, other_point)
            other_polygon_contains_point = other_polygon.contains(point)
            if other_point_is_external:
                return other_polygon_contains_point
            return not other_polygon_contains_point
        intersection_count = 0
        for i, segment in enumerate(self.segments):
            try:
                point = ray.intersection(segment)
            except SameGreatCircle:
                continue
            if point is None:
                continue
            if not segment.is_endpoint(point):
                intersection_count += 1
                continue
            # The ray crosses and endpoint of the segment.  We need to check
            # whether the polygon actually crosses the ray or bounces off.
            # At least two segments are involved: the incoming, possibly
            # some intermediate, and an outgoing.  We will consider this
            # crossing at the incoming segment and ignore the others.
            if segment.point2 != point:
                continue
            incoming_direction = sign(scalar_triple_product(ray.point1,
                                    ray.point2, segment.point1))
            j = (i + 1) % len(self.segments)
            while True:
                next_segment = self.segments[j]
                if not ray.contains(next_segment.point2):
                    outgoing_direction = sign(scalar_triple_product(ray.point1,
                                            ray.point2, next_segment.point2))
                    break
                j = (i + 1) % len(self.segments)
            # If the outgoing and incoming directions are the same, then
            # this was not a proper crossing.  Ignore it.
            if outgoing_direction == incoming_direction:
                continue
            intersection_count += 1
        return intersection_count % 2 != 0

    def find_internal_point(self):
        """ Tries to find an internal point. """
        # TODO is there a more efficient method?
        for i in xrange(len(self.vertices)):
            j = (i + 1) % len(self.vertices)
            k = (i + 2) % len(self.vertices)
            center = Segment(Segment(self.vertices[i],
                                self.vertices[j]).get_center(),
                    Segment(self.vertices[j],
                            self.vertices[k]).get_center()).get_center()
            for point in (center, -center):
                if self.border_contains(point):
                    continue
                if self.contains(point, assume_not_on_border=True):
                    return point
        raise DegeneratePolygon

    def complement(self):
        """ Returns the complement of the polygon. """
        internal_point = self.find_internal_point()
        return Polygon(self.vertices, internal_point)

    def open_split(self):
        """ Split this polygon into two polygons with overlap such that
            the union of the interior of the two polygons equals the
            interior of this polyogn. """
        # TODO be smarter about splitting: preferably we find the split
        # such that the area is roughly split in two.
        lengths = [(i, self.segments[i].get_tunnel_length()) for i
                    in xrange(len(self.segments))]
        lengths.sort(key=lambda x: -x[1])
        i = lengths[0][0]
        #try:
        return self._try_open_split_at(i, (i + len(self.segments)/2)
                                                    % len(self.segments))
        #except CouldNotSplitThere:
        #    pass
        for i in xrange(len(self.segments)):
            try:
                return self._try_open_split_at(i, (i+1) % len(self.segments))
            except CouldNotSplitThere:
                pass

    def _try_open_split_at(self, i1, i2):
        if i2 < i1:
            i2, i1 = i1, i2
        segment1 = self.segments[i1]
        segment2 = self.segments[i2]
        p1 = segment1.point_in_between(49, 100)
        p2 = segment1.point_in_between(51, 100)
        p3 = segment2.point_in_between(49, 100)
        p4 = segment2.point_in_between(51, 100)
        new_segment1 = Segment(p1, p4)
        new_segment2 = Segment(p2, p3)
        # Now check whether the split is ok.
        # TODO increase performance
        if new_segment1.intersection(new_segment2) is not None:
            raise CouldNotSplitThere
        for i in xrange(len(self.segments)):
            if i in (i1, i2):
                continue
            if new_segment1.intersection(self.segments[i]) is not None:
                raise CouldNotSplitThere
            if new_segment2.intersection(self.segments[i]) is not None:
                raise CouldNotSplitThere
        return [Polygon(self.vertices[i1:i2] + [p3, p2],
                            self.external_point),
                Polygon(self.vertices[:i1] + [p1, p4] + self.vertices[i2:],
                            self.external_point)]

    def __repr__(self):
        return "<sphere.Polygon %s without %s>" % (
                    self.vertices, self.external_point)
    def __eq__(self, other):
        for point in self.vertices:
            if not other.contains(point):
                return False
        for point in other.vertices:
            if not self.contains(point):
                return False
        return True
    def __ne__(self, other):
        for point in self.vertices:
            if not other.contains(point):
                return True
        for point in other.vertices:
            if not self.contains(point):
                return True
        return False

class Segment:
    """ A segment of a Great Circle with arc less than pi radians.  """
    def __init__(self, point1, point2):
        """ Creates a segment by its endpoints. """
        self.point1 = point1
        self.point2 = point2

    def get_tunnel_length(self):
        """ Returns the distance between the endpoints. """
        return self.point1.distance_to(self.point2)

    def get_arc(self):
        """ Returns the arc in radians covered by the segment. """
        return 2 * math.asin(.5 * self.point1.distance_to(self.point2))
    def get_greatcircle(self):
        """ Returns the great circle of which this arc is a part. """
        return GreatCircle.through(self.point1, self.point2)
    def is_endpoint(self, point):
        return point == self.point1 or point == self.point2
    def other_endpoint(self, endpoint):
        """ Given an endpoint of the segment, returns the other one. """
        if endpoint == self.point1:
            return self.point2
        assert endpoint == self.point2
        return self.point1
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
        return (self.point1.squared_distance_to(point) +
                self.point2.squared_distance_to(point)
                    <= self.point1.squared_distance_to(self.point2))
    def intersection(self, other):
        """ Returns intersection-point of two segments or None if
            the segments do not intersect. """
        c1 = self.get_greatcircle()
        c2 = other.get_greatcircle()
        p = c1.intersection(c2)
        if self.contains(p) and other.contains(p):
            return p
        p = -p
        if self.contains(p) and other.contains(p):
            return p
        return None

    def points_in_between(self, step):
        circle = self.get_greatcircle()
        arc1 = circle.angle_of(self.point1)
        arc2 = circle.angle_of(self.point2)
        diff = (arc2 - arc1) % (2 * math.pi)
        opposite = False
        if diff > math.pi:
            diff = 2 * math.pi - diff
            opposite = True
        steps = int(math.floor(diff / step))
        arc = arc1
        ret = []
        for i in xrange(steps):
            if opposite:
                arc -= step
            else:
                arc += step
            ret.append(circle.point_at(arc))
        return ret

    def point_in_between(self, num, denom):
        """ Gets the point that is num/denom way between the two endpoints. """
        x1, y1, z1 = self.point1._get_normalized_tuple()
        x2, y2, z2 = self.point2._get_normalized_tuple()
        return Point(num * x1 + (denom - num) * x2,
                     num * y1 + (denom - num) * y2,
                     num * z1 + (denom - num) * z2)

    def get_center(self):
        """ Returns the center point of the segment. """
        x1, y1, z1 = self.point1._get_normalized_tuple()
        x2, y2, z2 = self.point2._get_normalized_tuple()
        return Point(x1 + x2, y1 + y2, z1 + z2)

    def northern_orthocomplement(self):
        """ Returns a list of polygons containing exactly the points on the
            northern hemisphere orthogonal to any point of the segment.
            Assumes the segment is contained in the northern hemisphere. """
        external_point = self.get_center()
        # TODO drop assumption
        # We distinguish three cases:
        #   I) One of the endpoints is the north pole.
        #  II) Both endpoints are not the north pole and the
        #      orthocomplement (circles) of the two points intersect on
        #      the equator.
        # III) Both endpoints are not the borth pole and the orthocomplement
        #      (circles) of the two points do not intersect on the equator.
        if self.point1 == north_pole or self.point2 == north_pole:
            # (Case I)
            point = self.point2 if self.point1 == north_pole else self.point1
            point_perp = GreatCircle(point) # points orth. to point
            p = point_perp.intersection(equator)
            c = GreatCircle.through(point, north_pole)
            q = c.intersection(point_perp)
            if not q.on_northern_hemisphere:
                q = -q
            r = c.intersection(equator)
            if not Segment(r, point).contains(north_pole,
                    assume_on_greatcircle=True):
                r = -r
            return [Polygon([q, r, p], external_point),
                    Polygon([q, r, -p], external_point)]
        point1_perp = GreatCircle(self.point1)
        point2_perp = GreatCircle(self.point2)
        r = point1_perp.intersection(point2_perp)
        if r.on_equator:
            # (Case II)
            c = GreatCircle.through(self.point1, self.point2)
            q1 = c.intersection(point1_perp)
            q2 = c.intersection(point2_perp)
            if not q1.on_northern_hemisphere:
                q1 = -q1
            if not q2.on_northern_hemisphere:
                q2 = -q2
            return [Polygon([q1, q2, r], external_point),
                    Polygon([q1, q2, -r], external_point)]
        # (Case III)
        q = point1_perp.intersection(equator)
        p = point2_perp.intersection(equator)
        t = GreatCircle(Segment(p, q).get_center()).intersection(
                    GreatCircle.through(self.point1, self.point2))
        if (self.contains(t, assume_on_greatcircle=True) or
                self.contains(-t, assume_on_greatcircle=True)):
            return [Polygon([p, q, r], external_point),
                    Polygon([-p, -q, r], external_point)]
        return [Polygon([p, -q, r], external_point),
                Polygon([-p, q, r], external_point)]

    def __eq__(self, other):
        return ((self.point1, self.point2) == (other.point1, other.point2) or
                (self.point1, self.point2) == (other.point2, other.point1))
    def __ne__(self, other):
        return ((self.point1, self.point2) != (other.point1, other.point2) and
                (self.point1, self.point2) != (other.point2, other.point1))
    def __repr__(self):
        return "<sphere.Segment %s %s>" % (self.point1, self.point2)

class GreatCircle:
    """ A circle that splits the sphere in two.  Equivalently: the
        intersection of a plane with the sphere. """
    def __init__(self, normal):
        """ Creates a great circle given by the normal point of the
            plane it determines. """
        self.normal = normal

        # Linear transformations from greatcircle to the sphere and back.
        # Created by _generate_matrices.
        self._from_circle = None
        self._to_circle = None
    @staticmethod
    def through(point1, point2):
        """ Returns the great circle that goes through the given
            points. """
        return GreatCircle(cross_product(point1, point2))
    def contains(self, point):
        """ Checks whether point is contained in this great circle. """
        return self.normal.orthogonal_to(point)
    def intersection(self, other):
        """ Returns one of the two intersection point of this and the
            other great circle. """
        if self == other:
            raise SameGreatCircle
        return cross_product(self.normal, other.normal)
    def _generate_matrices(self):
        if self._from_circle is not None:
            return
        # First, normalize the great circle normal
        normal = self.normal
        if normal.z < 0:
            normal = -normal
        if normal.z == 0:
            if normal.x < 0:
                normal = -normal
        # Now, find the "canonical" orthogonal basis
        basis = orthonormal_basis_for(normal)
        self._from_circle = numpy.matrix([
                v._get_normalized_tuple() for v in basis]).T
        self._to_circle = numpy.linalg.inv(self._from_circle)
    def angle_of(self, point):
        """ Gets the angle of a point on this Great Circle. """
        self._generate_matrices()
        vec = self._to_circle * numpy.matrix(point._get_normalized_tuple()).T
        return math.atan2(vec[1], vec[2])
    def point_at(self, angle):
        """ Return the point on this great circle at the given angle. """
        self._generate_matrices()
        return Point(*map(float, self._from_circle * numpy.matrix(
                            [0, math.sin(angle), math.cos(angle)]).T))

    def get_quadrants(self):
        """ Returns four segments that make up this great circle. """
        return [Segment(self.point_at(0), self.point_at(math.pi / 2)),
                Segment(self.point_at(math.pi / 2), self.point_at(math.pi)),
                Segment(self.point_at(math.pi), self.point_at(3.0/2 * math.pi)),
                Segment(self.point_at(3.0/2 * math.pi),
                        self.point_at(2 * math.pi))]
                    

    def __repr__(self):
        return "<sphere.GreatCircle %s>" % self.normal

    def __eq__(self, other):
        return self.normal == other.normal or -self.normal == other.normal
    def __ne__(self, other):
        return self.normal != other.normal and -self.normal != other.normal

class Point:
    """ A point on the unit sphere. """
    def __init__(self, x, y, z):
        """ Creates a point by its unnormalized 3D coordinates """
        assert x != 0 or y != 0 or z != 0
        self.x, self.y, self.z = x, y, z

    def __repr__(self):
        return "<sphere.Point long %.1f lat %.1f>" % (
                        self.longitude / math.pi * 180,
                        self.latitude / math.pi * 180)

    def _get_norm(self):
        """ Returns the norm of the internal representative.
            Note: the actual norm is always 1. """
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def _get_normalized_tuple(self):
        """ Gets cartesian coordinates of the Point. """
        norm = self._get_norm()
        # TODO Is there a way to do this without using isinstance?
        if isinstance(self.x, fractions.Fraction):
            norm = fractions.Fraction(norm)
        return (self.x / norm, self.y / norm, self.z / norm)

    def squared_distance_to(self, other):
        """ Returns the square of the distance from this point to
            the given other point. """
        norm_self = self._get_norm()
        norm_other = other._get_norm()
        return ((self.x / norm_self - other.x / norm_other)**2 +
                (self.y / norm_self - other.y / norm_other)**2 +
                (self.z / norm_self - other.z / norm_other)**2)
    def distance_to(self, other):
        """ Returns the distance to another point. """
        return math.sqrt(self.squared_distance_to(other))
    def orthogonal_to(self, other):
        """ Returns whether this point is orthogonal to other. """
        scalar = scalar_product(self, other)
        # TODO Is there a way to do this without isinstance?
        if isinstance(scalar, float):
            return abs(scalar) < 1e-10
        return scalar == 0
    def collinear(self, other):
        """ Checks whether this point is on the same line through the
            origin as the other point.  Or equivalently: whether this
            point is equal or antipodal to the other. """
        return (self.x * other.y == self.y * other.x and
                self.y * other.z == self.z * other.y and
                self.x * other.z == self.z * other.x)
    @property
    def longitude(self):
        """ Returns the longitude of the point. """
        if self.x == 0:
            if self.y > 0:
                return .5 * math.pi
            if self.y == 0:
                return 0.0
            return -.5 * math.pi
        return math.atan(float(self.y) / self.x)
    @property
    def latitude(self):
        """ Return the latitude of the point. """
        x, y, z = self._get_normalized_tuple()
        return math.asin(z)
    @property
    def on_northern_hemisphere(self):
        """ Returns whether the point is on the northern hemisphere. """
        return self.z >= 0
    @property
    def on_equator(self):
        """ Returns whether the point is on the equator. """
        return self.z == 0
    def __eq__(self, other):
        if not isinstance(other, Point):
            return False
        return (self.x * other.y == self.y * other.x and
                self.y * other.z == self.z * other.y and
                self.x * other.z == self.z * other.x and
                (self.x * other.x > 0 or
                 self.y * other.y > 0 or
                 self.z * other.z > 0))
    def __ne__(self, other):
        if not isinstance(other, Point):
            return True
        return ((self.x * other.y != self.y * other.x or
                 self.y * other.z != self.z * other.y or
                 self.x * other.z != self.z * other.x) or
                (self.x * other.x < 0 or
                 self.y * other.y < 0 or
                 self.z * other.z < 0))
    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)
    def __hash__(self):
        return hash(self._get_normalized_tuple())

def scalar_product(point1, point2):
    """ Returns the (unnormalized) scalar product of two points. """
    return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z

def cross_product(point1, point2):
    """ Returns the cross product of point1 and point2. """
    return Point(point1.y * point2.z - point1.z * point2.y,
                 point1.z * point2.x - point1.x * point2.z,
                 point1.x * point2.y - point1.y * point2.x)

def scalar_triple_product(point1, point2, point3):
    """ Returns the scalar triple product of the three points. """
    return scalar_product(cross_product(point1, point2),  point3)

def orthonormal_basis_for(point):
    """ Deterministicly returns a orthonormal basis containing the point. """
    # First, find a basis.
    if not c_xy.contains(point):
        basis = [point, e_x, e_y]
    elif not c_yz.contains(point):
        basis = [point, e_y, e_z]
    else:
        basis = [point, e_x, e_z]
    # Now, orthogonalize
    return grahm_schmidt(basis)

def grahm_schmidt(basis):
    """ Returns basis into an orthonormal basis. """
    ret = []
    for v in basis:
        x, y, z = v.x, v.y, v.z
        for orth_vector in ret:
            factor = (scalar_product(v, orth_vector) * scalar_product(v, v)
                            / scalar_product(orth_vector, orth_vector))
            x -= factor * orth_vector.x
            y -= factor * orth_vector.y
            z -= factor * orth_vector.z
            assert Point(x, y, z).orthogonal_to(orth_vector)
        ret.append(Point(x, y, z))
    return ret

e_x = Point(1, 0, 0)
e_y = Point(0, 1, 0)
e_z = Point(0, 0, 1)

c_xy = GreatCircle.through(e_x, e_y)
c_yz = GreatCircle.through(e_y, e_z)
c_xz = GreatCircle.through(e_x, e_z)

north_pole = Point(0, 0, 1)
south_pole = Point(0, 0, -1)
equator = GreatCircle(north_pole)
eq0 = Point(1, 0, 0)
eq90 = Point(0, 1, 0)
eq180 = Point(-1, 0, 0)
eq270 = Point(0, -1, 0)

quadrant_n1 = Polygon([north_pole, eq0, eq90], south_pole)
quadrant_n2 = Polygon([north_pole, eq90, eq180], south_pole)
quadrant_n3 = Polygon([north_pole, eq180, eq270], south_pole)
quadrant_n4 = Polygon([north_pole, eq270, eq0], south_pole)
quadrant_s1 = Polygon([south_pole, eq0, eq90], north_pole)
quadrant_s2 = Polygon([south_pole, eq90, eq180], north_pole)
quadrant_s3 = Polygon([south_pole, eq180, eq270], north_pole)
quadrant_s4 = Polygon([south_pole, eq270, eq0], north_pole)

northern_hemisphere = Region([quadrant_n1, quadrant_n2,
                              quadrant_n3, quadrant_n4])
southern_hemisphere = Region([quadrant_s1, quadrant_s2,
                              quadrant_s3, quadrant_s4])
full_sphere = Region([quadrant_n1, quadrant_n2, quadrant_n3, quadrant_n4,
                      quadrant_s1, quadrant_s2, quadrant_s3, quadrant_s4])
