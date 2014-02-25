""" 3D viewer for objects on the sphere. """

import sphere

import sys
import math
import time
import colorsys

from OpenGL import GL, GLU, GLUT

COLOR_ANGLE = 0.294607050403

def view(objects):
    viewer = Viewer()
    for obj in objects:
        if isinstance(obj, sphere.Segment):
            viewer.add_segment(obj)
        elif isinstance(obj, sphere.Point):
            viewer.add_point(obj)
        elif isinstance(obj, sphere.GreatCircle):
            viewer.add_circle(obj)
        elif isinstance(obj, sphere.Polygon):
            viewer.add_polygon(obj)
        else:
            raise ValueError("I do not how to draw %r" % obj)
    viewer.run()

class Viewer:
    def __init__(self):
        self.points = []
        self.segments = []
        self.circles = []
        self.polygons = []

        self.running = False

        self.zoom = 4.0
        self.cam_long = 0.0
        self.cam_lat = .5 * math.pi

        self.n_colors = 0

    def next_color(self):
        """ Returns the next free color. """
        color_index = self.n_colors
        self.n_colors += 1
        return colorsys.hls_to_rgb((color_index * COLOR_ANGLE) % 1.0, 0.5, 1.0)

    def add_point(self, point):
        self.points.append((point, self.next_color()))
    def add_circle(self, circle):
        self.circles.append((circle, self.next_color()))
    def add_polygon(self, polygon):
        self.polygons.append((polygon, self.next_color()))
    def add_segment(self, segment):
        self.segments.append((segment, self.next_color()))

    def _display_point(self, point, color):
        GL.glBegin(GL.GL_POINTS)
        GL.glColor4f(color[0], color[1], color[2], 1)
        GL.glVertex3f(*point._get_normalized_tuple())
        GL.glEnd()

    def _display_segment(self, segment, color):
        GL.glColor4f(color[0], color[1], color[2], 1)
        GL.glBegin(GL.GL_POINTS)
        GL.glVertex3f(*segment.point1._get_normalized_tuple())
        GL.glVertex3f(*segment.point2._get_normalized_tuple())
        GL.glEnd()
        GL.glBegin(GL.GL_LINE_STRIP)
        GL.glVertex3f(*segment.point1._get_normalized_tuple())
        for point in segment.points_in_between(0.1):
            GL.glVertex3f(*point._get_normalized_tuple())
        GL.glVertex3f(*segment.point2._get_normalized_tuple())
        GL.glEnd()

    def display_segments(self):
        GL.glPointSize(4);
        GL.glLineWidth(2);
        GL.glDisable(GL.GL_LIGHTING)
        for segment, color in self.segments:
            self._display_segment(segment, color)
        GL.glEnable(GL.GL_LIGHTING)

    def display_points(self):
        GL.glPointSize(4);
        GL.glDisable(GL.GL_LIGHTING)
        for point, color in self.points:
            self._display_point(point, color)
        GL.glEnable(GL.GL_LIGHTING)

    def display_sphere(self):
        GL.glPointSize(1);
        GL.glLineWidth(0.5);
        GL.glMaterialfv(GL.GL_FRONT, GL.GL_DIFFUSE, [1, 1, 1, 0.5])
        GLUT.glutWireSphere(0.97, 20, 20)
        GLUT.glutSolidSphere(0.98, 20, 20)

    def display_box(self):
        GL.glLineWidth(2);
        GL.glDisable(GL.GL_LIGHTING)
        GL.glColor4f(1, 0, 0, 1)
        GL.glBegin(GL.GL_LINE_STRIP)
        GL.glVertex3f(-1, -1, -1)
        GL.glVertex3f(1, -1, -1)
        GL.glEnd()
        GL.glColor4f(0, 1, 0, 1)
        GL.glBegin(GL.GL_LINE_STRIP)
        GL.glVertex3f(-1, -1, -1)
        GL.glVertex3f(-1, 1, -1)
        GL.glEnd()
        GL.glColor4f(0, 0, 1, 1)
        GL.glBegin(GL.GL_LINE_STRIP)
        GL.glVertex3f(-1, -1, -1)
        GL.glVertex3f(-1, -1, 1)
        GL.glEnd()
        GL.glEnable(GL.GL_LIGHTING)

    def keyboard(self, key, x, y):
        redraw = False
        if key == 'q':
            self.close()
        elif key == 'w':
            self.cam_lat += 0.1
            redraw = True
        elif key == 's':
            self.cam_lat -= 0.1
            redraw = True
        elif key == 'a':
            self.cam_long -= 0.1
            redraw = True
        elif key == 'd':
            self.cam_long += 0.1
            redraw = True
        elif key == 'r':
            self.zoom += 0.1
            redraw = True
        elif key == 'f':
            self.zoom -= 0.1
            redraw = True
        else:
            print 'unknown key', key
        if redraw:
            GLUT.glutPostRedisplay()

    def display(self):
        """ Render the scene. """
        start = time.time()

        GL.glClear(GL.GL_COLOR_BUFFER_BIT |
                   GL.GL_DEPTH_BUFFER_BIT)
        GL.glPushMatrix()

        cam_x = self.zoom * math.sin(self.cam_lat) * math.cos(self.cam_long)
        cam_y = self.zoom * math.sin(self.cam_lat) * math.sin(self.cam_long)
        cam_z = self.zoom * math.cos(self.cam_lat)
        GLU.gluLookAt(cam_x, cam_y, cam_z, 0, 0, 0, 0, 0, 2)
        print cam_x, cam_y, cam_z

        self.display_box()
        self.display_points()
        self.display_segments()
        self.display_sphere()

        GL.glPopMatrix()
        GLUT.glutSwapBuffers()

        print 'display %.4f' % (time.time() - start)

    def close(self):
        # TODO GLUT.glutLeaveMainLoop does not always exist.  Is there a way
        # to quit the mainloop without exiting the whole program?
        sys.exit()

    def run(self):
        """ Set up the viewer. """
        if self.running:
            raise RuntimeError("Already running")
        self.running = True

        GLUT.glutInit([])
        GLUT.glutInitDisplayMode(GLUT.GLUT_DOUBLE |
                                 GLUT.GLUT_RGBA |
                                 GLUT.GLUT_ALPHA |
                                 GLUT.GLUT_DEPTH)
        GLUT.glutInitWindowSize(400, 400)
        GLUT.glutCreateWindow('sphere.viewer')
        
        GL.glClearColor(1, 1, 1, 1)
        GL.glShadeModel(GL.GL_SMOOTH)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)

        GL.glEnable(GL.GL_LIGHTING)
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, [5, 10, -10, 1])
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_DIFFUSE, [1, 1, 1, 1])
        GL.glLightf(GL.GL_LIGHT0, GL.GL_CONSTANT_ATTENUATION, 0.1)
        GL.glLightf(GL.GL_LIGHT0, GL.GL_LINEAR_ATTENUATION, 0.05)
        GL.glEnable(GL.GL_LIGHT0)

        GL.glMatrixMode(GL.GL_PROJECTION)
        GLU.gluPerspective(40, 1, 1, 40)
        GL.glMatrixMode(GL.GL_MODELVIEW)

        GLUT.glutKeyboardFunc(self.keyboard)
        GLUT.glutDisplayFunc(self.display)
        GLUT.glutWMCloseFunc(self.close)

        GLUT.glutMainLoop()

if __name__ == '__main__':
    view([
        sphere.Segment(sphere.Point(1,0,0),
                        sphere.Point(0,0,1))
        ])
