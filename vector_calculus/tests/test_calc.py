from vector_calculus.containers import Tensor, Vector
from vector_calculus.operators import *
from sympy import symbol, S, sin, cos
import unittest


class TestOperatorCalculus(unittest.TestCase):
    '''UnitTest of operators/calculus functionality.'''

    def test_dx(self):
        x, y, z = symbols('x, y, z')
        f = x**2 + y
        foo = Dx(f, x) == 2*x
        bar = Dx(f, y) == S(1)
        self.assertTrue(foo and bar)

        v = Vector([x, 2*y, 3*z])
        for i in range(1, 4):
            foo = Dx(v, symbols('x, y, z')[i-1])
            bar = Vector([0]*(i-1) + [i] + [0]*(3-i))
            self.assertEqual(foo, bar)

    def test_grad(self):
        # Grad of scalar from definition
        x, y, z = symbols('x, y, z')
        f = 2*x + 3*y**2 - sin(z)
        v = grad(f)
        v_ = Vector([2, 6*y, -cos(z)])
        self.assertEqual(v, v_)

        # Grad of vector from definition
        x, y, z = symbols('x, y, z')
        v = Vector([x**2*y, 5*x + sin(y)])
        g = grad(v)
        g_ = Tensor([[2*x*y, x**2], [5, cos(y)]])
        self.assertEqual(g, g_)

    def test_curl(self):
        # 3d identity curl(grad) = 0
        x, y, z = symbols('x, y, z')
        f = x*y*z
        v = grad(f)
        self.assertEqual(curl(v), Vector([0, 0, 0]))

        # Same 2d identity
        f = x*y
        v = grad(f)
        self.assertEqual(curl(v), S(0))

        # From defition
        v = curl(Vector([x**2*y, x*y*z, -x**2*y**2]))
        v_ = Vector([-2*x**2*y-x*y, 2*x*y**2, y*z - x**2])
        self.assertEqual(v, v_)

    def test_rot(self):
        x, y = symbols('x, y')
        # Rot from definition
        f = x**2 + y**2
        v_ = Vector([-2*y, 2*x])
        v = rot(f)
        self.assertEqual(v, v_)

    def test_div(self):
        x, y, z = symbols('x, y, z')
        # Div curl u = 0
        u = Vector([z-x**4, 2*x+y, x+y+z])
        v = curl(u)
        self.assertEqual(div(v), S(0))

        # Div vector from definition
        v = Vector([x**2*y, x*y*z, -x**2*y**2])
        d = div(v)
        d_ = 2*x*y + x*z
        self.assertEqual(d, d_)

        # Div tensor from defition
        v = Vector([x**2*y, x*y*z, -x**2*y**2])
        self.assertEqual(div(v), tr(grad(v, dim=3)))

        # curl curl = grad div - div grad
        u = Vector([x**2*y*z, x*y*z**2, -x**2*y**2*z])
        v = curl(curl(u)) - grad(div(u), dim=3) + div(grad(u, dim=3))
        v_ = Vector([0, 0, 0])
        self.assertEqual(v, v_)
