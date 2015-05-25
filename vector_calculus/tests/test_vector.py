from vector_calculus.containers import Vector
from sympy import symbols, S
import unittest


class TestVector(unittest.TestCase):
    '''UnitTest of Vector class.'''

    def test_len(self):
        for i in range(2, 4):
            self.assertEqual(len(Vector([1]*i)), i)

    def test_add(self):
        x = symbols('x')
        u = Vector([S(1), S(2) + x])
        v = Vector([S(1), S(2) + x])
        self.assertEqual(u+v, 2*u)

    def test_sub(self):
        x = symbols('x')
        u = Vector([S(1), S(2) + x])
        v = Vector([S(1), S(2) + x])
        self.assertEqual(u-v, Vector([S(0), S(0)]))

    def test_mul(self):
        x, y, z = symbols('x, y, z')
        u = Vector([x, y, z])
        self.assertEqual(-1*u, -u)
        self.assertEqual(u*2, 2*u)
        self.assertEqual(2*u, Vector([2*x, 2*y, 2*z]))

    def test_subs(self):
        x, y, z = symbols('x, y, z')
        u = Vector([x, y, z])
        
        s, t, r = symbols('s, t, r')
        v = Vector([s, t, r])
        self.assertEqual(v.subs({s: x, t: y, r: z}), u)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()
