from vector_calculus.containers import Tensor, Vector
from vector_calculus.operators import *
from sympy import symbol
import unittest
import numpy as np


class TestOperatorLinalg(unittest.TestCase):
    '''UnitTest of operators/linalg functionality.'''

    def test_tr(self):
        A = Tensor([[1, 2], [3, 4]])
        trace = tr(A)
        self.assertAlmostEqual(abs(trace), 5)

    def test_transpose(self):
        A_ = np.array([[1, 2], [3, 4]])
        B_ = A_.T
        A = Tensor(A_)
        self.assertEqual(transpose(A), Tensor(B_))

    def test_sym(self):
        A = sym(Tensor([[1, 2], [3, 4]]))
        v = Vector([3, 4])
        foo = inner(v, dot(A, v)) 
        bar = inner(v, dot(v, A))
        self.assertEqual(foo, bar)

    def test_skew(self):
        A = skew(Tensor([[1, 2], [3, 4]]))
        v = Vector([3, 4])
        foo = inner(v, dot(A, v)) 
        self.assertAlmostEqual(foo, 0)

    def test_deviatoric(self):
        A = Tensor([[1, 2], [3, 4]])
        Ad = deviatoric(A)
        trace = tr(Ad)
        self.assertAlmostEqual(abs(trace), 0)

    def test_cross(self):
        u_ = np.array([1, 2, 3])
        v_ = np.array([2, 0, -1])
        w_ = np.cross(u_, v_)

        u = Vector(u_)
        v = Vector(v_)
        w = cross(u, v)
        self.assertEqual(w, Vector(w_))

    def test_outer(self):
        u_ = np.array([1, 2, 3])
        v_ = np.array([2, 0, -1])
        w_ = np.outer(u_, v_)

        u = Vector(u_)
        v = Vector(v_)
        w = outer(u, v)
        self.assertEqual(w, Tensor(w_))

    def test_inner(self):
        # Vector
        u_ = np.array([1, 2, 3])
        v_ = np.array([2, 0, -1])
        i_ = sum(ui*vi for ui, vi in zip(u_, v_))

        u = Vector(u_)
        v = Vector(v_)
        i = inner(u, v)
        self.assertEqual(i, i_)

        # Matrix
        A_ = np.outer(u_, v_)
        i_ = sum(A_.T.dot(A_)[i, i] for i in range(3))

        A = Tensor(A_)
        i = inner(A, A)
        self.assertEqual(i, i_)

    def test_dof(self):
        A_ = np.array([[1, 2], [3, 4]])
        v_ = np.array([3, 4])
        Av_ = A_.dot(v_)
        vA_ = A_.T.dot(v_)

        A = Tensor(A_)
        v = Vector(v_)
        Av = dot(A, v)
        vA = dot(v, A)

        self.assertEqual(Av, Vector(Av_))
        self.assertEqual(vA, Vector(vA_))
