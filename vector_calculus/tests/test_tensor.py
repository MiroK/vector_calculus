from vector_calculus.containers import Tensor
from sympy import symbols, S
from numpy import eye, array
import unittest


class TestTensor(unittest.TestCase):
    '''UnitTest of Tensor class.'''

    def test_len(self):
        for i in range(2, 4):
            self.assertEqual(len(Tensor(eye(i))), i)

    def test_add(self):
        A = Tensor([[1, 2], [3, 4]])
        B = Tensor([[1, 0], [0, 1]])
        C = Tensor([[2, 2], [3, 5]])
        self.assertEqual(A+B, C)

    def test_sub(self):
        A = Tensor([[1, 2], [3, 4]])
        B = Tensor([[1, 0], [0, 1]])
        C = Tensor([[0, 2], [3, 3]])
        self.assertEqual(A-B, C)

    def test_div(self):
        A = Tensor([[2, 0], [0, 4]])
        B = Tensor([[1, 0], [0, 2]])
        self.assertEqual(A/2, B)

    def test_mul(self):
        A_ = array([[1, 2], [3, 4]])
        B_ = array([[10, -2], [-3, 34]])
        C_ = A_.dot(B_)

        A = Tensor(A_)
        B = Tensor(B_)
        C = Tensor(C_)
        self.assertEqual(A*B, C)

    def test_subs(self):
        try:
            A = Tensor([[1, 2], [2, 3, 4]])
        except AssertionError:
            self.assertTrue(True)

        x = symbols('x')
        A = Tensor([[x, 0], [0, 0]])
        B = Tensor([[1, 0], [0, 0]])
        self.assertEqual(A.subs({x: 1}), B)


# -----------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()
