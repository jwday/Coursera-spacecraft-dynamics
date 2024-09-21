import sympy as sm
from sympy.physics.vector import *
import numpy as np

N = ReferenceFrame('N')
B = ReferenceFrame('B')
F = ReferenceFrame('F')

b1 = sm.Rational(1,3)*sm.Array([1, 2, -2])
b2 = (1/sm.sqrt(2))*sm.Array([0, 1, 1])
b3 = (1/(3*sm.sqrt(2)))*sm.Array([4, -1, 1])

BN = sm.Matrix([b1, b2, b3])

f1 = sm.Rational(1,4)*sm.Array([3, -2, sm.sqrt(3)])
f2 = sm.Rational(1,2)*sm.Array([-1, 0, sm.sqrt(3)])
f3 = sm.Rational(-1,4)*sm.Array([sm.sqrt(3), 2*sm.sqrt(3), 1])

FN = sm.Matrix([f1, f2, f3])

B.orient_dcm(N, BN)
F.orient_dcm(N, FN)

'''
Notes

It is important to know what form of the direction cosine matrix is returned.
If B.dcm(A) is called, it means the “direction cosine matrix of B rotated relative to A”.
This is the matrix BA shown in the following relationship:

[b1, b2, b3] = BA * [a1, a2, a3]

In other words, B.dcm(A) returns the matrix used to transform a vector to B from A.
Put another way, BA is the matrix that expresses the B unit vectors in terms of the A unit vectors.
'''

print(B.dcm(N))     # DCM to B from N (BN)
print(B.dcm(F))     # DCM to B from F (BF)