import numpy as np
import sympy as sm
from sympy.physics.vector import *

N = ReferenceFrame('N')
B = ReferenceFrame('B')

NB = sm.Matrix([[-0.87097, -0.19355, 0.45161],
                [0.45161, -0.67742, 0.58065],
                [0.19355, 0.70968, 0.67742]])
BN = NB.transpose()
wBN = 0.1*B.x + 0.2*B.y + 0.3*B.z

B.orient_dcm(N, BN)
B.set_ang_vel(N, wBN)

def matcpop(x, R):
    x1 = x.to_matrix(R)[0]
    x2 = x.to_matrix(R)[1]
    x3 = x.to_matrix(R)[2]
    return sm.Matrix([[0, -x3, x2],
                      [x3, 0, -x1],
                      [-x2, x1, 0]])

dcm_rate = -matcpop(wBN, B)*B.dcm(N)