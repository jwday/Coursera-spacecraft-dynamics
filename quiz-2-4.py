import numpy as np
import sympy as sm

# adjust the return matrix values as needed
def result():
    BN = np.array([[0.33, 0.66, -0.66],
                    [0, 0.707, 0.707],
                    [0.943, -0.236, 0.236]])
    FN = np.array([[0.75, -0.5, 0.433],
                    [-0.5, 0, 0.866],
                    [-0.433, -0.866, -0.25]])
    BF = np.matmul(BN,np.transpose(FN))
    return BF

def result2():
    b1 = sm.Rational(1,3)*sm.Array([1, 2, -2])
    b2 = (1/sm.sqrt(2))*sm.Array([0, 1, 1])
    b3 = (1/(3*sm.sqrt(2)))*sm.Array([4, -1, 1])

    f1 = sm.Rational(1,4)*sm.Array([3, -2, sm.sqrt(3)])
    f2 = sm.Rational(1,2)*sm.Array([-1, 0, sm.sqrt(3)])
    f3 = sm.Rational(-1,4)*sm.Array([sm.sqrt(3), 2*sm.sqrt(3), 1])

    BN = sm.Matrix([b1, b2, b3])
    FN = sm.Matrix([f1, f2, f3])
    BF = BN*FN.transpose()
    return BF

print(np.array(result2()).astype(np.float64))
print(result())

# BF = [-0.372, -0.744, -0.555,
#         -0.047, 0.612, -0.789,
#         0.927, -0.267, -0.263]

# [[-0.37200847 -0.74401694 -0.55502117]
#  [-0.04736717  0.61237244 -0.78914913]
#  [ 0.92701998 -0.26728038 -0.26304971]]