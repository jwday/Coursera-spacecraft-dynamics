# A collection of all the functions I wrote while working on this course material in no particular order

import sympy as sm
import numpy as np
from numpy import sin, cos, tan, arcsin, arccos, arctan
from sympy.physics.mechanics import *
import itertools

# Coordinate system transformation functions
def findEquivalentRotationAngles(angles, first_rot_order, second_rot_order, test=False):
    t1, t2, t3 = sm.symbols('theta_1, theta_2, theta_3')   # Rotation angles for (first_rot_order)

    N = ReferenceFrame('N')
    B = ReferenceFrame('B')

    B.orient_body_fixed(N, (t1, t2, t3), first_rot_order)
    C = B.dcm(N)    # DCM from N to B

    syms, sym_eqns = constructEquivAngleFormulas(C, second_rot_order)
    alpha_eqn, beta_eqn, gamma_eqn, beta_check_eqn = sym_eqns

    a, b, g = syms

    t_vals = {t1: np.radians(angles[0]),
              t2: np.radians(angles[1]),
              t3: np.radians(angles[2])}
    equiv_a_angle = sm.solve(alpha_eqn.subs(t_vals), a, dict=True)[0]
    equiv_g_angle = sm.solve(gamma_eqn.subs(t_vals), g, dict=True)[0]
    equiv_b_angle = sm.solve(beta_eqn.subs(t_vals), b, dict=True)[0]
    equiv_bc_angle = sm.solve(beta_check_eqn.subs(t_vals).subs(equiv_a_angle).subs(equiv_g_angle), b, dict=True)[0]

    equiv_angles = equiv_a_angle | equiv_bc_angle | equiv_g_angle

    if not test:
        return equiv_angles
    else:
        # Test to see if the resultant rotation matrices are identical (to within 1E-10 precision)
        Bp = ReferenceFrame('Bp')
        Bp.orient_body_fixed(N, (a, b, g), second_rot_order)
        Cp = Bp.dcm(N)
        return (np.array(C.subs(t_vals)).astype(np.float64).round(10) == np.array(Cp.subs(equiv_angles)).astype(np.float64).round(10)).all()


 
def constructEquivAngleFormulas(C, newRotOrder):
    a, b, g = sm.symbols('alpha beta gamma')             # Rotation angles for (second_rot_order) set

    N = ReferenceFrame('N')
    Bp = ReferenceFrame('Bp')   # 'Dummy' B frame used to examine DCM and extract inverse transformation formulas

    Bp.orient_body_fixed(N, (a, b, g), newRotOrder)
    Cp = Bp.dcm(N)  # DCM from N to Bp

    # Construct the proper Euler angle formulas using the rotation order to inform which indicies of the DCM to call
    amounts, rot_order, rot_matrices = Bp._parse_consecutive_rotations((a, b, g), newRotOrder)
    if rot_order[0] == rot_order[2]:
        # Proper Euler angles
        idx1 = rot_order[1]-1
        idx2 = rot_order[0]-1
        idx3 = 5 - (rot_order[0] + rot_order[1])
        
        a_idxs = [(idx2, idx1), (idx2, idx3)]
        b_idxs = (idx2, idx2)
        g_idxs = [(idx1, idx2), (idx3, idx2)]

    else:
        # Tait-Bryan angles
        idx1 = rot_order[2]-1
        idx2 = rot_order[0]-1
        idx3 = 5 - (rot_order[0] + rot_order[2])

        a_idxs = [(idx1, idx3), (idx1, idx1)]
        b_idxs = (idx1, idx2)
        g_idxs = [(idx3, idx2), (idx2, idx2)]

    alpha_eqn = sm.Eq((C[a_idxs[0]]/C[a_idxs[1]]), (Cp[a_idxs[0]]/Cp[a_idxs[1]]).trigsimp())
    beta_eqn = sm.Eq(C[b_idxs], Cp[b_idxs])
    gamma_eqn = sm.Eq(C[g_idxs[0]]/C[g_idxs[1]], (Cp[g_idxs[0]]/Cp[g_idxs[1]]).trigsimp())
    beta_check_eqn = sm.Eq((C[a_idxs[1]]/C[b_idxs]), (Cp[a_idxs[1]]/Cp[b_idxs]).trigsimp())

    syms = [a, b, g]
    sym_eqns = [alpha_eqn, beta_eqn, gamma_eqn, beta_check_eqn]
    return syms, sym_eqns


def test_all_rotation_combinations():
    legal_rots = []

    for i in itertools.product([1,2,3], [1,2,3], [1,2,3]):
        if i[0] != i[1] and i[1] != i[2]:
            rot_str = '{}{}{}'.format(i[0],i[1],i[2])
            legal_rots.append(rot_str)

    for i in itertools.product(legal_rots, legal_rots):
        if i[0] != i[1]:
            print('{} to {}: {}'.format(i[0], i[1], findEquivalentRotationAngles((20,10,-10), i[0], i[1], test=True)))



# Spherical trigonometric functions
def sphlos(x,y,z):
    '''
    Spherical Law of Sines
    Useable for any general combination of either:
        QTY 2 sides and QTY 1 angle to find another angle, or...
        QTY 1 side and QTY 2 angles and to find another side
        x and z must be the same type of measure (either sides or angles)
        
        ex. If you are seeking angle A, then input sphlos(a,B,b)
        ex. If you are seeking arc a, then input sphlos(A,b,B)
    '''
    return arcsin(sin(x)*sin(y)/sin(z))

def sphloc_angle(x,y,z):
    '''
    Spherical Law of Cosines, formulated to return an angle
        ex. If you are seeking angle A opposite side a, then input sphloc_angle(a,b,c)
    '''
    return arccos(cos(x) - cos(y)*cos(z)/(sin(y)*sin(z)))

def sphloc_arc(x,y,z):
    '''
    Spherical Law of Cosines, formulated to return a side
        ex. If you are seeking side a opposite Angle A, then input sphloc_angle(A,b,c)
    '''
    return arccos(-cos(y)*cos(z) + sin(y)*sin(z)*cos(x))


def addRots(R1,R2):
    '''
    Add two consecutive symmetric (X-Y-X) rotations into one single (X-Y-X) rotation
    First rotation is R1 = (a,b,c)
    Second rotation is R2 = (d,e,f)
    Output is equivalent final rotation RF
    '''
    a, b, c = R1[0], R1[1], R1[2]
    d, e, f = R2[0], R2[1], R2[2]
    psi2 = arccos(cos(b)*cos(e)-sin(b)*sin(e)*cos*(c+d))
    psi1 = a + arctan(sin(b)*sin(e)*sin(c+d)/(cos(e)-cos(b)*cos(psi2)))
    psi3 = f + arctan(sin(b)*sin(e)*sin(c+d)/(cos(b)-cos(e)*cos(psi2)))
    RF = (psi1, psi2, psi3)
    return RF

def subRots(R1,RF):
    '''
    Subtract one of two consecutive symmetric (X-Y-X) rotations from one single (X-Y-X) rotation
    First rotation is R1 = (a,b,c)
    Equivalent final rotation is RF = (d,e,f)
    Output is second rotation R2
    '''
    a, b, c = R1[0], R1[1], R1[2]
    d, e, f = RF[0], RF[1], RF[2]
    phi2 = arccos(cos(b)*cos(e)+sin(b)*sin(e)*cos*(d-a))
    phi1 = -c + arctan(sin(b)*sin(e)*sin(d-a)/(cos(b)*cos(phi2)-cos(e)))
    phi3 = f - arctan(sin(b)*sin(e)*sin(d-a)/(cos(b)-cos(phi2)*cos(e)))
    R2 = (phi1, phi2, phi3)
    return R2


# DCMs and Principle Rotation Vectors
def getPRVfromDCM(C):
    '''
    Return the Principle Rotation Vector (PRV) in terms of rotation angle (phi) and unit vector (e) from a DCM input
    It should provide rotation angles between -pi and +pi and exclusively about the + rotation axis
    It should also be capable of handling singular rotations (i.e. rotations of 0, +/-90, or +/-180 about one of the principle axes)
    '''
    foo = (1/2)*(np.trace(C) - 1)
    # foo = 1 corresponds to a 0 deg rotation about any of the three principle axes
    # foo = 0 corresponds to either a 90 or -90 deg rotation about any of the three principle axes
    # foo = -1 corresponds to either a 180 or -180 deg rotation about any of the three principle axes
    u   = np.array([C[1,2]-C[2,1],
                    C[2,0]-C[0,2],
                    C[0,1]-C[1,0]])
    with np.errstate(all='raise'):
        try:
            n = abs(u/np.linalg.norm(u))
        except:
            # If it's a singular case, handle it by finding the eigenvectors and looking for the one that corresponds to an eigenvalue of 1. That's your rotation vector.
            eigs = np.linalg.eig(C)
            unique, counts = np.unique(eigs.eigenvalues, return_counts=True)
            d = dict(zip(unique,counts))
            if 1.0 in d and d[1.0] == 1:
                idx = np.where(eigs.eigenvalues == 1.0)[0][0]
                n = eigs.eigenvectors[idx]
            else:
                raise ValueError("Error: The QTY of eigenvectors with eigenvalues equal to unity is not exactly 1. Is your rotation matrix Identity or not right-handed?")
            
    Kn = np.array([[0, -n[2], n[1]],
                   [n[2], 0, -n[0]],
                   [-n[1], n[0], 0]])
    
    phis = np.arcsin(-np.trace(np.matmul(Kn,C))/2)
    phic = np.arccos(foo)
    
    if phis >= 0:
        phi = phic
    else:
        phi = -phic

    e = n
    return e, phi


def getPRVfromEuler(angles, rotOrder):
    '''
    Return the Principle Rotation Vector (PRV) in terms of rotation angle (phi) and unit vector (e) from Euler angle inputs (including rotation order)
    This function leverages the sympy.physics.mechanics package to shortcut some of the steps
    '''
    a, b, c = angles
    psi, theta, phi = sm.symbols('psi, theta, phi')
    N = ReferenceFrame('N', indices=('1', '2', '3'))
    B = ReferenceFrame('B', indices=('1', '2', '3'))
    B.orient_body_fixed(N, (psi, theta, phi), rotOrder)
    C = B.dcm(N).subs({psi: np.radians(a),
                       theta: np.radians(b),
                       phi: np.radians(c)})
    return getPRVfromDCM(np.array(C).astype(np.float64))


def addPRVs(phi1,e1,phi2,e2):
    '''
    Add two PRVs together
    '''
    sp1 = np.sin(phi1/2)
    sp2 = np.sin(phi2/2)
    cp1 = np.cos(phi1/2)
    cp2 = np.cos(phi2/2)
    phi = 2*np.arccos(cp1*cp2 - sp1*sp2*np.dot(e1,e2))
    e   = (cp2*sp1*e1 + cp1*sp2*e2 + sp1*sp2*np.cross(e1,e2))/np.sin(phi/2)
    return phi,e

def subPRVs(phi,e,phi1,e1):
    sp1  = np.sin(phi1/2)
    sp   = np.sin(phi/2)
    cp1  = np.cos(phi1/2)
    cp   = np.cos(phi/2)
    phi2 = 2*np.arccos(cp*cp1 + sp*sp1*np.dot(e,e1))
    e2   = (cp1*sp*e - cp*sp1*e1 + sp*sp1*np.cross(e,e1))/np.sin(phi2/2)
    return phi2,e2

def DCMfromQuaternion(b):
    # Build four matrices then add them together
    # 1. Base diagonal (diagonal of b0**2)
    # 2. Extension diagonal (b1, b2, b3 components of alternating sign)
    # 3. Base off-diagonal (2*br*bc components)
    # 4. Extension off-diagonal (2*b0*bd components)

    # Build the diagonal components
    base_diag = np.identity(3)*(b[0]**2)
    ext_diag = np.zeros((3,3))
    ext_diag_row_idx, ext_diag_col_idx = np.diag_indices(3)
    for i in zip(ext_diag_row_idx, ext_diag_col_idx):
        row = i[0]
        col = i[1]
        alt_b_val = b[row+1]**2
        b_group = -(b[1]**2)-(b[2]**2)-(b[3]**2)
        val = b_group + 2*alt_b_val
        ext_diag[row,col] = val

    # Bulld the off-diagonal components
    base_off_diag = np.zeros((3,3))
    ext_off_diag = np.zeros((3,3))
    triu_row_idx, triu_col_idx = np.triu_indices(3, k=1)

    for i in zip(triu_row_idx, triu_col_idx):
        row = i[0]
        col = i[1]
        val1 = 2*b[row+1]*b[col+1]
        val2 = 2*b[0]*b[4-(row+col)]
        base_off_diag[row,col] = val1
        ext_off_diag[row,col] = val2
    base_off_diag = base_off_diag + base_off_diag.T
    ext_off_diag = ext_off_diag + ext_off_diag.T

    sign_array = np.ones((3,3), dtype=int)
    sign_array[::2,::2] = -1
    sign_array[1::2,1::2] = -1
    ext_off_diag = np.multiply(ext_off_diag, sign_array)

    # Add the components together
    return base_diag + base_off_diag + ext_diag + ext_off_diag

def quaternionFromDCM(C):
    '''
    Use Sheppard's method to compute the quaternion parameters from a DCM
    '''
    b0_2 = (1/4)*(1+np.trace(C))
    b1_2 = (1/4)*(1+2*C[0,0]-np.trace(C))
    b2_2 = (1/4)*(1+2*C[1,1]-np.trace(C))
    b3_2 = (1/4)*(1+2*C[2,2]-np.trace(C))

    for i 