import math
import numpy
import scipy.integrate
import scipy.constants
import re
import os.path

# Module to:
#
#  1. Evaluate the cavity transit time factors from the multipole field maps.
#
#  2. Evaluate the transit time factors from the polynomial fits.
#
#  3. To propagate through a cavity.
#
# Author: Johan Bengtsson.


# Constants:
ps_dim = 6
x_ = 0; px_ = 1; y_ = 2; py_ = 3; ct_ = 4; delta_ = 5; orbit_ = 6


# scipy.constants.c    = 3.141592653589793
# scipy.constants.mu_0 = 1.256637061435917e-06
AU = 931.49432 # MeV/u

cav_homs = (
    'EFocus1', 'EDipole', 'HDipole', 'HMono', 'HQuad', 'EQuad', 'EFocus2'
    )


# ------------------------------------------------------------------------------
# Polynomial fits for the cavity transit time factors [T, S] for the multipole
# modes vs. 2pi/beta*lambda.
# ------------------------------------------------------------------------------

def Cav_TS_EFocus1_41(arg):
    p_T = numpy.array([
            -7.316972e+14, 2.613462e+14, -4.112154e+13, 3.739846e+12,
            -2.165867e+11, 8.280687e+09, -2.089452e+08, 3.354464e+06,
            -3.108322e+04, 1.256386e+02])
    p_S = numpy.array([
            -6.079177e+14,  2.229446e+14, -3.605780e+13,  3.374738e+12,
            -2.013750e+11,  7.942886e+09, -2.070369e+08,  3.438044e+06,
            -3.299673e+04,  1.394183e+02])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def Cav_TS_EFocus2_41(arg):
    p_T = numpy.array([
            -1.499544e+11, 5.612073e+10, -9.246033e+09, 8.799404e+08,
            -5.330725e+07, 2.132552e+06, -5.619149e+04, 8.943931e+02,
            -9.121320e+00, 1.038803e+00])
    p_S = numpy.array([
            -1.983302e+10, 8.570757e+09, -1.604935e+09, 1.714580e+08,
            -1.154148e+07, 5.095765e+05, -1.488249e+04, 2.696971e+02,
            -2.585211e+00, 1.305154e-02])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def Cav_TS_EDipole_41(arg):
    p_T = numpy.array([
            4.758398e+10, -1.656906e+10, 2.535541e+09, -2.237287e+08,
            1.255841e+07, -4.669147e+05, 1.125013e+04, -1.047651e+02,
            1.526489e+00, -1.005885e+00])
    p_S = numpy.array([
            1.155597e+11, -4.227114e+10, 6.817810e+09, -6.361033e+08,
            3.782592e+07, -1.488484e+06, 3.888964e+04, -6.407538e+02,
            5.884367e+00, -2.586200e-02])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def Cav_TS_EQuad_41(arg):
    p_T = numpy.array([
            -1.578312e+11, 5.896915e+10, -9.691159e+09, 9.192347e+08,
            -5.544764e+07, 2.206120e+06, -5.779110e+04, 9.127945e+02,
            -9.238897e+00, 1.038941e+00])
    p_S = numpy.array([
            -5.496217e+11, 2.008767e+11, -3.239241e+10, 3.024208e+09,
            -1.801113e+08, 7.094882e+06, -1.848380e+05, 3.069331e+03,
            -2.923507e+01, 1.248096e-01])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def Cav_TS_HMono_41(arg):
    p_T =numpy.array([
            -7.604796e+11, 4.632851e+11, -1.014721e+11, 1.165760e+10,
            -8.043084e+08, 3.518178e+07, -9.843253e+05, 1.697657e+04,
            -1.671357e+02, 1.703336e+00])
    p_S = numpy.array([
            -5.930241e+13, 2.189668e+13, -3.565836e+12, 3.360597e+11,
            -2.019481e+10, 8.022856e+08, -2.106663e+07, 3.524921e+05,
            -3.409550e+03, 1.452657e+01])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def Cav_TS_HDipole_41(arg):
    p_T = numpy.array([
            7.869717e+11, -3.116216e+11, 5.414689e+10, -5.420826e+09,
            3.446369e+08, -1.442888e+07, 3.985674e+05, -7.117391e+03,
            7.075414e+01, 6.853803e-01])
    p_S = numpy.array([
            -4.941947e+12, 1.791634e+12, -2.864139e+11, 2.649289e+10,
            -1.562284e+09, 6.090118e+07, -1.569273e+06, 2.575274e+04,
            -2.441117e+02, 1.021102e+00])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def Cav_TS_HQuad_41(arg):
    p_T = numpy.array([
            5.600545e+12, -2.005326e+12, 3.163675e+11, -2.885455e+10,
            1.676173e+09, -6.429625e+07, 1.627837e+06, -2.613724e+04,
            2.439177e+02, -1.997432e+00])
    p_S = numpy.array([
            1.131390e+13, -4.119861e+12, 6.617859e+11, -6.153570e+10,
            3.649414e+09, -1.431267e+08, 3.711527e+06, -6.135071e+04,
            5.862902e+02, -2.470704e+00])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

cav_transit_times_41 = {
    'EFocus1' : Cav_TS_EFocus1_41,
    'EFocus2' : Cav_TS_EFocus2_41,
    'EDipole' : Cav_TS_EDipole_41,
    'EQuad'   : Cav_TS_EQuad_41,
    'HMono'   : Cav_TS_HMono_41,
    'HDipole' : Cav_TS_HDipole_41,
    'HQuad'   : Cav_TS_HQuad_41
    }

def get_cav_41(cav_hom, f, beta):
    beta_rng = [0.025, 0.08]
    if beta_rng[0] <= beta and beta <= beta_rng[1]:
        lambda_ = scipy.constants.c/f
        [T, S] = cav_transit_times_41[cav_hom](2.0*math.pi/(beta*1e3*lambda_))
    else:
        print 'beta out of range: %5.3f [%5.3f, %5.3f]' % \
            (beta, beta_min, beta_max)
        exit(1)
    return T, S

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Cavity transverse thin lens model.
# ------------------------------------------------------------------------------

def M_EFocus1_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[px_, x_] = Z*V0/(beta**2*gamma*m)*(T*math.cos(phi)-S*math.sin(phi))/aper
    M[py_, y_] = M[px_, x_]
    return M

def M_EFocus2_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[px_, x_] = Z*V0/(beta**2*gamma*m)*(T*math.cos(phi)-S*math.sin(phi))/aper
    M[py_, y_] = M[px_, x_]
    return M

def M_EDipole_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[py_, orbit_] = Z*V0/(beta**2*gamma*m)*(T*math.cos(phi)-S*math.sin(phi))
    return M

def M_EQuad_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[px_, x_] = Z*V0/(beta**2*gamma*m)*(T*math.cos(phi)-S*math.sin(phi))/aper;
    M[py_, y_] = -M[px_, x_]

def M_HMono_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[px_, x_] = \
        -scipy.constants.mu_0*scipy.constants.c*Z*V0/(beta*gamma*m) \
        *(T*math.cos(phi+math.pi/2.0)-S*math.sin(phi+math.pi/2.0))/aper
    M[py_, y_] = M[px_, x_]
    return M

def M_HDipole_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[py_, orbit_] = \
        -scipy.constants.mu_0*scipy.constants.c*Z*V0/(beta*gamma*m) \
        *(T*math.cos(phi+math.pi/2.0)-S*math.sin(phi+math.pi/2.0))
    return M

def M_HQuad_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[px_, x_] = \
        -scipy.constants.mu_0*scipy.constants.c*Z*V0/(beta*gamma*m) \
        *(T*math.cos(phi+math.pi/2.0)-S*math.sin(phi+math.pi/2.0))/aper
    M[py_, y_] = -M[px_, x_]
    return M

def M_AccGap_41(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    return M

M_cav_41 = {
    'EFocus1' : M_EFocus1_41,
    'EFocus2' : M_EFocus2_41,
    'EDipole' : M_EDipole_41,
    'EQuad'   : M_EQuad_41,
    'HMono'   : M_HMono_41,
    'HDipole' : M_HDipole_41,
    'HQuad'   : M_HQuad_41,
    'AccGap'  : M_AccGap_41
    }

def get_cav_M_hom(cav_hom, L, m, Z, phi, aper):
    M = M_cav_41(L, m, Z, V0, T, S, phi, aper)
    arg = 2.0*math.pi/(beta*lambda_)
    M[ct_, delta_]= -2*math.pi/(lamda*beta**3*gamma*ionEs*L)
    M[delta_, ct_]= -Z*V0_1*(T_1*math.sin(phi+arg*L)+math.cos(phi+arg*L))


# ------------------------------------------------------------------------------


def rd_hom(file_name):
    xy = numpy.loadtxt(file_name)
    # z-axis is in [mm].
    return 1e-3*xy[:, 0], xy[:, 1]


def trapezoidal(x, y):
    # Trapezoid rule for numerical integration.
    return ((y[1:]+y[:-1])/2.0*(x[1:]-x[:-1])).sum()


def get_EML(z, EM):
    # Method: 'trapezoidal', 'trapz', 'simps', and 'romb'.
    method = 'trapezoidal'
    if method == 'trapezoidal':
        integ = trapezoidal(z, EM)
    if method == 'trapz':
        integ = numpy.trapz(EM, z)
    elif method == 'simps':
        # Simpson's method requires an odd number of samples.
        n = len(z)
        h = z[1] - z[0]
        if n % 2 == 0:
            z = numpy.append(z, z[n-1]+h)
            EM = numpy.append(EM, 0.0)
        integ = scipy.integrate.simps(EM, z)
    elif method == 'romb':
        # Romberg's method requires 2^k + 1 samples.
        n = len(z)
        k = int(numpy.floor(numpy.log2(n))) + 1
        npad = 2**k + 1 - n
        h = z[1] - z[0]
        for i in range(n, n+npad):
            z = numpy.append(z, z[n-1]+h)
            EM = numpy.append(EM, 0.0)
        integ = scipy.integrate.romb(EM, dx=h)
    return integ


def get_EM_center(z, EM):
    # Compute e-m center along z-axis.
    eml = ((EM[1:]+EM[:-1])/2.0*(z[1:]-z[:-1])).sum()
    em_mom = ((z[1:]+z[:-1])/2.0*(EM[1:]+EM[:-1])/2.0*(z[1:]-z[:-1])).sum()
    em_center = em_mom/eml
    return eml, em_center


def get_transit_time_factors(z, EM, beta, lambda_):
    # Compute transit time factors: [T, T', S, S'].
    [EML, em_center] = get_EM_center(z, numpy.absolute(EM))
    z -= em_center
    coef = 2.0*math.pi/(beta*lambda_)
    T = ((EM[1:]+EM[:-1])/2.0*numpy.cos(coef*(z[1:]+z[:-1])/2.0)
         *(z[1:]-z[:-1])).sum()/EML
    Tp = -((z[1:]+z[:-1])/2.0*(EM[1:]+EM[:-1])/2.0
          *numpy.sin(coef*(z[1:]+z[:-1])/2.0)*(z[1:]-z[:-1])).sum()/EML
    S = ((EM[1:]+EM[:-1])/2.0*numpy.sin(coef*(z[1:]+z[:-1])/2.0)
         *(z[1:]-z[:-1])).sum()/EML
    Sp = ((z[1:]+z[:-1])/2.0*(EM[1:]+EM[:-1])/2.0
          *numpy.cos(coef*(z[1:]+z[:-1])/2.0)*(z[1:]-z[:-1])).sum()/EML
    return T, Tp, S, Sp, EML, em_center


def get_cav(file_name, f, beta):
    # Compute: e-m center, transit time factors [T, T', S, S'], and integrated
    # field.
    [z, EM] = rd_hom(file_name)
    lambda_ = scipy.constants.c/f
    [T, Tp, S, Sp, EML, em_center] = \
        get_transit_time_factors(z, EM, beta, lambda_)
    return [em_center, T, Tp, S, Sp, EML]


def prt_transit_times(file_name, n, f, beta_min, beta_max):
    dbeta = (beta_max-beta_min)/(n-1)
    cav_hom = file_name.split('_')[1]
    outf = open(cav_hom+'.dat', 'w')
    [z, EM] = rd_hom(file_name)
    lambda_ = scipy.constants.c/f
    for k in range(n):
        beta = beta_min + k*dbeta
        [T, Tp, S, Sp, EML, em_center] = \
            get_transit_time_factors(z, EM, beta, lambda_)
        [T_pol, S_pol] = get_cav_41(cav_hom, f, beta)
        outf.write('%8.5f %8.5f %8.5f %8.5f %8.5f\n' % \
                       (beta, T, T_pol, S, S_pol))
    outf.close()


def prt_cav_tlm(file_name, outf, s, f, beta, sgn):
    cav_hom = file_name.split('_')[1]
    [em_center, T, Tp, S, Sp, EML] = get_cav(file_name, f_QWR, beta)
    em_center = sgn*em_center
    L = math.fabs(s-em_center)
    outf.write('%18.15f %-8s %18.15f %18.15f\n' % (s, 'drift', L, 0.0))
    s = em_center
    outf.write('%18.15f %-8s %18.15f %18.15f\n' % (s, cav_hom, 0.0, 1e-6*EML))
    return s


def prt_cav_tlm_41(home_dir, beta):
    # Generate QWR cavity thin lens model.
    # Thin lens kicks are located at the e-m centers.
    f_QWR = 80.5e6
    outf = open('tlm_cav_41.dat', 'w')
    s = -0.120
    outf.write('%18.15f %-8s %18.15f %18.15f\n' % (s, 'marker', 0.0, 0.0))
    # 1st gap.
    s = prt_cav_tlm(home_dir+'CaviMlp_EFocus2_41.txt', outf, s, f_QWR, beta, -1)
    s = prt_cav_tlm(home_dir+'CaviMlp_EDipole_41.txt', outf, s, f_QWR, beta, -1)
    s = prt_cav_tlm(home_dir+'CaviMlp_HDipole_41.txt', outf, s, f_QWR, beta, -1)
    s = prt_cav_tlm(home_dir+'CaviMlp_HMono_41.txt',   outf, s, f_QWR, beta, -1)
    s = prt_cav_tlm(home_dir+'CaviMlp_HQuad_41.txt',   outf, s, f_QWR, beta, -1)
    s = prt_cav_tlm(home_dir+'CaviMlp_EQuad_41.txt',   outf, s, f_QWR, beta, -1)
    s = prt_cav_tlm(home_dir+'CaviMlp_EFocus1_41.txt', outf, s, f_QWR, beta, -1)
    # 2nd gap.
    s = prt_cav_tlm(home_dir+'CaviMlp_EFocus1_41.txt', outf, s, f_QWR, beta, 1)
    s = prt_cav_tlm(home_dir+'CaviMlp_EQuad_41.txt',   outf, s, f_QWR, beta, 1)
    s = prt_cav_tlm(home_dir+'CaviMlp_HQuad_41.txt',   outf, s, f_QWR, beta, 1)
    s = prt_cav_tlm(home_dir+'CaviMlp_HMono_41.txt',   outf, s, f_QWR, beta, 1)
    s = prt_cav_tlm(home_dir+'CaviMlp_HDipole_41.txt', outf, s, f_QWR, beta, 1)
    s = prt_cav_tlm(home_dir+'CaviMlp_EDipole_41.txt', outf, s, f_QWR, beta, 1)
    s = prt_cav_tlm(home_dir+'CaviMlp_EFocus2_41.txt', outf, s, f_QWR, beta, 1)
    L = 0.120 - s
    outf.write('%18.15f %-8s %18.15f %18.15f\n' % (s, 'drift', L, 0.0))
    s += L
    outf.write('%18.15f %-8s %18.15f %18.15f\n' % (s, 'marker', 0.0, 0.0))
    outf.close()


def prt_get_cav(file_name, f, beta):
    cav_hom = file_name.split('_')[1]
    arg = 2.0*math.pi/(beta*lambda_)
    [em_center, T, Tp, S, Sp, EML] = get_cav(file_name, f_QWR, beta)
    print '%18.15f %18.15f %18.15f %18.15f %18.15f %18.15f %18.15f %s' % \
        (1e-3*arg, 1e3*em_center, T, S, Tp, Sp, 1e-6*EML, cav_hom)


def rd_tst_data(file_name):
    inf = open(file_name, 'r')
    print \
        '         k              EM Center             T                   S' \
        '                  T\'                 S\'                 EML'
    # Skip first line.
    line = inf.readline()
    line = inf.readline().strip('\r\n')
    while line:
        tokens = re.split(r'\s*', line)
        [arg, EM_center, T, S, Tp, Sp, EML, cav_hom] = \
            [float(tokens[0]), float(tokens[1]), float(tokens[2]),
             float(tokens[3]), float(tokens[4]), float(tokens[5]),
             float(tokens[6]), tokens[7]]
        print '%18.15f %18.15f %18.15f %18.15f %18.15f %18.15f %18.15f %s' % \
            (arg, EM_center, T, S, 1e-3*Tp, 1e-3*Sp, EML, cav_hom)
        line = inf.readline().strip('\r\n')


def get_cav_param(home_dir, cav_hom, beta):
    file_name = home_dir+'CaviMlp_'
    [em_center, T, Tp, S, Sp, EML] = get_cav(file_name, f_QWR, beta)
    cav.hom[cav_hom].param['em_center'] = em_center
    cav.hom[cav_hom].param['T']         = T
    cav.hom[cav_hom].param['Tp']        = S
    cav.hom[cav_hom].param['S']         = Tp
    cav.hom[cav_hom].param['Sp']        = Sp
    cav.hom[cav_hom].param['E0']        = EML
    print cav.hom[cav_hom].param
    print cav.hom['EDipole'].param
    return cav

def get_cav1():
    pass


def M_drift(L, m, Z, V0, T, S, phi, aper):
    M = numpy.identity(ps_dim)
    M[x_, px_] = L
    M[y_, py_] = L
    return M


def get_cav_M(L, m, Z, phi, aper):
    M = numpy.identity(ps_dim)
    M = numpy.dot(M, M1)
    return M


def rd_cav_tlm(file_name):
    inf = open(file_name, 'r')
    line = inf.readline().strip('\r\n')
    # Loop until blank line.
    while line:
        if line.startswith('%'):
            # Comment.
            pass
        else:
            tokens = re.split(r'\s*', line)
            [s, type, L, aper, ELM] = \
                [1e-3*float(tokens[0]), tokens[1], 1e-3*float(tokens[3]),
                 float(tokens[4]), float(tokens[5])]
            print '%18.15f %8s %18.15f %18.15f %5.3f' % (s, type, L, ELM, aper)
        line = inf.readline().strip('\r\n')


home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'

# HWR cavity.
f_QWR    = 80.5e6
# Cavity aperture.
aper_QWR = 17e-3
beta     =  0.041

# HWR cavity.
#f_HWR    = 322e6
#aper_HWR = 20e-3
#beta     =  0.29

gamma = 1.0/math.sqrt(1-beta**2)


#cav41 = cavity()
#for cav_hom in cav41.hom:
#    print cav_hom
#get_cav_param(home_dir, 'EFocus1', beta)
#exit(0)

#rd_cav_tlm(home_dir+'thinlenlon_41.txt')

if False:
    # Transit times from polynomial interpolation.
    print
    for cav_hom in cav_homs:
        get_cav_41(cav_hom, f_QWR, beta)

if False:
    # Cross check.
    print
    rd_tst_data(home_dir+'cross_check_41.dat')

    lambda_ = scipy.constants.c/f_QWR
    beta = 2.0*math.pi/(0.050887809949826*1e3*lambda_)
    print '\nbeta = %18.15f' % (beta)

    print
    for cav_hom in cav_homs:
        print home_dir+'CaviMlp_'+cav_hom+'_41.txt'
        prt_get_cav(home_dir+'CaviMlp_'+cav_hom+'_41.txt', f_QWR, beta)

prt_cav_tlm_41(home_dir, beta)

if True:
    for cav_hom in cav_homs:
        prt_transit_times(home_dir+'CaviMlp_'+cav_hom+'_41.txt',
                          25, f_QWR, 0.025, 0.08)
