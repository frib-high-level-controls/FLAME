import math
import numpy
import scipy.integrate
import scipy.constants
import os.path

# Module to evaluate the cavity transit time factors the multipole field maps.
#
# Author: Johan Bengtsson.


# Cavity transit time factors [T, S] for the multipole modes vs.
# 2pi/beta*lambda.


def CaviMlp_EFocus1_41(arg):
    p_T = numpy.array([
            -7.316972e+14, 2.613462e+14, -4.112154e+13, 3.739846e+12,
            -2.165867e+11, 8.280687e+09, -2.089452e+08, 3.354464e+06,
            -3.108322e+04, 1.256386e+02])
    p_S = numpy.array([
            -6.079177e+14,  2.229446e+14, -3.605780e+13,  3.374738e+12,
            -2.013750e+11,  7.942886e+09, -2.070369e+08,  3.438044e+06,
            -3.299673e+04,  1.394183e+02])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def CaviMlp_EFocus2_41(arg):
    p_T = numpy.array([
            -1.499544e+11, 5.612073e+10, -9.246033e+09, 8.799404e+08,
            -5.330725e+07, 2.132552e+06, -5.619149e+04, 8.943931e+02,
            -9.121320e+00, 1.038803e+00])
    p_S = numpy.array([
            -1.983302e+10, 8.570757e+09, -1.604935e+09, 1.714580e+08,
            -1.154148e+07, 5.095765e+05, -1.488249e+04, 2.696971e+02,
            -2.585211e+00, 1.305154e-02])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def CaviMlp_EDipole_41(arg):
    p_T = numpy.array([
            4.758398e+10, -1.656906e+10, 2.535541e+09, -2.237287e+08,
            1.255841e+07, -4.669147e+05, 1.125013e+04, -1.047651e+02,
            1.526489e+00, -1.005885e+00])
    p_S = numpy.array([
            1.155597e+11, -4.227114e+10, 6.817810e+09, -6.361033e+08,
            3.782592e+07, -1.488484e+06, 3.888964e+04, -6.407538e+02,
            5.884367e+00, -2.586200e-02])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def CaviMlp_EQuad_41(arg):
    p_T = numpy.array([
            -1.578312e+11, 5.896915e+10, -9.691159e+09, 9.192347e+08,
            -5.544764e+07, 2.206120e+06, -5.779110e+04, 9.127945e+02,
            -9.238897e+00, 1.038941e+00])
    p_S = numpy.array([
            -5.496217e+11, 2.008767e+11, -3.239241e+10, 3.024208e+09,
            -1.801113e+08, 7.094882e+06, -1.848380e+05, 3.069331e+03,
            -2.923507e+01, 1.248096e-01])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def CaviMlp_HMono_41(arg):
    p_T =numpy.array([
            -7.604796e+11, 4.632851e+11, -1.014721e+11, 1.165760e+10,
            -8.043084e+08, 3.518178e+07, -9.843253e+05, 1.697657e+04,
            -1.671357e+02, 1.703336e+00])
    p_S = numpy.array([
            -5.930241e+13, 2.189668e+13, -3.565836e+12, 3.360597e+11,
            -2.019481e+10, 8.022856e+08, -2.106663e+07, 3.524921e+05,
            -3.409550e+03, 1.452657e+01])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def CaviMlp_HDipole_41(arg):
    p_T = numpy.array([
            7.869717e+11, -3.116216e+11, 5.414689e+10, -5.420826e+09,
            3.446369e+08, -1.442888e+07, 3.985674e+05, -7.117391e+03,
            7.075414e+01, 6.853803e-01])
    p_S = numpy.array([
            -4.941947e+12, 1.791634e+12, -2.864139e+11, 2.649289e+10,
            -1.562284e+09, 6.090118e+07, -1.569273e+06, 2.575274e+04,
            -2.441117e+02, 1.021102e+00])
    return numpy.polyval(p_T, arg), numpy.polyval(p_S, arg)

def CaviMlp_HQuad_41(arg):
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
    'EFocus1' : CaviMlp_EFocus1_41,
    'EFocus2' : CaviMlp_EFocus2_41,
    'EDipole' : CaviMlp_EDipole_41,
    'EQuad'   : CaviMlp_EQuad_41,
    'HMono'   : CaviMlp_HMono_41,
    'HDipole' : CaviMlp_HDipole_41,
    'HQuad'   : CaviMlp_HQuad_41
    }


def get_cav(cav_hom, f, beta):
    lambda_ = scipy.constants.c/f
    arg = 2.0*math.pi/(beta*1e3*lambda_)
    beta_min = 0.025
    beta_max = 0.08
    if beta_min < beta and beta < beta_max:
        [T, S] = cav_transit_times_41[cav_hom](arg)
        print '                                T = %18.15f, S = %18.15f %s' % \
            (T, S, cav_hom)
    else:
        print 'beta out of range: %5.3f [%5.3f, %5.3f]' % \
            (beta, beta_min, beta_max)
        exit(1)
    return T, S


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


def get_cavity(file_name, f, beta):
    # Compute: e-m center, transit time factors [T, T', S, S'], ???.
    [z, EM] = rd_hom(file_name)
    lambda_ = scipy.constants.c/f
#    print '\nlambda %7.5f [m] = ' % (lambda_)
    [T, Tp, S, Sp, EML, em_center] = \
        get_transit_time_factors(z, EM, beta, lambda_)
    name = file_name.split('_')[1]
#    print 'EM center = %18.15f, T = %18.15f, S = %18.15f' \
#        ', EML = %18.15f %s' % \
#        (1e3*em_center, T, S, 1e-6*EML, name)
#    print 'EM center = %18.15f, T = %18.15f, S = %18.15f' \
#        ', Tp = %18.15f, Sp = %18.15f, EML = %18.15f %s' % \
#        (1e3*em_center, T, S, Tp, Sp, 1e-6*EML, name)
    return [em_center, T, Tp, S, Sp, EML]


home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'

# QWR cavity.
f_QWR = 80.5e6
beta  =  0.041
# HWR cavity.
f_HWR = 322e6

print
get_cav('EFocus1', f_QWR, beta)
get_cav('EDipole', f_QWR, beta)
get_cav('HDipole', f_QWR, beta)
get_cav('HMono',   f_QWR, beta)
get_cav('HQuad',   f_QWR, beta)
get_cav('EQuad',   f_QWR, beta)
get_cav('EFocus2', f_QWR, beta)

print
get_cavity(home_dir+'CaviMlp_EFocus1_41.txt', f_QWR, beta)
get_cavity(home_dir+'CaviMlp_EDipole_41.txt', f_QWR, beta)
get_cavity(home_dir+'CaviMlp_HDipole_41.txt', f_QWR, beta)
get_cavity(home_dir+'CaviMlp_HMono_41.txt',   f_QWR, beta)
get_cavity(home_dir+'CaviMlp_HQuad_41.txt',   f_QWR, beta)
get_cavity(home_dir+'CaviMlp_EQuad_41.txt',   f_QWR, beta)
get_cavity(home_dir+'CaviMlp_EFocus2_41.txt', f_QWR, beta)

outf = open('transit_times.dat', 'w')
beta_min = 0.025
beta_max = 0.08
n = 25
dbeta = (beta_max-beta_min)/(n-1)
for k in range(n):
    beta = beta_min + k*dbeta
    [em_center, T, Tp, S, Sp, EML] = \
        get_cavity(home_dir+'CaviMlp_EFocus1_41.txt', f_QWR, beta)
    outf.write('%7.5f %7.5f %7.5f\n' % (beta, T, S))
outf.close()
