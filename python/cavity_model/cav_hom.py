import math
import numpy
import scipy.integrate
import scipy.constants

# Module to evaluate the cavity transit time factors the multipole field maps.
#
# Author: Johan Bengtsson.


def rd_hom(file_name):
    xy = numpy.loadtxt(file_name)
    # z-axis is in [mm].
    xy[:,0] *= 1e-3
    return xy[:, 0], xy[:, 1]


def trapezoidal(x, y):
    # Trapezoid rule for numerical integration.
    integ = 0.0
    for k in range(0, len(x)-1):
        integ += (y[k+1]+y[k])/2.0*(x[k+1]-x[k])
    return integ


def get_EML(z, EM):
    # Method: 'trapezoidal', 'trapz', 'simps', and 'romb'.
    method = 'trapz'

    if method == 'trapezoidal':
        integ = trapezoidal(z, EM)
    if method == 'trapz':
        integ = numpy.trapz(EM, z)
    elif method == 'simps':
        # Simpson's method requires an odd number of samples.
        n = len(z); h = z[1] - z[0]
        if n % 2 == 0:
            z = numpy.append(z, z[n-1]+h); EM = numpy.append(EM, 0.0)
        integ = scipy.integrate.simps(EM, z)
    elif method == 'romb':
        # Romberg's method requires 2^k + 1 samples.
        n = len(z); k = int(numpy.floor(numpy.log2(n))) + 1
        npad = 2**k + 1 - n; h = z[1] - z[0]
        for i in range(n, n+npad):
            z = numpy.append(z, z[n-1]+h); EM = numpy.append(EM, 0.0)
        integ = scipy.integrate.romb(EM, dx=h)
    return integ


def get_EM_center(z, EM):
    # Compute e-m center along z-axis.
    eml = 0.0; em_mom = 0.0
    for k in range(0, len(z)-1):
        h = z[k+1] - z[k]
        z_avg = (z[k]+z[k+1])/2.0; EM_avg = (EM[k]+EM[k+1])/2.0
        eml += EM_avg*h; em_mom += EM_avg*z_avg*h
    em_center = em_mom/eml
    return eml, em_center


def get_transit_time_factors(z, EM, beta, lambda_):
    # Compute transit time factors: [T, T', S, S'].
    [EML, em_center] = get_EM_center(z, numpy.absolute(EM))
    z -= em_center
    T = 0.0; Tp = 0.0; S = 0.0; Sp = 0.0
    for k in range(0, len(z)-1):
        h = z[k+1] - z[k];
        z_avg = (z[k]+z[k+1])/2.0; EM_avg = (EM[k]+EM[k+1])/2.0;
        arg = 2.0*math.pi*z_avg/(beta*lambda_)
        T += EM_avg*math.cos(arg)*h; Tp -= EM_avg*z_avg*math.sin(arg)*h
        S += EM_avg*math.sin(arg)*h; Sp += EM_avg*z_avg*math.cos(arg)*h
    T /= EML; Tp /= EML; S /= EML; Sp /= EML
    return T, Tp, S, Sp, EML, em_center


def get_cavity(file_name, f, beta):
    # Compute: e-m center, transit time factors [T, T', S, S'], ???.
    [z, EM] = rd_hom(file_name)
    lambda_ = scipy.constants.c/f
#    print '\nlambda %7.5f [m] = ' % (lambda_)
    [T, Tp, S, Sp, EML, em_center] = \
        get_transit_time_factors(z, EM, beta, lambda_)
    name = file_name.split('_')[1]
    print 'EM center = %18.15f, T = %18.15f, S = %18.15f, EML = %18.15f %s' % \
        (1e3*em_center, T, S, 1e-6*EML, name)
    return [em_center, T, Tp, S, Sp, EML]


def get_cav(cav, f, beta):
    lambda_ = scipy.constants.c/f
    arg = 2.0*math.pi/(beta*1e3*lambda_)
    if cav == 'CaviMlp_EFocus2':
        T = -1.499544e+11*math.pow(arg,9)+5.612073e+10*math.pow(arg,8) \
            +-9.246033e+09*math.pow(arg,7)+8.799404e+08*math.pow(arg,6) \
            +-5.330725e+07*math.pow(arg,5)+2.132552e+06*math.pow(arg,4) \
            +-5.619149e+04*math.pow(arg,3)+8.943931e+02*math.pow(arg,2) \
            +-9.121320e+00*math.pow(arg,1)+1.038803e+00;
        S = -1.983302e+10*math.pow(arg,9)+8.570757e+09*math.pow(arg,8) \
            +-1.604935e+09*math.pow(arg,7)+1.714580e+08*math.pow(arg,6) \
            +-1.154148e+07*math.pow(arg,5)+5.095765e+05*math.pow(arg,4) \
            +-1.488249e+04*math.pow(arg,3)+2.696971e+02*math.pow(arg,2) \
            +-2.585211e+00*math.pow(arg,1)+1.305154e-02;
    elif cav == 'CaviMlp_EDipole':
        T = 4.758398e+10*math.pow(arg,9)+-1.656906e+10*math.pow(arg,8) \
            +2.535541e+09*math.pow(arg,7)+-2.237287e+08*math.pow(arg,6) \
            +1.255841e+07*math.pow(arg,5)+-4.669147e+05*math.pow(arg,4) \
            +1.125013e+04*math.pow(arg,3)+-1.047651e+02*math.pow(arg,2) \
            +1.526489e+00*math.pow(arg,1)+-1.005885e+00;
        S = 1.155597e+11*math.pow(arg,9)+-4.227114e+10*math.pow(arg,8) \
            +6.817810e+09*math.pow(arg,7)+-6.361033e+08*math.pow(arg,6) \
            +3.782592e+07*math.pow(arg,5)+-1.488484e+06*math.pow(arg,4) \
            +3.888964e+04*math.pow(arg,3)+-6.407538e+02*math.pow(arg,2) \
            +5.884367e+00*math.pow(arg,1)+-2.586200e-02;
    elif cav == 'CaviMlp_HDipole':
        T = 7.869717e+11*math.pow(arg,9)+-3.116216e+11*math.pow(arg,8) \
            +5.414689e+10*math.pow(arg,7)+-5.420826e+09*math.pow(arg,6) \
            +3.446369e+08*math.pow(arg,5)+-1.442888e+07*math.pow(arg,4) \
            +3.985674e+05*math.pow(arg,3)+-7.117391e+03*math.pow(arg,2) \
            +7.075414e+01*math.pow(arg,1)+6.853803e-01;
        S = -4.941947e+12*math.pow(arg,9)+1.791634e+12*math.pow(arg,8) \
            +-2.864139e+11*math.pow(arg,7)+2.649289e+10*math.pow(arg,6) \
            +-1.562284e+09*math.pow(arg,5)+6.090118e+07*math.pow(arg,4) \
            +-1.569273e+06*math.pow(arg,3)+2.575274e+04*math.pow(arg,2) \
            +-2.441117e+02*math.pow(arg,1)+1.021102e+00;
    elif cav == 'CaviMlp_HMono':
        T = -7.604796e+11*math.pow(arg,9)+4.632851e+11*math.pow(arg,8) \
            +-1.014721e+11*math.pow(arg,7)+1.165760e+10*math.pow(arg,6) \
            +-8.043084e+08*math.pow(arg,5)+3.518178e+07*math.pow(arg,4) \
            +-9.843253e+05*math.pow(arg,3)+1.697657e+04*math.pow(arg,2) \
            +-1.671357e+02*math.pow(arg,1)+1.703336e+00;
        S = -5.930241e+13*math.pow(arg,9)+2.189668e+13*math.pow(arg,8)\
        +-3.565836e+12*math.pow(arg,7)+3.360597e+11*math.pow(arg,6) \
            +-2.019481e+10*math.pow(arg,5)+8.022856e+08*math.pow(arg,4) \
            +-2.106663e+07*math.pow(arg,3)+3.524921e+05*math.pow(arg,2) \
            +-3.409550e+03*math.pow(arg,1)+1.452657e+01;
    elif cav == 'CaviMlp_HQuad':
        T = 5.600545e+12*math.pow(arg,9)+-2.005326e+12*math.pow(arg,8) \
            +3.163675e+11*math.pow(arg,7)+-2.885455e+10*math.pow(arg,6) \
            +1.676173e+09*math.pow(arg,5)+-6.429625e+07*math.pow(arg,4) \
            +1.627837e+06*math.pow(arg,3)+-2.613724e+04*math.pow(arg,2) \
            +2.439177e+02*math.pow(arg,1)+-1.997432e+00;
        S = 1.131390e+13*math.pow(arg,9)+-4.119861e+12*math.pow(arg,8) \
            +6.617859e+11*math.pow(arg,7)+-6.153570e+10*math.pow(arg,6) \
            +3.649414e+09*math.pow(arg,5)+-1.431267e+08*math.pow(arg,4) \
            +3.711527e+06*math.pow(arg,3)+-6.135071e+04*math.pow(arg,2) \
            +5.862902e+02*math.pow(arg,1)+-2.470704e+00;
    elif cav == 'CaviMlp_EQuad':
        T = -1.578312e+11*math.pow(arg,9)+5.896915e+10*math.pow(arg,8) \
            +-9.691159e+09*math.pow(arg,7)+9.192347e+08*math.pow(arg,6) \
            +-5.544764e+07*math.pow(arg,5)+2.206120e+06*math.pow(arg,4) \
            +-5.779110e+04*math.pow(arg,3)+9.127945e+02*math.pow(arg,2) \
            +-9.238897e+00*math.pow(arg,1)+1.038941e+00;
        S = -5.496217e+11*math.pow(arg,9)+2.008767e+11*math.pow(arg,8) \
            +-3.239241e+10*math.pow(arg,7)+3.024208e+09*math.pow(arg,6) \
            +-1.801113e+08*math.pow(arg,5)+7.094882e+06*math.pow(arg,4) \
            +-1.848380e+05*math.pow(arg,3)+3.069331e+03*math.pow(arg,2) \
            +-2.923507e+01*math.pow(arg,1)+1.248096e-01;
    elif cav == 'CaviMlp_EFocus1':
        T = -7.316972e14*math.pow(arg,9)+2.613462e14*math.pow(arg,8) \
            -4.112154e+13*math.pow(arg,7)+3.739846e+12*math.pow(arg,6) \
            +-2.165867e+11*math.pow(arg,5)+8.280687e+09*math.pow(arg,4) \
            +-2.089452e+08*math.pow(arg,3)+3.354464e+06*math.pow(arg,2) \
            +-3.108322e+04*math.pow(arg,1)+1.256386e+02;
        S = -6.079177e+14*math.pow(arg,9)+2.229446e+14*math.pow(arg,8) \
            +-3.605780e+13*math.pow(arg,7)+3.374738e+12*math.pow(arg,6) \
            +-2.013750e+11*math.pow(arg,5)+7.942886e+09*math.pow(arg,4) \
            +-2.070369e+08*math.pow(arg,3)+3.438044e+06*math.pow(arg,2) \
            +-3.299673e+04*math.pow(arg,1)+1.394183e+02;
    name = cav.split('_')[1]
    print '                                T = %18.15f, S = %18.15f %s' % \
        (T, S, name)
    return S, T


home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'

f_cav = 80.5e6; beta = 0.041

print
get_cav('CaviMlp_EFocus1', f_cav, beta)
get_cav('CaviMlp_EDipole', f_cav, beta)
get_cav('CaviMlp_HDipole', f_cav, beta)
get_cav('CaviMlp_HMono',   f_cav, beta)
get_cav('CaviMlp_HQuad',   f_cav, beta)
get_cav('CaviMlp_EQuad',   f_cav, beta)
get_cav('CaviMlp_EFocus2', f_cav, beta)

print
get_cavity(home_dir+'CaviMlp_EFocus1_41.txt', f_cav, beta)
get_cavity(home_dir+'CaviMlp_EDipole_41.txt', f_cav, beta)
get_cavity(home_dir+'CaviMlp_HDipole_41.txt', f_cav, beta)
get_cavity(home_dir+'CaviMlp_HMono_41.txt',   f_cav, beta)
get_cavity(home_dir+'CaviMlp_HQuad_41.txt',   f_cav, beta)
get_cavity(home_dir+'CaviMlp_EQuad_41.txt',   f_cav, beta)
get_cavity(home_dir+'CaviMlp_EFocus2_41.txt', f_cav, beta)
