import math
import numpy
import scipy.integrate

home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'


def rd_hom(file_name):

    xy = numpy.loadtxt(file_name)
    # z-axis is in [mm].
    xy[:,0] *= 1e-3
    return xy


def integrate_hom(file_name):
    xy = rd_hom(file_name)

    method = 'simps'
    if method == 'simps':
        integ = scipy.integrate.simps(xy[:,1], xy[:,0], even='avg')
    elif method == 'trapz':
        integ = numpy.trapz(xy[:,1], xy[:,0])
    elif method == 'trapz':
        # Romberg's method requires 2^k + 1 samples.
        integ = scipy.integrate.romb(xy[:,1], xy[:,0])
    return integ


print
Edip = integrate_hom(home_dir+'CaviMlp_EDipole_41.txt')
Equad = integrate_hom(home_dir+'CaviMlp_EQuad_41.txt')
Efocus1 = integrate_hom(home_dir+'CaviMlp_EFocus1_41.txt')
Efocus2 = integrate_hom(home_dir+'CaviMlp_EFocus2_41.txt')

Hdip = integrate_hom(home_dir+'CaviMlp_HDipole_41.txt')
Hmono = integrate_hom(home_dir+'CaviMlp_HMono_41.txt')
Hquad = integrate_hom(home_dir+'CaviMlp_HQuad_41.txt')
print
print 'Efocus1 = %13.10f [MV/m]' % (1e-6*Efocus1)
print 'Edip    = %13.10f [MV]' % (1e-6*Edip)
print 'Hdip    = %13.10f [MA]' % (1e-6*Hdip)
print 'Hmono   = %13.10f [MA/m]' % (1e-6*Hmono)
print 'Hquad   = %13.10f [MA/m]' % (1e-6*Hquad)
print 'Equad   = %13.10f [MV/m]' % (1e-6*Equad)
print 'Efocus2 = %13.10f [MV/m]' % (1e-6*Efocus2)
