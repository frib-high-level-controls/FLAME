import math
import numpy
import scipy.integrate

home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'


def rd_hom(file_name):

    xy = numpy.loadtxt(file_name)
    xy = numpy.array(xy)
    # z-axis is in [mm].
    xy[:,0] *= 1e-3
    return xy


def integrate_hom(file_name):
    xy = rd_hom(file_name)
    return {
        'simps': scipy.integrate.simps(xy[:,1], xy[:,0], even='avg'),
        'trapz': numpy.trapz(xy[:,1], xy[:,0])
        }['simps']


Edip = integrate_hom(home_dir+'CaviMlp_EDipole_41.txt')
print 'Edip    = %11.4e [MV]' % (1e-6*Edip)

Equad = integrate_hom(home_dir+'CaviMlp_EQuad_41.txt')
print 'Equad   = %11.4e [MV/m]' % (1e-6*Equad)

Efocus1 = integrate_hom(home_dir+'CaviMlp_EFocus1_41.txt')
print 'Efocus1 = %11.4e [MV/m]' % (1e-6*Efocus1)

Efocus2 = integrate_hom(home_dir+'CaviMlp_EFocus2_41.txt')
print 'Efocus2 = %11.4e [MV/m]' % (1e-6*Efocus2)
