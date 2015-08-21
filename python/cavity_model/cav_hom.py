import math
import numpy
import scipy.integrate

# Module to obtained the integrated field from the cavity field maps.

home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'


def rd_hom(file_name):

    xy = numpy.loadtxt(file_name)
    # z-axis is in [mm].
    xy[:,0] *= 1e-3
    return xy


def integrate_hom(file_name):
    xy = rd_hom(file_name)
    print len(xy[:, 0])

    method = 'trapz'
    if method == 'trapz':
        integ = numpy.trapz(xy[:, 1], xy[:, 0])
    elif method == 'simps':
        # Simpson's method requires an odd number of samples.
        n = len(xy[:, 0])
        h = xy[1, 0] - xy[0, 0]
        if n % 2 == 0:
            xy = numpy.append(xy, [[xy[n-1, 0]+h, 0.0]], axis=0)
        integ = scipy.integrate.simps(xy[:, 1], xy[:, 0])
    elif method == 'romb':
        # Romberg's method requires 2^k + 1 samples.
        n = len(xy[:, 0])
        k = int(numpy.floor(numpy.log2(n))) + 1
        npad = 2**k + 1 - n
        h = xy[1, 0] - xy[0, 0]
        for i in range(n, n+npad):
            xy = numpy.append(xy, [[xy[i-1, 0]+h, 0.0]], axis=0)
#        print '%3d %3d %3d %1d %3d %7.5f' % \
            (n+npad, len(xy[:, 0]), n, k, npad, h)
        integ = scipy.integrate.romb(xy[:, 1], dx=h)
    return integ


Edip    = integrate_hom(home_dir+'CaviMlp_EDipole_41.txt')
Equad   = integrate_hom(home_dir+'CaviMlp_EQuad_41.txt')
Efocus1 = integrate_hom(home_dir+'CaviMlp_EFocus1_41.txt')
Efocus2 = integrate_hom(home_dir+'CaviMlp_EFocus2_41.txt')

Hdip    = integrate_hom(home_dir+'CaviMlp_HDipole_41.txt')
Hmono   = integrate_hom(home_dir+'CaviMlp_HMono_41.txt')
Hquad   = integrate_hom(home_dir+'CaviMlp_HQuad_41.txt')
print
print 'Efocus2 = %18.15f [MV/m]' % (1e-6*Efocus2)
print 'Edip    = %18.15f [MV]'   % (1e-6*Edip)
print 'Hdip    = %18.15f [MA]'   % (1e-6*Hdip)
print 'Hmono   = %18.15f [MA/m]' % (1e-6*Hmono)
print 'Hquad   = %18.15f [MA/m]' % (1e-6*Hquad)
print 'Equad   = %18.15f [MV/m]' % (1e-6*Equad)
print 'Efocus1 = %18.15f [MV/m]' % (1e-6*Efocus1)
