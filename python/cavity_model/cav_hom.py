import math
import numpy
import scipy.integrate

# Module to evaluate the integrated field from the cavity multipole field maps.
#
# Author: Johan Bengtsson.


def rd_hom(file_name):
    xy = numpy.loadtxt(file_name)
    # z-axis is in [mm].
    xy[:,0] *= 1e-3
    return xy


def integrate_hom(file_name):
    # Method: 'trapz', 'simps', and 'romb'.
    method = 'trapz'

    xy = rd_hom(home_dir+file_name)

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


home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'

Edip    = integrate_hom('CaviMlp_EDipole_41.txt')
Equad   = integrate_hom('CaviMlp_EQuad_41.txt')
Efocus1 = integrate_hom('CaviMlp_EFocus1_41.txt')
Efocus2 = integrate_hom('CaviMlp_EFocus2_41.txt')

Hdip    = integrate_hom('CaviMlp_HDipole_41.txt')
Hmono   = integrate_hom('CaviMlp_HMono_41.txt')
Hquad   = integrate_hom('CaviMlp_HQuad_41.txt')

print
print 'Efocus2 = %18.15f [MV/m]' % (1e-6*Efocus2)
print 'Edip    = %18.15f [MV]'   % (1e-6*Edip)
print 'Hdip    = %18.15f [MA]'   % (1e-6*Hdip)
print 'Hmono   = %18.15f [MA/m]' % (1e-6*Hmono)
print 'Hquad   = %18.15f [MA/m]' % (1e-6*Hquad)
print 'Equad   = %18.15f [MV/m]' % (1e-6*Equad)
print 'Efocus1 = %18.15f [MV/m]' % (1e-6*Efocus1)
