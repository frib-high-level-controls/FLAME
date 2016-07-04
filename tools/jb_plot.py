
from __future__ import print_function

import sys
import numpy

from matplotlib.pylab import *


def plt_moment0(s, moment0):
    for i, L in zip(range(6), ('x', 'p_x', 'y', 'p_y', 'z', 'p_z')):
        subplot(3, 2, i+1)
        for j in range(5):
            plot(s, moment0[:, j, i], '-b')
        plot(s, moment0[:, 5, i], '-r')
        xlabel('s'); ylabel(L)


def plt_moment1(s, moment1):
    for i, L in zip(range(6), ('s_x', 's_p_x', 's_y', 's_p_y', 's_z', 's_p_z')):
        subplot(3, 2, i+1)
        for j in range(5):
            plot(s, moment1[:, j, i], '-b')
        plot(s, moment1[:, 5, i], '-r')
        xlabel('s'); ylabel(L)


def rd_data(file_name):
    file = open(file_name, 'r')

    s = []; moment0 = []; padding = numpy.array([NaN, NaN, NaN, NaN, NaN, NaN])
    first = True
    for line in file:
        fields = line.strip().split()
        s = numpy.append(s, float(fields[0]))
        if len(fields) == 19:
            term = numpy.array([[fields[1:7], fields[7:13],
                                 padding, padding, padding,
                                 fields[13:19]]], dtype=float)
        else:
            term = numpy.array([[fields[1:7], fields[7:13], fields[13:19],
                                 fields[19:25], fields[25:31],
                                 fields[31:37]]], dtype=float)
        if first:
            moment0 = term
            first = False
        else:
            moment0 = numpy.append(moment0, term, 0)

    return [s, moment0]

[file_name1, file_name2] = [sys.argv[1], sys.argv[2]]

[s, moment0] = rd_data(file_name1)
[s, moment1] = rd_data(file_name2)

plt.rcParams['savefig.dpi'] = 600 # For png.

fig1 = figure(1)
plt_moment0(s, moment0)
fig2 = figure(2)
plt_moment1(s, moment1)

#fig.savefig('fig.png')
fig1.savefig('fig1.ps', orientation='landscape')
fig2.savefig('fig2.ps', orientation='landscape')

#show()
ion(); show(); ioff()
raw_input('<ret> to continue>')
