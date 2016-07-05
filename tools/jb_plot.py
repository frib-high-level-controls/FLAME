
from __future__ import print_function

import sys
import numpy

from matplotlib.pylab import *


def plt_moment0(s, moment0):
    for i, L in zip(range(6),
                    ('x [mm]', 'p_x [mrad]', 'y [mm]', 'p_y [mrad]',
                     'z [rad]', 'p_z')):
        subplot(3, 2, i+1)
        for j in range(5):
            plot(s, moment0[:, j, i], '-b')
        plot(s, moment0[:, 5, i], '-r')
        xlabel('s [m]'); ylabel(L)


def plt_moment1(s, moment1):
    for i, L in zip(range(6),
                    ('s_x [mm]', 's_p_x [mrad]', 's_y [mm]', 's_p_y [mrad]',
                     's_z [rad]', 's_p_z')):
        subplot(3, 2, i+1)
        for j in range(5):
            plot(s, moment1[:, j, i], '-b')
        plot(s, moment1[:, 5, i], '-r')
        xlabel('s [m]'); ylabel(L)


def plt_moment0_diff(s, moment0):
    for i, L in zip(range(6),
                    ('x [mm]', 'p_x [mrad]', 'y [mm]', 'p_y [mrad]',
                     'z [rad]', 'p_z')):
        subplot(3, 2, i+1)
        for j in range(5):
            plot(s, moment0[:, j, i], '-b')
        plot(s, moment0[:, 5, i], '-r')
        xlabel('s [m]'); ylabel(L)


def plt_moment1_diff(s, moment1):
    for i, L in zip(range(6),
                    ('s_x [mm]', 's_p_x [mrad]', 's_y [mm]', 's_p_y [mrad]',
                     's_z [rad]', 's_p_z')):
        subplot(3, 2, i+1)
        for j in range(5):
            plot(s, moment1[:, j, i], '-b')
        plot(s, moment1[:, 5, i], '-r')
        xlabel('s [m]'); ylabel(L)


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

#[file_name1, file_name2] = [sys.argv[1], sys.argv[2]]

file_name1 = '/home/johan/git_repos/flame/build/src/moment0.txt'
file_name2 = '/home/johan/git_repos/flame/build/src/moment1.txt'
file_name3 = '/home/johan/tlm_workspace/TLM_JB/moment0_TLM.txt'
file_name4 = '/home/johan/tlm_workspace/TLM_JB/moment1_TLM.txt'

[s, moment0]         = rd_data(file_name1)
[s, moment1]         = rd_data(file_name2)
[s_TLM, moment0_TLM] = rd_data(file_name3)
[s_TLM, moment1_TLM] = rd_data(file_name4)

moment0_diff = moment0 - moment0_TLM
moment1_diff = moment1 - moment1_TLM

plt.rcParams['savefig.dpi'] = 600 # For png.

subplots_adjust(hspace=0.6) # default is 0.2.

fig1 = figure(1)
suptitle('Orbit for Corrected TLM')
plt_moment0(s, moment0_TLM)
fig2 = figure(2)
suptitle('RMS Beam Size for Corrected TLM')
plt_moment1(s, moment1_TLM)
fig3 = figure(3)
suptitle('Orbit for FLAME')
plt_moment0(s, moment0)
fig4 = figure(4)
suptitle('RMS Beam Size for FLAME')
plt_moment1(s, moment1)
fig5 = figure(5)
suptitle('Orbit Difference Between FLAME and Corrected TLM')
plt_moment0(s, moment0_diff)
fig6 = figure(6)
suptitle('RMS Beam Size Difference Between FLAME and Corrected TLM')
plt_moment1(s, moment1_diff)

fig1.savefig('fig1_FE.eps', orientation='landscape')
fig2.savefig('fig2_FE.eps', orientation='landscape')
fig3.savefig('fig3_FE.eps', orientation='landscape')
fig4.savefig('fig4_FE.eps', orientation='landscape')
fig5.savefig('fig5_FE.eps', orientation='landscape')
fig6.savefig('fig6_FE.eps', orientation='landscape')

ion(); show(); ioff()
raw_input('<ret> to continue>')
