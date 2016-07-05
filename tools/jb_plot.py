
from __future__ import print_function

import sys
import numpy

from matplotlib.pylab import *


def plt_ref_orbit(s, lng):
    subplot(1, 2, 1)
    plot(s, lng[:, 0], '-b')
    xlabel('s [m]'); ylabel('phase [rad]')
    subplot(1, 2, 2)
    plot(s, lng[:, 1], '-r')
    xlabel('s [m]'); ylabel('E_k [MeV]')


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


def rd_long(file_name):
    file = open(file_name, 'r')

    s = []; lng = []
    first = True
    for line in file:
        fields = line.strip().split()
        s = numpy.append(s, float(fields[0]))
        term = numpy.array([[fields[1], fields[2]]], dtype=float)
        if first:
            lng = term
            first = False
        else:
            lng = numpy.append(lng, term, 0)

    return [s, lng]


def rd_long_TLM(file_name):
    file = open(file_name, 'r')

    s = []; lng = []
    first = True
    for line in file:
        fields = line.strip().split()
        s = numpy.append(s, float(fields[0]))
        term = numpy.array([[fields[2], fields[1]]], dtype=float)
        if first:
            lng = term
            first = False
        else:
            lng = numpy.append(lng, term, 0)

    return [s, lng]


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


def plt_long(fig_no, title, s, lng):
    fig = figure(fig_no)
    subplots_adjust(hspace=0.6) # default is 0.2.
    subplots_adjust(wspace=0.4) # default is 0.2.
    suptitle(title)
    plt_ref_orbit(s, lng)
    return fig


def plt0(fig_no, title, s, moment0):
    fig = figure(fig_no)
    subplots_adjust(hspace=0.6) # default is 0.2.
    subplots_adjust(wspace=0.4) # default is 0.2.
    suptitle(title)
    plt_moment0(s, moment0)
    return fig


def plt1(fig_no, title, s, moment1):
    fig = figure(fig_no)
    subplots_adjust(hspace=0.6) # default is 0.2.
    subplots_adjust(wspace=0.4) # default is 0.2.
    suptitle(title)
    plt_moment0(s, moment1)
    return fig


#[file_name1, file_name2] = [sys.argv[1], sys.argv[2]]

file_name1 = '/home/johan/git_repos/flame/build/src/ref_orbit.txt'
file_name2 = '/home/johan/tlm_workspace/TLM_JB/tab_jb.txt'

file_name3 = '/home/johan/git_repos/flame/build/src/moment0.txt'
file_name4 = '/home/johan/git_repos/flame/build/src/moment1.txt'
file_name5 = '/home/johan/tlm_workspace/TLM_JB/moment0_TLM.txt'
file_name6 = '/home/johan/tlm_workspace/TLM_JB/moment1_TLM.txt'

[s, lng]             = rd_long(file_name1)
lng[:, 1] /= 1e6
[s, lng_TLM]         = rd_long_TLM(file_name2)

[s, moment0]         = rd_data(file_name3)
[s, moment1]         = rd_data(file_name4)
[s_TLM, moment0_TLM] = rd_data(file_name5)
[s_TLM, moment1_TLM] = rd_data(file_name6)

lng_diff     = lng - lng_TLM
moment0_diff = moment0 - moment0_TLM
moment1_diff = moment1 - moment1_TLM

fig1 = plt_long(1, 'Ref Orbit for Corrected TLM', s, lng_TLM)
fig2 = plt_long(2, 'Ref Orbit for FLAME', s, lng)
fig3 = plt_long(3, 'Ref Orbit Difference Between FLAME and Corrected TLM',
                s, lng_diff)

fig4 = plt0(4, 'Orbit for Corrected TLM', s, moment0_TLM)
fig5 = plt0(5, 'Orbit for FLAME', s, moment0)
fig6 = plt0(6, 'Orbit Difference Between FLAME and Corrected TLM',
            s, moment0_diff)
fig7 = plt1(7, 'RMS Beam Size for Corrected TLM', s, moment1_TLM)
fig8 = plt1(8, 'RMS Beam Size for FLAME', s, moment1)
fig9 = plt1(9, 'RMS Beam Size Difference Between FLAME and Corrected TLM',
            s, moment1_diff)

plt.rcParams['savefig.dpi'] = 600 # For png.

fig1.savefig('fig1_FE.eps', orientation='landscape')
fig2.savefig('fig2_FE.eps', orientation='landscape')
fig3.savefig('fig3_FE.eps', orientation='landscape')
fig4.savefig('fig4_FE.eps', orientation='landscape')
fig5.savefig('fig5_FE.eps', orientation='landscape')
fig6.savefig('fig6_FE.eps', orientation='landscape')
fig7.savefig('fig7_FE.eps', orientation='landscape')
fig8.savefig('fig8_FE.eps', orientation='landscape')
fig9.savefig('fig9_FE.eps', orientation='landscape')

ion(); show(); ioff()
raw_input('<ret> to continue>')
