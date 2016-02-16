import math
import numpy


home_dir1 = "/home/johan/git_repos/jmbgsddb/build/src/"
home_dir2 = "/home/johan/tlm_workspace/TLM_JB/"

def rd_data_1(file_name):

    inf = open(file_name, 'r')
    s = []; x = []; xp = []; y = []; yp = []; z = []; zp = [];
    for line in inf:
        if line[0] != '#':
            [s1, x1, xp1, y1, yp1, z1, zp1] = line.split()
            s.append(float(s1))
            x.append(float(x1))
            xp.append(float(xp1))
            y.append(float(y1))
            yp.append(float(xp1))
            z.append(float(z1))
            zp.append(float(zp1))
    inf.close()

    return s, x, xp, y, yp, z, zp


def rd_data_2(file_name):

    inf = open(file_name, 'r')
    s = []; E = []; phi = [];
    for line in inf:
        if line[0] != '#':
            [elem_type, name, s1, E1, phi1, beta, gamma] = line.split()
            s.append(float(s1))
            E.append(float(E1))
            phi.append(float(phi1))
    inf.close()

    return s, E, phi


def rd_data_3(file_name):

    inf = open(file_name, 'r')
    s = []; E = []; phi = [];
    for line in inf:
        if line[0] != '#':
            [s1, E1, phi1] = line.split()
            s.append(float(s1))
            E.append(float(E1))
            phi.append(float(phi1))
    inf.close()

    return s, E, phi


def analyze_data_1(file_name_1, file_name_2, file_name_3):

    ps_dim = 6

    [s, x, xp, y, yp, z, zp] = rd_data_1(file_name_1)
    n = numpy.array([0, 0])
    n[0] = len(s)
    data = numpy.zeros((2, ps_dim+1, n[0]))
    data[0, 0, :] = s;
    data[0, 1, :] = x; data[0, 2,:] = xp;
    data[0, 3, :] = y; data[0, 4,:] = yp;
    data[0, 5, :] = z; data[0, 6,:] = zp;

    [s, x, xp, y, yp, z, zp] = rd_data_1(file_name_2)
    n[1] = len(s)
    data[1, 0, :] = s;
    data[1, 1, :] = x; data[1, 2,:] = xp;
    data[1, 3, :] = y; data[1, 4,:] = yp;
    data[1, 5, :] = z; data[1, 6,:] = zp;

    delta = numpy.zeros((ps_dim+1, n[0]))

    delta[0, :] = data[0, 0, :]

    for k in range(1, ps_dim+1):
        delta[k] = data[1, k, :] - data[0, k, :]

    outf = open(file_name_3, 'w')
    for j in range(n[0]):
        for k in range(ps_dim+1):
            outf.write('%23.15e' % delta[k, j])
        outf.write('\n')
    outf.close()


def analyze_data_2(file_name_1, file_name_2, file_name_3):
    n_data = 3

    [s, E, phi] = rd_data_2(file_name_1)
    n = numpy.array([0, 0])
    n[0] = len(s)
    data = numpy.zeros((2, n_data, n[0]))
    data[0, 0, :] = s; data[0, 1, :] = E; data[0, 2,:] = phi;

    [s, E, phi] = rd_data_3(file_name_2)
    n[1] = len(s)
    data[1, 0, :] = s; data[1, 1, :] = E; data[1, 2,:] = phi;

    delta = numpy.zeros((n_data, n[0]))

    delta[0, :] = data[0, 0, :]

    for k in range(1, n_data):
        delta[k] = data[1, k, :] - data[0, k, :]

    outf = open(file_name_3, 'w')
    for j in range(n[0]):
        for k in range(n_data):
            outf.write('%23.15e' % delta[k, j])
        outf.write('\n')
    outf.close()


analyze_data_1(home_dir1+'CenofChg.out', home_dir2+'MCSModelCenVec.txt', \
               'analyze_res_1.dat')

analyze_data_1(home_dir1+'BeamRMS.out', home_dir2+'MCSModelRmsVec.txt', \
               'analyze_res_2.dat')

analyze_data_2(home_dir1+'long_tab.out', home_dir2+'MCSModelRf.txt', \
               'analyze_res_3.dat')
