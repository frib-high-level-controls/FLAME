import math
import numpy
import re
import StringIO
import sys

# Script to translate from TLM flat file to U-SCSI lattice.
#
#  type            name          length [m]  type [-18: Stripper, -21: Corrector, -28: BPM, -30: Matrix]
#  mark     LS1_CA01:BPM_D1129    0.000000                             -28.000000 
#
#  type            name          length [m]  aper [m]
#  drift           DRIFT          0.072000   0.020000 
#
#   type           name          length [m]  aper [m]   B [T]
# solenoid  LS1_CA01:SOL1_D1131   0.100000   0.020000  5.340000 
#
#   type           name          length [m]  aper [m]  phi [deg]  beta*gamma    type     phi1 [deg]  phi2 [deg]
#  dipole     FS1_CSS:DH_D2163    0.060000   0.020000  -1.000000   0.190370  400.000000   0.000000   -5.000000
#
#   type           name          length [m]  aper [m]  B2 [T/m]
# quadpole    FS1_CSS:QH_D2194    0.250000   0.025000  3.459800 
#
#   type           name          length [m]  aper [m]     f [MHz]  scl fact   phi [deg]
# rfcavity  LS1_CA01:CAV1_D1127   0.240000   0.017000  80.500000   0.640000  -35.000000
#
#   type           name          length [m]  aper [m]  phi [deg]     beta     cyl(0)/spher(1)  hor(0)/ver(1)...
#  ebend            EB1             1.0        0.02      90.0     0.00506953       1                1...
#   X fringe Y fringe  voltage asym fact
#     0.0      0.0           0.0
#
#   type           name          length [m]  aper [m]  voltage [V]  electrode radius [m]
#  equad           QE1H            0.1034      0.02     -358.574          0.0746
#
# Mis-alignment:
#
#   dx [m]  dy [m]  pitch [rad]  yaw [rad]  tilt [rad]


def get_misalign(tokens, ind):
    return '  dx = %s, dy = %s, pitch = %s, yaw = %s, tilt = %s;' \
        % (tokens[ind], tokens[ind+1], tokens[ind+2], tokens[ind+3], tokens[ind+4]) 


def get_index(tokens, token):
    return tokens.index(token) if token in tokens else None


def marker(line, tokens):
    global beam_line, n_marker, add_ind
    if float(tokens[2]) != 0:
        print '*** marker with non zero length: '
        exit(1)
    if add_ind: n_marker += 1; tokens[1] += '_%d' % (n_marker)
    beam_line.append(tokens[1])
    return '%s: marker;' % (tokens[1])

def drift(line, tokens):
    global beam_line, n_drift
    if add_ind: n_drift += 1; tokens[1] += '_%d' % (n_drift)
    beam_line.append(tokens[1])
    return '%s: drift, L = %s, aper = %s;' % (tokens[1], tokens[2], tokens[3])

def sbend(line, tokens):
    global beam_line, n_sbend, add_ind
    if add_ind: n_sbend += 1; tokens[1] += '_%d' % (n_sbend)
    beam_line.append(tokens[1])
    str = '%s: sbend, L = %s, phi = %s, phi1 = %s, phi2 = %s, bg = %s, aper = %s;' \
          % (tokens[1], tokens[2], tokens[4], tokens[7], tokens[8], tokens[5], tokens[3])
#    str += get_misalign(tokens, 9)
    return str

def solenoid(line, tokens):
    global beam_line, n_solenoid, add_ind
    if add_ind: n_solenoid += 1; tokens[1] += '_%d' % (n_solenoid)
    beam_line.append(tokens[1])
    str = '%s: solenoid, L = %s, B = %s, aper = %s,\n' \
        % (tokens[1], tokens[2], tokens[4], tokens[3])
    str += get_misalign(tokens, 5)
    return str

def quadrupole(line, tokens):
    global beam_line, n_quad, add_ind
    if add_ind: n_quad += 1; tokens[1] += '_%d' % (n_quad)
    beam_line.append(tokens[1])
    str = '%s: quadrupole, L = %s, B2 = %s, aper = %s,\n' \
        % (tokens[1], tokens[2], tokens[4], tokens[3])
    str += get_misalign(tokens, 5)
    return str

def rfcavity(line, tokens):
    global beam_line, n_cavity, add_ind
    if add_ind: n_cavity += 1; tokens[1] += '_%d' % (n_cavity)
    beam_line.append(tokens[1])
    str = '%s: rfcavity, cavtype = \"0.041QWR\", L = %s, f = %se6, phi = %s,\n  scl_fac = %s,' \
          ' aper = %s,\n' \
          % (tokens[1], tokens[2], tokens[4], tokens[6], tokens[5], tokens[3])
    str += get_misalign(tokens, 7)
    return str

def edipole(line, tokens):
    global beam_line, n_edipole
    n_edipole += 1; tokens[1] += '_%d' % (n_edipole)
    if add_ind: n_edipole += 1; tokens[1] += '_%d' % (n_ecavity)
    beam_line.append(tokens[1])
    return '%s: edipole, L = %s, phi = %s, x_frng = %s, y_frng = %s, beta = %s,' \
        ' spher = %s, asym_fac = %s, aper = %s;' \
        % (tokens[1], tokens[2], tokens[4], tokens[7], tokens[8],
           tokens[5], tokens[6], tokens[9], tokens[3])
#    str += get_misalign(tokens, 10)
    return str


def equad(line, tokens):
    global beam_line, n_equad
    n_equad += 1; tokens[1] += '_%d' % (n_equad)
    if add_ind: n_equad += 1; tokens[1] += '_%d' % (n_equad)
    beam_line.append(tokens[1])
    return '%s: equad, L = %s, V = %s, radius = %s, aper = %s;' \
        % (tokens[1], tokens[2], tokens[4], tokens[5], tokens[3])
#    str += get_misalign(tokens, 7)
    return str


tlm_dict = {
    'mark'     : marker,
    'drift'    : drift,
    'solenoid' : solenoid,
    'dipole'   : sbend,
    'quadpole' : quadrupole,
    'combquad' : quadrupole,
    'rfcavity' : rfcavity,
    'ebend'    : edipole,
    'equad'    : equad,
    }


def parse_definition(line, tokens):
    n_elem = 10; # No of elements per line.

    for k in range(len(tokens)):
        # Remove white space; unless a string.
        if not tokens[k].startswith('"'):
            tokens[k] = re.sub('[\s]', '', tokens[k])
    try:
        str = tlm_dict[tokens[0]](line, tokens)
    except KeyError:
        print '\n*** undefined token: ', tokens[0]
        print line
        print tokens
        exit(1)
    return str


def parse_line(line, outf):
    line_lc = line.lower()
    if not line_lc.rstrip():
        # Blank line.
        outf.write('\n')
    elif line_lc.startswith('#'):
        # Comment.
        outf.write('{ %s }\n' % (line.strip('!')))
    else:
        # Definition.
#        tokens = re.split(r'[ ]', line_lc)
        tokens = re.split(r'\s+', line_lc)
        # Replace ':' with '_' in name.
        tokens[1] = tokens[1].replace(':', '_', 1)
        outf.write('%s\n' % (parse_definition(line_lc, tokens)))


def prt_decl(outf):
    outf.write('# Beam envelope simulation.\n')
    outf.write('\nsim_type = "MomentMatrix2";\n\n')


def transl_file(file_name):
    global beam_line
    str = file_name.split('.')[0]+'.lat'
    inf = open(file_name, 'r')
    outf = open(str, 'w')
    prt_decl(outf)
    line = inf.readline()
    while line:
        line = line.strip('\r\n')
        while line.endswith('&'):
            # Line
            line = line.strip('&')
            line += (inf.readline()).strip('\r\n')
        parse_line(line, outf)
        line = inf.readline()
    outf.write('\ncell: LINE = (\n  S,\n')
    n = len(beam_line); n_max = 8;
    outf.write('  ')
    for k in range(n):
        outf.write(beam_line[k]);
        if (k+1 != n): outf.write(', ')
        if (k+1) % n_max == 0: outf.write('\n'); outf.write('  ')
    if n % n_max != n_max: outf.write('\n')
    outf.write(');\n')
    outf.write('\nUSE: cell;\n')


home_dir = ''

n_marker = 0; n_drift = 0; n_sbend = 0; n_solenoid = 0
n_quad = 0; n_cavity = 0; n_ecavity = 0; n_edipole = 0; n_equad = 0

add_ind = True

beam_line = [];

transl_file(home_dir+sys.argv[1])
