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
#   type           name          length [m]  aper [m]  phi [deg]  beta*gamma    type     phi1 [deg]  phi2 [deg]  b2 [1/m^2]
#  dipole     FS1_CSS:DH_D2163    0.060000   0.020000  -1.000000   0.190370  400.000000   0.000000   -5.000000
#
#   type           name          length [m]  aper [m]  B2 [T/m]
# quadpole    FS1_CSS:QH_D2194    0.250000   0.025000  3.459800 
#
#   type           name          length [m]  aper [m]     f [MHz]  scl fact   phi [deg]
# rfcavity  LS1_CA01:CAV1_D1127   0.240000   0.017000  80.500000   0.640000  -35.000000


def get_index(tokens, token):
    return tokens.index(token) if token in tokens else None


def marker(line, tokens):
    global beam_line, n_marker, add_ind
    if float(tokens[3]) != 0:
        print '*** marker with non zero length: '
        exit(1)
    if add_ind: n_marker += 1; tokens[2] += '_%d' % (n_marker)
    beam_line.append(tokens[2])
    return '%s: marker;' % (tokens[2])

def drift(line, tokens):
    global beam_line, n_drift
    if add_ind: n_drift += 1; tokens[2] += '_%d' % (n_drift)
    beam_line.append(tokens[2])
    return '%s: drift, L = %s, aper = %s;' % (tokens[2], tokens[3], tokens[4])

def sbend(line, tokens):
    global beam_line, n_sbend, add_ind
    if add_ind: n_sbend += 1; tokens[2] += '_%d' % (n_sbend)
    beam_line.append(tokens[2])
    str = '%s: sbend, L = %s, phi = %s, phi1 = %s, phi2 = %s, bg = %s, aper = %s;' \
          % (tokens[2], tokens[3], tokens[5], tokens[8], tokens[9], tokens[6], tokens[4])
    return str

def solenoid(line, tokens):
    global beam_line, n_solenoid, add_ind
    if add_ind: n_solenoid += 1; tokens[2] += '_%d' % (n_solenoid)
    beam_line.append(tokens[2])
    return '%s: solenoid, L = %s, B = %s, aper = %s;' \
        % (tokens[2], tokens[3], tokens[5], tokens[4])

def quadrupole(line, tokens):
    global beam_line
    beam_line.append(tokens[2])
    return '%s: quadrupole, L = %s, B2 = %s, aper = %s;' \
        % (tokens[2], tokens[3], tokens[5], tokens[4])

def rfcavity(line, tokens):
    global beam_line
    beam_line.append(tokens[2])
    str = '%s: rfcavity, L = %s, f = %se6, phi = %s, scl_fac = %s,' \
          ' aper = %s;' \
          % (tokens[2], tokens[3], tokens[5], tokens[7], tokens[6], tokens[4])
    return str

def e_dipole(line, tokens):
    global beam_line, n_e_dipole
    # n_e_dipole += 1; tokens[2] += '_%d' % (n_e_dipole)
    beam_line.append(tokens[2])
    return '%s: edipole, L = %s, phi = %s, phi1 = %s, phi2 = %s, E= %s, scl_fac = %s,' \
        ' aper = %s;' \
        % (tokens[2], tokens[4], tokens[7], tokens[9], tokens[10], tokens[8], tokens[6],
           tokens[5])


# TLM -> Tracy-2,3 dictionary.
tlm2tracy = {
    'mark'     : marker,
    'drift'    : drift,
    'solenoid' : solenoid,
    'dipole'   : sbend,
    'quadpole' : quadrupole,
    'combquad' : quadrupole,
    'rfcavity' : rfcavity,
    'ebend'    : e_dipole,
    }


def parse_definition(line, tokens):
    n_elem = 10; # No of elements per line.

    for k in range(len(tokens)):
        # Remove white space; unless a string.
        if not tokens[k].startswith('"'):
            tokens[k] = re.sub('[\s]', '', tokens[k])
    try:
        str = tlm2tracy[tokens[1]](line, tokens)
    except KeyError:
        print '\n*** undefined token!'
        print line
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
        tokens = re.split(r'[ ]', line_lc)
        # Replace ':' with '_' in name.
        tokens[2] = tokens[2].replace(':', '_', 1)
        outf.write('%s\n' % (parse_definition(line_lc, tokens)))


def prt_decl(outf):
    outf.write('# Beam envelope simulation.\n')
    outf.write('\nsim_type = "MomentMatrix";\n\n')


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
    outf.write('\ncell: LINE = (\n')
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

n_marker = 0; n_drift = 0; n_sbend = 0; n_solenoid = 0; n_e_dipole = 0

add_ind = True

beam_line = [];

transl_file(home_dir+sys.argv[1])
