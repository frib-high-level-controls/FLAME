import math
import numpy
import sys
import time

home_dir = '/home/bengtsson/FRIB/Cavity Model/Multipole41/'


def rd_hom(file_name):
    lines = [line.rstrip('\n') for line in open(file_name)]
    for line in lines:
        print '%s' % (line)


#def integrate_hom():


file_name = 'CaviMlp_EDipole_41.txt'

rd_hom(home_dir+file_name)
