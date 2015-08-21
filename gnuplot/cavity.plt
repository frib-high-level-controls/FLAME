# Script to plot the cavity multipole field maps.
#
# Author: Johan Bengtsson.

home_dir = '/home/bengtsson/FRIB/Cavity Model/';

# Terminal type is: 0 - X11, 1 - PS, 2 - PNG, 3 - JPG
term_type = 0;

# On-axis long. E-field [V/m] vs. z [mm].
file1 = home_dir."axisData_41.txt";

# Transverse higher order E-field [V/m] ([V] for dipole) vs. z [mm].
# Dipole.
file2 = home_dir."Multipole41/CaviMlp_EDipole_41.txt";
# First half of defocusing component.
file3 = home_dir."Multipole41/CaviMlp_EFocus1_41.txt";
# Second half of defocusing component.
file4 = home_dir."Multipole41/CaviMlp_EFocus2_41.txt";
# Quadrupole.
file5 = home_dir."Multipole41/CaviMlp_EQuad_41.txt";

# Tranverse higher order H-field [A/m] ([A] for dipole) vs. z [mm]
file6 = home_dir."Multipole41/CaviMlp_HDipole_41.txt";
file7 = home_dir."Multipole41/CaviMlp_HMono_41.txt";
file8 = home_dir."Multipole41/CaviMlp_HQuad_41.txt";
file9 = home_dir."Multipole41/thinlenlon_41.txt";

if (term_type == 0) set terminal x11;
if (term_type == 1) set terminal postscript landscape enhanced color solid \
  lw 2 "Times-Roman" 15;
#if (term_type == 2) set terminal pngcairo size 350, 262 enhanced \
#  font 'Verdana,10';
if (term_type == 2) set terminal pngcairo enhanced font "Times-Roman";
if (term_type == 3) set terminal jpeg enhanced;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "cyan";
set style line 3 lt 1 lw 1 lc rgb "green";
set style line 4 lt 1 lw 1 lc rgb "red";
set style line 5 lt 1 lw 1 lc rgb "dark-green";
set style line 6 lt 1 lw 1 lc rgb "dark-orange";
set style line 7 lt 1 lw 1 lc rgb "dark-violet";
set style line 8 lt 1 lw 1 lc rgb "orange-red";

if (term_type == 1) set output "cavity_1.ps"
if (term_type == 2) set output "cavity_1.png"
if (term_type == 3) set output "cavity_1.jpg"

set title "Cavity Fundamental Mode: axisData";
set xlabel "z [cm]"; set ylabel "E_{/Symbol \174\174} [MV/m]";
plot file1 using (1e-1*$1):(1e-6*$2) notitle with lines ls 1;

if (term_type == 0) pause mouse "click on graph to cont.\n";

if (term_type == 1) set output "cavity_2.ps"
if (term_type == 2) set output "cavity_2.png"
if (term_type == 3) set output "cavity_2.jpg"

set multiplot;

set size 0.5, 0.5; set origin 0.0, 0.5;
set title "Cavity Mode: EDipole";
set xlabel "z [cm]"; set ylabel "E_{/Symbol \\136} [MV]";
plot file2 using (1e-1*$1):(1e-6*$2) notitle with lines ls 3;

set origin 0.5, 0.5;
set title "Cavity Mode: EQuad";
set xlabel "z [cm]"; set ylabel "E_{/Symbol \\136} [MV/m]";
plot file5 using (1e-1*$1):(1e-6*$2) notitle with lines ls 3;

set origin 0.0, 0.0;
set title "Cavity: EFocus1 & EFocus2";
set xlabel "z [cm]"; set ylabel "E_{/Symbol \\136} [MV/m]";
plot file3 using (1e-1*$1):(1e-6*$2) notitle with lines ls 3, \
     file4 using (1e-1*$1):(1e-6*$2) notitle with lines ls 3;

unset multiplot;
if (term_type == 0) pause mouse "click on graph to cont.\n";

if (term_type == 1) set output "cavity_3.ps"
if (term_type == 2) set output "cavity_3.png"
if (term_type == 3) set output "cavity_3.jpg"

set multiplot;

set size 0.5, 0.5; set origin 0.0, 0.5;
set title "Cavity Mode: HDipole";
set xlabel "z [cm]"; set ylabel "H_{/Symbol \\136} [kA]";
plot file6 using (1e-1*$1):(1e-3*$2) notitle with lines ls 3;

set origin 0.5, 0.5;
set title "Cavity Mode: HMono";
set xlabel "z [cm]"; set ylabel "H_{/Symbol \\136} [kA/m]";
plot file7 using (1e-1*$1):(1e-3*$2) notitle with lines ls 3;

set origin 0.0, 0.0;
set title "Cavity Mode: HQuad";
set xlabel "z [cm]"; set ylabel "H_{/Symbol \\136} [kA/m]";
plot file8 using (1e-1*$1):(1e-3*$2) notitle with lines ls 3;

unset multiplot;
if (term_type == 0) pause mouse "click on graph to cont.\n";

if (term_type == 1) set output "cavity_4.ps"
if (term_type == 2) set output "cavity_4.png"
if (term_type == 3) set output "cavity_4.jpg"

set size 1.0, 1.0; set origin 0.0, 0.0;
set title "Cavity Thin Lens Model: {/Symbol \\362}Eds,  {/Symbol \\362}Hds";
set xlabel "z [cm]";
set ylabel "[MV, MV/m, MA, MA/m]";
set yrange [-0.01:];
plot file9 using (1e-1*($1+$4/2.0)):6:(1e-1*$4) notitle with boxes ls 3;
if (term_type == 0) pause mouse "click on graph to cont.\n";
