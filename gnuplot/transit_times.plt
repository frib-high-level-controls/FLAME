# Script to plot the cavity transit times vs. 2pi/beta*lambda.
#
# Author: Johan Bengtsson.

home_dir = "/home/bengtsson/git_repos/jmbgsddb/python/cavity_model/";

# Terminal type is: 0 - X11, 1 - PS, 2 - PNG, 3 - JPG
term_type = 0;

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


T(arg) = -7.316972e+14*arg**9 + 2.613462e+14*arg**8 \
         -4.112154e+13*arg**7 + 3.739846e+12*arg**6 \
         -2.165867e+11*arg**5 + 8.280687e+09*arg**4 \
         -2.089452e+08*arg**3 + 3.354464e+06*arg**2 \
         -3.108322e+04*arg**1 + 1.256386e+02;

S(arg) = -6.079177e+14*arg**9 + 2.229446e+14*arg**8 \
         -3.605780e+13*arg**7 + 3.374738e+12*arg**6 \
         -2.013750e+11*arg**5 + 7.942886e+09*arg**4 \
         -2.070369e+08*arg**3 + 3.438044e+06*arg**2 \
         -3.299673e+04*arg**1 + 1.394183e+02;

c0 = 299792458e0;

#
f_QWR = 80.5e6;
lambda = c0/f_QWR;
# arg = 2.0*pi/(beta*1e3*lambda)

if (term_type == 1) set output "cavity_1.ps"
if (term_type == 2) set output "cavity_1.png"
if (term_type == 3) set output "cavity_1.jpg"

set title "QWR Transit Times: EFocus1";
set xlabel "2 {/Symbol \p}/{/Symbol \b \l}"; set ylabel "[sec]";
set xrange [0.025:0.08]
plot T(2.0*pi/(x*1e3*lambda)) title "T" with lines ls 1, \
     S(2.0*pi/(x*1e3*lambda)) title "S" with lines ls 3;

if (term_type == 0) pause mouse "click on graph to cont.\n";

if (term_type == 1) set output "cavity_2.ps"
if (term_type == 2) set output "cavity_2.png"
if (term_type == 3) set output "cavity_2.jpg"

set title "QWR Transit Times: EFocus1";
set xlabel "2 {/Symbol \p}/{/Symbol \b \l}"; set ylabel "[sec]";
plot home_dir."transit_times.dat" using 1:2 title "T" with lines ls 1, \
     home_dir."transit_times.dat" using 1:3 title "S" with lines ls 3;

if (term_type == 0) pause mouse "click on graph to cont.\n";
