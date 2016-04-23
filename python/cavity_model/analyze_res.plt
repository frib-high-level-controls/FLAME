# terminal 0, PostScript 1, PDF 2.
term_type = 2

home_dir_1 = "/home/johan/tlm_workspace/TLM_JB/"
home_dir_2 = "/home/johan/git_repos/jmbgsddb/build/src/"
home_dir_3 = "/home/johan/git_repos/jmbgsddb/python/cavity_model/"

file_name_1 = home_dir_1."MCSModelCenVec.txt"
file_name_2 = home_dir_1."MCSModelRf.txt"
file_name_3 = home_dir_1."MCSModelRmsVec.txt"

file_name_4 = home_dir_2."CenofChg.out"
file_name_5 = home_dir_2."long_tab.out"
file_name_6 = home_dir_2."BeamRMS.out"

file_name_7 = home_dir_3."analyze_res_1.dat"
file_name_8 = home_dir_3."analyze_res_3.dat"
file_name_9 = home_dir_3."analyze_res_2.dat"

if (term_type == 0) set terminal x11
if (term_type == 1) set terminal postscript landscape enhanced color solid \
  lw 2 "Times-Roman" 15
if (term_type == 2) set terminal pdfcairo enhanced color solid

set grid

# Left align legend.
set key Left

set style line 1 lw 1 lc rgb "red"
set style line 2 lw 1 lc rgb "dark-orange"
set style line 3 lw 1 lc rgb "blue"
set style line 4 lw 1 lc rgb "dark-green"
set style line 5 lw 1 lc rgb "purple"
set style line 6 lw 1 lc rgb "cyan"
set style line 7 lw 1 lc rgb "green"

if (term_type != 0) set output "Trans_Orbit_TLM.pdf"
set title "Transverse Orbit (TLM)"
set xlabel "s [m]"; set ylabel ""
plot file_name_1 using 1:2 title "x [mm]" with lines ls 3, \
     file_name_1 using 1:(1e3*$3) title "x' [mrad]" with lines ls 6, \
     file_name_1 using 1:4 title "y [mm]" with lines ls 1, \
     file_name_1 using 1:(1e3*$5) title "y' [mrad]" with lines ls 2
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Long_Orbit_TLM.pdf"
set title "Longitudinal Orbit (TLM)"
set xlabel "s [m]"; set ylabel ""
plot file_name_1 using 1:6 title "z [rad]" with lines ls 4, \
     file_name_1 using 1:7 title "z' [MeV/u]" with lines ls 5
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "E_Tot_TLM.pdf"
set title "Total Energy and Phase (TLM)"
set xlabel "s [m]"; set ylabel "[MeV/u]"; set y2label "[rad]"
set ytics nomirror; set y2tics;
plot file_name_2 using 1:2 title "E_tot [MeV/u]" with lines ls 3, \
     file_name_2 using 1:3 axis x1y2 title "{/Symbol j} [rad]" \
     with lines ls 1
if (term_type == 0) pause mouse "click on graph to cont.\n"

set ytics nomirror; unset y2tics; unset y2label

if (term_type != 0) set output "Trans_Beam_Size_TLM.pdf"
set title "Transverse RMS Beam Size (TLM)"
set xlabel "s [m]"; set ylabel ""
plot file_name_3 using 1:2 title "x RMS [mm]" with lines ls 3, \
     file_name_3 using 1:(1e3*$3) title "x' RMS [rad]" with lines ls 6, \
     file_name_3 using 1:4 title "y RMS [mm]" with lines ls 1, \
     file_name_3 using 1:(1e3*$5) title "y' RMS [rad]" with lines ls 2
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Long_Beam_Size_TLM.pdf"
set title "Longitudinal RMS Beam Size (TLM)"
set xlabel "s [m]"; set ylabel ""
plot file_name_3 using 1:6 title "z RMS [rad]" with lines ls 4, \
     file_name_3 using 1:7 title "z' RMS [MeV/u]" with lines ls 5
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Trans_Orbit_New_Code.pdf"
set title "Transverse Orbit (New Code)"
set xlabel "s [m]"; set ylabel ""
plot file_name_4 using 1:2 title "x [mm]" with lines ls 3, \
     file_name_4 using 1:(1e3*$3) title "x' [mrad]" with lines ls 6, \
     file_name_4 using 1:4 title "y [mm]" with lines ls 1, \
     file_name_4 using 1:(1e3*$5) title "y' [mrad]" with lines ls 2
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Long_Orbit_New_Code.pdf"
set title "Longitudinal Orbit (New Code)"
set xlabel "s [m]"; set ylabel ""
plot file_name_4 using 1:6 title "z [rad]" with lines ls 4, \
     file_name_4 using 1:7 title "z' [MeV/u]" with lines ls 5
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "E_Tot_New_Code.pdf"
set title "Total Energy and Phase (New Code)"
set xlabel "s [m]"; set ylabel "[MeV/u]"; set y2label "[rad]"
set ytics nomirror; set y2tics;
plot file_name_5 using 3:4 title "E_tot [MeV/u]" with lines ls 3, \
     file_name_5 using 3:5 axis x1y2 title "{/Symbol j} [rad]" \
     with lines ls 1
if (term_type == 0) pause mouse "click on graph to cont.\n"

set ytics nomirror; unset y2tics; unset y2label

if (term_type != 0) set output "Trans_Beam_Size_New_Code.pdf"
set title "Transverse RMS Beam Size (New Code)"
set xlabel "s [m]"; set ylabel ""
plot file_name_6 using 1:2 title "x RMS [mm]" with lines ls 3, \
     file_name_6 using 1:(1e3*$3) title "x' RMS [rad]" with lines ls 6, \
     file_name_6 using 1:4 title "y RMS [mm]" with lines ls 1, \
     file_name_6 using 1:(1e3*$5) title "y' RMS [rad]" with lines ls 2
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Long_Beam_Size_New_Code.pdf"
set title "Longitudinal RMS Beam Size (New Code)"
set xlabel "s [m]"; set ylabel ""
plot file_name_6 using 1:6 title "z RMS [rad]" with lines ls 4, \
     file_name_6 using 1:7 title "z' RMS [MeV/u]" with lines ls 5
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Trans_Orbit_Diff.pdf"
set title "Transverse Orbit Difference Between New Code and TLM"
set xlabel "s [m]"; set ylabel ""
plot file_name_7 using 1:2 title "x [mm]" with lines ls 3, \
     file_name_7 using 1:(1e3*$3) title "x' [mrad]" with lines ls 6, \
     file_name_7 using 1:4 title "y [mm]" with lines ls 1, \
     file_name_7 using 1:(1e3*$5) title "y' [mrad]" with lines ls 2
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Long_Orbit_Diff.pdf"
set title "Longitudinal Orbit Difference Between New Code and TLM"
set xlabel "s [m]"; set ylabel ""
plot file_name_7 using 1:6 title "z [rad]" with lines ls 4, \
     file_name_7 using 1:7 title "z' [MeV/u]" with lines ls 5
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "E_Tot_Diff.pdf"
set title "Total Energy and Phase Difference Between New Code and TLM"
set xlabel "s [m]"; set ylabel "[MeV/u]"; set y2label "[rad]";
set ytics nomirror; set y2tics;
plot file_name_8 using 1:2 title "E_tot" with lines ls 3, \
     file_name_8 using 1:3 axis x1y2 title "{/Symbol j}" with lines ls 1
if (term_type == 0) pause mouse "click on graph to cont.\n"

set ytics nomirror; unset y2tics; unset y2label

if (term_type != 0) set output "Trans_Beam_Size_Diff.pdf"
set title "Transverse RMS Beam Size Difference Between New Code and TLM"
set xlabel "s [m]"; set ylabel ""
plot file_name_9 using 1:2 title "x RMS [mm]" with lines ls 3, \
     file_name_9 using 1:(1e3*$3) title "x' RMS [mrad]" with lines ls 6, \
     file_name_9 using 1:4 title "y RMS [mm]" with lines ls 1, \
     file_name_9 using 1:(1e3*$5) title "y' RMS [mrad]" with lines ls 2
if (term_type == 0) pause mouse "click on graph to cont.\n"

if (term_type != 0) set output "Long_Beam_Size_Diff.pdf"
set title "Longitudinal RMS Beam Size Difference Between New Code and TLM"
set xlabel "s [m]"; set ylabel ""
plot file_name_9 using 1:6 title "z RMS [rad]" with lines ls 4, \
     file_name_9 using 1:7 title "z' RMS [MeV/u]" with lines ls 5
if (term_type == 0) pause mouse "click on graph to cont.\n"
