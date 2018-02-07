#!/usr/bin/gnuplot

set terminal pdf enhanced font 'times new roman, 32' size 11 in, 4 in
set output 'rate_cep180.pdf'
set multiplot layout 1, 3

set grid
set border lw 5
set style line 1 lt 1 lw 5 lc 'black'
set style line 2 lt 2 lw 7 lc 'red'
set key left top Left reverse
#---------
set macros
POS = 'graph 0.1, 0.85
LMAR = 'set lmargin at screen 0.10; set rmargin at screen 0.34'
MMAR = 'set lmargin at screen 0.42; set rmargin at screen 0.66'
RMAR = 'set lmargin at screen 0.74; set rmargin at screen 0.98'
#---------

dirsae  = '../output-sae-test-cep180/'
dirsaep = '../output-saep-test-cep180/'
file    = 'rate.dat'
w0  = 0.057
tau = 2*pi/w0

#--- graph 1
@LMAR
set label '(a) θ = 0^0' at graph 0.1, 0.65
set xlabel 'Time (optical cycles)'
set ylabel 'Ionization rate (x 10^{-2})' offset 1

set xrange [0.5:2]
set yrange [0:0.8]
set format y '%.2f'
set xtics 0.25, 0.5
set ytics 0.2

plot dirsae.'0/'.file  u ($1/tau):($2*100) w l ls 1 t 'SAE', \
	 dirsaep.'0/'.file u ($1/tau):($2*100) w l ls 2 t 'SAE+P', \

#--- graph 2
@MMAR
unset label
unset ylabel
unset key
set label '(b) θ = 90^0 (x3)' at @POS
plot dirsae.'90/'.file  u ($1/tau):($2*300) w l ls 1 t 'SAE', \
	 dirsaep.'90/'.file u ($1/tau):($2*300) w l ls 2 t 'SAE+P', \

#--- graph 3
@RMAR
unset label
unset ylabel
set label '(c) θ = 180^0' at @POS
plot dirsae.'180/'.file  u ($1/tau):($2*100) w l ls 1 t 'SAE', \
	 dirsaep.'180/'.file u ($1/tau):($2*100) w l ls 2 t 'SAE+P', \
