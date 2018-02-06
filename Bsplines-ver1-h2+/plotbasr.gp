set term pdf enhanced font 'cmr, 26' size 11 in, 8 in
set output 'bspl.pdf'

set title 'B-splines order k = 9'

set xlabel 'r (a.u.)'
set ylabel 'B_i(r)'
set xrange [0:40]

plot for [i=1:30] 'basisr.dat' u 1:i+1 w l t 'B'.i
