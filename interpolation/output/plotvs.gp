set term pdf enhanced font 'cmr, 26' size 10 in, 5.5 in
set output "rV-iter.pdf"
set multiplot

set xrange [0:309] 
set yrange [0:177]

unset tics
unset border

set lmargin at screen 0.175
set rmargin at screen 0.9
set bmargin at screen 0.15
set tmargin at screen 0.9
plot "GraphCO2.jpg" binary filetype=jpg w rgbimage t ''

set xrange [0:20]
set yrange [-3:-0.5]
set xtics 0, 2, 20
set ytics -3, 0.5, -0.5

set xlabel "r (a.u.)"
set ylabel "r x V(r)"
set key cent

plot \
	 "RPOT.DAT" u 1:2 w l lw 4 lc 'green' title "one iteration"
