set terminal pdf enhanced font 'cmr, 26' size 11.69 in, 8.27 in
set output 'pot.pdf'
set title 'SAE potential - TISE wavefunctions'

set xlabel 'r (a.u.)'
set ylabel 'cos{/Symbol q}'
set zlabel 'V(r,{/Symbol q})'

set xrange [0:60]
set zrange [-3:0]
set xtics 0,10,60
set ytics -1,0.4,1

set view 75,345

sp \
	"POT.DAT" u 1:2:3 w l t ''
