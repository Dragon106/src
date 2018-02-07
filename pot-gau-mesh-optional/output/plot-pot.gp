set terminal pdf enhanced font 'cmr, 26' size 11.69 in, 8.27 in
set output 'pot.pdf'
set title 'SAE potential - wf from TISE'

set xlabel 'r (a.u.)'
set ylabel 'cos{/Symbol q}'
set zlabel 'V(r,{/Symbol q})'

set xrange [0:20]
set zrange [-5:0]
set xtics 0,2,20
set ytics -1,0.4,1

set view 75,345
sp \
	"POT.DAT" u 1:2:3 w l t ''
