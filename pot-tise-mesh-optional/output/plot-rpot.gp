set terminal pdf enhanced font 'cmr, 26' size 11.69 in, 8.27 in
set output 'rpot.pdf'
set title 'effective charge - wf from TISE'

set xlabel 'r (a.u.)'
set ylabel 'r x V(r)'

set xrange [0:20]
set yrange [-5:0]
set xtics 0,2,20
set ytics -5, 0.5, 0


plot \
	"RPOT.DAT" u 1:2 t 'cos {/Symbol q} = -1' w l lw 3, \
	"RPOT.DAT" u 1:3 t 'cos {/Symbol q} =  0' w l lw 3
