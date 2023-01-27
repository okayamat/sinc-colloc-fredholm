set term postscript eps enhanced "Times-Roman" 20
set output "example3_t.eps"
set size 0.82
set logscale y
set key spacing 3
set xrange [0:1]
set yrange [1e-16:1]
set xlabel "{/Times-Roman=24 computation time [s]}"
set ylabel "{/Times-Roman=24 maximum error}"
plot "SE_orig_ex3.dat" using 3:2 w lp title "Original SE-Sinc" lt 3 pt 2 ps 2, "SE_new_ex3.dat" using 3:2 w lp title "New SE-Sinc" lt 3 pt 8 ps 2, "DE_orig_ex3.dat" using 3:2 w lp title "Original DE-Sinc" lt 4 pt 1 ps 2, "DE_new_ex3.dat" using 3:2 w lp title "New DE-Sinc" lt 4 pt 6 ps 2
