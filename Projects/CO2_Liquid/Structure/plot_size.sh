set style fill solid
set xrange [4:96]

set multiplot layout 4,1
p "9.4/3000K/size_proba_1.8stride1.dat" u 1:2 with boxes title "42GPa"
p "9.2/3000K/size_proba_1.8stride1.dat" u 1:2 with boxes title  "50GPa"
p "9.0/3000K/size_proba_1.75.dat" u 1:2 with boxes title  "55GPa"
set xlabel "Size (atoms)"
set key left
p "8.82/3000K/size_proba_1.75.dat" u 1:2 with boxes title  "65GPa"

