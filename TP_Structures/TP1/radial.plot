# commandes gnuplot pour tracer les fonctions radiales
# ce fichier de commandes fait tous les tracés demandés à la suite, mais
# on peut les faire un à un en interactif
# il peut etre invoque dans gnuplot par: load 'radial.plot'
# les '#' indiquent les commentaires...
set xzeroaxis	# tracer l'axe des x
set xl 'r'   	# label de l'axe des x
set title 'Fonctions radiales'
plot 'radial.001' ti 'n=1', 'radial.002' ti 'n=2, l=0', 'radial.002' u 1:3 ti 'n=2, l=1'	# on peut en faire autant qu'on veut
pause -1	# tapez sur 'Return' pour continuer
set title 'Fonctions radiales au carre'
plot 'radial.001' u 1:($2**2) ti 'n=1', 'radial.002' u 1:($2**2) ti 'n=2, l=0', 'radial.002' u 1:($3**2) ti 'n=2, l=1' 
pause -1	# tapez sur 'Return' pour continuer
set title 'Fonctions radiales au carre'
set yr[0:0.5]	# echelle de l'axe des y
plot 'radial.001' u 1:($2**2) ti 'n=1', 'radial.002' u 1:($2**2) ti 'n=2, l=0', 'radial.002' u 1:($3**2) ti 'n=2, l=1' 
pause -1	# tapez sur 'Return' pour continuer
set title 'Probabilités radiales'
set yr[0:0.6]	# echelle de l'axe des y
plot 'radial.001' u 1:($2**2*$1**2) ti 'n=1', 'radial.002' u 1:($2**2*$1**2) ti 'n=2, l=0', 'radial.002' u 1:($3**2*$1**2) ti 'n=2, l=1' 
