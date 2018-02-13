set pm3d map
unset key
set xl 'x'
set yl 'z'
set xr [-15:15]
set yr [0:15]
set title 'n=1'
splot 'hydro.1_0' u ($1*sin($2)):($1*cos($2)):3
pause -1
set title 'n=2, l=0'
splot 'hydro.2_0' u ($1*sin($2)):($1*cos($2)):3
pause -1
set title 'n=2, l=1'
splot 'hydro.2_1' u ($1*sin($2)):($1*cos($2)):3
pause -1

