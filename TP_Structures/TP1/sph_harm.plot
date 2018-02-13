set mapping spherical ; set hidden3d nooffset
unset key
se xr [-.5:.5]
se yr [-.5:.5]
se zr [-.5:.5]
set title 'l=0, m=0'
spl 'sph_harm_0' u 2:(pi/2-$1):(abs($3))
pause -1
set title 'l=1, m=0'
spl 'sph_harm_1' u 2:(pi/2-$1):(abs($3))
pause -1
set title 'l=1, m=1'
spl 'sph_harm_1' u 2:(pi/2-$1):(abs($4))
unset multiplot
