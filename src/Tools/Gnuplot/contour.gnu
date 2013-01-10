#===============================================
# data section to be changed by the user
#===============================================
xmin=-4.
xmax=6.
ymin=-5.
ymax=5.
zmin=-0.015
zmax=0.3
ncontour=30
infile="example.dat"
#
#===============================================
# define line styles to be used with ls in splot
#===============================================
# define a linestyle for splot (used with ls)
set style line 1 lt 1 lc rgb "black" lw 1
set palette rgbformula 22,13,-31
#
#==================================
# data related statements
#===================================
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set zrange [zmin:zmax]
# sample input data onto a 60x60 grid
set dgrid3d  60,60,1
set cntrparam levels auto ncontour
#
#==================================
# surface plot
#===================================
set surface   
# remove axes
unset border
#  remove tics from the axes
unset xtics
unset ytics
unset ztics
# no title written
set key off
set contour both
set view map
set size square
splot  infile with pm3d 
