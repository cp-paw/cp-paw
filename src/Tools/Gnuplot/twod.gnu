#
#===============================================
# data section to be changed by the user
#===============================================
xmin=-4.
xmax=4.
ymin=-4.
ymax=4.
zmin=-5.
zmax=0.
rot_x=30
rot_z=20
scale=1.0
scale_z=2.5
infile="out.dat"
#
#===============================================
# define line styles to be used with ls in splot
#===============================================
# define a linestyle for splot (used with ls)
set style line 1 lt 1 lc rgb "black" lw 1
# map hight values to colors
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
#
#==================================
# surface plot
#===================================
set surface   
set hidden3d
set data style lines
# places the zero of the z-axis into the xy plane
set xyplane at 0.
# remove axes
unset border
#  remove tics from the axes
unset xtics
unset ytics
unset ztics
# no title written
set key off
# place contours onto the surface
set pm3d explicit hidden3d 1
#
#  angle, angle, overall scale, scale z-axis
set view rot_x,rot_z,scale,scale_z
#
splot  infile with pm3d
