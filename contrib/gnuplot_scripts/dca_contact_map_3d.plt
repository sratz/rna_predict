# visualize quality of dca contacts
# data files are assumed to be in the following format:
# <nuc1> <nuc2> <dca_score> <native_rmsd>
# native rmsd is plotted on the z axis.
# dca score is visualized using color.

# maximum number of dca predictions to plot
limit = 100

# which dca input file to plot
input = 'AVrnaDCA.txt'

# full contact map
full = 'all.txt'


set size square
set view equal xy
set xrange [0:]
set yrange [0:]

set xlabel "nuc 1"
set ylabel "nuc 2"
set zlabel "native rmsd"

splot full using 1:2:4:(0) with points title "full contact map" lc rgb 'black' lt 1, \
      full using 2:1:4:(0) with points title "full contact map" lc rgb 'black' lt 1, \
      input every ::0::limit using 2:1:4:3 with points palette title input." (max ".limit.")"

pause -1
