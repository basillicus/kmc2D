set terminal postscript enhanced color "Helvetica,18"


set key bottom right
set format x "%5.3f"

set output "statistics-1.eps"
# plot "statistics.kmc" u 1:3 w l t 'n Cis'     ,\
#      "statistics.kmc" u 1:4 w l t 'n L trans' ,\
#      "statistics.kmc" u 1:5 w l t 'n D trans' 

set output "statistics-moves.eps"
plot "statistics.kmc" u 1:($7/1000)  w l t 'Free Diff/100 '  ,\
     "statistics.kmc" u 1:($8/10)    w l t 'Hbond Diff/100' ,\
     "statistics.kmc" u 1:9          w l t 'Free Isom'   ,\
     "statistics.kmc" u 1:10         w l t 'Hbond Isom' 

set output "statistics-interactions.eps"
plot "statistics.kmc" u 1:11 w l t 'N H-bond '  ,\
     "statistics.kmc" u 1:12 w l t 'N Neighb.' 
