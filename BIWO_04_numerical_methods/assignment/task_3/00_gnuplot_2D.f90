! Program for producing two dimensional graphs with gnuplot
!**************************************************************************!
! Prerequisite: Downloading of gnuplot from www.gnuplot.info and copying
!               the folder "binary" in the working directory
!**************************************************************************!
!**************************************************************************!

subroutine gnuplot_2D(x_range,y_range)

implicit none

integer x_range
real y_range

! Opening, writing and closing the gnuplot-input file for total runtime
!**************************************************************************!
open(10,file='plot_2D.plt')
write(10,*) 'set title "Gaussianx10^-4 and Cholesky runtime comparison"'
write(10,*) 'set xrange [0:',x_range,']'
write(10,*) 'set yrange [0:',y_range,']'
write(10,*) 'set xlabel "n size of system"'
write(10,*) 'set ylabel "t in seconds"'
write(10,*) 'plot "data_cholesky.txt" with lines,"data_gaussian.txt" with lines, "preconditions_cholesky.txt" with lines, &
"triangular_matrix_cholesky.txt" with lines, "triangular_matrix_gaussian.txt" with lines, "solution_cholesky.txt" with lines,&
 "solution_gaussian.txt" with lines'
write(10,*) 'pause -1 "Hit return to continue"'
close(10)

! Calling the wgnuplot.exe and handing over the gnuplot-input file
!**************************************************************************!
call system('binary\wgnuplot plot_2D.plt') !here we can apply DOS Commands as we're using windows OS

end subroutine