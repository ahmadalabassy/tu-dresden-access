# Numerical Methods

Implement a program in the the programming language Fortran for:
- Solving a **quadratic system of linear equations** using different methods.
- Analysing the efficiency of the methods.

## Task 1
Implement a program for reading the elements of a quadratic system of linear equations of any given size from a text file, and to check whether the preconditions for the application of Cholesky decomposition are fulfilled using **Laplace expansion** for the determinant criterion.

## Task 2
Extend the program for an efficient solution of the system using **Cholesky decomposition** and **Gaussian elemination**.

## Task 3
**Analyse** and **compare** the required **computing times (CPU times)** for the individual steps of the two methods depending on the **size** of the coefficient matrix. **Visualize the results** by plotting the computing times against the size of the coefficient matrix.

### Note
- To generate a system automatically, edit the system size at line 18 in `system_generator/00_system_generator.f90`. The default size is preset to 10, open `./system_generator/compile.bat` in terminal to create an executable, open `a.exe`, copy the generated `system.txt` to working directory.
- Alternatively, you can test a system of your own choosing, in a file named `system.txt`, that should be added to working directory. First line is a positive integer (`n`) that sets the system size, followed by `n` lines (rows), each line (row) containing `n` columns of coefficients, each coefficient is seperated from the coefficient prior by whitespace. Coefficient matrix rows are to be followed by n lines, each has a numerical value representing that of the right-hand-side vector. For further guidance, check the auto-generated format of `system.txt`.
- For the third task, system size can be changed for Gaussian elimination at `line 20` in `Task_3/01_quad_system.f90`. Size is capped to `10` for Cholesky decomposition; to maintain tolerable processing time at testing, but could be changed from the subsequent line, `line 21` of the same script. To plot computing times against the size of the coefficient matrix, add gnuplot binary to `./task_3/binary/`. Binaries can be downloaded from: [http://gnuplot.info/download.html](http://gnuplot.info/download.html).