Problem #2:

to compile the program, run
gcc -std=c99 -lm -o 3D 3D.c nrutil.c nrutil.h fgausselim.c

In the loop with n running up to nsteps, both
FTCS and Crank-Nicolson are implemented; comment
one out to run the other.  Crank-Nicolson is NOT
implemented with a sparse matrix and thus takes
forever to run with gaussian elimination; 
FTCS is more sprightly.

Problem #3:

- code sets Gaussian initial condition
- code has direchlet and/or periodic boundary conditions
- code has "source term"

Problem #5:

Using Crank-Nicolson the way it is implemented currently took
58750 ms for a 20x20x20 tensor; this was sufficiently bad
performance that I did not test further.

Please see included graph for FTCS running time. X axis is
"problem size" n=nx=ny=nz.
