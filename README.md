# Molekel
## General Information
Molekel is a simple, slow, and inefficient program to calculate molecular energies
(currently only RHF). The integrals are evaluated using the McMurchie-Davidson scheme.
The code is inspired by the work of Joshua Goings (https://joshuagoings.com/2017/04/28/integrals/)
, but is written in Fortran. Other sources are the books Molecular Electronic-Structure
Theory by Helgaker *et. al.* and Modern Quantum Chemistry by Attila Szabo and Neil S.
Ostlund. The cmake build is based on the work of Leonard Reuter (https://github.com/Leonard-Reuter/cmake-fortran-template)
 and some modules of the code were copied from Amolqc (https://github.com/luechow-group/Amolqc).
The program is intended to serve for educational purposes, however, feel free to do what you
want with it!
## Build Instructions
In order for this code to work, you need the gfortran compiler, and the lapack and blas 
libraries.


Clone the repository to your desired location, then use the following commands to build
the program.
```
export MOLEKEL=<path/to/molekel>
cd $MOLEKEL
mkdir build
cd build
cmake ..
make molekel
export PATH=$PATH:MOLEKEL/build/bin
```
And you are all set! You can try your build with the example.in in the test folder. 
You have to pipe the output in a file, else it will be printed in your terminal.
```
molekel $MOLEKEL/test/example.in > <output>
```
## Documentation
The documentation can be found in the docs directory.
## Remarks
This program is work in progress! More features will be added from time to time! 
Take the results with a grain of salt, since I have only done little testing.
If you encounter any problems, feel free to report them. 

