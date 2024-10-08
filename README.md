# Generic L-function calculator

A generic L-function calculator for motivic L-functions, originally written by Dave J. Platt.

For details of the methods we refer to:
 - *Computing motivic L-functions*, (under preparation) by Jonathan W. Bober, Andrew R. Booker, Edgar Costa, Min Lee, David J. Platt, and Andrew V. Sutherland.


## Dependencies
It majorly depends on:
 - [flint - A C Library for doing number theory.](https://flintlib.org)
 - [primesieve - Fast prime number generator](https://github.com/kimwalisch/primesieve)
 
which depend on:

 - [GMP: GNU Multiple Precision Arithmetic Library](https://gmplib.org/) (for FLINT)
 - [MPFR: GNU Multiple Precision Floating-Point Reliably](http://www.mpfr.org/) (for FLINT)


## Installation

1. Download `glfunc`
```
git clone https://github.com/djplatt/glfunc.git
```

2. Change your working directory and run the configure file
```
cd glfunc
./configure <options>
```

3. If any of these libraries are installed in some other location than the default path `/usr/local`, pass `--with-primesieve=...`, `--with-gmp=...`, `--with-mpfr=...`, or `--with-flint=...` with the correct path to configure (type ./configure --help to show more options).

4. Compile everything by doing
```
make
```

5. You can use it as shared library, the inteface is defined in 
interface is defined in `include/glfunc.h` and the shared library is provided in
`build/liblfun.so`




