# Introduction


A set of C/C++ utilities to solve geodesics using boost or GSL in metric spaces.

# Boost

The more developed geodesic solver tool uses Boost.  The main classes in use are **metric** and **geodesicSolver**.  The **metric** class takes in C++ function objects which represent the metric, inverse metric, and inverse metric derivative (see the doxygen documentation for their use).  The **geodesicSolver** solves the geodesic system given a metric class and initial velocities and positions.

The functions the metric class uses are of the form

```c++
std::function<double(int , int,std::vector<double>)>
```
Where the first two integers are the standard metric indeces `g_ij` and `std::vector<double>` is the coordinates in increasing order.


## Optimal accuracy and runtime

In order to get the most optimal runtime, supplying the inverse metric and the derivative of the inverse metrix `d g^ij / d x^k` as a function of the form `std::function<double(int i, int j, int k std::vector<double> coords)>` is needed.

Supplying only the metric to the constructor will requre the computation of numerical derivatives via finite differences.

# Gsl

This library has a collection of metric files and a solver which can be used with GSL to solve geodesic equations.  

The file solve_geodesics.cc implements the rkp8d (dormand prince method) to solve the geodesic equations for 
an infalling observer.  Feel free to take or modify this code for whichever purposes you require.

To compile you need to have gsl installed and type the command

gcc -static solve_geodesics.cc -lgsl -lgslcblas -lm

or

gcc solve_geodesics.cc -lgsl -lgslcblas -lm

depending on your gsl installation.