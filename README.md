# DMC
An implementation of the Diffusion Monte Carlo method (DMC) for calculating the ground-state energy and
wave function of a quantum system. This code is based on the algorithm described in Ref. [1].

## Compilation

This program relies on the C++ functionality only present in C++11 or newer. Notably it makes use of
the STL thread library. The code can be compiled with a command of the form
```
g++ -O3 -pthreads dmc.cpp -o dmc
```
or the equivalent for your compiler of choice.

## Running

DMC calculations may be run by passing appropriate arguments to the template function
```c++
template< int maxN, int dp >
void runSimulation( const int pNum, const int N0, const int nb, const double xmin, const double xmax,
                    const double timeStep, const double relaxTime, const double alpha,
                    std::array<double, dp> &x0, std::function<double (std::array<double, dp> &)> &V,
                    const int numThreads = 1, const bool printWF = false, const bool printDist = false )
```
The parameters are:

-`maxN`: The maximum number of points in the distribution.

-`dp`: The effective dimension of the problem, the number of particles multiplied by the spatial
dimension.

-`pNum`: The number of particles.

-`N0`: The initial number of points in the distribution.

-`nb`: The number of buckets used for determining the wave function.

-`xmin` and `xmax`: Defines the support of the wave function [`xmin`,`xmax`]^`dp`.

-`timeStep`: The time step in atomic units.

-`alpha`: Update parameter for the self consistent energy calculation.

-`x0`: The initial position of all points in the distribution.

-`V`: The potential.

-`numThreads`: The number of concurrent simulations run.

-`printWF`: A boolean flag that controls output of the wave function

-`printDist`: A boolean flag that controls output of the distribution of points.

Two other functions that may be of interest to the user are
```c++
double uniform( const double a, const double b)
```
```c++
double gaussian( const double mean, const double stdDev)
```
which can be used to generate random numbers drawn from uniform/Gaussian distributions. These may be
useful for generating initial values for `x0`.

The sample parameters provided in `dmc.cpp` will calculate the ground state of the quantum harmonic
oscillator.

## References
[1] I. Kosztin, B. Faber, and K. Schulten, [Introduction to the Diffusion Monte Carlo
Method](https://arxiv.org/abs/physics/9702023). arXiv:physics/9702023 [physics.comp-ph].

