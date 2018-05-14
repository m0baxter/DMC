
#include <array>
#include <functional>
#include <cmath>
#include "dmc.hpp"


double V( std::array<double, 1> &x ) {

   return 0.5 * x[0] * x[0];

}

int main() {

   const int d = 1;                   //Spatial dimension
   const int s = 1;                   //Number of particles
   const int dp = d*s;                //Effective dimension.
   int N0 = 1000;                     //Intial numberof points.
   const int maxN = 2000;            //Maximum number of points.
   int nb = 200;                        //Number of bins for wave function histogram.
   double xmin = -5;                 //Lower bound of grid
   double xmax = 5;                  //Upper bound of grid
   double timeStep = 0.001;           //Time step.
   double relaxTime = 50;            //Time to wait for convergence.
   double alpha = 5.0;               //Update parameter for energy.

   std::array<double, 1> x0 = { uniform(-1, 1) };  //Initial position of particles.
   std::function<double (std::array<double, dp> &)> potential = V; //The Potential.

   runSimulation<maxN, dp>( s, N0, nb, xmin, xmax, timeStep, relaxTime, alpha, x0, potential, 1, true, true );

   return 0;
}

