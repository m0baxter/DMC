
#include <array>
#include <functional>
#include <cmath>
#include "dmc.hpp"


double V( std::array<double, 2> &x ) {

   return 0.5*( x[0]*x[0] + x[1]*x[1] );
}


int main() {

   const int d = 2;                   //Spatial dimension
   const int s = 1;                   //Number of particles
   const int dp = d*s;                //Effective dimension.
   int N0 = 4000;                     //Intial numberof points.
   const int maxN = 10000;            //Maximum number of points.
   int nb = 2;                        //Number of bins for wave function histogram.
   double xmin = -10;                 //Lower bound of grid
   double xmax = 10;                  //Upper bound of grid
   double timeStep = 0.001;           //Time step.
   double relaxTime = 200;            //Time to wait for convergence.
   double alpha = 0.95;               //Update parameter for energy.
   std::array<double, 2> x0 = {0.0, 0.0 };  //Initial position of particles.
   std::function<double (std::array<double, dp> &)> potential = helium; //The Potential.

   runSimulation<maxN, dp>( s, N0, nb, xmin, xmax, timeStep, relaxTime, alpha, x0, potential, 8 );

   return 0;
}

