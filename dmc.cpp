
#include <array>
#include <functional>
#include <cmath>
#include "dmc.hpp"


double helium( std::array<double, 6> &x ) {

   return -2/std::sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] )
          -2/std::sqrt( x[3]*x[3] + x[4]*x[4] + x[5]*x[5] )
          + 1/std::sqrt( (x[0] - x[3])*(x[0] - x[3])
                       + (x[1] - x[4])*(x[1] - x[4])
                       + (x[2] - x[5])*(x[2] - x[5]));
}


double V( std::array<double, 1> &x ) {

   return x[0]*x[0] * 0.5;
}


int main() {

   const int d = 3;                   //Spatial dimension
   const int s = 2;                   //Number of particles
   const int dp = d*s;                //Effective dimension.
   int N0 = 4000;                     //Intial numberof points.
   const int maxN = 10000;            //Maximum number of points.
   int nb = 2;                        //Number of bins for wave function histogram.
   double xmin = -10;                 //Lower bound of grid
   double xmax = 10;                  //Upper bound of grid
   double timeStep = 0.001;           //Time step.
   double relaxTime = 200;            //Time to wait for convergence.
   double alpha = 0.95;               //Update parameter for energy.
   std::array<double, 6> x0 = {0.0, 0.0, 1.0, 0.0, 0.0, -1.0};  //Initial position of particles.
   std::function<double (std::array<double, dp> &)> potential = helium; //The Potential.

   runSimulation<maxN, dp>( N0, nb, xmin, xmax, timeStep, relaxTime, alpha, x0, potential, 8 );

   return 0;
}

