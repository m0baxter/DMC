
#ifndef DMC_HPP
#define DMC_HPP

#include <iostream>
#include <array>
#include <utility>
#include <functional>
#include <cmath>
#include <algorithm>
#include <map>
#include <thread>
#include <fstream>
#include <sstream>
#include <random>


double gaussian( const double mean, const double stdDev ) {

   thread_local std::mt19937 engine( std::random_device{}() );
   std::normal_distribution<double> dist( mean, stdDev );
   return dist(engine);
}

double uniform( const double a, const double b ) {

   thread_local std::mt19937 engine( std::random_device{}() );
   std::uniform_real_distribution<double> dist( a, std::nextafter(b, std::numeric_limits<double>::max()) );
   return dist(engine);
}


template< int maxN, int dp > //dp = d * num particles
std::pair< std::array<int, maxN>, std::array< std::array<double, dp>, maxN> > initialize( const int N0, std::array<double, dp> &x0 ) {
   /*Initializes the array of flags and the array of points representing the ground state wave function
    * distribution.*/

   std::array<int, maxN> flags;
   std::array<std::array<double, dp>, maxN> points;

   for ( int i = 0; i < N0; ++i ) {

      flags[i]  = 1;
      points[i] = x0;
   }

   for ( int i = N0; i < maxN; ++i ) {

      flags[i] = 0;
      points[i].fill(0.0);
   }

   return std::make_pair( flags, points );
}

template< int maxN, int dp >
void walk( const double timeStep, std::array<int, maxN> &flags, std::array< std::array<double, dp>, maxN> &points ) {

   for ( int i = 0; i < maxN; ++i ) {

      if ( flags[i] ) {
         for ( int j = 0; j < dp; ++j ) {
            points[i][j] += timeStep * gaussian(0.0, 1.0);
         }
      }
   }
}

template< int dp >
int spawnNumber( std::array<double, dp> &x, const double E, const double u, double dt,
                 std::function<double (std::array<double, dp> &)> &V ) {

  return std::min( (int)std::floor(std::exp( (E - V(x))*dt ) + u), 3 );
}


template< int maxN, int dp >
int branch( const int N, const double dt, const double alpha, std::array<int, maxN> &flags,
            std::array< std::array<double, dp>, maxN> &points, std::vector<double> &Es,
            std::function<double (std::array<double, dp> &)> &V ) {

   int nextN = N;

   for (int i = 0; i < maxN; ++i ) {

      if ( flags[i] == 1 ) {

         int m = spawnNumber<dp>( points[i], Es.back(), uniform(0.0, 1.0), dt, V );

         if ( m ) {
            for (int j = 0; j < m - 1; ++j ) {
               int newPoint = std::distance( flags.begin(), std::find(flags.begin(), flags.end(), 0) );

               if ( newPoint < maxN ) {
                  points[newPoint] = points[i];
                  ++nextN;

                  if ( newPoint < i ) {
                     flags[newPoint] = 1;
                  }
                  else {
                     flags[newPoint] = 2;
                  }
               }
            }
         }
         else {
            flags[i] = 0;
            --nextN;
         }

      }
      else if ( flags[i] == 2 ){
         flags[i] = 1;
      }
   }

   Es.push_back( Es.back() + alpha*(1 - (1.0 * nextN)/N ) );

   return nextN;
}

int bucketNumber( const int nb, const double a, const double b, const double x ) {
   /*Returns the index of the bucket that x falls in given that the interval
    * [a,b] is divided into nb + 1 buckets.*/

   if ( x < a ) {
      return 0;
   }
   else if ( x > b ) {
      return nb;
   }
   else {
      return std::floor( (x - a) * nb / (b -  a) );
   }
}

template< int maxN, int dp >
void count( const int nb, const double a, const double b, std::map< std::array<int, dp>, int> &counts,
            std::array<int, maxN> &flags, std::array< std::array<double, dp>, maxN> &points ) {
   /*Counts the number of points in each bucket, updates the current counts.*/

   for ( int i = 0; i < maxN; ++i ) {

      if ( flags[i] ) {

         std::array<int, dp> key;

         for ( int d = 0; d < dp; ++d ) {
            key[d] = bucketNumber( nb, a, b, points[i][d] );
         }

         if ( counts.count(key) ) {
            counts[key] += 1;
         }
         else {
            counts[key] = 1;
         }
      }
   }
}

template< int dp >
void writeWavefunction( const double a, const double b, const int nb,
                        std::map< std::array<int, dp>, int> &counts ) {

   std::stringstream ss;
   ss << "./wavefunction/wf-" << std::this_thread::get_id() << ".txt";
   std::ofstream writeFile( ss.str().c_str() );

   int total = 0;
   double s = (b - a)/nb;

   for ( auto &p : counts ) {
      total += p.second;
   }

   for ( auto &p : counts ) {
      for ( int i = 0; i < dp; ++i ) {
         writeFile << a + p.first[i] * s << "   ";
      }

      writeFile << (1.0*p.second)/total << "\n";
   }

   writeFile.close();
}

void writeEnergy( const std::vector<double> &Es ) {

   std::stringstream ss;
   ss << "./energy/energy-" << std::this_thread::get_id() << ".txt";
   std::ofstream writeFile( ss.str().c_str() );

   for ( auto &E : Es ) {
      writeFile << E << "\n";
   }

   writeFile.close();
}

template< int maxN, int dp >
void writeDistribution( const int pNum, std::array<int, maxN> &flags,
                        std::array< std::array<double, dp>, maxN> &points ) {

   std::stringstream ss;
   ss << "./distros/dist-" << std::this_thread::get_id() << ".txt";
   std::ofstream writeFile( ss.str().c_str() );

   int d = dp/pNum;

   for ( int i = 0; i < maxN; ++i ) {

      if ( flags[i] ) {
         for ( int p = 0; p < pNum; ++p ) {
            for ( int j = 0; j < d; ++j ) {
               writeFile << points[i][ p * d + j ] << "  ";
            }

            writeFile << "\n";
         }
      }
   }

   writeFile.close();
}

template< int maxN, int dp >
void diffusionMC( const int pNum, const int N0, const int nb, const double xmin, const double xmax,
                  const double timeStep, const double relaxTime, const double alpha,
                  std::array<double, dp> x0, std::function<double (std::array<double, dp> &)> &V,
                  const bool printWF, const bool printDist ) {
   /*Executes the main DMC loop.*/

   if ( N0 > maxN ) {
      std::cout << "Initial number of points must be less than maximum." << std::endl;
      return;
   }

   if ( xmin >= xmax ) {
      std::cout << "Invalid interval: xmin >= xmax." << std::endl;
      return;
   }

   std::vector<double> Es = { V(x0) };
   std::map< std::array<int, dp>, int> counts;
   int currentN = N0;
   double tau = 0;
   double dt = std::sqrt( timeStep );

   auto p = initialize<maxN, dp>( N0, x0 );
   std::array<int, maxN> flags = p.first;
   std::array< std::array<double, dp>, maxN> points = p.second;

   while ( tau < 2 * relaxTime ) {

      walk<maxN, dp>( dt, flags, points );

      currentN = branch<maxN, dp>( currentN, timeStep, alpha, flags, points, Es, V );

      if ( printWF and tau > relaxTime ) {
         count<maxN, dp>( nb, xmin, xmax, counts, flags, points );
      }

      tau += timeStep;
   }

   writeEnergy(Es);

   if ( printWF ) {
      writeWavefunction<dp>( xmin, xmax, nb, counts );
   }

   if ( printDist ) {
      writeDistribution<maxN, dp>( pNum, flags, points );
   }
}

template< int maxN, int dp >
void runSimulation( const int pNum, const int N0, const int nb, const double xmin, const double xmax,
                    const double timeStep, const double relaxTime, const double alpha,
                    std::array<double, dp> x0, std::function<double (std::array<double, dp> &)> &V,
                    const int numThreads = 1, const bool printWF = false, const bool printDist = false ) {
   /*Simultaneously runs numThreads diffusion Monte Carlo simulations.*/

   std::vector<std::thread> threads;

   for ( int i = 0; i < numThreads; ++i ) {
      auto f = std::bind( diffusionMC<maxN, dp>, pNum, N0, nb, xmin, xmax, timeStep, relaxTime, alpha,
                          x0, V, printWF, printDist );
      threads.push_back( std::thread(f) );
   }

   for ( auto &t : threads ) {
      t.join();
   }

   //std::vector<double> avgEs;
   //int numEs = results.back().first.size();

   //for ( int i = 0; i < numEs; ++i ) {
   //   std::vector<double> temp;

   //   for ( auto &p : results ) {
   //      temp.push_back( p.first.at(i) );
   //   }
   //   avgEs.push_back( std::accumulate(temp.begin(), temp.end(), 0.0) / temp.size() );
   //}

   //std::cout << avgEs.size() << "  numEs: " << numEs << std::endl;

   //writeEnergy(avgEs);
}

#endif

