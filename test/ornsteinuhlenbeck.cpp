/***********************************************************************
* File        : ornsteinuhlenbeck.cpp
* Author      : Nikos E Kouavris
* Date        : 07/12/2016
* Update      : 12/12/2016
* Copyright   : GPL
* Description : Ornstein-Uhlenbeck process
*
* build : g++ -std=c++11 -Ofast ornsteinuhlenbeck.cpp -o test.run
* run   : ./test.run  
***********************************************************************/

#include <iostream>
#include <cassert>
#include <vector>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>
#include <ctime>
#include "/home/nikos/nkcpplibs/sodeint/stochastic_euler.h"
#include "ornsteinuhlenbeck_class.h"


using namespace std;
using namespace boost::numeric::odeint;


int main (int argc, char* argv[]) {  
    cout << "\n";    
    size_t seed = 35074;            // seed for random numbers
    
    // initialize state array
    state_type v0;        
    boost::random::mt19937 rng;     // random number generator
    // rng.seed(seed);                
    rng.seed(static_cast<unsigned int>(std::time(0)));     
    boost::random::uniform_real_distribution<double> uni_rand(0,1);

    for (size_t i=0; i<oup.N; ++i) {
        v0[i] = uni_rand(rng);
    }
    
    cout << "Solve " << v0.size() << " equations in a system of " << oup.N << " nodes...\n";
    
    // System's parameters
    oup.noise_intensity = 0.01;
    oup.tau = 1.0;

    // Construct an object as a pair of the classes *_det and *_stoch
    auto ornsteinuhlenbeck = make_pair(ornsteinuhlenbeck_det(),ornsteinuhlenbeck_stoch());

    // Define the integrator stepper 
    StochasticEuler<v0.size()> integrator(seed);


    const double dt = 0.1;     
    const double tmax = 1.0;          // real time
    size_t int_steps = tmax / dt;      // integration steps
    double t = 0.0;

    // Integration loop
    for (size_t it=0; it<(int_steps); ++it) {
        t += dt;
        integrator.do_step(ornsteinuhlenbeck,v0,t,dt);  
        cout << v0[0] << "\n";
    }

    return 0;
}
