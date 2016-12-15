/***********************************************************************
* File        : fhn_ring.cpp
* Author      : Nikos Kouavris
* Date        : 02/12/2016
* Update      : 15/12/2016
* Copyright   : GPL
* Description : Kuramoto model for general coupling architecture 
*
* build : g++ -std=c++11 -Ofast kuramoto_ring.cpp -o ringKuramoto.run
* run   : ./ringKuramoto.run  
***********************************************************************/

#include <iostream>
#include <cassert>
#include <vector>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

#include "/home/nikos/nkcpplibs/sodeint/stochastic_heun.h"
#include "kuramoto_class.h"


using namespace std;
using namespace boost::numeric::odeint;


// Create a ring of N nodes.
// Each node is joined to the 2 nearest neighbors: 1 to the left and 1 to the right.
typedef std::vector<std::vector<double>> adjacency_matrix;
adjacency_matrix ring_network(const size_t N) {
    adjacency_matrix adj(N, std::vector<double>(N));
    for (int i=0; i<N; ++i) {
        for (int j=i-1; j<i+1; ++j) {
            if (j<0)    j = N-1;
            if (j>N-1)  j = 0;
            if (i!=j)   adj[i][j] = adj[j][i] = 1.0; 
            else        adj[i][j] = adj[j][i] = 0.0;             
        }
    }
    return adj;
}  




int main (int argc, char* argv[]) {  
    cout << "\n";    
    size_t seed = 35074;            // seed for random numbers
    
    // initialize state array
    state_type v0;        
    boost::random::mt19937 rng;     // random number generator
    rng.seed(seed);                     
    boost::random::normal_distribution<double> normal_rand(0.0,1.0);
    boost::random::uniform_real_distribution<double> uni_rand(0,1);

    for (size_t i=0; i<kurp.N; ++i) v0[i] = uni_rand(rng);
  
    cout << "Solve " << v0.size() << " equations in a system of " << kurp.N << " nodes...\n";
    
    // System's parameters
    kurp.noise_intensity = 0.01; 
    kurp.coupling = 0.08;
    kurp.omega.fill(0.5);
    kurp.adj_mat = ring_network(kurp.N);

    // folder and files to work with
    string folder = "./";
    string prefix = str(boost::format("ring_N%d_D%.4f_coupling%.5f")
                                 %kurp.N%kurp.noise_intensity%kurp.coupling);
    string subfix = ".out";   
    ofstream fout_states(folder + prefix + "_states" + subfix);
    ofstream fout_times(folder + prefix + "_times" + subfix);

    cout << "The files will have the prefix: " << prefix << "\n\n";
    cout << "start the integration...\n";

    // Construct an object as a pair of the classes kuramoto_det and kuramoto_stoch
    auto kuramoto = make_pair(kuramoto_det(),kuramoto_stoch());

    // Define the integrator stepper 
    StochasticHeun<v0.size()> integrator(seed);


    vector<double> times;
    vector<state_type> states;

    const double dt = 0.1;     
    const double tmax = 150.0;      // real time
    size_t total_time = tmax / dt;  // integration steps
    size_t shift_time = 0.0;//total_time/2.0;
    double stream_step = 1.0;
    double t = 0.0;

    // Integration loop
    for (size_t it=0; it<(total_time); ++it) {
        t += dt;
        integrator.do_step(kuramoto,v0,t,dt);  
        
        if ((it>shift_time) && (fmod(it,(stream_step))==0.0)) {
            times.push_back(t);
            states.push_back(v0);
        }
    }

    cout << "real integration time: " << t << ",  integration steps: " << total_time << "\n";

    // save times
    for (double tt : times) {
        fout_times << tt << '\n';
    }

    // save states spacetime
    for (state_type t_state : states) {
        for (double i_state : t_state) {
            fout_states << '\t' << fmod(i_state,6.283185307179586);
        }
         fout_states << '\n';
    }

    return 0;
}
