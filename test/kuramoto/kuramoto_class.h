/***********************************************************************
 *  File        : kuramoto_class.h
 *  Author      : Nikos E Kouvaris
 *  Date        : 02/12/2016
 *  Update      : 15/12/2016
 *  Copyright   : GPL
 *  Description : Kuramoto model for general coupling architecture 
 ***********************************************************************/

#ifndef _KURAMOTO_H_
#define _KURAMOTO_H_

#include <iostream>
#include <vector>
#include <array>
#include <cmath> // sin

// Network size is a Global Static Variable 
// In the class kuramoto_params I pass its value to kurp.N
// Throughout the code I use it as kurp.N 
const static size_t SYSTEM_SIZE_USED_ONLY_HERE = 100;

// state_type 
// It stores the system's dynamical variables in 1D array as [θ0,θ1,...,θN]
typedef std::array<double,SYSTEM_SIZE_USED_ONLY_HERE> state_type;                              


class kuramoto_params {
public:    
    kuramoto_params () {}     // constructor
    ~kuramoto_params () {}    // destructor
    double coupling;
    double noise_intensity;
    state_type omega;
    std::vector<std::vector<double>> adj_mat;
    const size_t N = SYSTEM_SIZE_USED_ONLY_HERE;
} kurp;


// Deterministic part
class kuramoto_det {
public:
    kuramoto_det () {}   // constructor     
    ~kuramoto_det () {}  // destructor
    
    // overload operator() for use with odeint steppers
    void operator () (const state_type& v0, state_type& dvdt) const {
        for (size_t i=0; i<kurp.N; ++i) {                       
            // coupling terms: Modify here for implementing the coupling scheme required for each problem
            double sum_phi = 0.0;
            for(int j=0; j<kurp.N; ++j) {
                sum_phi += kurp.adj_mat[i][j] * sin(v0[j] - v0[i]);
            } 
            // coupled odes
            dvdt[i] = kurp.omega[i] + kurp.coupling * sum_phi / (1.0*kurp.N);
        }
    }
};//:~kuramoto_det
   

// Stochastic part 
class kuramoto_stoch {
public:
    kuramoto_stoch() {}       // constructor that does nothing
    ~kuramoto_stoch() {}      // destructor that does nothing
    void operator() (const state_type& v0, state_type& dvdt) const {
        for (size_t i=0; i<kurp.N; ++i) {    
            dvdt[i] = sqrt(2.0*kurp.noise_intensity);
        }
    }

};//:~kuramoto_stoch


#endif
