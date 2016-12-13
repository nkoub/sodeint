/***********************************************************************
 *  File        : ornsteinuhlenbeck_class.h
 *  Author      : Nikos E Kouvaris
 *  Date        : 07/12/2016
 *  Update      : 12/12/2016
 *  Copyright   : GPL
 *  Description : Ornstein-Uhlenbeck model 
 ***********************************************************************/

#ifndef _ORNSTEIN_UHLENBECK_H_
#define _ORNSTEIN_UHLENBECK_H_

#include <iostream>
#include <vector>
#include <array>

// Number of EQuations
// In the class ornsteinuhlenbeck_params I pass its value to oup.N
// Throughout the code I use it as oup.N 
const static size_t NEQ = 1;

// state_type 
// It stores the system's dynamical variables in 1D array as [x0,x1,...,xN]
typedef std::array<double,NEQ> state_type;                              

// System parameters
class ornsteinuhlenbeck_params {
public:    
    ornsteinuhlenbeck_params () {}     // constructor
    ~ornsteinuhlenbeck_params () {}    // destructor
    double tau;
    double noise_intensity;
    const size_t N = NEQ;
} oup;


// Deterministic part
class ornsteinuhlenbeck_det {
public:
    ornsteinuhlenbeck_det () {}   // constructor     
    ~ornsteinuhlenbeck_det () {}  // destructor
    
    // overload operator() for use with odeint steppers
    void operator () (const state_type& v0, state_type& dvdt) const {
        // odes
        dvdt[0] = -v0[0] / oup.tau;
    }
};//:~ornsteinuhlenbeck_det
   

// Stochastic part 
class ornsteinuhlenbeck_stoch {
public:
    ornsteinuhlenbeck_stoch() {}       // constructor that does nothing
    ~ornsteinuhlenbeck_stoch() {}      // destructor that does nothing
    void operator() (const state_type& v0, state_type& dvdt) const {
        dvdt[0] = sqrt(2.0*oup.noise_intensity);
    }

};//:~ornsteinuhlenbeck_stoch


#endif
