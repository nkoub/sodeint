/*********************************************************************************
* File        : stochastic_heun.h
* Author      : Nikos E Kouavris
* Date        : 17/11/2016
* Update      : 21/11/2016
* Copyright   : GPL
* Description : Heun integrator for Stochastic ODEs with additive white noise.
*               (also called as Improved Euler, Euler-Heun, Heun, or Runge-Kutta 2nd order)
*********************************************************************************/
#ifndef _STOCHASTIC_HEUN_H_
#define _STOCHASTIC_HEUN_H_

#include <iostream>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

boost::random::mt19937 RNG_HEUN; // Random Number Generator for Heun



/************************************************************
* Euler-Heun scheme. Explicit order 0.5 strong Taylor scheme
************************************************************/
// NEQ: number of Equations
// e.g., NEQ=20 for a network of 10 nodes with dynamics in each node given by 2 odes 
template <size_t NEQ>
class StochasticHeun {
private:
	// boost::random::mt19937 RNG_HEUN; 
public:
	typedef std::array <double,NEQ> state_type;
	typedef boost::numeric::odeint::stepper_tag stepper_category;
	StochasticHeun (int seed) { RNG_HEUN.seed(seed); } 				// seed for RNG_HEUN   	
	~StochasticHeun () {}   	

	template <class S>
	void do_step(S system, state_type& v0, double t, const double dt) const {
		// define Gaussian noise generator  
    	boost::random::normal_distribution<double> gaussian_noise(0.0,1.0);	 
    	// and generate white noise for every node i
    	state_type noise;
		for (size_t i=0; i<NEQ; ++i) {
			noise[i] = gaussian_noise(RNG_HEUN);
		}
	
		// pass values from v0 to det and stoch
		state_type det, stoch;
		system.first(v0,det);
		system.second(v0,stoch);				

		// do one integration step (only with the stochastic part)		
		state_type v0_hat;
		for (size_t i=0; i<NEQ; ++i) {
			v0_hat[i] += stoch[i] * noise[i] * sqrt(dt);
		}

		// pass values from y_hat to det_hat and stoch_hat
		state_type det_hat, stoch_hat;
		system.first(v0_hat,det_hat);
		system.second(v0_hat,stoch_hat);

		// do one integration step
		for (size_t i=0; i<NEQ; ++i) {
			v0[i] += det[i] * dt + 0.5 * (stoch[i] + stoch_hat[i]) * noise[i] * sqrt(dt);
		}
	}
};//:~StochasticHeun
#endif
