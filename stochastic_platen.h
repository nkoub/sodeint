/************************************************************
* File        : stochastic_platen.h
* Author      : Nikos E Kouvaris
* E-mail 	  : nkouba@gmail.com
* Date        : 17/11/2016
* Update      : 15/12/2016
* Description : Explicit Platen integrator for Stochastic ODEs 
*				to be used within boost::odeint
************************************************************/

#ifndef _STOCHASTIC_PLATEN_H_
#define _STOCHASTIC_PLATEN_H_

#include <iostream>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

boost::random::mt19937 RNG_PLATEN; // Random Number Generator for Euler


/********************************************************
* Platen scheme. Explicit order 1.0 strong Taylor scheme
* see Kloeden's book pp.374--375
********************************************************/
// NEQ: number of Equations
// e.g., NEQ=20 for a network of 10 nodes with dynamics in each node given by 2 odes 
template <size_t NEQ>
class StochasticPlaten {
private:
	// boost::random::mt19937 RNG_PLATEN; // Does not work. Why ???
public:
	typedef std::array <double,NEQ> state_type;
	typedef boost::numeric::odeint::stepper_tag stepper_category;
	StochasticPlaten (int seed) { RNG_PLATEN.seed(seed); } 	// seed for RNG_PLATEN   
	~StochasticPlaten () {}    
	template <class S>
	void do_step(S system, state_type& v0, double t, const double dt) const {
		// define Gaussian noise generator  
    	boost::random::normal_distribution<double> gaussian_noise(0.0,1.0);	 
    	// and generate white noise for every node i
    	state_type noise;
		for (size_t i=0; i<NEQ; ++i) {
			noise[i] = gaussian_noise(RNG_PLATEN);
		}		
		
		// pass values of v0 to det and stoch
		state_type det, stoch;
		system.first(v0,det);
		system.second(v0,stoch);

		// intermadiate step
		state_type v0_hat;
		for (size_t i=0; i<NEQ; ++i) {
			v0_hat[i] = v0[i] + det[i] * dt + stoch[i] * sqrt(dt);
		}
		
		// pass values of v0_hat to det_hat and stoch_hat
		state_type det_hat, stoch_hat;
		system.first(v0_hat,det_hat);
		system.second(v0_hat,stoch_hat);
		
		// final step
		for (size_t i=0; i<NEQ; ++i) {
			v0[i] += det[i] * dt 
			      + stoch[i] * noise[i] * sqrt(dt)
				  + (stoch_hat[i]-stoch[i]) * sqrt(dt) * (noise[i] * noise[i] - dt) / (2.0 * sqrt(dt));
		}
	}
};//:~StochasticPlaten
#endif






