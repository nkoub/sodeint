/*********************************************************************************
* File        : stochastic_euler.h
* Author      : Nikos E Kouavris
* Date        : 17/11/2016
* Update      : 21/11/2016
* Copyright   : GPL
* Description : Euler and Euler-Predictor-Corrector integrators 
*				for Stochastic ODEs with additive white noise.
*********************************************************************************/
#ifndef _STOCHASTIC_EULER_H_
#define _STOCHASTIC_EULER_H_

#include <iostream>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

boost::random::mt19937 RNG_EULER; // Random Number Generator for Euler


/*******************************************************
* Euler scheme. Explicit order 0.5 strong Taylor scheme
*******************************************************/
// NEQ: number of Equations
// e.g., NEQ=20 for a network of 10 nodes with dynamics in each node given by 2 odes 
template <size_t NEQ>
class StochasticEuler {
private:
	// boost::random::mt19937 RNG_EULER;
public:
	typedef std::array <double,NEQ> state_type;
	typedef boost::numeric::odeint::stepper_tag stepper_category;
	StochasticEuler (int seed) { RNG_EULER.seed(seed); } 			// seed for RNG_EULER   	
	~StochasticEuler () {}
    
	template <class S>
	void do_step(S system, state_type& v0, double t, const double dt) const {
		// define Gaussian noise generator  
    	boost::random::normal_distribution<double> gaussian_noise(0.0,1.0);	 
    	// and generate white noise for every node i
    	state_type noise;
		for (size_t i=0; i<NEQ; ++i) {
			noise[i] = gaussian_noise(RNG_EULER);
		}
	
		// pass values of v0 to det and stoch
		state_type det, stoch;
		system.first(v0,det);
		system.second(v0,stoch);		
  
		// do one integration step
		for (size_t i=0; i<NEQ; ++i) {
			v0[i] += det[i] * dt + stoch[i] * noise[i] * sqrt(dt);
		}
	}
}; //:~StochasticEuler


/*************************************************************************
* Euler Predictor-Corretor scheme. Explicit order 0.5 strong Taylor scheme
*************************************************************************/
// NEQ: number of Equations
// e.g., NEQ=20 for a network of 10 nodes with dynamics in each node given by 2 odes 
template <size_t NEQ>
class StochasticEulerPC {
private:
	// boost::random::mt19937 RNG_EULER;
public:
	typedef std::array <double,NEQ> state_type;
	typedef boost::numeric::odeint::stepper_tag stepper_category;
	StochasticEulerPC (int seed) { RNG_EULER.seed(seed); } 	// seed for RNG_EULER   	
	~StochasticEulerPC () {}   	

	template <class S>
	void do_step(S system, state_type& v0, double t, const double dt) const {
		// define Gaussian noise generator  
    	boost::random::normal_distribution<double> gaussian_noise(0.0,1.0);	 
    	// and generate white noise for every node i
    	state_type noise;
		for (size_t i=0; i<NEQ; ++i) {
			noise[i] = gaussian_noise(RNG_EULER);
		}

		// pass values of v0 to det and stoch
		state_type det, stoch;
		system.first(v0,det);
		system.second(v0,stoch);

		// Predictor: do one integration step
		state_type v0_hat;
		for (size_t i=0; i<NEQ; ++i) {
			v0_hat[i] = v0[i] + det[i] * dt + stoch[i] * noise[i] * sqrt(dt);
		}
		
		// pass values of v0_hat to det_hat and stoch_hat
		state_type det_hat;
		// state_type stoch_hat;
		system.first(v0_hat,det_hat);
		// system.second(v0_hat,stoch_hat);
		
		// Corrector: do one integration step
		for (size_t i=0; i<NEQ; ++i) {
			v0[i] += 0.5 * (det_hat[i] + det[i]) * dt + stoch[i] * noise[i] * sqrt(dt);
		}
	}
}; //:~StochasticEulerPC


#endif





/*
//Euler scheme, explicit order 0.5 strong 
template <size_t N>
class StochasticEuler1 {
public:
	typedef std::array <double,N> state_type;
	typedef boost::numeric::odeint::stepper_tag stepper_category;
	  
	template <class S>
	void do_step(S system, state_type& y0, double t, const double dt) const {
		state_type det, stoch, noise;
		
		// pass values of y0 to det,stoch,noise
		std::get<0>(system)(y0,det);
		std::get<1>(system)(y0,stoch);
		std::get<2>(system)(y0,noise);
		
		for (size_t i=0; i<N; ++i) {
			y0[i] += det[i] * dt + stoch[i] * noise[i] * sqrt(dt) ;
		}
	}
};
//:~StochasticEuler
*/
