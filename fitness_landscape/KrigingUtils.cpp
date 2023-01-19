#include "KrigingUtils.hpp"
#include "./libIntegrate/libIntegrate/Integrate.hpp"
#include <math.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <random>
using namespace std;

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

double agingGompertzMakaham(double age, double a, double b, double c){
    return c + a*exp(b*age);
}

double learning(double age, double Lmax, double k_learning, double n){
    return (Lmax/(1 + exp(n*(age-k_learning)))) - Lmax;
}

double growth(double age, double Gmax, double growth_rate){
    if (age < 30){
        return (Gmax/(1 + pow(age,growth_rate)) - Gmax);
    }
    else {
        return (Gmax/(1 + pow(30, growth_rate)) - Gmax);
    }
}

double mortalityWithLearning(double age, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate){
    double res = agingGompertzMakaham(age, a, b, c) + learning(age, Lmax, k_learning, n) + growth(age, Gmax, growth_rate);
    if (res > 0 and res < 1){
        return res;
    }
    else if (res > 1){
        return 1.0;
    }
    else{
        return 0;
    }
}

double brassPolynomial(double age, double c, double d, double w){
    if (age > d and age < (d+w)){
        return c*(age-d)*(pow((d+w-age),2));
    }
    else {
        return 0.0;
    }
}

double l(double age, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate){
    _1D::SimpsonRule<double> integrate;
    double I = integrate([a, b, c, Lmax, k_learning, n, Gmax, growth_rate](double x){return mortalityWithLearning(x, a, b, c, Lmax, k_learning, n, Gmax, growth_rate);}, 0, age, 1000);
    return exp(-I);
}

double eulerLotka(double fin, double r, double c_f, double d, double w, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate){
    _1D::SimpsonRule<double> integrate;
    double I = integrate([r, c_f, d, w, a, b, c, Lmax, k_learning, n, Gmax, growth_rate](double x){return l(x, a, b, c, Lmax, k_learning, n, Gmax, growth_rate) * brassPolynomial(x, c_f, d, w) * exp(-r*x);}, 0, fin, 1000);
    return I - 1;
}

double brents_fun(std::function<double (double)> f, double lower, double upper, double tol, unsigned int max_iter)
{
	double a = lower;
	double b = upper;
	double fa = f(a);	// calculated now to save function calls
	double fb = f(b);	// calculated now to save function calls
	double fs = 0;		// initialize 
 
	if (!(fa * fb < 0))
	{
		std::cout << "Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
		return -11;
	}
 
	if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
	{
		std::swap(a,b);
		std::swap(fa,fb);
	}
 
	double c = a;			// c now equals the largest magnitude of the lower and upper bounds
	double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
	bool mflag = true;		// boolean flag used to evaluate if statement later on
	double s = 0;			// Our Root that will be returned
	double d = 0;			// Only used if mflag is unset (mflag == false)
 
	for (unsigned int iter = 1; iter < max_iter; ++iter)
	{
		// stop if converged on root or error is less than tolerance
		if (std::abs(b-a) < tol)
		{
			std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
			return s;
		} // end if
 
		if (fa != fc && fb != fc)
		{
			// use inverse quadratic interopolation
			s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
				+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
				+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );
		}
		else
		{
			// secant method
			s = b - fb * (b - a) / (fb - fa);
		}
 
			// checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
		if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
				( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
				( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
				( mflag && (std::abs(b-c) < tol) ) ||
				( !mflag && (std::abs(c-d) < tol))	)
		{
			// bisection method
			s = (a+b)*0.5;
 
			mflag = true;
		}
		else
		{
			mflag = false;
		}
 
		fs = f(s);	// calculate fs
		d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
		c = b;		// set c equal to upper bound
		fc = fb;	// set f(c) = f(b)
 
		if ( fa * fs < 0)	// fa and fs have opposite signs
		{
			b = s;
			fb = fs;	// set f(b) = f(s)
		}
		else
		{
			a = s;
			fa = fs;	// set f(a) = f(s)
		}
 
		if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
		{
			std::swap(a,b);		// swap a and b
			std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
		}
 
	} // end for
    
	std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
}

double fitness(double fin, double c_f, double d, double w, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate){
    auto f = [fin, c_f, d, w, a, b, c, Lmax, k_learning, n, Gmax, growth_rate](double r){return eulerLotka(fin, r, c_f, d, w, a, b, c, Lmax, k_learning, n, Gmax, growth_rate);};
    double fitness = brents_fun(f, -1, 1, 0.00000001, 100);
    if (fitness < 0){
        fitness = 0;
    }
    return fitness;
}

double fertilityPeter(double age, double b_0){
	return (b_0*(1 + (age*0.05)));
	// return (b_0/(1 + (age*0.05)));
	// return b_0;
}

double eulerLotkaPeter(double fin, double r, double b_0, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate){
    _1D::SimpsonRule<double> integrate;
    double I = integrate([r, b_0, a, b, c, Lmax, k_learning, n, Gmax, growth_rate](double x){return l(x, a, b, c, Lmax, k_learning, n, Gmax, growth_rate) * fertilityPeter(x, b_0) * exp(-r*x);}, 0, fin, 2000);
    return I - 1;
}

double fitnessPeter(double fin, double b_0, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate){
	auto f = [fin, b_0, a, b, c, Lmax, k_learning, n, Gmax, growth_rate](double r){return eulerLotkaPeter(fin, r, b_0, a, b, c, Lmax, k_learning, n, Gmax, growth_rate);};
    double fitness = brents_fun(f, -1, 1, 0.000001, 100);
    if (fitness < 0){
        fitness = 0;
    }
    return fitness;
}