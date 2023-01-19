#include "utils.hpp"
#include <math.h>
#include <cmath>
#include <random>
using namespace std;

// Calculates the value of the Gompertz-Makeham mortality function for a given age
// age: The age of the individual
// a: Parameter of the Gompertz-Makeham function
// b: Parameter of the Gompertz-Makeham function
// c: Parameter of the Gompertz-Makeham function
double agingGompertzMakaham(float age, double a, double b, double c){
    return c + a*exp(b*age);
}

// Calculates the effect of learning on mortality
// age: The age of the individual
// Lmax: The maximum effect of learning on mortality
// k_learning: The age at which the effect of learning on mortality is half of Lmax
// n: The steepness of the sigmoid function that describes the effect of learning on mortality
double learning(float age, double Lmax, double k_learning, double n){
    return (Lmax/(1 + exp(n*(age-k_learning)))) - Lmax;
}

// Calculates the effect of growth on mortality
// age: The age of the individual
// Gmax: The maximum effect of growth on mortality
// growth_rate: The steepness of the function that describes the effect of growth on mortality
double growth(float age, double Gmax, double growth_rate){
    if (age < 30){
        return (Gmax/(1 + pow(age,growth_rate)) - Gmax);
    }
    else {
        return (Gmax/(1 + pow(30, growth_rate)) - Gmax);
    }
}

// Calculates the hazard rate for agents with learning
// age: The age of the individual
// a: Parameter of the Gompertz-Makeham function
// b: Parameter of the Gompertz-Makeham function
// c: Parameter of the Gompertz-Makeham function
// Lmax: The maximum effect of learning on mortality
// k_learning: The age at which the effect of learning on mortality is half of Lmax
// n: The steepness of the sigmoid function that describes the effect of learning on mortality
// Gmax: The maximum effect of growth on mortality
// growth_rate: The steepness of the function that describes the effect of growth on mortality
double mortalityWithLearning(float age, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate){
    double res = agingGompertzMakaham(age, a, b, c) + learning(age, Lmax, k_learning, n) + growth(age, Gmax, growth_rate);
    if (res > 1e-05 and res < 1){
        return res;
    }
    else if (res > 1){
        return 1.0;
    }
    else{
        return 1e-05;
    }
}

// Calculates the hazard rate for agents with learning
// age: The age of the individual
// a: Parameter of the Gompertz-Makeham function
// b: Parameter of the Gompertz-Makeham function
// c: Parameter of the Gompertz-Makeham function
// Lmax: The maximum effect of learning on mortality
// k_learning: The age at which the effect of learning on mortality is half of Lmax
// n: The steepness of the sigmoid function that describes the effect of learning on mortality
// Gmax: The maximum effect of growth on mortality
// growth_rate: The steepness of the function that describes the effect of growth on mortality
double mortalityWithoutLearning(float age, double a, double b, double c, double Gmax, double growth_rate){
    double res = agingGompertzMakaham(age, a, b, c) + growth(age, Gmax, growth_rate);
    if (res > 1e-05 and res < 1){
        return res;
    }
    else if (res > 1){
        return 1.0;
    }
    else{
        return 1e-05;
    }
}

double brassPolynomial(float age, double c, double d, double w){
    if (age > d and age < (d+w)){
        return c*(age-d)*(pow((d+w-age),2));
    }
    else {
        return 0.0;
    }
}

double normalizedBrassPolynomial(float age, double c, double d, double w, double max_fertility){
    return brassPolynomial(age, c, d, w)/max_fertility;
}

bool randomSex(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<int> distribution(1, 0.5);

    int sex = distribution(gen);
    return sex;
}