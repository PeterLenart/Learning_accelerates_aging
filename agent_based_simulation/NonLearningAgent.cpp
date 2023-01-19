#include "NonLearningAgent.hpp"
#include "utils.hpp"
#include <iostream>
#include <random>
using namespace std;

//constructor
NonLearningAgent::NonLearningAgent(bool sex, float age, double b){
    this->sex = sex;
    this->age = age;
    this->b = b;
}

//function used to increment the agent's age
void NonLearningAgent::incrementAge(double time_step){
    this->age = this->age + time_step;
}

//function used to give the newborn agent a chance to mutate
void NonLearningAgent::mutate(double mutation_rate){
    std::random_device rd;
    std::mt19937 gen(rd());

    //binomial distribution used to check if the agent mutates or not
    std::binomial_distribution<int> distribution(1, mutation_rate);

    //mutation normal distribution
    std::normal_distribution<double> n_distribution(this->b, 0.006);

    int mutation = distribution(gen);

    //if the binomial trial succeeds, mutate the agent's aging rate
    if (mutation == 1){
        this->b = n_distribution(gen);
    }
}

//function used to kill the agent if the binomial trial succeeds
bool NonLearningAgent::death(double a, double c, double Gmax){
    double mortality = mortalityWithoutLearning(this->age, a, this->b, c, Gmax);
    std::random_device rd;
    std::mt19937 gen(rd());

    //binomial trial where the chance of sucess is the chance the agent has of dying during the year
    std::binomial_distribution<int> distribution(1, mortality);
    int death = distribution(gen);
    if(death == 1){
        return true;
    }
    else{
        return false;
    }
}

//checks if an agent can reproduce
bool NonLearningAgent::fertilityCheck(double c, double d, double w, double max_fertility){
    //computes normalized fertility (chance of reproducing this year)
    double fertility = normalizedBrassPolynomial(this->age, c, d, w, max_fertility);
    std::random_device rd;
    std::mt19937 gen(rd());

    //binomial trial
    std::binomial_distribution<int> distribution(1, fertility);
    int fert = distribution(gen);
    if(fert == 1){
        return true;
    }
    else{
        return false;
    }
}