#ifndef UTILS_HPP
#define UTILS_HPP

double agingGompertzMakaham(float age, double a, double b, double c);
double learning(float age, double Lmax, double k_learning, double n);
double growth(float age, double Gmax, double growth_rate = 0.1);
double mortalityWithLearning(float age, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate = 0.1);
double mortalityWithoutLearning(float age, double a, double b, double c, double Gmax, double growth_rate = 0.1);
double brassPolynomial(float age, double c, double d, double w);
double maleFertilityFunction(float age, double c, double d, double w, float age_half);
double normalizedBrassPolynomial(float age, double c, double d, double w, double max_fertility);
double normalizedMaleFertilityFunction(float age, double c, double d, double w, float age_half, double max_fertility);
bool randomSex();

#endif