#ifndef KRIGINGUTILS_HPP
#define KRIGINGUTILS_HPP

double agingGompertzMakaham(double age, double a, double b, double c);
double learning(double age, double Lmax, double k_learning, double n);
double growth(double age, double Gmax, double growth_rate = 0.1);
double mortalityWithLearning(double age, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate = 0.1);
double l(double age, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate = 0.1);
double eulerLotka(double fin, double r, double c_f, double d, double w, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate);
double fitness(double fin, double c_f, double d, double w, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate);
double brassPolynomial(double age, double c, double d, double w);

double fertilityPeter(double age, double b_0);
double eulerLotkaPeter(double fin, double r, double b_0, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate);
double fitnessPeter(double fin, double b_0, double a, double b, double c, double Lmax, double k_learning, double n, double Gmax, double growth_rate);

#endif