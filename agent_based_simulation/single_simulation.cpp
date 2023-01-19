#include "LearningAgent.hpp"
#include "utils.hpp"
#include <iostream>
#include <random>
#include <vector>
#include <bits/stdc++.h>
#include <omp.h>
#include <fstream>
using namespace std;

int main(int argc, char *argv[]){

    //initial parameters fitted on population data from France from 1950

    double a = 1.15106610e-02; 
    // double initial_b = 2.71975789e-02;
    double initial_b = 0.07; //initial aging rate of the populations
    double c = 4.26805906e-02;
    // double Lmax =  2.93238557e-02;
    double Lmax = 0.1; //initial Lmax of the learning population
    // double Lmax = 0.075;
    double k_learning = 3.94526937e+01;
    double n =  9.41817943e-02;
    double Gmax = 1.00399978e-01;
    double growth_rate = 9.23916941e-02;

    //female fertility parameters
    double c_f_female = 2.445e-5;
    double d_female = 14.8;
    double w_female = 32.836;
    double female_max_fertility = 0.12824068527053395;

    //male fertility parameters
    double c_f_male = c_f_female/2.5;
    double d_male = 14.8;
    double w_male = w_female+15;
    double male_max_fertility = 0.15924727639749914;

    //population size at time = 0
    int initial_population_size = 20000; 
    //max population size
    int const population_threshold = 100000;

    //number of simulation
    int simulation_number = 10;
    //number of years in one simulation
    double const simulation_time = 2000;
    double time_step = 1;

    //aging rate and learning mutation rate
    double mutation_rate = 0.02;

    //initialization of the learning male and female populations
    std::vector<LearningAgent> malePopulation;
    std::vector<LearningAgent> femalePopulation;

    double average_b_male, average_b_female, average_b_total, average_Lmax;
    double sum_b_male, sum_b_female, sum_Lmax_male, sum_Lmax_female;

    //initilization of the averages
    average_b_female = 0;
    average_b_male = 0; 
    average_b_total = 0;

    sum_b_female = 0;
    sum_b_male = 0;

    average_Lmax = 0;

    sum_Lmax_female = 0;
    sum_Lmax_male = 0;

    //parallelization parameters
    int nthreads = 6;
    bool parallel_on = true;
    omp_set_num_threads(nthreads);
    
    double time;
    double start;
    double end;
    //clock start
    start = omp_get_wtime();

    //initialization of the deadMales and deadFemales vectors
    std::vector<int> deadMales;
    std::vector<int> deadFemales;

    int malePopulationSize;
    int femalePopulationSize;
    int maleBabyCounter;
    int femaleBabyCounter;


    std::vector<double> b_vector_learning;
    std::vector<double> Lmax_vector_learning;
    //simulation loop
    for(int s = 0; s < simulation_number; s++){

        cout << s << endl;
        //clear all the vectors from the previous simulation
        b_vector_learning.clear();
        Lmax_vector_learning.clear();

        malePopulation.clear();
        femalePopulation.clear();

        //creation of the initial male agents (both learning and non learning)
        for(int i = 0; i < initial_population_size/2; i++){
            std::random_device rd;
            std::mt19937 gen(rd());

            //the agent's aging rate is distributed following a normal distribution around the initial aging rate
            std::normal_distribution<double> b_distribution(initial_b, 0.001);

            //random age
            std::normal_distribution<double> age_distribution(20, 10);
            double age = age_distribution(gen);
            while (age < 0){
                age = age_distribution(gen);
            }

            double b = b_distribution(gen);

            //the male agent is added to the male population vector
            malePopulation.push_back(LearningAgent(false, round(age), b, Lmax));
        }

        //creation of the initial female agents (both learning and non learning)
        for(int i = 0; i < initial_population_size/2; i++){
            std::random_device rd;
            std::mt19937 gen(rd());

            //the agent's aging rate is distributed following a normal distribution around the initial aging rate
            std::normal_distribution<double> b_distribution(initial_b, 0.001);

            //random age
            std::normal_distribution<double> age_distribution(20, 10);
            double age = age_distribution(gen);
            while (age < 0){
                age = age_distribution(gen);
            }
            double b = b_distribution(gen);

            //the female agent is added to the female population vector
            femalePopulation.push_back(LearningAgent(true, round(age), b, Lmax));
        }

        cout << "Simulation with learning ..." << endl;

        //simulation for the learning population
        time = 0;
        while (time < simulation_time){
            //baby counting counters reset
            maleBabyCounter = 0;
            femaleBabyCounter = 0;

            //compute the female and male population sizes
            malePopulationSize = malePopulation.size();
            femalePopulationSize = femalePopulation.size();

            //clearing of the dead agents vectors
            deadMales.clear();
            deadFemales.clear();

            //loop over the male agents to handle deaths
            #pragma omp parallel for if (parallel_on)
            for(int i = 0; i < malePopulationSize; i++){
                //if the agent can't reproduce anymore, the agent is as good as dead
                //we remove it to increase our effective population size
                if(malePopulation[i].getAge() > d_male + w_male){
                    deadMales.push_back(i);
                }
                //if the agent's age is below this threshold, check if it's going to die during the year
                else if(malePopulation[i].death(a, c, k_learning, n, Gmax)){
                    //if yes, add it to the dead male agents vector
                    deadMales.push_back(i);
                }
                //if the agent does not die, increment its age
                else{
                    malePopulation[i].incrementAge(time_step);
                }
            }

            //sorting of the dead males indexes to facilitate their removal
            std::sort(deadMales.begin(), deadMales.end(), greater<int>());

            //removal of all the dead male agents
            for (int i = 0; i < deadMales.size(); i++){
                malePopulation[deadMales[i]] = malePopulation.back();
                malePopulation.pop_back();
            }

            //loop over the female population to handle deaths
            #pragma omp parallel for if (parallel_on)
            for(int i = 0; i < femalePopulationSize; i++){
                //if the agent can't reproduce anymore, the agent is as good as dead
                //we remove it to increase our effective population size
                if(femalePopulation[i].getAge() > d_female + w_female){
                    deadFemales.push_back(i);
                }
                //if the agent's age is below this threshold, check if it's going to die during the year
                else if(femalePopulation[i].death(a, c, k_learning, n, Gmax)){
                    //if yes, add it to the dead female agents vector
                    deadFemales.push_back(i);
                }
                //if the agent does not die, increment its age
                else{
                    femalePopulation[i].incrementAge(time_step);
                }
            }

            //sorting of the dead females indexes to facilitate their removal
            std::sort(deadFemales.begin(), deadFemales.end(), greater<int>());
            //removal of all the dead female agents
            for (int i = 0; i < deadFemales.size(); i++){
                femalePopulation[deadFemales[i]] = femalePopulation.back();
                femalePopulation.pop_back();
            }

            //create random couples
            std::vector<std::pair<LearningAgent, LearningAgent>> couples;
            if (femalePopulation.size() > malePopulation.size()){
                for(int i = 0; i < malePopulation.size(); i++){
                    couples.push_back(std::make_pair(malePopulation[i], femalePopulation[i]));
                }
            }
            else {
                for(int i = 0; i < femalePopulation.size(); i++){
                    couples.push_back(std::make_pair(malePopulation[i], femalePopulation[i]));
                }
            }

            //loop over all couples
            for (auto couple : couples){
                //if the current population is bellow the maximum population threshold
                if(femalePopulation.size() + malePopulation.size() < population_threshold){
                    //if the male and the female can reproduce this year
                    if (couple.first.fertilityCheck(c_f_male, d_male, w_male, male_max_fertility) and couple.second.fertilityCheck(c_f_female, d_female, w_female, female_max_fertility)){
                        
                        //create a baby of a random sex
                        bool sex = randomSex();
                        std::random_device rd;
                        std::mt19937 gen(rd());

                        //the baby's aging rate is the mean of its parents' aging rate.

                        double baby_b = (couple.first.getB() + couple.second.getB())/2;
                        double baby_Lmax = (couple.first.getLmax() + couple.second.getLmax())/2;

                        //if the baby is a female, add it to the female population
                        if (sex){

                            //create a new LearningAgent object
                            LearningAgent baby = LearningAgent(sex, 0.0, baby_b, baby_Lmax);

                            //give it a chance to mutate 
                            baby.mutate(mutation_rate);
                            // uncomment this line to allow Lmax to mutate and evolve
                            // baby.mutate_learning(mutation_rate);
                            femalePopulation.push_back(baby);
                            femaleBabyCounter++;
                        }
                        //if the baby is a male, add it to the male population
                        else {

                            //create a new LearningAgent object
                            LearningAgent baby = LearningAgent(sex, 0.0, baby_b, baby_Lmax);

                            //give it a chance to mutate
                            baby.mutate(mutation_rate);
                            // uncomment this line to allow Lmax to mutate and evolve
                            // baby.mutate_learning(mutation_rate);
                            malePopulation.push_back(baby);
                            maleBabyCounter++;
                        }
                    }
                }
            }


            //compute the average aging rate and maximum benefit of learning in both populations
            average_b_female = 0;
            average_b_male = 0; 
            average_b_total = 0;

            sum_b_female = 0;
            sum_b_male = 0;

            average_Lmax = 0;

            sum_Lmax_female = 0;
            sum_Lmax_male = 0;

            for (auto male : malePopulation){
                sum_b_male += male.getB();
                sum_Lmax_male += male.getLmax();
            }

            for (auto female : femalePopulation){
                sum_b_female += female.getB();
                sum_Lmax_female += female.getLmax();
            }

            average_b_male = sum_b_male/malePopulation.size();
            average_b_female = sum_b_female/femalePopulation.size();
            average_b_total = (average_b_male*malePopulation.size() + average_b_female*femalePopulation.size())/(malePopulation.size() + femalePopulation.size());

            average_Lmax = (sum_Lmax_male + sum_Lmax_female)/(malePopulation.size() + femalePopulation.size());
            // cout << "new babies : " << maleBabyCounter + femaleBabyCounter << endl;

            //add those averages to the vectors
            b_vector_learning.push_back(average_b_total);
            Lmax_vector_learning.push_back(average_Lmax);
            //increment time
            time += time_step;
        }

        end = omp_get_wtime();
        double time_taken = end - start;
        std::cout << "TIME TAKEN BY THE PROGRAM : " << time_taken << endl;

        //save the vectors in text files for latter analysis
        std::ofstream outFileB("/simulation_results.txt", std::ios::app);
        for (const auto &e : b_vector_learning) outFileB << e << " ";
        outFileB << "\n";
        outFileB.close();

        // std::ofstream outFileLmax("./simulation_with_learning_evolution_Lmax.txt", std::ios::app);
        // for (const auto &e : Lmax_vector_learning) outFileLmax << e << " ";
        // outFileLmax << "\n";
        // outFileLmax.close();
    }
    return 0;
}