#include "LearningAgent.hpp"
#include "NonLearningAgent.hpp"
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

    //fitted on cohort data from the United Kingdom from 1950.

    // double a = 1.59253264e-02;
    // double initial_b = 1.39462044e-02;
    // double c = 1.55190395e-02;
    // double Lmax =  1.14012524e-02;
    // double k_learning = 3.33767043e+01;
    // double n =  1.07676768e-01;
    // double Gmax = 5.81282414e-02;
    // double growth_rate = 1.07539252e-01;

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
    int initial_population_size = 5000; 
    //max population size
    int const population_threshold = 10000;

    //number of simulation
    int simulation_number = 10;
    //number of years in one simulation
    double const simulation_time = 40000;
    double time_step = 1;

    //aging rate and learning mutation rate
    double mutation_rate = 0.02;

    //initialization of the learning male and female populations
    std::vector<LearningAgent> malePopulation;
    std::vector<LearningAgent> femalePopulation;

    double average_b_male, average_b_female, average_b_total, average_Lmax;
    double sum_b_male, sum_b_female, sum_Lmax_male, sum_Lmax_female;

    //initialization of the non learning male and female populations
    std::vector<NonLearningAgent> malePopulationWithoutLearning;
    std::vector<NonLearningAgent> femalePopulationWithoutLearning;

    //initilization of the averages
    average_b_female = 0;
    average_b_male = 0; 
    average_b_total = 0;

    sum_b_female = 0;
    sum_b_male = 0;

    //parallelization parameters
    int nthreads = 8;
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
    std::vector<double> b_vector_without_learning;

    //simulation loop
    for(int s = 0; s < simulation_number; s++){

        cout << s << endl;
        //clear all the vectors from the previous simulation
        b_vector_learning.clear();
        b_vector_without_learning.clear();

        malePopulation.clear();
        femalePopulation.clear();
        malePopulationWithoutLearning.clear();
        femalePopulationWithoutLearning.clear();

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
            //an identical agent, without learning, is added to the non learning male population vector
            malePopulationWithoutLearning.push_back(NonLearningAgent(false, round(age), b));
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
            //an identical agent, without learning, is added to the non learning female population vector
            femalePopulationWithoutLearning.push_back(NonLearningAgent(true, round(age), b));
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

                        //if the baby is a female, add it to the female population
                        if (sex){

                            //create a new LearningAgent object
                            LearningAgent baby = LearningAgent(sex, 0.0, baby_b, Lmax);

                            //give it a chance to mutate 
                            baby.mutate(mutation_rate);
                            femalePopulation.push_back(baby);
                            femaleBabyCounter++;
                        }
                        //if the baby is a male, add it to the male population
                        else {

                            //create a new LearningAgent object
                            LearningAgent baby = LearningAgent(sex, 0.0, baby_b, Lmax);

                            //give it a chance to mutate
                            baby.mutate(mutation_rate);
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

            for (auto male : malePopulation){
                sum_b_male += male.getB();
            }

            for (auto female : femalePopulation){
                sum_b_female += female.getB();
            }

            average_b_male = sum_b_male/malePopulation.size();
            average_b_female = sum_b_female/femalePopulation.size();
            average_b_total = (average_b_male*malePopulation.size() + average_b_female*femalePopulation.size())/(malePopulation.size() + femalePopulation.size());

            // cout << average_b_total << endl;
            // cout << time << endl;
            // cout << "male population : " << malePopulation.size() << endl;
            // cout << "female population : " << femalePopulation.size() << endl;
            cout << "new babies : " << maleBabyCounter + femaleBabyCounter << endl;
            // cout << "new female babies : " << femaleBabyCounter << endl;

            //add those averages to the vectors
            b_vector_learning.push_back(average_b_total);

            //increment time
            time += time_step;
        }

        cout << "Simulation without learning ..." << endl;
        //simulation for the non learning population
        time = 0;
        while (time < simulation_time){
            maleBabyCounter = 0;
            femaleBabyCounter = 0;
            malePopulationSize = malePopulationWithoutLearning.size();
            femalePopulationSize = femalePopulationWithoutLearning.size();
            deadMales.clear();
            deadFemales.clear();

            #pragma omp parallel for if (parallel_on)
            for(int i = 0; i < malePopulationSize; i++){
                if(malePopulationWithoutLearning[i].getAge() > d_male + w_male){
                    deadMales.push_back(i);
                }
                else if(malePopulationWithoutLearning[i].death(a, c, Gmax)){
                    deadMales.push_back(i);
                }
                else{
                    malePopulationWithoutLearning[i].incrementAge(time_step);
                }
            }

            std::sort(deadMales.begin(), deadMales.end(), greater<int>());
            for (int i = 0; i < deadMales.size(); i++){
                malePopulationWithoutLearning[deadMales[i]] = malePopulationWithoutLearning.back();
                malePopulationWithoutLearning.pop_back();
            }

            // cout << femalePopulation.size() << endl;

            #pragma omp parallel for if (parallel_on)
            for(int i = 0; i < femalePopulationSize; i++){
                if(femalePopulationWithoutLearning[i].getAge() > d_female + w_female){
                    deadFemales.push_back(i);
                }
                else if(femalePopulationWithoutLearning[i].death(a, c, Gmax)){
                    // malePopulation.erase(malePopulation.begin() + i);
                    deadFemales.push_back(i);
                }
                else{
                    femalePopulationWithoutLearning[i].incrementAge(time_step);
                }
            }

            std::sort(deadFemales.begin(), deadFemales.end(), greater<int>());
            for (int i = 0; i < deadFemales.size(); i++){
                femalePopulationWithoutLearning[deadFemales[i]] = femalePopulationWithoutLearning.back();
                femalePopulationWithoutLearning.pop_back();
            }

            std::vector<std::pair<NonLearningAgent, NonLearningAgent>> couples;
            if (femalePopulationWithoutLearning.size() > malePopulationWithoutLearning.size()){
                for(int i = 0; i < malePopulationWithoutLearning.size(); i++){
                    couples.push_back(std::make_pair(malePopulationWithoutLearning[i], femalePopulationWithoutLearning[i]));
                }
            }
            else {
                for(int i = 0; i < femalePopulationWithoutLearning.size(); i++){
                    couples.push_back(std::make_pair(malePopulationWithoutLearning[i], femalePopulationWithoutLearning[i]));
                }
            }

            // #pragma omp parallel for if (parallel_on)
            for (auto couple : couples){
                if(femalePopulationWithoutLearning.size() + malePopulationWithoutLearning.size() < population_threshold){
                    if (couple.first.fertilityCheck(c_f_male, d_male, w_male, male_max_fertility) and couple.second.fertilityCheck(c_f_female, d_female, w_female, female_max_fertility)){
                        bool sex = randomSex();
                        std::random_device rd;
                        std::mt19937 gen(rd());
                        double baby_b = (couple.first.getB() + couple.second.getB())/2;

                        if (sex){
                            NonLearningAgent baby = NonLearningAgent(sex, 0.0, baby_b);
                            baby.mutate(mutation_rate);
                            femalePopulationWithoutLearning.push_back(baby);
                            femaleBabyCounter++;
                        }
                        else {
                            NonLearningAgent baby = NonLearningAgent(sex, 0.0, baby_b);
                            baby.mutate(mutation_rate);
                            malePopulationWithoutLearning.push_back(baby);
                            maleBabyCounter++;
                        }
                    }
                }
            }

            average_b_female = 0;
            average_b_male = 0; 
            average_b_total = 0;

            sum_b_female = 0;
            sum_b_male = 0;

            for (auto male : malePopulationWithoutLearning){
                sum_b_male += male.getB();
            }

            for (auto female : femalePopulationWithoutLearning){
                sum_b_female += female.getB();
            }

            average_b_male = sum_b_male/malePopulationWithoutLearning.size();
            average_b_female = sum_b_female/femalePopulationWithoutLearning.size();
            average_b_total = (average_b_male*malePopulationWithoutLearning.size() + average_b_female*femalePopulationWithoutLearning.size())/(malePopulationWithoutLearning.size() + femalePopulationWithoutLearning.size());

            // cout << average_b_male << endl;
            // cout << time << endl;
            // cout << "male population : " << malePopulationWithoutLearning.size() << endl;
            // cout << "female population : " << femalePopulationWithoutLearning.size() << endl;
            cout << "new babies : " << maleBabyCounter + femaleBabyCounter << endl;
            // cout << "new female babies : " << femaleBabyCounter << endl;

            b_vector_without_learning.push_back(average_b_total);

            time += time_step;
        }

        end = omp_get_wtime();
        double time_taken = end - start;
        std::cout << "TIME TAKEN BY THE PROGRAM : " << time_taken << endl;

        //save the vectors in text files for latter analysis
        std::ofstream outFile("./base_results_for_paper.txt", std::ios::app);
        for (const auto &e : b_vector_learning) outFile << e << " ";
        outFile << "\n";
        for (const auto &e : b_vector_without_learning) outFile << e << " ";
        outFile << "\n";
        outFile.close();
    }
    return 0;
}