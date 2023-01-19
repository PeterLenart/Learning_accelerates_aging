#ifndef NONLEARNINGAGENT_HPP
#define NONLEARNINGAGENT_HPP

class NonLearningAgent
{

private:
    bool sex;
    float age;
    double b;

public:
    NonLearningAgent(bool sex, float age, double b);

    void incrementAge(double time_step);
    void mutate(double mutation_rate);
    bool death(double a, double c, double Gmax);
    bool fertilityCheck(double c, double d, double w, double max_fertility);

    float getAge(){
        return this->age;
    }

    double getB(){
        return this->b;
    }
};

#endif