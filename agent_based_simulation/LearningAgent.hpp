#ifndef LEARNINGAGENT_HPP
#define LEARNINGAGENT_HPP

class LearningAgent
{

private:
    bool sex;
    float age;
    double b;
    double Lmax;

public:
    LearningAgent(bool sex, float age, double b, double Lmax);
    void incrementAge(double time_step);
    void mutate(double mutation_rate);
    void mutate_learning(double mutation_rate);
    bool death(double a, double c, double k_learning, double n, double Gmax);
    bool fertilityCheck(double c, double d, double w, double max_fertility);

    float getAge(){
        return this->age;
    }

    double getB(){
        return this->b;
    }

    double getLmax(){
        return this->Lmax;
    }
};

#endif