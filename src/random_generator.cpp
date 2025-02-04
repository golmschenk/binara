#include "random_generator.h"

#include <random>

class RandomGenerator
{
    std::mt19937 generator;
    std::uniform_real_distribution<> uniformDistribution;
    std::normal_distribution<> normalDistribution;

public:
    explicit RandomGenerator(const unsigned int seed = 0)
        : generator(seed), uniformDistribution(0.0, 1.0), normalDistribution(0.0, 1.0)
    {
    }

    double get_uniform_random()
    {
        return uniformDistribution(generator);
    }

    double get_normal_random()
    {
        return normalDistribution(generator);
    }
};

extern "C" {
RandomGenerator* create_random_generator(const unsigned int seed)
{
    return new RandomGenerator(seed);
}

double get_uniform_random_value(RandomGenerator* generator)
{
    return generator->get_uniform_random();
}

double get_normal_random_value(RandomGenerator* generator)
{
    return generator->get_normal_random();
}

void destroy_random_generator(RandomGenerator* generator)
{
    delete generator;
}
}
