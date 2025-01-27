#include "random_generator.h"

#include <random>

class RandomGenerator
{
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

public:
    explicit RandomGenerator(const unsigned int seed = 0)
        : gen(seed), dis(0.0, 1.0)
    {
    }

    double get_random()
    {
        return dis(gen);
    }
};

extern "C" {
RandomGenerator* create_random_generator(const unsigned int seed)
{
    return new RandomGenerator(seed);
}

double get_random_value(RandomGenerator* generator)
{
    return generator->get_random();
}

void destroy_random_generator(RandomGenerator* generator)
{
    delete generator;
}
}
