#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

class RandomGenerator;

RandomGenerator* create_random_generator(unsigned int seed);

double get_uniform_random_value(RandomGenerator* generator);

double get_normal_random_value(RandomGenerator* generator);

void destroy_random_generator(RandomGenerator* generator);

#endif // RANDOM_GENERATOR_H
