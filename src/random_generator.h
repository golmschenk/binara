#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct RandomGenerator RandomGenerator;

RandomGenerator* create_random_generator(unsigned int seed);

double get_random_value(RandomGenerator* generator);

void destroy_random_generator(RandomGenerator* generator);

#ifdef __cplusplus
}
#endif

#endif // RANDOM_GENERATOR_H
