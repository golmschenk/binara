#include <nanobind/nanobind.h>

#include "mcmc_wrapper.h"

NB_MODULE(binara_ext, m)
{
    m.def("run_mcmc", &Run_MCMC);
}
