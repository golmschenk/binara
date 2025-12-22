#include <nanobind/nanobind.h>

#include "mcmc_wrapper.h"

NB_MODULE(binara_ext, module)
{
    module.def("internal_run_mcmc", &Run_MCMC);
    module.def("internal_calculate_light_curve", &Calculate_Lightcurve);
    module.def("internal_calculate_log_likelihood", &Log_Likelihood);
}
