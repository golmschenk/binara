#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <iostream>

#include "mcmc_wrapper.h"
#include "likelihood.h"


namespace nb = nanobind;

nb::ndarray<double, nb::numpy, nb::ndim<1>> internal_calculate_light_curve(
    const nb::ndarray<double>& times,
    const nb::ndarray<double>& parameters
)
{
    const size_t times_size = times.size();
    auto* model_fluxes = new double[times_size];
    // Create the Python ownership object.
    nb::capsule owner(model_fluxes, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    std::cout << "Process: " << getprogname() << std::endl;
    Calculate_Lightcurve(times.data(), times_size, parameters.data(), model_fluxes);
    const size_t model_fluxes_shape[1] = {times_size};
    return {model_fluxes, 1, model_fluxes_shape, owner};
}

NB_MODULE(binara_ext, module)
{
    module.def("internal_run_mcmc", &Run_MCMC);
    module.def("internal_calculate_light_curve", &internal_calculate_light_curve,
               nb::arg("times"), nb::arg("parameters"));
    module.def("internal_calculate_log_likelihood", &Log_Likelihood);
}
