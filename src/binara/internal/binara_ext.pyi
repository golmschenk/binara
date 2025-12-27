from typing import overload

import numpy
from numpy.typing import NDArray


def internal_run_mcmc(arg0: int, arg1: int, arg2: int, arg3: int, arg4: int, arg5: int, /) -> None: ...

@overload
def internal_calculate_light_curve(arg0: float, arg1: int, arg2: float, arg3: float, /) -> None: ...

@overload
def internal_calculate_light_curve(times: NDArray[numpy.float64], parameters: NDArray[numpy.float64], model_fluxes: NDArray[numpy.float64]) -> None: ...

def internal_calculate_log_likelihood(arg0: float, arg1: float, arg2: float, arg3: int, arg4: int, arg5: float, arg6: float, arg7: float, arg8: int, arg9: int, arg10: int, /) -> float: ...
