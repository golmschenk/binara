from io import StringIO
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd
import xarray


def load_states_data_frame_from_states_file(file_handle: Path | StringIO) -> pd.DataFrame:
    """
    Loads the states output text file as a Pandas data frame.

    :param file_handle: The handle to the output file.
    :return: The data frame.
    """
    data_frame = pd.read_csv(file_handle, index_col=False, header=None, sep=r'\s+')
    column_names = ['iteration', 'log_likelihood'] + [
        f'parameter{index}' for index in range(len(data_frame.columns) - 2)]
    data_frame.columns = column_names
    return data_frame


def convert_from_states_pandas_data_frame_to_xarray_dataset(data_frame: pd.DataFrame) -> xarray.Dataset:
    iterations = data_frame['iteration'].to_numpy(dtype=np.int32)
    log_likelihoods = data_frame['log_likelihood'].to_numpy(dtype=np.float32)
    parameters = data_frame.drop(['iteration', 'log_likelihood'], axis='columns').to_numpy(dtype=np.float32)
    dataset = xarray.Dataset(
        data_vars=dict(
            parameter=(["iteration", "parameter_index"], parameters),
            log_likelihood=(["iteration"], log_likelihoods),
        ),
        coords=dict(
            iterations=iterations,
            parameter_index=np.arange(parameters.shape[1], dtype=np.int32),
        ),
    )
    return dataset


def load_log_likelihoods_from_states_file(file_handle: Path | StringIO) -> npt.NDArray:
    """
    Loads the log likelihood from the output text file.

    :param file_handle: The handle to the output file.
    :return: The log likelihood array.
    """
    data_frame = load_states_data_frame_from_states_file(file_handle)
    return data_frame['log_likelihood'].to_numpy()


def load_iterations_from_states_file(file_handle: Path | StringIO) -> npt.NDArray:
    """
    Loads the log likelihood from the output text file.

    :param file_handle: The handle to the output file.
    :return: The log likelihood array.
    """
    data_frame = load_states_data_frame_from_states_file(file_handle)
    return data_frame['iteration'].to_numpy()
