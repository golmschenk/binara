import typing
from pathlib import Path

import pandas as pd


def load_states_data_frame_from_file(file_handle: Path | typing.TextIO) -> pd.DataFrame:
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
