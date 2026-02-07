import os
import shutil
from pathlib import Path

from binara.init_data import write_mcmc_data

from tests.data_directory_manipulations import verify_directories_match
from tests.working_directory_context import use_working_directory

# Python and the extension library each seem to have their own copy of OpenMP. This allows that during tests.
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

def test_basic_init_data():
    working_directory = Path(__file__).parent
    with use_working_directory(working_directory):
        input_data_directory = working_directory.joinpath('input_data')
        if input_data_directory.exists():
            shutil.rmtree(input_data_directory)
        shutil.copytree(working_directory.joinpath('template_starting_input_data_directory'), input_data_directory)
        write_mcmc_data(tic_id=220052771, sector=6)
        verify_directories_match(input_data_directory, working_directory.joinpath('expected_resulting_input_data_directory'))
        if input_data_directory.exists():
            try:
                shutil.rmtree(input_data_directory)
            except PermissionError:  # In the Windows tests, this can sporadically fail, but it's fine to skip it.
                pass
