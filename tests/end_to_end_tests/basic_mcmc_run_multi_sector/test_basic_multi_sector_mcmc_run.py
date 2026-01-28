import shutil
from pathlib import Path

import binara

from tests.data_directory_manipulations import verify_directories_match
from tests.working_directory_context import use_working_directory


def test_basic_multi_sector_mcmc_run():
    working_directory = Path(__file__).parent
    with use_working_directory(working_directory):
        data_directory = working_directory.joinpath('data')
        if data_directory.exists():
            shutil.rmtree(data_directory)
        shutil.copytree(working_directory.joinpath('template_starting_data_directory'), data_directory)
        binara.enforce_expected_data_directory_tree()
        binara.internal_run_mcmc(150284425, -1, 1, 1, 0, 0)
        verify_directories_match(data_directory, working_directory.joinpath('expected_resulting_data_directory'))
        if data_directory.exists():
            try:
                shutil.rmtree(data_directory)
            except PermissionError:  # In the Windows tests, this can sporadically fail, but it's fine to skip it.
                pass
