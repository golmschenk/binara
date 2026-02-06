import shutil
from pathlib import Path

import pytest

import binara

from tests.data_directory_manipulations import verify_directories_match
from tests.working_directory_context import use_working_directory


def test_can_run_mcmc_multi_sector():
    working_directory = Path(__file__).parent
    with use_working_directory(working_directory):
        data_directory = working_directory.joinpath('data')
        if data_directory.exists():
            shutil.rmtree(data_directory)
        sessions_directory = working_directory.joinpath('sessions')
        if sessions_directory.exists():
            shutil.rmtree(sessions_directory)
        shutil.copytree(working_directory.joinpath('template_starting_data_directory'), data_directory)
        binara.internal_run_mcmc(150284425, -1, 1, 1, 0, 0)
        verify_directories_match(sessions_directory,
                                 working_directory.joinpath('expected_resulting_sessions_directory'))
        if data_directory.exists():
            try:
                shutil.rmtree(data_directory)
            except PermissionError:  # In the Windows tests, this can sporadically fail, but it's fine to skip it.
                pass
        if sessions_directory.exists():
            try:
                shutil.rmtree(sessions_directory)
            except PermissionError:  # In the Windows tests, this can sporadically fail, but it's fine to skip it.
                pass
