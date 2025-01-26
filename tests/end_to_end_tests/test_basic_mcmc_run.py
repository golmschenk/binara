from pathlib import Path

import binara

from tests.working_directory_context import use_working_directory


def test_can_run_mcmc():
    with use_working_directory(Path(__file__).parent):
        binara.enforce_expected_data_directory_tree()
        binara.run_mcmc(220052771, 6, 1, 0, 0, 0)
