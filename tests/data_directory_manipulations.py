import re
from pathlib import Path

import pytest


def verify_directories_match(data_directory, expected_resulting_run_directory):
    expected_path_list = list(expected_resulting_run_directory.glob('*'))
    run_path_list = list(data_directory.glob('*'))
    verify_file_lists_match(data_directory, run_path_list, expected_path_list)
    for expected_path in expected_path_list:
        run_path = data_directory.joinpath(expected_path.name)
        assert run_path.exists()
        if Path(expected_path).is_dir():
            verify_directories_match(run_path, expected_path)
        else:
            verify_run_files_match(run_path, expected_path)


def verify_file_lists_match(run_directory, run_path_list, expected_path_list):
    assert set([path.name for path in run_path_list]) == set([path.name for path in expected_path_list])


def verify_run_files_match(run_path, expected_run_path):
    with expected_run_path.open() as expected_file, run_path.open() as run_file:
        line_number = 1
        expected_list = re.split(r'(\s+)', expected_file.read())
        actual_list = re.split(r'(\s+)', run_file.read())
        for (expected_item, actual_item) in zip(expected_list, actual_list):
            expected_white_space_match = re.fullmatch(r'\s+', expected_item)
            if expected_white_space_match is not None:
                actual_white_space_match = re.fullmatch(r'\s+', actual_item)
                assert actual_white_space_match is not None, f'''
                    When comparing the expected {expected_run_path} and the actual {run_path}
                    on line {line_number}, expected a segment of white space and found {actual_item}.  
                '''
                if '\n' in expected_item:
                    line_number += 1
                continue
            try:
                expected_number = float(expected_item)
                actual_number = float(actual_item)
                relative_tolerance = 0.01
                assert actual_number == pytest.approx(expected_number, rel=relative_tolerance), f'''
                    When comparing the expected {expected_run_path} and the actual {run_path}
                    on line {line_number}, the number {expected_number} was expected and the actual was {actual_number}
                    which does not match to a relative tolerance of {relative_tolerance}.  
                '''
                continue
            except ValueError:
                assert actual_item == expected_item, f'''
                    When comparing the expected {expected_run_path} and the actual {run_path}
                    on line {line_number}, the string {expected_item} was expected and the actual was {actual_item}.  
                '''
                continue
