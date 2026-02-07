import contextlib
import os
from pathlib import Path


@contextlib.contextmanager
def use_working_directory(working_directory: Path):
    original_working_directory = Path(os.getcwd())
    os.chdir(working_directory)
    try:
        yield
    finally:
        os.chdir(original_working_directory)
