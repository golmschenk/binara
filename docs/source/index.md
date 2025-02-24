% binara documentation master file, created by
% sphinx-quickstart on Mon Jan  1 22:48:08 2024.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# Welcome to binara's documentation!

## Installation

To install `binara` use

```shell
pip install binara
```

Although not required, as is generally good practice for any development project, we highly recommend creating a separate virtual environment for each distinct project. For example, via Conda, creating a virtual environment for a project using `binara` might look like

```
conda create -n binara_env python=3.11
```

Then before working, be sure to activate your environment with

```shell
conda activate binara_env
```

Then install `binara` within this environment.

```{toctree}
:caption: 'Contents:'
:maxdepth: 2

user_guides/index
developer_guides/index
```

# Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
