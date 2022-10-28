# winterdrp

[![CI](https://github.com/winter-telescope/winterdrp/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/winter-telescope/winterdrp/actions/workflows/continous_integration.yml)
[![Coverage Status](https://coveralls.io/repos/github/winter-telescope/winterdrp/badge.svg?branch=main)](https://coveralls.io/github/winter-telescope/winterdrp?branch=main)

Open-source modular python package for astronomy image reduction.

in addition to python requirements, winterdrp requires:

* Source-extractor/sextractor
* Swarp
* Scamp

Instructions for .env
1. Copy the `.env.example` file to the root of the project and update the environment variables.
2. Name this file `.env`

Install Instructions

winterdrp uses poetry for dependency management. You can install the code in the following way:

```
git clone git@github.com:winter-telescope/winterdrp.git
cd winterdrp
poetry install
```

This will install the combination of python dependencies which we have tested. You can then use winterdrp from the command line:

```python -m winterdrp ......```
