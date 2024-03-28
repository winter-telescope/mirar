Since `mirar` is an open-source software project, *anyone* can contribute it. You simply need to create a fork of the repository, commit your changes, and then make a pull request.

We have a few general guidelines that are helpful for keeping things organised:

* **Use Github Pull Requests** We like to make sure that the code stays working. So if you are developing something, create a fork of the repo and open a branch. Develop away, and when you are ready, open a pull request. We can then review the code and approve the PR to merge your changes into the main codebase.
* **Use Github Issues** to coordinate your development. Whether you found a bug, you want to request an enhancement, or you're actively developing a new feature, Github Issues is a great place to keep everyone informed about what you're working on. Click on the label button to provide more info about your topic. Every time you make a relevant PR, remember to tag the issue (e.g `git commit -m 'progress on #12'`), and when you finish and issue you can close it with a commit too! (e.g `git commit -m 'Close #12`').
* Keep Github Actions Happy! Github Actions runs all of our unit tests, to make sure you didn't break anything with your commit. You can see if the CI is happy by checking on the github page (look for a tick/cross next to your commit). If your commit failed, be sure to check the Github Actions website logs, to see exactly what went wrong. You can also check the badge:

   .. image:: https://github.com/winter-telescope/mirar/actions/workflows/continuous_integration.yml/badge.svg?branch=main
      :target: https://github.com/winter-telescope/mirar/actions/workflows/continuous_integration.yml?branch=main
      :align: center

* Our unit tests are set up to run on forks and thus feel free to enable Github Actions on your fork. These tests are designed to clearly fail if you are missing any required tokens or authentication.

* **Keep Github Actions Busy!** Github Actions will only run unit tests if we make the unit tests first. When you add a new feature, you also need to add some unit tests so that we can ensure this feature continues to work in the future. Your tests should be saved in the `tests/` directory, and you can find plenty of examples there to copy. Coveralls.io checks how much of the code is covered by tests, and helps you see which lines still need to be covered. If your commit adds a lot of new code but does not add unit tests, your commit will be tagged on github with a red cross to let you know that the code coverage is decreasing. If you want to know more about how to design unit tests, you can check out a guide `here <https://medium.com/swlh/introduction-to-unit-testing-in-python-using-unittest-framework-6faa06cc3ee1>`_. You can see all of this on the website: https://coveralls.io/repos/github/winter-telescope/mirar or the summary badge:

     .. image:: https://coveralls.io/repos/github/winter-telescope/mirar/badge.svg?branch=main
        :target: https://coveralls.io/github/winter-telescope/mirar?branch=main
        :align: center
* **Keep the code well-documented** When you write code, it is easier to understand 'what' than 'why'. People are not mind-readers, and this includes your future self. This is where documentation helps. If you add doctstrings following the `standard python style <https://peps.python.org/pep-0287/>`_, the code can be automatically converted to documentation.
* **Keep the code style clean** Everyone has their own coding style. But for collaborative projects, we need to ensure consistency. We use various automated tools to check and reformat code, and all new code committed should inform to these standards. See :ref:`Pre-commit Hooks`.

Updating the documentation
--------------------------

The documentation (generated primarily from docstrings) can be modified with the following command, **executed from the docs directory**:

.. code-block:: bash

    sphinx-apidoc -o source/ ../mirar --module-first --force


You can then build the documentation locally:

.. code-block:: bash

    sphinx-build source build


Checking the tests locally
--------------------------

You can run the tests with:

.. code-block:: bash

    TESTDATA_CHECK="True" python -m unittest discover tests/


This will check that the correct test data version is available, and then run all the tests.

You can also check the code contained within the docstrings/documentation:

.. code-block:: bash

    poetry run make -C docs/ doctest

Pre-commit Hooks
----------------

We use `pre-commit hooks <https://pre-commit.com/>`_ to ensure our coding style is consistent.

You can see the checks we follow in `.pre-commit-config.yaml`:

.. literalinclude:: ../../.pre-commit-config.yaml

Most of these checks will run automatically, and fix themselves automatically.
So you if run `git commit -m ...` and see something failed, you can just do `git add x` and add the now-corrected file, then commit again.

The only exception to this rule is the use of `pylint <https://pylint.pycqa.org/>`_.
Pylint ensures adherence to the PEP8 and other coding standards.
It outputs a list of all mistakes it identified, with lines and clear messages describing the problem.
Because it is so strict, we have the following policy:

**All new code should conform to the pylint standards for check categories of `FATAL` (F), `ERROR` (E), `WARNING` (W) and `REFACTOR` (R), when committed.**

**All new code should conform to the full pylint standards (including the additional `CONVENTION` (C) which covers documentation) at the time a PR is merged**

In the meantime, we are trying to upgrade older code to meet these standards.
