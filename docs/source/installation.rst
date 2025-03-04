Installation
============


Installing the package
----------------------

You need to install the package itself. The code is built using python.
We suggest creating a dedicated `conda <https://www.anaconda.com/products/distribution>`_ environment, but you could also use a virtual environment.

Prerequisite: Creating a conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In most cases, you can do the following:

.. code-block:: bash

    conda create -n mirar python=3.11
    conda activate mirar
    pip install --upgrade pip

However, if you are using a new Mac with an arm chip, you might run into trouble.
Instead, we suggest:

.. code-block:: bash

    conda create -n mirar
    conda activate mirar
    conda config --env --set subdir osx-64
    conda install python=3.11
    pip install --upgrade pip

Option 1: Installing via pip
----------------------------

The easiest way to install mirar is via pip:

.. code-block:: bash

    pip install mirar

This will install the latest version of mirar from `PyPI <https://pypi.org/project/mirar/>`_.

This method is recommended if you just want to use mirar without making changes.

Option 2: Installing via git
----------------------------

Downloading
^^^^^^^^^^^

You can alternatively grab the latest version of the code from github:

.. code-block:: bash

    git clone git@github.com:winter-telescope/mirar.git
    cd mirar

This method is recommended if you want to contribute to the code.

Installing python dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next you need to actually install mirar. We use `poetry <https://python-poetry.org/>`_ to manage dependencies.
Firstly ensure you have poetry installed, in an isolated environment:

.. code-block:: bash

    conda deactivate
    conda create -n pipx python=3.11
    conda activate pipx
    pip install pipx
    pipx install poetry
    pipx ensurepath


Then exit the terminal and open a new one. You should now be able to run poetry commands.
You can now install the dependencies. Navigate to the root of the mirar directory,  and run:

.. code-block:: bash
    conda activate mirar
    poetry install

Lastly, you need to install the `pre-commit hooks <https://pre-commit.com/>`_ (see :doc:`contributing-guide` for more info about what these do):

.. code-block:: bash

    pre-commit install

Now you should have installed mirar. You can check it worked by opening up python and trying to import it:

.. code-block:: bash

    python

.. doctest::

    >>> from mirar.paths import PACKAGE_NAME
    >>> print(f"This is the {PACKAGE_NAME} package")
    This is the mirar package

Non-python dependencies
-----------------------

Finally you meed to install any optional dependencies that you might want to use.
We again recommend using conda. Whether you need these dependencies depends on your intended usage of mirar.

Dependencies include:

* `source-extractor <https://www.astromatic.net/software/sextractor/>`_ (a.k.a sextractor)
* `scamp <https://www.astromatic.net/software/scamp/>`_
* `swarp <https://www.astromatic.net/software/swarp/>`_
* `psfex <https://www.astromatic.net/software/psfex/>`_
* `astrometry.net <https://nova.astrometry.net/>`_
* `postgreSQL <https://www.postgresql.org/download/>`_

PostgreSQL is relatively straightforward to install via the `official website <https://www.postgresql.org/download/>`_.
The other packages might be more complicated, and will depend on your platform.

In general, you can install these packages in any way you like. We provide you with a few tips below, but if the packages are already available on your system, you should not need to install them again.

Astrometry.net
^^^^^^^^^^^^^^

To run astrometry solutions with Astrometry.net (the default for SEDMv2), you'll need to download Astrometry.net
locally, as outlined `here <http://astrometry.net/use.html>`_. Once you have a local version, there should be an
astrometry-net folder somewhere on your machine. If you used Homebrew, it should be here:

.. code-block:: bash

    /opt/homebrew/Cellar/astrometry-net/

Then, make sure to also grab index files from
`this directory <https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/>`_

Once you have downloaded the index files, you can specify the path to the astrometry.net folder and the index files via envirnoment variable:

.. code-block:: bash

    export ANET_INDEX_DIR=/path/to/astrometry-net

or specify this via the .env file in the root of the repository.

PostgreSQL
^^^^^^^^^^

Database management is done through PostgreSQL. You can install it via the `official website <https://www.postgresql.org/download/>`_.

Some pipelines require a database to store the results. If you want to use this functionality, you will need to install PostgreSQL.
These pipelines also typically require q3c, which is a PostgreSQL extension. You can install it via the `official website <https://github.com/segasai/q3c>`_.



astromatic software with apt-get (Linux only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use apt-get if you are running Debian-based Linux:

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install -y sextractor scamp swarp psfex
    sudo ln -s /usr/bin/source-extractor /usr/bin/sex
    sudo ln -s /usr/bin/SWarp /usr/bin/swarp

The latter two lines are to ensure that source-extractor/swarp can be called from the command line in the way expected by mirar.

astromatic software with conda (Linux, Mac or Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can try installing things via conda:

.. code-block:: bash

    conda install -c conda-forge astromatic-source-extractor astromatic-scamp astromatic-swarp astromatic-psfex astrometry gsl
