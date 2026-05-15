Installation
============

You can either install mirar on your machine (using the steps below)
or you can alternatively follow the separate instructions for docker :doc:`docker`.


Cloning the repository
----------------------

Firstly, you need to clone the repository. You can do this with:

.. code-block:: bash

    git clone git@github.com:winter-telescope/mirar.git

Installing the package
----------------------

You need to install the package itself. The code is built using python.
It also uses the astromatic software suite (e.g. source-extractor, scamp, swarp, psfex) and astrometry.net, which are not python packages. You will need to install these separately (see below).


Option 1: Installing via install script
---------------------------------------

We provide an install script that will install mirar and all the dependencies for you.
This is the recommended way to install mirar, as it is the easiest and most likely to work.
It should run on linux, and on mac (new M1 chips only).

You will need to have conda installed on your machine to use this method. You can download it from the `official website <https://www.anaconda.com/products/distribution>`_.

Then, once you have cloned the repository, you can run the install script:

.. code-block:: bash

    cd mirar
    bash install_mirar.sh

This script will run for several minutes, with a lot of output in the terminal.
It will create a conda environment called `mirar` and install all the necessary dependencies for you. It will also install mirar itself.
You can specify an alternative name for the conda environment by instead running:

.. code-block:: bash

    bash install_mirar.sh mycustomenv

Either way, you should now have mirar installed.
You can activate the conda environment with:

.. code-block:: bash

    conda activate mirar

Scroll down to the end of this page for Final Steps.

Option 2: Creating the environment and installing manually
---------------------------------------------------------

You also have the option to create the conda environment and install mirar manually.
This is recommended if you want more control over the installation process, or if you run into issues with the install script.

Step 2a: Creating a conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In most cases, you can do the following:

.. code-block:: bash

    conda create -n mirar python=3.11
    conda activate mirar
    pip install --upgrade pip

However, if you are using a new Mac with an M1 chip, you might run into trouble.
Instead, we suggest:

.. code-block:: bash

    conda create -n mirar
    conda activate mirar
    conda config --env --set subdir osx-64
    conda install python=3.11
    pip install --upgrade pip

Step 2b: Installing via pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to install mirar is via pip:

.. code-block:: bash

    pip install -e ".[dev]"

This will install the latest version of mirar from `PyPI <https://pypi.org/project/mirar/>`_.

This method is recommended if you just want to use mirar without making changes.

Step 2c: Installing pre-commit hooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Step 2d: Installing non-python dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally you meed to install any optional dependencies that you might want to use.
We again recommend using conda. Whether you need these dependencies depends on your intended usage of mirar.

Dependencies include:

* `source-extractor <https://www.astromatic.net/software/sextractor/>`_ (a.k.a sextractor)
* `scamp <https://www.astromatic.net/software/scamp/>`_
* `swarp <https://www.astromatic.net/software/swarp/>`_
* `psfex <https://www.astromatic.net/software/psfex/>`_
* `astrometry.net <https://nova.astrometry.net/>`_

In general, you can install these packages in any way you like. We provide you with a few tips below, but if the packages are already available on your system, you should not need to install them again.

Software with conda (Linux, some Mac support)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can try installing everything via conda:

.. code-block:: bash

    conda install -c conda-forge astromatic-source-extractor astromatic-scamp astromatic-swarp astromatic-psfex astrometry gsl

This is likely to only work for linux at present, because some of the packages are not available for mac via conda.
If you are on mac, we recommend installing as many packages as possible with conda, any then any missing packages via the respective software websites.

Astrometry.net via Homebrew (Mac only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To run astrometry solutions with Astrometry.net (the default for SEDMv2), you'll need to download Astrometry.net
locally, as outlined `here <http://astrometry.net/use.html>`_. Once you have a local version, there should be an
astrometry-net folder somewhere on your machine. If you used Homebrew, it should be here:

.. code-block:: bash

    /opt/homebrew/Cellar/astrometry-net/

astromatic software with apt-get (Linux only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use apt-get if you are running Debian-based Linux:

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install -y sextractor scamp swarp psfex
    sudo ln -s /usr/bin/source-extractor /usr/bin/sex
    sudo ln -s /usr/bin/SWarp /usr/bin/swarp

The latter two lines are to ensure that source-extractor/swarp can be called from the command line in the way expected by mirar.

Final Steps
-----------------

Astrometry.net
^^^^^^^^^^^^^^^

For either option 1 or 2, you will have installed astrometry.net. This solves astrometry for images, but requires 'index files' to work. These are files that contain information about the positions of stars in the sky, and are used to match the stars in your images to known star positions. Index files need to be downloaded separately, and are not included in the astrometry.net installation.

You can grab index files from
`this directory <https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/>`_

Once you have downloaded the index files, you can specify the path to the astrometry.net folder and the index files via envirnoment variable:

.. code-block:: bash

    export ANET_INDEX_DIR=/path/to/astrometry-net

or specify this via the .env file in the root of the mirar repository.

PostgreSQL
^^^^^^^^^^

Some (but not all) pipelines use database integration. In general surveys (e.g. WINTER) need postgresQL, while follow-up instruments (such as LMI) do not.

If you want to use a pipeline that requires database integration, you will need to install PostgreSQL. You can do this via the `official website <https://www.postgresql.org/download/>`_.
These same database-integrated pipelines also typically require q3c, which is a PostgreSQL extension. You can install it via the `official website <https://github.com/segasai/q3c>`_.

Using mirar
^^^^^^^^^^^

Check out the `usage guide <usage>`_ to get started with using mirar, and the `contributing guide <contributing-guide>`_ if you want to contribute to the codebase.
