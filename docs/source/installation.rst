Installation
============

Downloading
-----------

To use winterdrp, first you need to clone it:

.. code-block:: bash

    git clone git@github.com:winter-telescope/winterdrp.git


Installing the package
----------------------

You next need to install the package itself. The code is built using python.
We suggest creating a dedicated `conda <https://www.anaconda.com/products/distribution>`_ environment (), but you could also use a virtual environment.

.. code-block:: bash

    conda create -n winter_env python=3.10
    conda activate winter_env
    pip install --upgrade pip

Next you need to actually install winterdrp. We use `poetry <https://python-poetry.org/>`_ to manage dependencies:

.. code-block:: bash

    pip install poetry
    poetry install winterdrp

Lastly, you need to install the `pre-commit hooks <https://pre-commit.com/>`_ (see :doc:`contributing-guide` for more info about what these do):

.. code-block:: bash

    pre-commit install

Now you should have installed winterdrp. You can check it worked by opening up python and trying to import it:

.. code-block:: bash

    python

.. doctest::

    >>> from winterdrp.paths import package_name
    >>> print(f"This is the {package_name} package")
    This is the winterdrp package

Non-python dependencies
-----------------------

Finally you meed to install any optional dependencies that you might want to use.
We again recommend using conda. Whether you need these dependencies depends on your intended usage of winterdrp.

Dependencies include:

* `source-extractor <https://www.astromatic.net/software/sextractor/>`_ (a.k.a sextractor)
* `scamp <https://www.astromatic.net/software/scamp/>`_
* `swarp <https://www.astromatic.net/software/swarp/>`_
* `psfex <https://www.astromatic.net/software/psfex/>`_
* `postgreSQL <https://www.postgresql.org/download/>`_

PostgreSQL is relatively straightforward to install via the `official website <https://www.postgresql.org/download/>`_.
The other packages might be more complicated, and will depend on your platform.

astromatic software with apt-get (Linux only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use apt-get if you are running :

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install -y sextractor scamp swarp psfex
    sudo ln -s /usr/bin/source-extractor /usr/bin/sex
    sudo ln -s /usr/bin/SWarp /usr/bin/swarp

The latter two lines are to ensure that source-extractor/swarp can be called from the command line in the way expected by winterdrp.

astromatic software with conda (Linux, Mac or Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can try installing things via conda:

.. code-block:: bash

    conda install -c conda-forge astromatic-source-extractor astromatic-scamp astromatic-swarp astromatic-psfex
