Usage
=====

1 - Configuring variables
---------------------

mirar has several options for execution.
You can pass some of these as arguments, but some e.g tokens or passwords are best included as environment variables.

You can find a full list of variables in `env.example`:

.. literalinclude:: ../../env.example


Option 1a: Using a .env file
^^^^^^^^^^^^^^^^^^^^^^^^

We recommend that you create an .env file in the root of your repository, and include all the necessary variables there.

.. code-block:: bash

    cp env.example .env

Then edit the .env file to include the correct values for your machine.
mirar will try to automatically load the .env file whenever you use the code.


Option 1b: Setting variables in the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also set individual variables in the command line:

.. code-block:: bash

    export RAW_DATA_DIR=/home/astronomer/rawdata

Option 1c: Setting variables in the command line with a .env file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your .env file is located elsewhere, you can load all of these variables at once using the command line:

.. code-block:: bash

    set -o allexport
    source .env
    set +o allexport


2 - Configuring data
---------------------

mirar is designed to work with data organized in a specific way.

Let's assume you have a single directory you want to put data in, called RAW_DATA_DIR=/home/astronomer/rawdata.
You should set the RAW_DATA_DIR variable to this path (see above).

Then, within this directory, you have a subdirectory for each instrument, aa sub-subdirectory for each night of data, and then a sub-sub-subdirectory called 'raw' containing the raw exposures.

Let's say you have a single night of data from LMI, taken on 2023-07-01. You would organize your data as follows:

.. code-block:: text

    /home/astronomer/rawdata/
        lmi/
            20230701/
                raw/
                    lmi_20230701_0001.fits
                    lmi_20230701_0002.fits
                    ...

You might have multiple nights or multiple instruments.
An example of such a structure is shown below:

.. code-block:: text

    /home/astronomer/rawdata/
        lmi/
            20230701/
                raw/
                    lmi_20230701_0001.fits
                    lmi_20230701_0002.fits
                    ...
            20230702/
                raw/
                    lmi_20230702_0001.fits
                    lmi_20230702_0002.fits
                    ...
        gmos/
            20260701/
                raw/
                    gmos_20230701_0001.fits
                    gmos_20230701_0002.fits
                    ...
            20260702/
                raw/
                    gmos_20230702_0001.fits
                    gmos_20230702_0002.fits
                    ...

It doesn't matter what the exact names of the files are,
as long as they are in the correct directory structure.

3 - Deciding which pipeline and configuration to run
----------------------------------------------------

We have one mirar 'pipeline' per instrument (e.g. GMOS, LMI, WINTER, SPRING etc).
The pipeline can run in different configurations, for example one might fully reduce the data while another might just create a csv log of the data.

How can you know which instruments are available? You can check the documentation here: :doc:`autogen/pipelines`.
Alternatively, you can check via the code:

.. doctest::

    >>> from mirar.pipelines import Pipeline
        >>> print(sorted([x for x in Pipeline.pipelines.keys()]))
        ['git', 'gmos', 'lmi', 'sedmv2', 'summer', 'wasp', 'winter', 'wirc']
    >>> print(sorted([x for x in Pipeline.pipelines.keys()]))
    ['git', 'gmos', 'lmi', 'sedmv2', 'summer', 'wasp', 'winter', 'wirc']

The configurations are also listed here :doc:`autogen/pipelines`.

4 - Running the code
--------------------

Once you've chosen your pipeline and config, you can execute mirar via the command line:

.. code-block:: bash

    mirar-run -p name-of-pipeline -n night-to-reduce -c configuration-to-use

One example is the following:

.. code-block:: bash

    mirar-run -p lmi -n 20230701 -c log

In general, starting with a 'log' configuration is a good idea,
as it will quickly create a csv file with all the information about the data.
You can check that this looks correct before running a more
computationally expensive configuration.

Pipelines generally have a 'default' configuration,
which you can choose by not specifying any configuration or by choosing `-c default`.
For data reduction, that's usually what you want. For example:

.. code-block:: bash

    mirar-run -p lmi -n 20230701
