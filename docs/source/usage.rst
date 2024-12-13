Usage
=====

Configuring variables
---------------------

mirar has several options for execution.
You can pass some of these as arguments, but some e.g tokens or passwords are best included as environment variables.

You can find a full list of variables in `env.example`:

.. literalinclude:: ../../env.example

If you have created your own .env file, mirar will try to automatically load the .env file.

You can also set individual variables in the command line:

.. code-block:: bash

    export RAW_DATA_DIR=/home/astronomer/rawdata

If you have installed mirar via pip, or your .env file is located elsewhere, you can load all of these variables at once using the command line:

.. code-block:: bash

    set -o allexport
    source .env
    set +o allexport


Running the code
----------------

You can execute mirar via the command line:

.. code-block:: bash

    python -m mirar -p name-of-pipeline -n night-to-reduce

One example is the following:

.. code-block:: bash

    python -m mirar -p summer

How can you know which pipelines are available? You can check the documentation here: :doc:`mirar.pipelines`.
Alternatively, you can check via the code:

.. doctest::

    >>> from mirar.pipelines import Pipeline
        >>> print(sorted([x for x in Pipeline.pipelines.keys()]))
        ['git', 'gmos', 'sedmv2', 'summer', 'wasp', 'winter', 'wirc']
    >>> print(sorted([x for x in Pipeline.pipelines.keys()]))
    ['git', 'gmos', 'sedmv2', 'summer', 'wasp', 'winter', 'wirc']
