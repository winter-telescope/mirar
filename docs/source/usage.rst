Usage
=====

You can execute winterdrp via the command line:

.. code-block::
    python -m winterdrp -p name-of-pipeline -n night-to-reduce

One example is the following:

.. code-block::
    python -m winterdrp -p summer

How can you know which pipelines are available? You can check the documentation here: :doc:`winterdrp.pipelines`.
Alternatively, you can check via the code:

.. doctest::

    >>> from winterdrp.pipelines import Pipeline
    >>> print([x for x in Pipeline.pipelines.keys()])
    ['wirc', 'summer']
