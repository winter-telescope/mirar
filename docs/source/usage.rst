Usage
=====

.. _installation:

Installation
------------

To use Lumache, first install it using pip:

.. code-block:: console

   (.venv) $ pip install lumache

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``winterdrp.pipelines.get_pipeline()`` function:

.. autofunction:: winterdrp.pipelines.get_pipeline

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`winterdrp.pipelines.get_pipeline`
will raise an exception.

For example:

>>> from winterdrp.pipelines import Pipeline
>>> print(Pipeline.pipelines.keys())
