Usage
=====

.. _installation:

Installation
------------

To use winterdrp, first install it using pip:

.. code-block:: console

   (.venv) $ pip install lumache

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``winterdrp.pipelines.get_pipeline()`` function:

.. autofunction:: winterdrp.pipelines.get_pipeline

Otherwise, :py:func:`winterdrp.pipelines.get_pipeline`
will raise an exception.

For example:

.. testsetup::
   :skipif: pd is None

   data = pd.Series([42])

.. doctest::
   :skipif: pd is None

   >>> data.iloc[0]

.. doctest::

    >>>from winterdrp.pipelines import Pipeline
    >>>print(Pipeline.pipelines.keys())
