STALK: The file I/O
===================

This is an introduction to the basic principles of file I/O in STALK. The basic function of
the file I/O is to cache the results of expensive calculations to disk so that they can be
re-used later and inspected by the user. The file I/O also enables asynchronous modes of
operation, where the cost functions are evaluated separately. Furthermore, since the
algorithms have a random component, serialization of the results is a way to ensure
reproducibility.

No I/O
------

The default behavior of STALK is to not write any files to disk. Instead, the objects are
are only stored in memory. This is suitable for lightweight calculations with the perks of
minimal input configuration and no cluttering of the working directory with files.

For example, first part of the :ref:`orientation` has no file I/O. So, calculating a
line-search could look as follows:

.. code-block:: python

    p = ParameterSet([0.2])
    ls = LineSearch(p, d=0, offsets=np.linspace(-1.0, 1.0, 7), sigma=sigma)
    # No I/O
    ls.evaluate(pes, add_sigma=True)

The same principle applies to all other evaluation calls in STALK, including the evaluation
of individual parameter sets, calculation of Hessians, optimization of the parallel
line-searches, and so on. The following calls would produce no files on disk:

.. code-block:: python

    # p is a ParameterSet object
    p.evaluate(pes, path=None)
    # hessian is a ParameterHessian object
    hessian.compute_fdiff(pes, path=None)
    # surrogate is a Surrogate object
    surrogate.optimize(path=None, **optimize_kwargs)
    # lsi is a LineSearchIteration object
    lsi.propagate(pes, path=None)

where the path argument would be None by default.

With I/O
--------

In contrast, the file I/O is immediately enabled by passing a 'path' argument to evaluation
call. Let us see how 'path' works in different evaluation contexts, starting from the ground
up.

Parameter set and value caching
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest file I/O context is a single parameter set, a *structure*, whose file I/O is
triggered by the PES evaluation. Each instance is assigned a directory specified by a unique
label in a user-specified subpath. The generic file layout looks as follows:

.. code-block:: bash

    path/
    └── label/
        ├── params.in
        ├── sigma.in
        ├── error.out
        └── value.out

and it could be produced by the following code:

.. code-block:: python

    p1 = ParameterSet([1.0], label='label')
    pes(p1, path='path')

The files have the following meanings:

* params.in: the parameter values of the structure (array)
* sigma.in: (optional) the target uncertainty of the structure (float)
* error.out: (optional) the uncertainty of the structure upon evaluation (float)
* value.out: the energy/value of the structure upon evaluation (float)

Notes:
* The path could be './' or '', reducing the depth to one directory.
* If 'label' is not specified, a unique hash used instead based on the parameter values.
* In different modes of operation, the files may be created in different sequences. For
example, in the 'I/O' mode, only 'params.in' and 'sigma.in' are created, and
'value.out' is left to be supplied by the user after external evaluation.

To understand the caching of the results, consider the following examples:

.. code-block:: python

    # First, evaluate and write to disk
    pes(p1, path='path')
    # Subsequent evaluation calls will not trigger re-evaluation with the result in memory
    pes(p1, path=None)
    pes(p1, path='path')
    # Reset the value from memory
    p1.reset_value()
    # The result is recovered to memory from disk, no need to re-evaluate
    pes(p1, path='path')
    # The result will be re-evaluated, because the memory is flushed and a new path is assigned
    pes(p1, path='another_path', reset_value=True)

Caching of the results is crucial, because it avoids unnecessary re-evaluation of the PES,
which can be expensive. Moreover, if the results are stochastic, result serialization is a
way to ensure reproducibility. For example, tracing back and forth a sequence of
line-searches is only possible if the values are serialized.

Understanding the caching of results may be important for the correct troubleshooting of the
workflows.

Line-search
^^^^^^^^^^^

The next unit of the file I/O is a line-search, comprising a set of structures along a line
in the parameter space. It produces a file layout such as follows:

.. code-block:: bash

    ls/
    ├── d0_+0.2000/
    │   ├── params.in
    │   ├── sigma.in
    │   ├── error.out
    │   └── value.out
    ├── d0_+0.1000/
    ├── d0_-0.1000/
    ├── d0_-0.2000/
    ├── eqm/
    └── result/
        ├── fit.out
        ├── ls.out
        ├── x0.out
        └── y0.out

where the repeated .in/.out files for parameter evaluations outside 'd0\_+0.2000/' are
omitted. However, a 'result/' subdirectory is also featured containing the following files:

* fit.out: Fitting coefficients of the line-search (array)
* ls.out: The raw line-search data: offset/value/error (array)
* x0.out: The horizontal estimate of the minimum point and errorbar (float, float)
* y0.out: The vertical estimate of the minimum point and errorbar (float, float)

The 'results/' is shown here separately but the data can be written to any directory,
including the parent 'ls/', which is in fact the default. The above file layout can be
produced by the following code:

.. code-block:: python

    ls = LineSearch(structure=p, M=5, R=0.2, sigma=0.2, d=0)
    pes(ls, path='ls')
    # After evaluation, serialize the result to disk.
    ls.save_result('ls/result/')
    # Afterwards, the line-search could be reloaded with:
    ls_load = LineSearch(structure=p, M=5, R=0.2, sigma=0.2, d=0, path='ls/result/')

So, the line-search is centered around the supplied structure 'p' and then expanded along
the direction d=0 on a regular grid of M=5 points with a range of R=0.2. The equilibrium is
labeled as 'eqm', whereas the other structures are labeled according to their direction and
offsets up to 4 decimal places.


Parallel line-search
^^^^^^^^^^^^^^^^^^^^

The parallel line-search is a collection of line-searches along directions assumed to be
independent apart from the common equilibrium. Since the parallel line-search is also
considered iterative, each step produces a file layout such as follows:

.. code-block:: bash

    pls/
    ├── d0_+0.2000/
    │   ├── params.in
    │   ├── sigma.in
    │   ├── error.out
    │   └── value.out
    ├── d0_-0.2000/
    ├── d1_+0.1500/
    ├── d1_+0.0750/
    ├── d1_-0.0750/
    ├── d1_-0.1500/
    ...
    ├── eqm/
    ├── ls0/
    │   ├── fit.out
    │   ├── ls.out
    │   ├── x0.out
    │   └── y0.out
    ├── ls1/
    ...

so it is similar to a single line-search, but with more directions. The solved line-searches
have been serialized in 'ls0/' and 'ls1' subdirectories. It can be produced as follows,
assuming that the structure 'p' has at least 2 parameters:

.. code-block:: python

    hessian = ParameterHessian(structure=p, hessian=[[1.0, 0.0], [0.0, 1.0]])
    pls = ParallelLineSearch(hessian=hessian, windows=[0.3, 0.5], sigmas=[0.04, 0.04])
    pls.propagate(pls, path='pls')

Line-search iteration
^^^^^^^^^^^^^^^^^^^^^

Finally, a collection of parallel line-searches form an iteration sequence under a common
path, such as follows:

.. code-block:: bash

    lsi/
    ├── pls0/
    │   ├── d0_+0.2000/
    │   │   ├── params.in
    │   │   ├── sigma.in
    │   │   ├── error.out
    │   │   └── value.out
    │   ...
    │   └── eqm/
    ├── pls1/
    │   └── ...
    └── pls2/
        └── ...
        ...

Summary
-------

In summary, the file I/O in STALK is an optional feature but highly useful for reducing
heavy operations and reproducing series subject to random noise. The feature is enabled by
supplying 'path' keyword to operations, such as PES evaluation, line-search and
optimization of line-search parameters.

The file I/O creates a layout of files, which may be used to inspect and edit files, and
restore and branch extensive workflows.