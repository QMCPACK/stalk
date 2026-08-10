STALK: The orientation
======================

This orientations is an introduction to the basic concepts of STALK in terms of Python
implementation.

Potential energy surface
------------------------

What is the problem we want to solve using STALK? It can be one of a few, but really the
least common denominator is the existence of a potential energy surface (PES). The PES is a
mapping of a set of parameters `p` into a scalar-valued number `E`. Typically, we want to
find the special set of parameters `p^*` that minimizes `E`, denoted `E_0`.

Consider a simple 1-dimensional PES, a parabola:

.. code-block:: python

    # NB: this is not quite a STALK PES, but close enough to illustrate the idea.
    def pes(params: list[float]) -> float:
        c = [1.0, 2.0, 3.0]
        E = c[0] * params[0]**2 + c[1] * params[0] + c[2]
        return E
    # end def

Where params is a 1-item list of parameters. Whatever scalar is supplied in `params[0]` maps
into a scalar return value `E`. The coefficients `c` characterize the PES, determining all
values of `E(p)`, including the minimum `E_0 = E(p^*)`.

In STALK, the set of parameters is treated as a class :func:`stalk.params.ParameterSet`,
which contains the array of numerical parameters, which can be multidimensional. The
dimensionality (here: 1) should match that expected by the PES function. Therefore, a
consistent set of basic STALK objects could look as follows:

.. code-block:: python

    # A raw PES function with fixed coefficients
    def pes_func_fixed_c(params: ParameterSet):
        c = [1.0, 2.0, 3.0]
        E = c[0] * params[0]**2 + c[1] * params[0] + c[2]
        return E
    # end def
    # Define properly wrapped-up PES function
    pes = PesFunction(pes_func_fixed_c)
    # An instance of a set of parameters
    p = ParameterSet([1.0])
    Ef = pes(p)  # Ef = 6.0

The stalk.PesFunction wrapper allows for the PES function to take extra arguments for
convenience. Thus, the above definition could be generalized to any 1D parabola as follows:

.. code-block:: python

    # A generic PES function with kwargs coefficients
    def pes_func(params: ParameterSet, c: list = [1.0, 2.0, 3.0]) -> float:
        E = c[0] * params[0]**2 + c[1] * params[0] + c[2]
        return E
    # end def
    # The PES wrapper
    pes1 = PesFunction(pes_func, c=[1.0, 2.0, 3.0])
    pes2 = PesFunction(pes_func, c=[1.0, -2.0, 0.0])

The two PESs thus defined are different (just not in the parabolic shape), so they will
evaluate differently: 

.. code-block:: python

    p = ParameterSet([1.0])
    pes1(p)  # p.value = 6.0
    pes2(p)  # p.value = -1.0

Consequently, they should also have different solutions. This is a fundamental idea in STALK
and in surrogate acceleration, where one (cheap) PES is used to inform the solving of
another, which is typically noisy, more expensive but also more important to solve.

PES optimization
----------------

In optimization, we our task is to find the minimum point of the PES, that is, the set of
parameters that result in the lowest energy. Here, the terms *PES* and *energy* are borrowed
from the context of quantum chemistry, but they generalize to any minimization (or
maximization) problem, where the energy is some *cost*, making PES a *cost function*.

In any case, to succeed in the optimization, we should first decide *how to* optimize.

To being with a trivial example, let us assume that we already *know* the analytic shape
of the PES. For example, a upward-facing parabola looks like this:

.. math:: 
  E(p) = c_0 p^2 + c_1 p + c_2,

where :math:`c_0 > 0`. The minimum point is then given by the zero of the derivative,
namely,

.. math::
    p^* = -\frac{c_1}{2 c_0}.

Now, if we already know the coefficients :math:`c_0` and :math:`c_1`, we can just calculate
:math:`p^*`. If not, we can *deduce* them: We evaluate the PES at three points, and then
solve for the coefficients. Why three? Because a parabola is a second-order polynomial
(:math:`n=2`) with :math:`n+1=3` coefficients, and we need, generally, :math:`n` different
points to solve for :math:`n` coefficients.

What if we do not know the analytic shape of the PES (or the shape is too complicated to
solve)? Then, the least we must require to make our task meaningful is that the PES have the
following properties:
* The PES has a local minimum with the parameters.
* The PES can be evaluated at any point near the presumed minimum.

In principle, we could then evaluate the PES at all points and choose that with the lowest
energy. In practice, this in unfeasible, as the PES is typically multidimensional and
continuous. Indeed, we should further assume:
* The PES is (ideally) smooth.
* The PES is (ideally) continuous.
The meaning of *ideally* will be clarified later.

Even when we do not know the analytic shape, we can still use deduction: making educated
guesses about the shape of the PES based on a few evaluations. For example, if we keep
taking moderate steps toward the direction of decreasing energy (the gradient), we will
eventually arrivate at a local minimum. Likewise, if we assume that the PES has a minimum
between two points, we can bracket it by evaluating the PES at a few points and then fitting
a curve to those points, whose minimum we can take as our next best guess.

These are outlines for numerical optimization methods. There are many available (STALK
included), as this is a well-studied and central problem in many disciplines. The approaches
can be classified into *global* and *local* methods. We will only treat the local
optimization: this means, we will try to find the nearest minimum to a given starting point.
The global optimization is more complicated, and we will not treat it here. 

The local optimization approaches can further be classified into *gradient-based* and
*gradient-free* methods. Gradients mean the partial derivatives of the PES with respect to
the parameters. When available, they can be used to substantially speed up the optimization.
This is how the most performant optimization algorithms work. STALK does not provide new
gradient-based methods, but it provides an interface for those implemented in SciPy. This
happens through the PesFunction.relax method.

.. code-block:: python

    # The initial parameter set
    p1 = ParameterSet([1.0])
    # The PES
    pes1 = PesFunction(pes_func, c = [1.0, 2.0, 3.0])
    # Classical optimization using the BFGS method at 1e-6 tolerance
    res1 = pes1.relax(p, method='BFGS', tol=1e-6)  # res1.x0 = x0, res1.y0 = y0

Different algorithms may have different virtues, including speed, robustness, and ability to
handle noise. In STALK, we will focus on the latter, and for that reason, let us consider
the meaning of a *noisy PES*.

Noisy PES
---------

Difficulties arise when the PES is *noisy*. Loosely, this means that each evaluation of the
PES has an element of statistical uncertainty, random error. Then, the PES is no more a
smooth and continous function upon subsequent evaluations but it has the form:

.. math::
  E(p, \sigma) = E(p, \sigma=0) + \sigma

where :math:`\sigma` is a measure of the noise that always evaluates differently. We may
only know its approximate statistical distribution. For example, :math:`\sigma` could be the
standard deviation of a normal distribution, white noise.

Let us try it out with a simple example in STALK:

.. code-block:: python

    # Let us associate our parameter set with finite noise.
    sigma = 0.1
    p = ParameterSet([1.0])
    n = 1000
    # We use the original PES, as any STALK PES can be made noisy upon request.
    pes = PesFunction(pes_func, c=[1.0, 2.0, 3.0])
    # With add_sigma=True, STALK adds random (white) noise to the PES measurement,
    # because the analytic PES does not have noise intrinsically. This appears artificial,
    # because it is. But it is a good way to illustrate noise and also to simulate its effect 
    # on derived properties like the optimization.
    energies = []
    for _ in range(n):
        pes(p, sigma=sigma, add_sigma=True)
        energies.append(p.value)
    # end for
    # Calculate the apparent standard deviation of the noisy PES evaluations
    sigma_out = np.std(energies)
    # Calculate the aggregate mean of the noisy PES evaluations
    E_mean = np.mean(energies)
    # Calculate the exact PES for reference
    E_exact = pes(p, add_sigma=False)
    # The histogram should resemble a normal distribution with standard deviation sigma.
    plt.hist(energies, bins=20)
    plt.axvline(E_exact, color='r', linestyle='dashed', linewidth=1, label='Exact PES value')
    plt.axvline(E_mean, color='g', linestyle='dashed', linewidth=1, label=f'Mean of {n} PES evaluations')
    plt.show()

You should notice that the histogram resembles a normal distribution with standard deviation
:math:`\sigma`. The apparent standard deviation :math:`\sigma_{out}` is close to the input
:math:`\sigma` but not exactly. This serves to demonstrate that working with random noise
using finite samples (like 1000) always has an element of uncertainty, including not only
the first moment (the mean) but also the second moment (the standard deviation) and beyond.

Let us now consider the uncertainty of the PES evaluation. In the STALK convention, we use
the words *value*, *error* for each set of *params* (parameters) that were evaluated with a
given PES. Try it yourself:

.. code-block:: python

    # Try before and after evaluation
    print(f'The input parameters: {p.params}')
    print(f'The input uncertainty sigma={p.sigma}')
    print(f'The result mean value={p.value}')
    print(f'The result uncertainty sigma={p.error}')

You should notice that the error matches the input uncertainty :math:`\sigma`, by design: we
requested the error level of :math:`\sigma` and STALK has provided it. According to the
earlier example, we also know we got it.

Under the hood, we may consider *error* as an estimate made from a finite number of
statistical samples :math:`N`. When the uncertainty is normally distributed (which is often
a decent approximation), the standard error of the mean (SEM) is calculated as:

.. math:: \mathrm{SEM} = \frac{\sigma}{\sqrt{N}}
    :label: eqsem

where the subscript :math:`N` in :math:`\sigma_N` emphasizes that we use the apparent
standard deviation from the actual samples without knowing :math:`\sigma` precisely. So, in
the above example, the error of each individual evaluation was :math:`\sigma=0.1`, but the
error of the mean value after :math:`N=1000` samples would :math:`\mathrm{SEM} = 0.0032`.
That is, we know that about 68% of the time, the mean value of the PES evaluation will be
within :math:`\pm 0.0032` of the true value.

The finite uncertainty is bad, because it disturbs the evaluation and optimization. Eq.
:eq:`eqsem` shows that  we can control error with sampling: :math:`\sigma` decreases with
:math:`N` and vanishes at infinity. However, since each measurement has a *cost* (whether in
time, money or other resources), we cannot make :math:`N` infinite in practice. Therefore,
we must suffer a finite error and pay a cost of

.. math::

    \mathrm{cost} \propto N \propto \sigma^{-2}

So, the smaller the uncertainty, the more expensive the evaluation.

When the PES is noisy, so will be its gradients. Therefore, the traditional gradient-based
methods face complications in the presence of noise. We will not treat them here, but
instead turn focus on the gradient-free methods, and the virtues of STALK.

*Note: the "error" only represents the statistical uncertainty, not systematic bias that may
also be in effect in any PES of a real-life context. The bias of a snapshot PES evaluation
is considered an intrinsic property of the PES, so in some way not a bias at all. Whether
the PES is biased (e.g. how much a numerical simulation deviates from a physical
measurement) can only be judged externally.*

Line-search optimization
------------------------

The usual gradient-based optimization methods face problems when the gradients cannot be
depended upon. This can be due to the noise or the some reason why the gradients are beyond
grasp (e.g. the cost function is a not numerical at all but an experimental observable). In
either case, there are times to use gradient-free methods.

One of the most robust is the line-search method. It is simple: we sample a number of points
along a line in the parameter space, and then find the minimum point along that line. In
doing so, we presume that the minimum lie between the end points. If it does not, it is not
hard to extend the search until it does.

Furthermore, because every evaluation has a cost, we need a prescription for estimating the
minimum point among a finite set of values. A simple but crude approach would pick the
lowest-energy point among the samples. A more dynamic approach is to fit a parameterized
curve to the sampled values and then find the analytic minimum of that curve.

The latter way is how STALK operates, and it is easy to try out with the following example:

.. code-block:: python

    import numpy as np
    points = np.linspace(-3, 3, 7)
    values = [pes(p.copy(params=[point])).value for point in points]
    # 'pf3' for 3rd order polynomial fit
    lsb = LineSearchBase(points, fit_kind='pf3')
    lsb.values = values
    res1 = lsb.search()  # res1.x0 = -1, res1.y0 = 2

The part of extracting the values from the PES looks a bit hacky, because the
:func:`stalk.ls.LineSearchBase` class only treats points and values, not the parametric
context in the PES. Thus, a more human-readable way to accomplish this is to use the
:func:`stalk.ls.LineSearch` class:

.. code-block:: python

    offsets = np.linspace(-3, 3, 7)
    ls = LineSearch(p, d=0, offsets)
    ls.evaluate(pes)
    print(ls1.fit_res)  # res1.x0 = -2, res1.y0 = 2

Comparison shows that the two line-searches yield different results for the x-displacement
:math:`x_0` (the energy-minimizing parameter) but the same results for the y-displacement
:math:`y_0` (the minimum energy). Indeed, they both find the right minimum but express it in
different bases: The first one uses absolute parameters :math:`p` while the second uses
relative offsets from the starting point :math:`p = [1.0]`.

In 1D this may seem inconsequential, but in higher dimensions it will be much more
convenient to use the relative offsets (the latter kind). Its virtues include that by
assuming that the reference point is close to the minimum, we always anticipate a solution
close to zero. In the end, the two representations are equivalent and the difference is
merely a matter of human comprehension.


Noisy line-search
-----------------

Let us now consider how noise affects the line-search. The following example shows how the
:func:`stalk.ls.LineSearch`` is set up to be evaluated with the *input noise*
:math:`\sigma`.

.. code-block:: python

    sigma = 0.1  # noise level for the PES evaluation
    R = 3  # spatial extent of the grid
    M = 7  # number of grid points
    offsets = np.linspace(-R, R, M)
    # Supply finite sigma
    ls = LineSearch(p, d=0, offsets=offsets, N=500, sigma=sigma)
    ls.evaluate(pes, add_sigma=True)
    print(ls.x0, ls.x0_err)  # ls.x0 ~ -2, ls.x0_err ~ 0.07
    print(ls.y0, ls.y0_err)  # ls.y0 ~ 2, ls.y0_err ~ 0.14

We notice that finite noise has propagated into finite errorbars of the line-search results.
The errorbars have been estimated by STALK by resampling of the fitting :math:`N` times with
different instances of (simulated) random noise. This is nothing more but a simple
estimation of the statistical uncertainty, apart of which :math:`N` has no other purpose.
While the estimation of the errorbars has its own uncertainty, it should be reliable and
reproducible for sufficiently large :math:`N`.

It can be very instructive to see how the errorbar responds to different values of the
line-search settings, like the input noise :math:`\sigma`, the grid extent :math:`R`, the
number of grid points :math:`M` or the initial proximity to the minimum.

We will investigate this more thoroughly later, because it is a central question in STALK
and the noisy optimization in general: If the way how we set up the line-search indeed
affects its performance and accuracy, we should thrive to do it in an optimal way.
