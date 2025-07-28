================================================================
One-dimensional Discrete Models of Maximum Likelihood Degree One
================================================================

| This page contains the source code and explanations related to the computational
| results presented in the paper:
| Carlos Améndola, Viet Duc Nguyen Janike Oldekop: One-dimensional Discrete Models of Maximum Likelihood Degree One
| ARXIV: 

ABSTRACT: We settle a conjecture by Bik and Marigliano stating that the degree of a one-dimensional discrete model with rational maximum likelihood estimator is bounded above by a linear function in the number of states, therefore showing that there are only finitely many fundamental such models for any given number of states. We study these models from a combinatorial perspective with regard to their existence and enumeration. In particular, sharp models, those whose degree attains the maximal bound, enjoy special properties and have been studied as rational maps between unit spheres. In this way, we present a new link between Cauchy-Riemann geometry and algebraic statistics.





Required Software
~~~~~~~~~~~~~~~~~

- Python 3
- Cython
- NumPy

Fundamental Models
~~~~~~~~~~~~~~~~~~

A one-dimensional discrete model of maximum likelihood degree one is a subset :math:`\mathcal{M}` of the probability simplex :math:`\Delta_n` given by a parameterization of the form

.. math::
    p : [0,1] \to \Delta_n, \qquad t \mapsto (c_it^{\nu_i}(1-t)^{\mu_i})_{i=0}^n,

where :math:`\nu_i` and :math:`\mu_i` are non-negative exponents and :math:`c_i` are positive real coefficients. Such a model is fundamental if the exponent pairs :math:`(\nu_i,\mu_i)` are pairwise distinct, different from :math:`(0,0)` and, given :math:`(\nu_i,\mu_i)`, the scalings :math:`c_i` are uniquely determined by the constraint

.. math::
    c_0t^{\nu_0}(1-t)^{\mu_0} + c_1t^{\nu_1}(1-t)^{\mu_1} + \ldots + c_nt^{\nu_n}(1-t)^{\mu_n} = 1.

Its degree is defined as :math:`\textup{deg}(\mathcal{M}) = \max \{ \nu_i + \mu_i \mid i \in [n] \}`. According to Section 6 in the paper, fundamental models exist whenever :math:`n\le d\le 2n-1`. The number of fundamental models for small :math:`n,d\in\mathbb{N}^+` is presented in Table 1 and Figure 9 in the paper. The source code can be viewed here:

.. toctree::
   :maxdepth: 1
   :glob:
   :titlesonly:
   
   support_candidates
   fundamental_models

:download:`support_candidates.ipynb <support_candidates.ipynb>`

:download:`fundamental_models.ipynb <fundamental_models.ipynb>`

The notebooks make use of a `chipsplitting library`_ developed for this project.

A dump of candidate supports can be downloaded from the repository (see the ``data/`` directory). Note that cases for degree :math:`d \geq 10` have been excluded due to large file sizes (over 100mb). They can be recomputed by running the ``support_candidates.ipynb`` notebook, though this may take several hours.

The dataset of all fundamental models, computed by ``fundamental_models.ipynb``, is also available there (see the ``fundamental-models/`` directory).

.. _chipsplitting library: https://gitlab.mis.mpg.de/rok/mathrepo/-/tree/master/source/FundamentalModels

-------------------------------------------------------------------------

Project page created: 25/07/2025

Project contributors: Carlos Améndola, Viet Duc Nguyen, Janike Oldekop

Corresponding author of this page: 

Software used: The computations were performed using Python 3, with the libraries Cython and NumPy. The code was executed within a Jupyter Notebook environment.

System setup used: The analysis was run on a MacBook Pro (2023) equipped with an Apple M3 Pro processor and 36 GB of RAM, running macOS 15.5.

License for code of this project page: MIT License (https://spdx.org/licenses/MIT.html)

License for all other content of this project page: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

Last updated 25/07/2025.



