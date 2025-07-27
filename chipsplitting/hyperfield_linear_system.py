"""
This module contains everything related to solving linear systems 
of homogeneous linear forms in a hyperfield.
"""

from collections import deque

import numpy as np

from .hyperfield_linear_form import HyperfieldLinearForm
from .hyperfield_vector import HyperfieldVector
from . import solver_ext

class HyperfieldHomogeneousLinearSystem:
    """
    A system of homogeneous linear forms in a hyperfield.
    It contains a list of linear forms that we want to roots for.
    """

    def __init__(self, linear_forms: list[HyperfieldLinearForm]):
        self.linear_forms = linear_forms

        conditions = []
        for form in linear_forms:
            if form.support_pos[0]:
                cond = np.copy(form.support_pos)
                cond[0] = False
                conditions.append(form.support_pos)
            elif form.support_neg[0]:
                cond = np.copy(form.support_neg)
                cond[0] = False
                conditions.append(form.support_neg)
            elif np.any(form.support_pos) and np.any(form.support_neg):
                conditions.append(form.support_pos)
                conditions.append(form.support_neg)
        self.conditions = np.array(conditions)

    def __str__(self) -> str:
        return "\n".join(str(form) for form in self.linear_forms)

    def is_solved_by_valid(self, vector: HyperfieldVector) -> bool:
        """
        Returns true if the vector is a root of every linear form in the system.
        Vector must be valid, i.e. its negative support must be empty or may only contain zero.
        """
        assert vector.is_valid, "Vector must be valid."
        support_without_neg = np.array(vector.values, dtype=bool)
        support_without_neg[0] = False
        return bool(np.all(self.conditions.dot(support_without_neg)))

    def make_constraints(self) -> list[set[int]]:
        """
        Given the list of linear forms, compute a list constraints.
        Each constraint must be satisfied for a configuration to be considered valid.
        A configuration satisfies a constraint if and only if some component
        of the configuration is in contained in the constraint.
        """
        constraints = []
        # contains tuples of positive and negative support for each linear form
        supports = [
            (set(form.support_pos.nonzero()[0]), set(form.support_neg.nonzero()[0]))
            for form in self.linear_forms
        ]

        for pos, neg in supports:
            if 0 in pos:
                new_constr = {i for i in pos if i > 0}
                if new_constr not in constraints:
                    constraints.append(new_constr)
            elif 0 in neg:
                new_constr = {i for i in neg if i > 0}
                if new_constr not in constraints:
                    constraints.append(new_constr)
            elif len(pos) and len(neg):
                pos_constr = set(pos)
                neg_constr = set(neg)
                if pos_constr not in constraints:
                    constraints.append(pos_constr)
                if neg_constr not in constraints:
                    constraints.append(neg_constr)
            else:
                raise UnsolvableSystemException("Unsolvable system of linear forms")

        to_remove = []
        for c in constraints:
            for d in constraints:
                if c == d:
                    continue
                intersected = c.intersection(d)
                if len(intersected) == len(c):
                    to_remove.append(d)
                elif len(intersected) == len(d):
                    to_remove.append(c)

        return [list(x) for x in constraints if x not in to_remove]

    def quick_solve_loop_fast(self, support_size: int):
        constraints = self.make_constraints()
        return solver_ext.quick_solve_loop_cython_int16(constraints, support_size)


class UnsolvableSystemException(Exception):
    """
    The system is not solvable.
    """
