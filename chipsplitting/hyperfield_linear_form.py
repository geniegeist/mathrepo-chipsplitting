"""
Module contains linear forms for hyperfields.
"""

import numpy as np
from numpy._typing import NDArray
from .base_linear_form import BaseLinearForm
from . import utils
from .hyperfield_vector import HyperfieldVector


class HyperfieldLinearForm(BaseLinearForm):
    """
    Linear form for hyperfield
    """

    def __init__(self, support_pos: NDArray[np.bool_], support_neg: NDArray[np.bool_]):
        """
        A linear form in a hyperfield is just a sum of x_ij whose coefficients are either 1 or -1.

        :param support_pos: A list of boolean values that represent the positive support.
        :param support_neg: A list of boolean values that represent the negative support.
        """

        self._support_pos = np.array(support_pos)
        self._support_neg = np.array(support_neg)

    @property
    def support_neg(self):
        return self._support_neg

    @property
    def support_pos(self):
        return self._support_pos

    @property
    def degree(self):
        """
        The degree of the linear form. It is defined as argmax i+j,
        where (i,j) is contained in the support
        """
        return int(-1.5 + np.sqrt(0.25 + 2 * self.support_pos.size))

    def __repr__(self) -> str:
        return (
            f"HyperfieldLinearForm(positive support: {self.support_pos}, "
            f"negative support: {self.support_neg})"
        )

    def __str__(self) -> str:
        txt = ""
        for len_row, row_index in enumerate(range(self.degree, -1, -1)):
            for col_index in range(len_row + 1):
                pos_val = self.support_pos[utils.get_array_index(col_index, row_index)]
                neg_val = self.support_neg[utils.get_array_index(col_index, row_index)]
                val = "+" if pos_val else "-" if neg_val else "."
                txt += f"{''.join([' '] * (4 - len(val))) + val} "
            txt += "\n"
        return txt

    def __call__(self, v: HyperfieldVector):
        for x in v[self.support_pos | self.support_neg]:
            if np.isnan(x):
                return np.nan

        has_pos, has_neg = False, False

        for x in v[self.support_pos]:
            if x > 0:
                has_pos = True
            elif x < 0:
                has_neg = True

            if has_neg and has_pos:
                return np.nan

        for x in v[self.support_neg]:
            if x > 0:
                has_neg = True
            elif x < 0:
                has_pos = True

            if has_neg and has_pos:
                return np.nan

        if has_neg and has_pos:
            return np.nan

        if has_neg:
            return -1

        if has_pos:
            return 1

        return 0

    def __neg__(self):
        return HyperfieldLinearForm(self.support_neg, self.support_pos)

    def __eq__(self, other):
        return np.all(self.support_pos == other.support_pos) and np.all(
            self.support_neg == other.support_neg
        )

    def __hash__(self):
        return hash(tuple(list(self.support_pos) + list(self.support_neg)))

    def contract(self, contraction_size: int):
        """
        Converts the hyperfield linear form to a hyperfield contraction form.
        """
        assert self.degree >= contraction_size * 3 - 1, "Degree must be at least 12"

        num_coordinates = contraction_size**2 * 3 + contraction_size * 4
        support_pos = np.full(num_coordinates, False)
        support_neg = np.full(num_coordinates, False)

        x_coordinates = [
            utils.get_array_index(i, j)
            for i in range(contraction_size)
            for j in range(contraction_size)
        ]
        y_coordinates = [
            utils.get_array_index(i, self.degree - i - (contraction_size - 1) + j)
            for i in range(contraction_size)
            for j in range(contraction_size)
        ]
        z_coordinates = [
            utils.get_array_index(self.degree - (contraction_size - 1) - j + i, j)
            for i in range(contraction_size)
            for j in range(contraction_size)
        ]
        b_coordinates = [
            utils.get_array_index(contraction_size, i) for i in range(contraction_size)
        ]
        c_coordinates = [
            utils.get_array_index(i, contraction_size) for i in range(contraction_size)
        ]
        d0_coordinates = [
            utils.get_array_index(contraction_size, self.degree - contraction_size - i)
            for i in range(contraction_size)
        ]
        d1_coordinates = [
            utils.get_array_index(
                contraction_size + 1, self.degree - contraction_size - 1 - i
            )
            for i in range(contraction_size)
        ]

        for c, i in enumerate(x_coordinates):
            index_to_update = (
                c // contraction_size
            ) * contraction_size + c % contraction_size
            if self.support_pos[i]:
                support_pos[index_to_update] = True
            elif self.support_neg[i]:
                support_neg[index_to_update] = True

        for c, i in enumerate(y_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2
            )
            if self.support_pos[i]:
                support_pos[index_to_update] = True
            elif self.support_neg[i]:
                support_neg[index_to_update] = True

        for c, i in enumerate(z_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 2
            )
            if self.support_pos[i]:
                support_pos[index_to_update] = True
            elif self.support_neg[i]:
                support_neg[index_to_update] = True

        for c, i in enumerate(b_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
            )
            if self.support_pos[i]:
                support_pos[index_to_update] = True
            elif self.support_neg[i]:
                support_neg[index_to_update] = True

        for c, i in enumerate(c_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
                + contraction_size
            )
            if self.support_pos[i]:
                support_pos[index_to_update] = True
            elif self.support_neg[i]:
                support_neg[index_to_update] = True

        for c, i in enumerate(d0_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
                + contraction_size * 2
            )
            if self.support_pos[i]:
                support_pos[index_to_update] = True
            elif self.support_neg[i]:
                support_neg[index_to_update] = True

        for c, i in enumerate(d1_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
                + contraction_size * 3
            )
            if self.support_pos[i]:
                support_pos[index_to_update] = True
            elif self.support_neg[i]:
                support_neg[index_to_update] = True

        return HyperfieldLinearForm(support_pos, support_neg)
