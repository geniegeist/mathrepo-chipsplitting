"""
Module contains general linear forms.
"""

import numpy as np
from numpy._typing import NDArray

from .base_linear_form import BaseLinearForm
from .hyperfield_linear_form import HyperfieldLinearForm

from .utils import gauss
from .utils.coordinate_transformation import get_array_index


class LinearForm(BaseLinearForm):
    """
    Class for general linear forms
    """

    def __init__(self, support_pos: NDArray[np.int_], support_neg: NDArray[np.int_]):
        """
        A linear form is just a sum of x_ij.

        """
        assert np.all(support_pos >= 0)
        assert np.all(support_neg >= 0)

        self._support_pos = np.array(support_pos)
        self._support_neg = np.array(support_neg)

    @classmethod
    def zero(cls, degree: int):
        """
        Returns a zero linear form of the given degree.
        """
        return cls(np.array([0] * gauss(degree + 1)), np.array([0] * gauss(degree + 1)))

    @property
    def support_neg(self) -> NDArray[np.int_]:
        return self._support_neg

    @property
    def support_pos(self) -> NDArray[np.int_]:
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
            f"LinearForm(positive support: {self.support_pos}, "
            f"negative support: {self.support_neg})"
        )

    def __str__(self) -> str:
        txt = ""
        for len_row, row_index in enumerate(range(self.degree, -1, -1)):
            for col_index in range(len_row + 1):
                pos_val = self.support_pos[get_array_index(col_index, row_index)]
                neg_val = self.support_neg[get_array_index(col_index, row_index)]
                val = (
                    str(pos_val)
                    if pos_val > 0
                    else f"-{str(neg_val)}" if neg_val > 0 else "."
                )
                txt += f"{''.join([' '] * (4 - len(val))) + val} "
            txt += "\n"
        return txt

    def __call__(self, v):
        pos = 0
        neg = 0

        for i in range(self.support_pos.size):
            pos += self.support_pos[i] * v[i]

        for i in range(self.support_neg.size):
            neg += self.support_neg[i] * v[i]

        return pos - neg

    def __hash__(self):
        return hash((tuple(self.support_pos), tuple(self.support_neg)))

    def __eq__(self, other):
        return np.all(self.support_pos == other.support_pos) and np.all(
            self.support_neg == other.support_neg
        )

    def __neg__(self):
        return LinearForm(self.support_neg, self.support_pos)

    def __add__(self, other):
        support_pos = np.array([0] * gauss(self.degree + 1))
        support_neg = np.array([0] * gauss(self.degree + 1))

        for i in range(support_pos.size):
            val = (
                self.support_pos[i]
                - self.support_neg[i]
                + other.support_pos[i]
                - other.support_neg[i]
            )
            if val > 0:
                support_pos[i] = val
            elif val < 0:
                support_neg[i] = -val

        return LinearForm(support_pos, support_neg)

    def __sub__(self, other):
        support_pos = np.array([0] * gauss(self.degree + 1))
        support_neg = np.array([0] * gauss(self.degree + 1))

        for i in range(support_pos.size):
            val = (
                self.support_pos[i]
                - self.support_neg[i]
                - other.support_pos[i]
                + other.support_neg[i]
            )
            if val > 0:
                support_pos[i] = val
            elif val < 0:
                support_neg[i] = -val

        return LinearForm(support_pos, support_neg)

    def to_hyperfield(self) -> HyperfieldLinearForm:
        """
        Converts the Pascal form to a hyperfield linear form.
        """
        support_pos = np.full(gauss(self.degree + 1), False)
        support_neg = np.full(gauss(self.degree + 1), False)

        for index, val in enumerate(self.support_pos):
            if val > 0:
                support_pos[index] = True

        for (
            index,
            val,
        ) in enumerate(self.support_neg):
            if val > 0:
                support_neg[index] = True

        return HyperfieldLinearForm(support_pos, support_neg)

    def get(self, contraction_size, key):
        if type(key) is not str:
            raise NotImplementedError(f"key {key} not implemented")

        letter, index = key.split("_")
        index = int(index)

        if index >= contraction_size:
            raise Exception("Invalid index")

        res = []

        if letter in ("d0", "d"):
            for col in range(
                contraction_size, self.degree - contraction_size - index + 1, 2
            ):
                res.append(
                    self.support_pos[get_array_index(col, self.degree - col - index)]
                    - self.support_neg[get_array_index(col, self.degree - col - index)]
                )
        elif letter in ("d1", "e"):
            for col in range(
                contraction_size + 1, self.degree - contraction_size - index + 1, 2
            ):
                res.append(
                    self.support_pos[get_array_index(col, self.degree - col - index)]
                    - self.support_neg[get_array_index(col, self.degree - col - index)]
                )
        elif letter == "b":
            for col in range(
                contraction_size, self.degree - contraction_size - index + 1
            ):
                res.append(
                    self.support_pos[get_array_index(col, index)]
                    - self.support_neg[get_array_index(col, index)]
                )
        elif letter == "c":
            for row in range(
                contraction_size, self.degree - contraction_size - index + 1
            ):
                res.append(
                    self.support_pos[get_array_index(index, row)]
                    - self.support_neg[get_array_index(index, row)]
                )
        else:
            raise NotImplementedError(f"key {key} not implemented")

        return np.array(res)

    def is_contractable(self, contraction_size):
        assert self.degree >= contraction_size * 3 - 1

        for index in range(contraction_size):
            for letter in ["b", "c", "d", "e"]:
                vector = self.get(contraction_size, f"{letter}_{index}")
                if not (
                    np.all(np.sign(vector) == [1] * len(vector))
                    or np.all(np.sign(vector) == [-1] * len(vector))
                    or np.all(np.sign(vector) == [0] * len(vector))
                ):
                    return False

        return True
