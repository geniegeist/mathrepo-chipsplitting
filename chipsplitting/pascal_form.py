"""
Module for Pascal forms
"""

import numpy as np

from .hyperfield_linear_form import HyperfieldLinearForm
from .utils.binomial import ncr

from .linear_form import LinearForm
from .utils import gauss, get_array_index


def make_supports_diag_pascal(degree: int, unit: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Generates the positive and negative support for a diagonal Pascal form.

    :param degree: The degree of the Pascal form.
    :param unit: The unit of the Pascal form.
    :return: A tuple containing the positive and negative support.
    """
    support_pos = np.full(gauss(degree + 1), 0)
    support_neg = np.full(gauss(degree + 1), 0)

    for c in range(degree + 1):
        for r in range(degree + 1):
            if c + r > degree:
                continue
            n = degree - c - r
            k = unit - c
            support_pos[get_array_index(c, r)] = 0 if k < 0 or n < k else ncr(n, k)

    return support_pos, support_neg


def make_supports_col_pascal(degree: int, unit: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Generates the support of a column pascal basis function.
    """
    support_pos = np.full(gauss(degree + 1), 0)
    support_neg = np.full(gauss(degree + 1), 0)
    for c in range(degree + 1):
        for r in range(degree + 1):
            if c + r > degree:
                continue
            sign = 1 if unit % 2 == 0 else -1
            inner_sign = 1 if r % 2 == 0 else -1
            n = c
            k = unit - r
            if sign * inner_sign == 1:
                support_pos[get_array_index(c, r)] = 0 if k < 0 or n < k else ncr(n, k)
            else:
                support_neg[get_array_index(c, r)] = 0 if k < 0 or n < k else ncr(n, k)
    return support_pos, support_neg


def make_supports_row_pascal(degree: int, unit: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Generates the support of a row pascal basis function
    """
    support_pos = np.full(gauss(degree + 1), 0)
    support_neg = np.full(gauss(degree + 1), 0)
    for c in range(degree + 1):
        for r in range(degree + 1):
            if c + r > degree:
                continue
            sign = 1 if unit % 2 == 0 else -1
            inner_sign = 1 if c % 2 == 0 else -1
            n = r
            k = unit - c
            if sign * inner_sign == 1:
                support_pos[get_array_index(c, r)] = 0 if k < 0 or n < k else ncr(n, k)
            else:
                support_neg[get_array_index(c, r)] = 0 if k < 0 or n < k else ncr(n, k)
    return support_pos, support_neg


class PascalForm(LinearForm):
    """
    Class for Pascal forms.
    """

    def __init__(self, degree: int, mode: str, unit: int):
        """
        Constructor for Pascal forms.

        :param degree: The degree of the Pascal form.
        :param mode: The mode of the Pascal form. Can be 'diagonal', 'column' or 'row'.
        :param unit: The unit of the Pascal form. Can be any integer from 0 to degree.
        """
        support_pos = np.full(gauss(degree + 1), 0)
        support_neg = np.full(gauss(degree + 1), 0)
        self._mode = mode

        if mode in ("diagonal", "diag"):
            support_pos, support_neg = make_supports_diag_pascal(degree, unit)
        elif mode in ("column", "col"):
            support_pos, support_neg = make_supports_col_pascal(degree, unit)
        elif mode == "row":
            support_pos, support_neg = make_supports_row_pascal(degree, unit)
        else:
            raise ValueError("Mode must be 'diagonal', 'column' or 'row'.")

        super().__init__(support_pos, support_neg)
