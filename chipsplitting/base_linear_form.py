"""
Module containing the abstract class for linear forms
"""

import abc

from numpy._typing import NDArray

class BaseLinearForm(abc.ABC):
    """
    Abstract class for linear forms.
    A linear form is a polynomial of degree one without any constant.
    """

    @property
    @abc.abstractmethod
    def support_neg(self) -> NDArray:
        """
        Returns the negative support of the linear form.
        """

    @property
    @abc.abstractmethod
    def support_pos(self) -> NDArray:
        """
        Returns the positive support of the linear form.
        """

    def __repr__(self) -> str:
        return f"LinearForm(positive support: {self.support_pos}, negative support: {self.support_neg})"
