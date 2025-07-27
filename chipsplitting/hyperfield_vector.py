import numpy as np

from . import utils

class HyperfieldVector:
    """
    A hyperfield vector is a vector that only contains -1, 0 and 1 values.

    :param values: A list of values that represent the vector. It can be a 2D triangle or a flat list.
    """

    def __init__(self, values: list):
        if len(values) == 0:
            self.values = np.array([])
        elif isinstance(values[0], list):
            flat = np.full(len(np.concatenate(values)), 0)
            for r, row in enumerate(values):
                for c, val in enumerate(row):
                    flat[utils.get_array_index(c, len(values) - 1 - r)] = val
                self.values = np.array(flat)
        else:
            self.values = np.array(values)

        for x in self.values:
            assert x in [
                -1,
                0,
                1,
            ], "Hyperfield vector may only contain -1, 0 and 1 values."

    @property
    def degree(self):
        return int(-1.5 + np.sqrt(0.25 + 2 * self.values.size))

    @property
    def is_valid(self):
        if self.values[0] > 0:
            return False
        for v in self.values[1:]:
            if v < 0:
                return False
        return True

    def __getitem__(self, key):
        return self.values[key]

    def __repr__(self) -> str:
        return f"HyperfieldVector({self.values})"

    def __str__(self) -> str:
        txt = ""
        for len_row, row_index in enumerate(range(self.degree, -1, -1)):
            for col_index in range(len_row + 1):
                val = str(self.values[utils.get_array_index(col_index, row_index)])
                txt += f"{''.join([' '] * (4 - len(val))) + val} "
            txt += "\n"
        return txt

    def __lor__(self, other):
        return self.values | other

    def contract(self, contraction_size: int):
        """
        Converts the hyperfield vector to a hyperfield contraction vector.
        """
        assert self.degree >= contraction_size * 3, "Degree must be at least 12"

        num_coordinates = contraction_size**2 * 3 + contraction_size * 4
        values = np.full(num_coordinates, 0)

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
            values[index_to_update] = self.values[i]

        for c, i in enumerate(y_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2
            )
            values[index_to_update] = self.values[i]

        for c, i in enumerate(z_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 2
            )
            values[index_to_update] = self.values[i]

        for c, i in enumerate(b_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
            )
            indexes = [
                (i, c)
                for i in range(contraction_size, self.degree - contraction_size + 1 - c)
            ]
            s = [self.values[utils.get_array_index(i, j)] for i, j in indexes]

            if -1 in s and 1 in s:
                values[index_to_update] = np.nan
            elif 1 in s:
                values[index_to_update] = 1
            elif -1 in s:
                values[index_to_update] = -1
            else:
                values[index_to_update] = 0

        for c, i in enumerate(c_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
                + contraction_size
            )
            indexes = [
                (c, i)
                for i in range(contraction_size, self.degree - contraction_size + 1 - c)
            ]
            s = [self.values[utils.get_array_index(i, j)] for i, j in indexes]

            if -1 in s and 1 in s:
                values[index_to_update] = np.nan
            elif 1 in s:
                values[index_to_update] = 1
            elif -1 in s:
                values[index_to_update] = -1
            else:
                values[index_to_update] = 0

        for c, i in enumerate(d0_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
                + contraction_size * 2
            )
            indexes = [
                (col, self.degree - col - c)
                for col in range(
                    contraction_size, self.degree - contraction_size + 1 - c, 2
                )
            ]
            s = [self.values[utils.get_array_index(i, j)] for i, j in indexes]

            if -1 in s and 1 in s:
                values[index_to_update] = np.nan
            elif 1 in s:
                values[index_to_update] = 1
            elif -1 in s:
                values[index_to_update] = -1
            else:
                values[index_to_update] = 0

        for c, i in enumerate(d1_coordinates):
            index_to_update = (
                (c // contraction_size) * contraction_size
                + c % contraction_size
                + contraction_size**2 * 3
                + contraction_size * 3
            )
            indexes = [
                (col, self.degree - col - c)
                for col in range(
                    contraction_size + 1, self.degree - contraction_size + 1 - c, 2
                )
            ]
            s = [self.values[utils.get_array_index(i, j)] for i, j in indexes]

            if -1 in s and 1 in s:
                values[index_to_update] = np.nan
            elif 1 in s:
                values[index_to_update] = 1
            elif -1 in s:
                values[index_to_update] = -1
            else:
                values[index_to_update] = 0

        return HyperfieldVector(values)
