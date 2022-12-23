"""
Module for DBQueryConstraints to carefully specify postgres query constraints
"""
from typing import Optional

import numpy as np

POSTGRES_ACCEPTED_COMPARISONS = ["=", "<", ">", "between"]


class DBQueryConstraints:
    """
    Object containing one or more postgres query constraints
    """

    def __init__(
        self,
        columns: str | list[str],
        accepted_values: str | int | float | list[str | float | int | list],
        comparison_types: Optional[str | list[str]] = None,
    ):
        self.columns = []
        self.accepted_values = []
        self.comparison_types = []

        if not isinstance(columns, list):
            columns = [columns]
        if not isinstance(accepted_values, list):
            accepted_values = [accepted_values]

        assert len(columns) == len(accepted_values)

        if comparison_types is None:
            comparison_types = ["="] * len(accepted_values)
        assert len(comparison_types) == len(accepted_values)

        for i, column in enumerate(columns):
            self.add_constraint(
                column=column,
                accepted_values=accepted_values[i],
                comparison_type=comparison_types[i],
            )

    def add_constraint(
        self,
        column: str,
        accepted_values: str | int | float | tuple[float, float] | tuple[int, int],
        comparison_type: str = "=",
    ):
        """
        Add a new constraint

        :param column: column
        :param accepted_values: accepted value for comparison
        :param comparison_type: type of comparison, e.g '='
        :return: None
        """
        assert comparison_type in POSTGRES_ACCEPTED_COMPARISONS

        if comparison_type == "between":
            assert np.logical_and(
                isinstance(accepted_values, tuple), len(accepted_values) == 2
            )

        self.columns.append(column)
        self.accepted_values.append(accepted_values)
        self.comparison_types.append(comparison_type)

    def __add__(self, other):
        new = self.__class__(self.columns, self.accepted_values, self.comparison_types)
        for args in other:
            new.add_constraint(*args)
        return new

    def __iadd__(self, other):
        for args in other:
            self.add_constraint(*args)
        return self

    def __len__(self):
        return self.columns.__len__()

    def __iter__(self):
        return iter(zip(self.columns, self.accepted_values, self.comparison_types))

    def parse_constraints(
        self,
    ) -> str:
        """
        Converts the list of constraints to sql

        :return: sql string
        """
        constraints = []
        for i, column in enumerate(self.columns):
            if self.comparison_types[i] == "between":
                constraints.append(
                    f"{column} between {self.accepted_values[i][0]} "
                    f"and {self.accepted_values[i][1]}"
                )
            else:
                constraints.append(
                    f"{column} {self.comparison_types[i]} {self.accepted_values[i]}"
                )

        return " AND ".join(constraints)
