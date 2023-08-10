"""
Module for DBQueryConstraints to carefully specify postgres query constraints
"""
import numpy as np

POSTGRES_ACCEPTED_COMPARISONS = ["=", "<", ">", "<=", ">=", "between", "<>", "!="]


class DBQueryConstraints:
    """
    Object containing one or more postgres query constraints
    """

    def __init__(
        self,
        columns: str | list[str] | None = None,
        accepted_values: str
        | int
        | float
        | list[str | float | int | list]
        | None = None,
        comparison_types: str | list[str] | None = None,
    ):
        self.columns = []
        self.accepted_values = []
        self.comparison_types = []

        self.q3c_query = None

        if columns is None:
            assert accepted_values is None

        else:
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

    def add_q3c_constraint(
        self,
        ra: float,
        dec: float,
        crossmatch_radius_arcsec: float,
        ra_field_name: str = "ra",
        dec_field_name: str = "dec",
    ):
        """
        Add a q3c constraint

        :param ra: ra of source
        :param dec: dec of source
        :param crossmatch_radius_arcsec: crossmatch radius in arcsec
        :param ra_field_name: ra field name in database
        :param dec_field_name: dec field name in database
        :return: None
        """

        crossmatch_radius_deg = crossmatch_radius_arcsec / 3600.0

        constraints = (
            f"q3c_radial_query({ra_field_name},{dec_field_name},"
            f"{ra},{dec},{crossmatch_radius_deg}) "
        )

        self.q3c_query = constraints

    def __add__(self, other):
        new = self.__class__(self.columns, self.accepted_values, self.comparison_types)
        new.q3c_query = self.q3c_query
        for args in other:
            new.add_constraint(*args)
            if other.q3c_query is not None:
                new.add_q3c_constraint(*other.q3c_query)
        return new

    def __iadd__(self, other):
        for args in other:
            self.add_constraint(*args)
            if other.q3c_query is not None:
                self.add_q3c_constraint(*other.q3c_query)
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

        if self.q3c_query is not None:
            constraints.append(self.q3c_query)

        for i, column in enumerate(self.columns):
            if self.comparison_types[i] == "between":
                constraints.append(
                    f"{column.lower()} between {self.accepted_values[i][0]} "
                    f"and {self.accepted_values[i][1]}"
                )
            else:
                constraints.append(
                    f"{column.lower()} {self.comparison_types[i]} "
                    f"'{self.accepted_values[i]}'"
                )

        return " AND ".join(constraints)
