"""
Models for the 'diff' table
"""
import os
from typing import ClassVar

import pandas as pd
from pydantic import Field, field_validator
from sqlalchemy import VARCHAR, Column, Float, ForeignKey, Integer, Sequence
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.database.base_model import BaseDB
from mirar.database.constraints import DBQueryConstraints
from mirar.database.transactions.update import _update_database_entry
from mirar.pipelines.winter.models._candidates import Candidate
from mirar.pipelines.winter.models.base_model import WinterBase


class DiffsTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Diffs table in database
    """

    __tablename__ = "diffs"
    __table_args__ = {"extend_existing": True}

    diffid = Column(
        Integer,
        Sequence(start=1, name="diff_udiffid_seq"),
        autoincrement=True,
        unique=True,
        primary_key=True,
    )

    stackid: Mapped[int] = mapped_column(ForeignKey("stacks.stackid"))
    stack_id: Mapped["StacksTable"] = relationship(back_populates="diff")

    scormean = Column(Float, nullable=True)
    scormed = Column(Float, nullable=True)
    scorstd = Column(Float, nullable=True)
    maglim = Column(Float, nullable=True)
    zp = Column(Float, nullable=True)
    savepath = Column(VARCHAR(255), unique=True)

    candidates = relationship("CandidatesTable", back_populates="diff_id")


class Diff(BaseDB):
    """
    A pydantic model for a diff database entry
    """

    sql_model: ClassVar = DiffsTable

    savepath: str = Field(min_length=1)
    stackid: int = Field(ge=0)
    scormean: float | None = Field()
    scormed: float | None = Field()
    scorstd: float | None = Field(ge=0)
    zp: float | None = Field()
    maglim: float | None = Field()

    @field_validator("savepath")
    @classmethod
    def validate_savepath(cls, savepath: str) -> str:
        """
        Ensure that path exists

        :param savepath: field value
        :return: field value
        """
        assert os.path.exists(savepath)
        return savepath

    def insert_entry(
        self, returning_key_names: str | list[str] | None = None
    ) -> pd.DataFrame:
        """
        Insert entry into database
        :param returning_key_names: names of keys to return
        :return: dataframe of inserted entries
        """
        dbconstraints = DBQueryConstraints()
        dbconstraints.add_constraint(
            column="stackid",
            accepted_values=self.stackid,
        )
        _update_database_entry(
            update_dict={"deprecated": True},
            sql_table=Candidate.sql_model,
            db_constraints=dbconstraints,
        )

        return self._insert_entry(returning_key_names=returning_key_names)
