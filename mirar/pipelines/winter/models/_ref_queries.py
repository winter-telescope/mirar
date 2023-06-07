"""
Module to make reference components table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Float, ForeignKey, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.pipelines.winter.models.base_model import WinterBase, dec_field, ra_field
from mirar.processors.sqldatabase.base_model import BaseDB


class RefQueriesTable(WinterBase):
    """
    Table for queries made to the WFAU server
    """

    __tablename__ = "refqueries"

    queryid = Column(Integer, primary_key=True)
    query_ra = Column(Float)
    query_dec = Column(Float)
    query_filt = Column(VARCHAR(10))
    compid: Mapped[int] = mapped_column(ForeignKey("refcomponents.compid"))
    components: Mapped["RefComponentsTable"] = relationship(back_populates="queries")


class RefQueries(BaseDB):
    """
    Pydantic model for Reference queries
    """

    sql_model: ClassVar = RefQueriesTable

    query_ra: float = ra_field
    query_dec: float = dec_field
    query_filt: str = Field(min_length=1)
    compid: int = Field(ge=0)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(
            values=[self.query_ra, self.query_dec, self.query_filt],
            keys=["query_ra", "query_dec", "query_filt"],
        )
