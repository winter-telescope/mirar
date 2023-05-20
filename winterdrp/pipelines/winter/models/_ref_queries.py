"""
Module to make reference components table
"""
from typing import ClassVar

from pydantic import Field
from sqlalchemy import Column, Float, ForeignKey, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship

from winterdrp.pipelines.winter.models.base_model import WinterBase, dec_field, ra_field
from winterdrp.processors.sqldatabase.base_model import BaseDB


class RefQueriesTable(WinterBase):
    """
    Table for queries made to the WFAU server
    """

    __tablename__ = "refqueries"

    queryid = Column(Integer, primary_key=True)
    query_ra = Column(Float)
    query_dec = Column(Float)
    compid: Mapped[int] = mapped_column(ForeignKey("refcomponents.compid"))
    components: Mapped["RefComponentsTable"] = relationship(back_populates="queries")


class RefQueries(BaseDB):
    """
    Pydantic model for Reference queries
    """

    sql_model: ClassVar = RefQueriesTable

    query_ra: float = ra_field
    query_dec: float = dec_field
    compid: int = Field(ge=0)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self.sql_model().exists(
            values=[self.query_ra, self.query_dec], keys=["query_ra", "query_dec"]
        )
