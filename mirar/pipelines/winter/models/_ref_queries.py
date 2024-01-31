"""
Module to make reference components table
"""

from typing import ClassVar

from pydantic import Field
from sqlalchemy import VARCHAR, Column, Float, ForeignKey, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.pipelines.winter.models.base_model import WinterBase


class RefQueriesTable(WinterBase):
    """
    Table for queries made to the WFAU server
    """

    __tablename__ = "refqueries"

    queryid = Column(Integer, primary_key=True)
    qry_ra = Column(Float)
    qry_dec = Column(Float)
    qry_filt = Column(VARCHAR(10))
    compid: Mapped[int] = mapped_column(ForeignKey("refcomponents.compid"))
    components: Mapped["RefComponentsTable"] = relationship(back_populates="queries")

    ra_column_name = "qry_ra"
    dec_column_name = "qry_dec"


class RefQuery(BaseDB):
    """
    Pydantic model for a reference query entry
    """

    sql_model: ClassVar = RefQueriesTable

    qry_ra: float = ra_field
    qry_dec: float = dec_field
    qry_filt: str = Field(min_length=1)
    compid: int = Field(ge=0)

    def exists(self) -> bool:
        """
        Checks if the pydantic-ified data exists the corresponding sql database

        :return: bool
        """
        return self._exists(
            values=[self.qry_ra, self.qry_dec, self.qry_filt],
            keys=["qry_ra", "qry_dec", "qry_filt"],
        )
