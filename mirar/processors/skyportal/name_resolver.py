"""
Processor to resolve target coordinates from a target name
"""

import logging
from typing import Optional

from mirar.data import ImageBatch
from mirar.paths import TARGET_KEY
from mirar.processors.base_processor import BaseImageProcessor
from mirar.processors.skyportal.client import SkyportalClient

logger = logging.getLogger(__name__)


class SkyportalNameResolver(BaseImageProcessor):
    """
    Processor to set target coordinates from the target name, using SkyPortal

    Headers generally record the telescope pointing, which is dithered and need
    not sit on the target at all, so the position is looked up by name instead.
    Images whose name cannot be resolved keep the coordinates they already have,
    since a missing position is no reason to abandon a night of data.
    """

    base_key = "skyportalnameresolver"

    def __init__(
        self,
        ra_key: str,
        dec_key: str,
        name_key: str = TARGET_KEY,
        skyportal_client: SkyportalClient = SkyportalClient(),
    ):
        super().__init__()
        self.ra_key = ra_key
        self.dec_key = dec_key
        self.name_key = name_key
        self.skyportal_client = skyportal_client
        self.resolved_positions: dict[str, Optional[tuple[float, float]]] = {}

    def description(self) -> str:
        return (
            f"Resolving '{self.name_key}' to '{self.ra_key}'&'{self.dec_key}' "
            f"via {self.skyportal_client.base_url}"
        )

    def resolve(self, name: str) -> Optional[tuple[float, float]]:
        """
        Get the position of one target by name, caching each lookup

        :param name: Name of the target, e.g. 'ZTF26abzrlze'
        :return: Right ascension and declination in degrees, or None if not found
        """
        if name not in self.resolved_positions:
            position = self.skyportal_client.get_source_position(name)

            if position is None:
                logger.warning(
                    f"Could not resolve '{name}', so its images keep the "
                    f"coordinates they already have"
                )

            self.resolved_positions[name] = position

        return self.resolved_positions[name]

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        for image in batch:
            position = self.resolve(image[self.name_key])

            if position is not None:
                image[self.ra_key], image[self.dec_key] = position

        return batch
