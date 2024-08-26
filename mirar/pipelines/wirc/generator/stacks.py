"""
Module to group images based on the target coordinates into planned stack groups
"""

from astropy import coordinates as coords
from astropy import units as u
from astropy.coordinates import Angle

from mirar.data import ImageBatch
from mirar.paths import OBSCLASS_KEY
from mirar.processors.utils.image_selector import split_images_into_batches

MAX_RADIUS_DEG = 5.0 / 60.0  # WIRC is 10 arc minutes each side


def label_stack_id(batch: ImageBatch) -> ImageBatch:
    """
    Label the stack id of the images in the batch
    :param batch: Original batch of images
    :return: Labeled batch of images
    """

    ras = []
    decs = []

    for image in batch:

        if image[OBSCLASS_KEY] != "science":
            image["targnum"] = -1
            continue

        target_ra = Angle(image["CRVAL1"], unit="hourangle").degree
        target_dec = Angle(image["CRVAL2"], unit="degree").degree

        position = coords.SkyCoord(target_ra, target_dec, unit="deg")

        match = None

        if len(ras) > 0:
            idx, d2d, _ = position.match_to_catalog_sky(
                coords.SkyCoord(ra=ras, dec=decs, unit="deg")
            )

            mask = d2d < (MAX_RADIUS_DEG * u.deg)
            if mask:
                match = idx

        if match is None:
            ras.append(target_ra)
            decs.append(target_dec)
            image["targnum"] = int(len(ras) - 1)
        else:
            image["targnum"] = int(match)

    new_batches = split_images_into_batches(
        batch, ["targnum", "filter", "telfocus", "exptime"]
    )

    combined_batch = ImageBatch()

    for i, split_batch in enumerate(new_batches):
        label = f"stack{i}"
        for image in split_batch:
            image["stackid"] = label
            combined_batch.append(image)

    return combined_batch
