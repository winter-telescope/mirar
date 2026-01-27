from pathlib import Path
from mirar.data import Image
from mirar.io import open_fits, open_raw_image
from mirar.paths import core_fields
# TODO: Add more imports as required


def load_raw_spring_fits(path: str | Path) -> tuple[np.array, astropy.io.fits.Header]:
    """
    Function to load a raw GIT image

    :param path: path of file
    :return: data and header of image

    # TODO : Implement the function to read SPRING FITS files, and append header so that
    # TODO: it has all required keyword values. The required keywords are listed in the variable core_fields (see above imports for path).

    """
    data, header = open_fits(path)

    ## INSERT CODE HERE ##

    for field in core_fields:
        if field not in header:
            raise KeyError(f"Core field {field} not found in header")


def load_raw_spring_image(path: str | Path) -> Image:
    """
    Function to load a raw GIT image

    :param path: Path to the raw image
    :return: Image object
    """
    return open_raw_image(path, load_raw_spring_fits)
