from astropy.io import fits
# from winterdrp.pipelines import reformat_raw_data


# def open_fits(path):
#     img = fits.open(path)
#     img = reformat_raw_data(img)
#     return img


def create_fits(data, header, history):
    proc_hdu = fits.PrimaryHDU(data)
    proc_hdu.header = header  # Copy over the header from the raw file
    proc_hdu.header.add_history(history)  # Add a note to the header
    return proc_hdu
