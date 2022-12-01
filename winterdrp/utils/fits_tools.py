def save_HDU_as_fits(HDU, path):
    HDU.writeto(path, overwrite=True)
