import numpy as np


# Median combine with scaled normalization of individual and final frames
def median_combine_normed(
        img_array: np.ma.masked_array
) -> np.ma.masked_array:
    """
    Median combine with normalization of individual and final frames
    Returns the combined image
    """
    norm_array = np.ma.zeros(img_array.shape)

    # Store normalized versions of each image
    for i in range(len(img_array[0, 0, :])):
        norm_array[:, :, i] = img_array[:, :, i] / np.ma.median(img_array[:, :, i])

        # median combination of normalized images
    combined_img = np.ma.median(norm_array, axis=2)

    # Normalize the median combination by itself
    combined_img = combined_img / np.ma.median(combined_img)

    return combined_img
