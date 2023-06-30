"""
Module containing helper functions for autoastrometry
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)


def median(float_list: list[float]) -> float:
    """
    Calculate median value of input list

    :param float_list: list of floats
    :return: median
    """
    float_array = np.array(float_list)
    return float(np.median(float_array))


def stdev(float_list: list[float]) -> float:
    """
    Calculate standard deviation value of input list

    :param float_list: list of floats
    :return: std deviation
    """
    float_array = np.array(float_list)
    return float(np.std(float_array))


def mode(float_list: list[float]) -> float:
    """
    Calculate mode value of input list

    :param float_list: list of floats
    :return: mode
    """
    if len(float_list) == 0:
        err = "Float list is empty, cannot calculate mode."
        logger.error(err)
        raise ValueError(err)

    sorted_array = np.array(sorted(float_list))
    diff = sorted_array[1:] - sorted_array[:-1]
    n_diffs = len(diff)
    if n_diffs >= 32:
        step = n_diffs / 16
    elif n_diffs >= 6:
        step = 2
    else:
        step = 1

    min_mean = diff.sum()
    i_mean = n_diffs / 2

    for i in range(n_diffs):
        step_index = [int(max(i - step, 0)), int(min(i + step, n_diffs))]
        step_mean = diff[step_index[0] : step_index[1]].mean()
        if step_mean < min_mean:
            min_mean = step_mean
            i_mean = i

    list_mode = sorted_array[int(i_mean)]  # + s[i_mean+1])/2

    return list_mode


def ra_str_2_deg(ra_str: str) -> float:
    """
    Convert ra string to degrees

    :param ra_str: hexadecismal ra
    :return: ra (degrees)
    """
    ra_str = str(ra_str).strip()
    ra = ra_str.split(":")
    if len(ra) == 1:
        return float(ra_str)
    return 15 * (float(ra[0]) + float(ra[1]) / 60.0 + float(ra[2]) / 3600.0)


def dec_str_2_deg(dec_str: str) -> float:
    """
    Convert dec string to degrees

    :param dec_str: hexadecismal ra
    :return: dec (degrees)
    """
    dec_str = str(dec_str).strip()
    dec = dec_str.split(":")
    if len(dec) == 1:
        return float(dec_str)
    sign = 1
    if dec_str[0] == "-":
        sign = -1
    return sign * (abs(float(dec[0])) + float(dec[1]) / 60.0 + float(dec[2]) / 3600.0)


def unique(input_list: list) -> list:
    """
    Returns unique values in a list

    :param input_list: list to check
    :return: unique values
    """
    lis = sorted(input_list[:])  # make a copy
    llen = len(lis)
    i = 0
    while i < llen - 1:
        if lis[i + 1] == lis[i]:
            del lis[i + 1]
            llen = llen - 1
        else:
            i = i + 1
    return lis
