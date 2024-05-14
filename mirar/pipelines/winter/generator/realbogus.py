"""
Functions to apply rbscore
"""

import numpy as np
import pandas as pd
import torch
from torch import nn
from winterrb.utils import make_triplet


def apply_rb_to_table(model: nn.Module, table: pd.DataFrame) -> pd.DataFrame:
    """
    Apply the realbogus score to a table of sources

    :param model: Pytorch model
    :param table: Table of sources
    :return: Table of sources with realbogus score
    """

    rb_scores = []

    for _, row in table.iterrows():
        triplet = make_triplet(row, normalize=True)
        triplet_reshaped = np.transpose(np.expand_dims(triplet, axis=0), (0, 3, 1, 2))
        with torch.no_grad():
            outputs = model(torch.from_numpy(triplet_reshaped))

        rb_scores.append(float(outputs[0]))

    table["rb"] = rb_scores

    return table
