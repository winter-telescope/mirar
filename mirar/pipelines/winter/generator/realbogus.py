"""
Functions to apply rbscore
"""

import numpy as np
import pandas as pd
import torch
from torch import nn
from winterrb.tree import get_numpy_from_df
from winterrb.utils import make_triplet
from xgboost import XGBClassifier


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
            logits = model(torch.from_numpy(triplet_reshaped))
            outputs = torch.sigmoid(logits).cpu().numpy().ravel()

        rb_scores.append(float(outputs[0]))

    table["rb"] = rb_scores

    return table


def apply_xrb_to_table(clf: XGBClassifier, table: pd.DataFrame) -> pd.DataFrame:
    """
    Apply the xrealbogus score to a table of sources

    :param clf: xgboost model
    :param table: DataFrame of sources
    :return: Table of sources with 'xrealbogus' score
    """

    numpy_array = get_numpy_from_df(table)

    scores = clf.predict_proba(numpy_array).T[1]

    table["xrb"] = scores

    return table
