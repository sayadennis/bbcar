import sys
from typing import Union

import numpy as np
import pandas as pd
from sklearn import metrics


class PriorPredictor:
    """
    This is a baseline model class.
    """

    def __init__(self, p: float, random_state: int = 38):
        self.p = p
        self.random_state = random_state

    def predict(self, X: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
        """
        Create a random binary prediction results with a given prior.
        """
        np.random.seed(self.random_state)

        # generate predictions
        pred = np.random.choice([0, 1], size=X.shape[0], p=[1 - self.p, self.p])

        return pred

    def predict_proba(self, X: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
        """
        Create a random continuous prediction results with a given prior.
        """
        # create empty array
        scores = np.empty(X.shape[0])

        # generate predictions
        pred = self.predict(X)

        # generate probabilities
        for i, y in enumerate(pred):
            if y == 1:
                scores[i] = np.random.uniform(0.5, 1)
            else:
                scores[i] = np.random.uniform(0, 0.5)

        return scores


if __name__ == "__main__":
    # load input and target
    X = pd.read_csv(sys.argv[1], index_col=0)
    y = pd.read_csv(sys.argv[2], index_col=0)
    outfn = sys.argv[3]

    ol = list(set(X.index).intersection(y.index))
    X = X.loc[ol, :]
    y = y.loc[ol, :]

    if isinstance(y, pd.DataFrame):
        y = y.to_numpy()

    # with a random model, CV score is the same as score on whole set
    random_states = [38, 39, 40, 41, 42]
    pf = pd.DataFrame(
        index=[f"seed {i}" for i in random_states],
        columns=["ROC", "Precision", "Recall", "F1"],
    )

    # iterage over seeds and record performance
    for i in random_states:
        basemodel = PriorPredictor(
            p=np.sum(y) / y.shape[0], random_state=i
        )  # or p = 0.5
        y_pred = basemodel.predict(X)
        y_score = basemodel.predict_proba(X)
        pf.loc[f"seed {i}", "ROC"] = metrics.roc_auc_score(y, y_score)
        pf.loc[f"seed {i}", "Precision"] = metrics.precision_score(y, y_pred)
        pf.loc[f"seed {i}", "Recall"] = metrics.recall_score(y, y_pred)
        pf.loc[f"seed {i}", "F1"] = metrics.f1_score(y, y_pred)

    # save performance table
    pf.to_csv(outfn)
