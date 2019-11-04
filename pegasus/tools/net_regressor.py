import time
import numpy as np

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import check_array
from sklearn.neural_network import MLPRegressor
import logging

logger = logging.getLogger("pegasus")


class MaxStdScaler(BaseEstimator, TransformerMixin):
    def __init__(self, copy=True, factor=1.0):
        self.factor = float(factor)
        self.copy = copy

    def fit(self, X):
        X = check_array(X, copy=self.copy, estimator=self, dtype=np.float64)
        self.scaler = np.max(np.std(X, axis=0, ddof=1)) / self.factor
        return self

    def transform(self, X):
        X = check_array(X, copy=self.copy, estimator=self, dtype=np.float64)
        X /= self.scaler

        return X

    def inverse_transform(self, X, copy=None):
        if copy is None:
            copy = self.copy

        X = check_array(X, copy=copy, estimator=self, dtype=np.float64)
        X *= self.scaler

        return X


def net_train_and_predict(X_train, y_train, X_pred, alpha, random_state, verbose=False):
    start_time = time.perf_counter()

    scaler_x = MaxStdScaler()
    X_train = scaler_x.fit_transform(X_train)
    scaler_y = MaxStdScaler(factor=15.0)
    y_train = scaler_y.fit_transform(y_train)

    regressor = MLPRegressor(
        hidden_layer_sizes=(100, 75, 50, 25),
        activation="relu",
        solver="sgd",
        learning_rate="adaptive",
        alpha=alpha,
        random_state=random_state,
    )
    regressor.fit(X_train, y_train)
    logger.info(regressor.loss_)

    y_pred = scaler_y.inverse_transform(
        regressor.predict(scaler_x.transform(X_pred)), copy=False
    )

    end_time = time.perf_counter()

    if verbose:
        logger.info(
            "Deep regressor traning and predicting finished. Time spent = {:.2f}s.".format(
                end_time - start_time
            )
        )

    return y_pred
