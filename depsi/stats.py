"""general estimation algorithms."""

import numpy as np
from scipy.stats import chi2


def blue(A, y, Qyy):
    """Calculate the Best Linear Unbiased Estimator (BLUE).

    Parameters
    ----------
    A : np.ndarray
        The A matrix
        Shape (m, n), where m is the number of observations and n is the number of unknown parameters.
    y : np.ndarray
        The vector of observations. Shape (m, 1).
    Qyy : np.ndarray
        The variance-covariance matrix of the observations. Shape (m, m).

    Returns
    -------
    x_hat : np.ndarray
        The vector of estimates for the unknown parameters. Shape (n, 1).
    Qx_hat : np.ndarray
        The variance-covariance matrix of the estimated unknown parameters. Shape (n, n).

    """
    x_hat = np.linalg.inv(A.T @ np.linalg.inv(Qyy) @ A) @ A.T @ np.linalg.inv(Qyy) @ y
    Qx_hat = np.linalg.inv(A.T @ np.linalg.inv(Qyy) @ A)
    return x_hat, Qx_hat


def blue_q_yy_inv(A, y, Qyy_inv):
    """Calculate the Best Linear Unbiased Estimator (BLUE) using the inverse of the variance-covariance matrix.

    This function calculates the Best Linear Unbiased Estimator (BLUE) for the unknown parameters using
    the inverse of the variance-covariance matrix of the observations.

    Parameters
    ----------
    A : np.ndarray
        The A matrix.
        Shape (m, n), where m is the number of observations and n is the number of unknown parameters.
    y : np.ndarray
        The vector of observations. Shape (m, 1).
    Qyy_inv : np.ndarray
        The inverse of the variance-covariance matrix of the observations. Shape (m, m).

    Returns
    -------
    x_hat : np.ndarray
        The vector of estimates for the unknown parameters. Shape (n, 1).
    Qx_hat : np.ndarray
        The variance-covariance matrix of the estimated unknown parameters. Shape (n, n).

    """
    try:
        Qx_hat = np.linalg.inv(A.T @ (Qyy_inv) @ A)
    except np.linalg.LinAlgError:
        Qx_hat = np.linalg.pinv(A.T @ (Qyy_inv) @ A)
        print("Warning: Normal inverse failed, therefore tried pseudo-inverse")

    x_hat = Qx_hat @ A.T @ Qyy_inv @ y
    return x_hat, Qx_hat


def overall_model_test(alpha, ehat, Qyy_inv, q, output):
    """Perform the overall model test and prints the outcome.

    This function computes the overall model test (OMT)
    It calculates the test statistic (Tq) and compares it with the critical value (k) from the Chi-squared distribution.
    If the statistic is lower than the critical value, the model is accepted; otherwise, it is rejected.

    Parameters
    ----------
    alpha : float
        The false alarm probability (significance level) for the Chi-squared test.
    ehat : numpy.ndarray
        The residuals from the network adjustment (the difference between observed and adjusted values).
    Qyy_inv : numpy.ndarray
        The inverse of the variance-covariance matrix (precision matrix) of the observations.
    q : int
        The degrees of freedom of the Chi-squared distribution, typically corresponding to the number of equations
        or observations minus the number of unknowns.
    output : int
        If set to 1, the function prints the result of the model test (whether it is accepted or rejected).

    Returns
    -------
    k : float
        The critical value (threshold) from the Chi-squared distribution for the given false alarm probability (alpha).
    Tq : float
        The test statistic calculated from the residuals and the inverse of the variance-covariance matrix.

    Example
    -------
    k, Tq = overall_model_test(alpha=0.05, ehat=ehat, Qyy_inv=Qyy_inv, q=q, output=1)
    """
    k = chi2.ppf(1 - alpha, q)
    Tq = (ehat.T @ Qyy_inv @ ehat).item()

    if output == 1:
        if Tq < k:
            print(f"(T = {Tq:.1f}) < (K = {k:.1f}), OMT is accepted.")
        else:
            print(f"(T = {Tq:.1f}) > (K = {k:.1f}), OMT is rejected.")

    return k, Tq
