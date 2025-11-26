"""utilities for statistical estimation."""

import logging

import numpy as np
from scipy.optimize import fsolve
from scipy.stats import chi2, ncx2, norm

logger = logging.getLogger(__name__)


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


def pretest(b, alpha0=0.001, gamma0=0.80):
    """Perform preparations for statistical hypothesis testing.

    This function calculates critical values and parameters needed for hypothesis testing,
    including the non-centrality parameter and significance levels for both one-dimensional
    and multi-dimensional tests. It is particularly useful for pre-testing when redundancy
    (degrees of freedom) is involved.

    Parameters
    ----------
    b : int
        Redundancy, representing the number of conditions or degrees of freedom.
    alpha0 : float, optional
        The level of significance for a 1-dimensional test (default is 0.001).
    gamma0 : float, optional
        The desired power of the test (default is 0.80).

    Returns
    -------
    lam0 : float
        The non-centrality parameter, calculated based on the desired test power and significance level.
    k1 : float
        The critical value for a 1-dimensional test (normal distribution).
    kb : float
        The critical value for a b-dimensional test (non-central chi-squared distribution).
        Returns NaN if redundancy b is zero.
    alphab : float
        The effective level of significance for the b-dimensional test.
        Returns NaN if redundancy b is zero.

    Notes
    -----
    - When b = 1, `k1` is derived from a chi-squared distribution and then converted to the normal domain.
    - When b > 0, the function calculates `kb` using the non-central chi-squared distribution with `b` degrees
      of freedom.
    - If b = 0, the function handles it gracefully by setting `kb` and `alphab` to NaN.
    - The function relies on `lambda0()` to compute the non-centrality parameter and uses the `chi2` and `ncx2`
      distributions from `scipy.stats`.
    """
    # Compute critical value for b=1 dof (chi-squared)
    k1 = chi2.ppf(1.0 - alpha0, 1)

    # Compute non-centrality parameter
    lam0 = lambda0(gamma0, 1, k1)

    # Compute critical value for b dof and level of significance for b-dimensional test
    if b > 0:
        kb = ncx2.ppf(1.0 - gamma0, b, lam0)
        alphab = 1.0 - chi2.cdf(kb, b)
    else:
        logger.warning("zero redundancy (b=0) encountered in pretest")
        kb = np.nan
        alphab = np.nan

    # Compute critical value for b=1 dof (normal)
    k1 = np.sqrt(k1)

    return lam0, k1, kb, alphab


def _lambda_approx(lambda_init, df, gam0, cv):
    """Approximation of the non-centrality parameter lambda.

    This function computes an approximation of the non-centrality parameter (λ)
    for hypothesis testing using an approach based on the FORTRAN routine CHILNC.F by F. Kenselaar.
    The method utilizes the normal distribution and follows Abramowitz and Stegun's formula 26.4.28.

    Parameters
    ----------
    lambda_init : float
        Approximate initial value of the non-centrality parameter λ.
    df : int
        Degrees of freedom for the test.
    gam0 : float
        Desired right-hand probability (test power).
    cv : float
        Critical value for the hypothesis test.

    Returns
    -------
    lam0ap : float
        The approximated non-centrality parameter λ.

    Notes
    -----
    - The function uses a normal approximation to the chi-squared distribution.
    - For small `xn` values, a safeguard is included to handle potential numerical stability issues.
    - The approximation involves transforming the critical value into the normal domain and computing the quantile `qn`.

    References
    ----------
    - Abramowitz, M., & Stegun, I. A. (1964). Handbook of Mathematical Functions.
      Dover Publications. (Formula 26.4.28)
    - Original FORTRAN routine: CHILNC.F by F. Kenselaar.
    - Script based on `_lambda_approx.m` function written in Matlab by Marcel Martens.
    - Translated to Python by Wietske Brouwer on 21-03-2024.
    """
    a = df + lambda_init
    b = lambda_init / a
    nn = a / (1.0 + b)
    q29n = 2.0 / (9.0 * nn)
    xn = ((cv / a) ** (1.0 / 3.0) - (1.0 - q29n)) / np.sqrt(q29n)

    if xn == 0.0:
        qn = 0.5
    else:
        x_squared_half = 0.5 * xn * xn
        qn = (1.0 - norm.cdf(x_squared_half)) / 2.0

    if xn <= 0.0:
        qn = 1.0 - qn

    lam0ap = qn - gam0

    return lam0ap


def _lam0_accurate(lambda_init, df, gam0, cv):
    """Compute the accurate value of the non-centrality parameter lambda.

    This function calculates the precise non-centrality parameter (λ) for statistical hypothesis testing.
    The approach is based on the original FORTRAN routine CHILNC.F by F. Kenselaar, and uses the
    non-central chi-squared distribution.

    Parameters
    ----------
    lambda_init : float
        Initial approximation of the non-centrality parameter λ.
    df : int
        Degrees of freedom for the test.
    gam0 : float
        Desired right-hand probability (test power).
    cv : float
        Critical value for the hypothesis test.

    Returns
    -------
    lam0ac : float
        The accurate non-centrality parameter λ.

    Notes
    -----
    - If the input `lambda_init` is negative, it is reset to 0 to maintain validity.
    - The function uses the cumulative distribution function (CDF) of the non-central chi-squared
      distribution to compute the difference between the observed and desired probabilities.

    References
    ----------
    - Original FORTRAN routine: CHILNC.F by F. Kenselaar.
    - Script based on `_lam0_accurate.m` function written in Matlab by Marcel Martens.
    - Translated to Python by Wietske Brouwer on 21-03-2024.
    """
    # If the input `lambda_init` is negative, reset to 0 to maintain validity.
    lambda_init = max(lambda_init, 0.0)
    return (1.0 - ncx2.cdf(cv, df, lambda_init)) - gam0


def lambda0(gam0, df, cv):
    """Compute of the non-centrality parameter lambda (λ).

    This function computes the non-centrality parameter λ for the non-central chi-squared
    distribution using a bisection iteration method. It is based on the original FORTRAN
    routine CHILNC.F by F. Kenselaar and relies on helper functions `_lam0_accurate` and `_lambda_approx`
    for approximation and accuracy.

    Parameters
    ----------
    gam0 : float
        Right-hand probability (test power), should be between 0 and 1.
    df : int
        Degrees of freedom for the test, must be >= 1.
    cv : float
        Critical value for the hypothesis test, must be >= 0.

    Returns
    -------
    lam0 : float
        The computed non-centrality parameter λ.

    Raises
    ------
    ValueError
        If degrees of freedom (`df`) is less than 1.
        If the critical value (`cv`) is less than 0.
        If `gam0` is not between 0 and 1.
        If the combination of `gam0` and `cv` leads to a negative λ.

    Notes
    -----
    - The function uses a bisection method to approximate the initial λ value and then refines
      it using a root-finding method (`fsolve`).
    - If the critical value (`cv`) is close to the value obtained from the chi-squared
      distribution (`cv1`), λ is set to 0 directly.
    - The maximum number of iterations (`MAX_ITER`) is set to 50, with an expansion factor (`EXPANSION_FACTOR`) of 1.6.

    References
    ----------
    - Original FORTRAN routine: CHILNC.F by F. Kenselaar.
    - Script based on `lambda0.m` function written in Matlab by Marcel Martens.
    - Translated to Python by Wietske Brouwer on 21-03-2024.
    """
    # Initializations
    MAX_ITER = 50
    EXPANSION_FACTOR = 1.6

    # Checks on input
    if df < 1:
        raise ValueError("LAMBDA0: degree of freedom smaller than 1.")
    if cv < 0:
        raise ValueError("LAMBDA0: critical value less than zero.")
    if gam0 <= 0 or gam0 > 1:
        raise ValueError("LAMBDA0: gamma0 not between 0 and 1")

    cv1 = chi2.ppf(1 - gam0, df)

    if cv < cv1:
        raise ValueError("LAMBDA0: this gamma0 and Critical Value lead to negative lambda")
    elif cv == cv1:
        lam0 = 0.0
    else:
        # Compute start value for lambda (lstrt)
        x1 = 0
        x2 = 5
        ii = 0
        f1 = _lambda_approx(x1, df, gam0, cv)
        f2 = _lambda_approx(x2, df, gam0, cv)
        step = x2 - x1

        while (ii <= MAX_ITER) and (f1 * f2 >= 0):
            ii += 1

            if abs(f1) < abs(f2):
                x1 -= step
                f1 = _lambda_approx(x1, df, gam0, cv)
            else:
                x2 += step
                f2 = _lambda_approx(x2, df, gam0, cv)
            step *= EXPANSION_FACTOR
        lstrt = (x1 + x2) / 2

        # Compute approximate value for lambda using an approximation formula for non-central chi-square distribution
        if cv - cv1 < 0.1:
            lambda_init = 0
        else:
            lambda_init = fsolve(lambda lambda_init: _lambda_approx(lambda_init, df, gam0, cv), lstrt)

        # Compute an accurate value for lambda
        lam0 = fsolve(lambda lambda_init: _lam0_accurate(lambda_init, df, gam0, cv), lambda_init)

        # If multiple solutions are found, raise an error
        if lam0.size > 1:
            raise ValueError("lambda0: multiple solutions found for lambda!")

    return lam0 if np.isscalar(lam0) else lam0[0]
