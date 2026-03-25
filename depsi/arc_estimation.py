"""arc estimation algorithms."""

from typing import Literal

import dask.array as da
import numpy as np
import xarray as xr
from scipy.optimize import curve_fit

import depsi.model_definition as md
import depsi.stats as est
from depsi.utils import get_distance, wrap_phase

# Constants
STOP_HEIGHT = 1e-4  # Stop search step for height [m]
STOP_VEL = 1e-7  # Stop search step for velocity [m/y]
MAX_COUNT = 10  # Maximum number of search iterations


def _compute_dd(sd_complex_i, sd_complex_j, sd_quality_i, sd_quality_j):
    """Compute the Double-Difference (DD) phase observation given complex time series of points i and j.

    Point i serves as the reference point and is subtracted from point j:
    - `sd_complex_conj_i = sd_complex_i.conj()`
    - `dd_arc = sd_complex_j * sd_complex_conj_i`

    The function also calculates the diagonal of the Variance-Covariance Matrix (VCM) of the DD phase
    (`Q_yy_diagonal_sigma`), which represents the standard deviation (sigma) of the DD phase quality.

    Parameters
    ----------
    sd_complex_i : np.ndarray
        Complex time series of the reference point i. Shape (n, ).
    sd_complex_j : np.ndarray
        Complex time series of point j. Shape (n, ).
    sd_quality_i : np.ndarray
        Quality values (sigma) of the single-difference (SD) phase for point i. Shape (n, ).
    sd_quality_j : np.ndarray
        Quality values (sigma) of the single-difference (SD) phase for point j. Shape (n, ).

    Returns
    -------
    tuple
        - dd_arc : np.ndarray
          Double-difference (DD) phase observations. Shape (n, ).
        - Q_yy_diagonal_sigma : np.ndarray
          Diagonal values of the Variance-Covariance Matrix (VCM) of the DD phases. Shape (n, ).

    Example
    -------
    >>> dd_arc, Q_yy_diagonal_sigma = _compute_dd(sd_complex_i, sd_complex_j,
                                                    sd_quality_i, sd_quality_j)
    """
    sd_complex_conj_i = sd_complex_i.conj()  # Compute the complex conjugate for the reference point i
    dd_arc = sd_complex_j * sd_complex_conj_i

    # Compute the diagonal of the VCM of the dd phases
    Q_yy_diagonal_sigma = np.sqrt((sd_quality_i) ** 2 + (sd_quality_j) ** 2)

    return dd_arc, Q_yy_diagonal_sigma


def _unwrap_phases_filter(filter_length, arc_dd, phase_arc, jump):
    """Filter and unwrap double-difference (DD) arc phase observations in the complex domain.

    This function filters the real and imaginary components of the DD arc observations, computes the filtered phase
    (angle), and unwraps the phase based on the filtered function. It simplifies the detection of 2π jumps and applies
    corrections to produce an unwrapped phase time series.

    Parameters
    ----------
    filter_length : int
        Length of the moving average filter used to smooth the real and imaginary components of the arc observation.
    arc_dd : np.ndarray
        Double-difference arc observations in the complex domain. Shape (n, ).
    phase_arc : np.ndarray
        Wrapped phase observations for the arc. Shape (n, ).
    jump : float
        Threshold as a fraction of 2π used to detect and correct phase jumps during unwrapping.

    Returns
    -------
    tuple
        - phase_arc_unwrap : np.ndarray
          Unwrapped phase observations for the arc. Shape (n, ).
        - pi_diff : np.ndarray
          Integer multiple of 2π differences between the unwrapped filtered phase and the original phase. Shape (n, ).
        - filtered_phase_wrap : np.ndarray
          Filtered and wrapped phase observations. Shape (n, ).
        - filter_real : np.ndarray
          Filtered real component of the arc observations. Shape (n, ).
        - filter_imag : np.ndarray
          Filtered imaginary component of the arc observations. Shape (n, ).

    Notes
    -----
    - The function applies a moving average filter to smooth the real and imaginary parts of the input arc observations.
    - The filtered phase is used to detect 2π phase jumps, and a correction is applied to unwrap the phase.
    - Corrections are based on the difference between the original phase and the unwrapped filtered phase.
    """
    ## Filter in Re and Im domain to detect ambiguity levels easily
    filter_real = np.convolve(arc_dd.real, np.ones(filter_length) / filter_length, mode="same")
    filter_imag = np.convolve(arc_dd.imag, np.ones(filter_length) / filter_length, mode="same")
    complex_filtered = filter_real + filter_imag * 1j
    filtered_phase_wrap = np.angle(complex_filtered)

    # Detect 2pi jumps in the filtered double difference and calculate shift of time series
    abs_diff = np.zeros(len(filtered_phase_wrap) - 1)
    for k in range(len(filtered_phase_wrap) - 1):
        abs_diff[k] = filtered_phase_wrap[k + 1] - filtered_phase_wrap[k]

    shift = np.zeros(len(arc_dd))
    for k in range(len(abs_diff) - 1):
        if abs_diff[k] < -jump * 2 * np.pi or abs_diff[k] > jump * 2 * np.pi:
            if abs_diff[k] < -jump * 2 * np.pi:
                shift[k + 1 :] = shift[k + 1 :] + 2 * np.pi
            if abs_diff[k] > jump * 2 * np.pi:
                shift[k + 1 :] = shift[k + 1 :] - 2 * np.pi
        else:
            shift[k + 1 : 0] = 0

    # Unwrap the filtered phase function
    filtered_phase_unwrap = filtered_phase_wrap + shift

    # Unwrap the DD phase observation for the arc based on the unwraped filtered function
    # Correct for integer pi values differences between the filtered unwraped function and the DD observations
    pi_diff = np.around((phase_arc - filtered_phase_unwrap) / (2 * np.pi), 0)
    phase_arc_unwrap = phase_arc - 2 * np.pi * pi_diff

    return phase_arc_unwrap, pi_diff, filtered_phase_wrap, filter_real, filter_imag


def _scipy_fit_partition_2nd_order_bounds(breakpoints, x_data, arc_obs, initial_guess, bounds, vcm, n_max_iter):
    """Estimate the parameters for an arc using a partitioned second-order polynomial fit.

    This function splits the time series of the arc observations into multiple partitions at the specified breakpoints.
    For each partition, a 2nd order polynomial is fit, and the function ensures that the time series is continuous
    at the breakpoints. The fitting process uses the `curve_fit` function with specified bounds and variance-covariance
    matrix (vcm).

    Parameters
    ----------
    breakpoints : list of int
        Indices where the time series is divided into partitions.
    x_data : np.ndarray
        The input data for the model, typically including time-related variables. Shape (m, n).
    arc_obs : np.ndarray
        The observed arc, including both real and imaginary parts. Shape (m,).
    initial_guess : np.ndarray
        Initial guess for the unknown parameters. Shape (n,).
    bounds : tuple of (lower_bounds, upper_bounds)
        The bounds for the parameters during fitting. Each bound is an array of length n.
    vcm : np.ndarray
        The variance-covariance matrix of the observations. Shape (m, m).
    n_max_iter : np.ndarray
        The maximum nr of iterations for non-linear lsq per arc

    Returns
    -------
    estimated_params : np.ndarray
        The estimated parameters after fitting. Shape (n,).
    pcov : np.ndarray
        The covariance matrix of the estimated parameters. Shape (n, n).

    """

    def _model_arc_2nd_order(x_data, *model_params, bkps=breakpoints):
        """Model for an arc using a second-order polynomial for each partition in the time series.

        The arc time series is divided into partitions at the specified breakpoints, and for each partition, a 2nd
        polynomial is fit. The model accounts for amplitude variations, displacement, cross range, and temperature.

        This function needs to be defined inside _scipy_fit_partition_2nd_order_bounds since the model
        uses breakpoints.

        Parameters
        ----------
        x_data : np.ndarray
            The input data for the model, including time (t), temperature (T), and CR (cr2ph). Shape (m, n).
        model_params : list of float
            The parameters for the model. These include amplitude, displacement model parameters, and cross range
            and temperature. The exact number and order of parameters depend on the number of breakpoints.
        bkps : list of int, optional
            Breakpoints at which the time series is divided into partitions. Default is breakpoints.

        Returns
        -------
        np.ndarray
            The modeled arc, which includes both the real and imaginary parts of the arc observations.
            Shape (2 * m,).
        """
        t, temp, cr2ph = x_data

        # Define the amplitude per partition
        aa = model_params[0 : len(bkps)]

        # Define the parameters for the displacement model (third order polynomial)
        intercept = model_params[len(bkps)]
        p1 = model_params[len(bkps) + 1 : 2 * len(bkps) + 1]
        p2 = model_params[2 * len(bkps) + 1 : 3 * len(bkps) + 1]

        # Define parameters for the cross range and temperature
        height = model_params[-2]
        expansion = model_params[-1]

        # Define the displacement phase and ampltiudes (they vary per partition)
        displ = np.zeros(len(t))
        ampl = np.zeros(len(t))

        # c is the start of a new partition
        c = 0
        for i in range(len(bkps)):
            # Define till what index the function should go (which is the end of the partition)
            idx = int(bkps[i]) + 1

            # Define the displacement values
            displ[c:idx] = intercept - (p1[i] * t[c] + p2[i] * t[c] ** 2) + (p1[i] * t[c:idx] + p2[i] * t[c:idx] ** 2)
            ampl[c:idx] = aa[i]

            # Define the 'intercept' of the new partition (that is the end of the next partition)
            intercept = displ[idx - 1]

            # c is the starting point of a new partition
            c = idx - 1

        # Define the Real complex observation
        real_part = ampl * np.cos(height * cr2ph + expansion * temp + displ)
        imag_part = ampl * np.sin(height * cr2ph + expansion * temp + displ)
        return np.append(real_part, imag_part)

    estimated_params, pcov, infodict, _, _ = curve_fit(
        f=_model_arc_2nd_order,
        xdata=x_data,
        ydata=arc_obs,
        p0=initial_guess,
        bounds=bounds,
        sigma=vcm,
        absolute_sigma=True,
        full_output=True,
        max_nfev=n_max_iter,
    )

    return estimated_params, pcov


def _scipy_fit_partition_2nd_order_bounds_derivative(
    breakpoints, x_data, arc_obs, initial_guess, bounds, vcm, n_max_iter
):
    """Estimate the parameters for an arc using a partitioned second-order polynomial fit.

    This function splits the time series of the arc observations into multiple partitions at the specified breakpoints.
    For each partition, a 2nd order polynomial is fit, and the function ensures that the time series is continuous
    at the breakpoints. The fitting process uses the `curve_fit` function with specified bounds and variance-covariance
    matrix (vcm).

    Parameters
    ----------
    breakpoints : list of int
        Indices where the time series is divided into partitions.
    x_data : np.ndarray
        The input data for the model, typically including time-related variables. Shape (m, n).
    arc_obs : np.ndarray
        The observed arc, including both real and imaginary parts. Shape (m,).
    initial_guess : np.ndarray
        Initial guess for the unknown parameters. Shape (n,).
    bounds : tuple of (lower_bounds, upper_bounds)
        The bounds for the parameters during fitting. Each bound is an array of length n.
    vcm : np.ndarray
        The variance-covariance matrix of the observations. Shape (m, m).
    n_max_iter : np.ndarray
        The maximum nr of iterations for non-linear lsq per arc

    Returns
    -------
    estimated_params : np.ndarray
        The estimated parameters after fitting. Shape (n,).
    pcov : np.ndarray
        The covariance matrix of the estimated parameters. Shape (n, n).

    """

    def _model_arc_2nd_order_derivative(x_data, *model_params, bkps=breakpoints):
        """Model for an arc using a second-order polynomial for each partition in the time series.

        The arc time series is divided into partitions at the specified breakpoints, and for each partition, a 2nd
        polynomial is fit. The model accounts for amplitude variations, displacement, cross range, and temperature.

        This function needs to be defined inside _scipy_fit_partition_2nd_order_bounds since the model
        uses breakpoints.

        Parameters
        ----------
        x_data : np.ndarray
            The input data for the model, including time (`t`), temperature (`T`), and CR (`cr2ph`). Shape (m, n).
        model_params : list of float
            The parameters for the model. These include amplitude, displacement model parameters, and cross range
            and temperature. The exact number and order of parameters depend on the number of breakpoints.
        bkps : list of int, optional
            Breakpoints at which the time series is divided into partitions. Default is `breakpoints`.

        Returns
        -------
        np.ndarray
            The modeled arc, which includes both the real and imaginary parts of the arc observations.
            Shape (2 * m,).
        """
        t, temp, cr2ph = x_data

        nr_bkps = len(bkps)

        # Define the amplitude per partition
        aa = model_params[0:nr_bkps]

        # Define the parameters for the displacement model (third order polynomial)
        intercept = model_params[nr_bkps]
        p1 = np.array(model_params[nr_bkps + 1 : 2 * nr_bkps + 1])
        p2 = np.array(model_params[2 * nr_bkps + 1 : 3 * nr_bkps + 1])

        # Define parameters for the cross range and temperature
        height = model_params[-2]
        expansion = model_params[-1]

        # Define the displacement phase and ampltiudes (they vary per partition)
        displ = np.zeros(len(t))
        ampl = np.zeros(len(t))

        # c is the start of a new partition
        c = 0
        prev_slope = 0

        for i in range(nr_bkps):
            # Define till what index the function should go (which is the end of the partition)
            idx = int(bkps[i]) + 1
            if i > 0:
                p1[i] = prev_slope  # Zorg dat de eerste afgeleide overeenkomt met de vorige
                intercept = displ[c]

            # Define the displacement values
            displ[c:idx] = intercept - (p1[i] * t[c] + p2[i] * t[c] ** 2) + (p1[i] * t[c:idx] + p2[i] * t[c:idx] ** 2)
            ampl[c:idx] = aa[i]

            # Define the 'intercept' of the new partition (that is the end of the next partition)
            intercept = displ[idx - 1]

            prev_slope = p1[i] + 2 * p2[i] * t[idx - 1]

            # c is the starting point of a new partition
            c = idx - 1

        # Define the Real complex observation
        real_part = ampl * np.cos(height * cr2ph + expansion * temp + displ)
        imag_part = ampl * np.sin(height * cr2ph + expansion * temp + displ)
        return np.append(real_part, imag_part)

    estimated_params, pcov, infodict, _, _ = curve_fit(
        f=_model_arc_2nd_order_derivative,
        xdata=x_data,
        ydata=arc_obs,
        p0=initial_guess,
        bounds=bounds,
        sigma=vcm,
        absolute_sigma=True,
        full_output=True,
        max_nfev=n_max_iter,
    )

    return estimated_params, pcov


def _model_arc_partitions_2nd_order_phases(x_data, model_params):
    """Calculate the forward model of the total, displacement, cross range, and thermal phase for an arc.

    The function divides the time series into different partitions. For each partition, a 2nd order polynomial is used.

    Parameters
    ----------
    x_data : tuple
        Contains the following elements:
        - bkps (list of int): Breakpoints that divide the time series into partitions.
        - t (np.ndarray): Time data for the arc.
        - temp (np.ndarray): Temperature data for the arc.
        - cr2ph (np.ndarray): Cross-range data for the arc.
    model_params : list of float
        Model parameters for the arc:
        - Amplitudes per partition (A).
        - Displacement model parameters (intercept, p1, p2).
        - Parameters for cross-range (H) and temperature (expansion).

    Returns
    -------
    tuple
        Contains the following elements:
        - phase_total (np.ndarray): The total phase (displacement + cross range + thermal).
        - phase_thermal (np.ndarray): The thermal phase.
        - phase_cross_range (np.ndarray): The cross-range phase.
        - phase_displacement (np.ndarray): The displacement phase.
        - real_part (np.ndarray): The real part of the arc observation.
        - imag_part (np.ndarray): The imaginary part of the arc observation.
    """
    bkps, t, temp, cr2ph = x_data

    # Define the amplitude per partition
    aa = model_params[0 : len(bkps)]

    # Define the parameters for the displacement model (third order polynomial)
    intercept = model_params[len(bkps)]
    p1 = model_params[len(bkps) + 1 : 2 * len(bkps) + 1]
    p2 = model_params[2 * len(bkps) + 1 : 3 * len(bkps) + 1]

    # Define parameters for the cross range and temperature
    height = model_params[-2]
    expansion = model_params[-1]

    # Define the displacement phase and ampltiudes (they vary per partition)
    displ = np.zeros(len(t))
    ampl = np.zeros(len(t))

    # c is the start of a new partition
    c = 0
    for i in range(len(bkps)):
        # Define till what index the function should go (which is the end of the partition)
        idx = int(bkps[i]) + 1

        # Define the displacement values
        displ[c:idx] = intercept - (p1[i] * t[c] + p2[i] * t[c] ** 2) + (p1[i] * t[c:idx] + p2[i] * t[c:idx] ** 2)
        ampl[c:idx] = aa[i]

        # Define the 'intercept' of the new partition (that is the end of the next partition)
        intercept = displ[idx - 1]

        # c is the starting point of a new partition
        c = idx - 1

    # Define the Real complex observation
    real_part = ampl * np.cos(height * cr2ph + expansion * temp + displ)
    imag_part = ampl * np.sin(height * cr2ph + expansion * temp + displ)

    # Define the phases
    phase_thermal = expansion * temp
    phase_cross_range = height * cr2ph
    phase_displacement = displ
    phase_total = phase_thermal + phase_cross_range + phase_displacement

    return phase_total, phase_thermal, phase_cross_range, phase_displacement, real_part, imag_part


def _unwrap_phases(observed_phase, estimated_phase):
    """Unwrap the observed phases based on modeled/estimated phase.

    This function corrects phase jumps by adjusting the observed phase values with respect to the estimated phase.

    Parameters
    ----------
    observed_phase : np.ndarray
        The observed phase values (typically in radians).
    estimated_phase : np.ndarray
        The estimated or modeled phase values based on the model parameters.

    Returns
    -------
    np.ndarray
        The unwrapped phase values, which are corrected to account for phase wrapping.
    """
    pi_diff = np.around((estimated_phase - observed_phase) / (2 * np.pi), 0)
    phase_unwrap = observed_phase + 2 * np.pi * pi_diff

    return phase_unwrap


def _compute_residuals_per_partition_stm(y_arc, y_est, Q_dd, bkps):
    """Compute the rmse and std of the residuals per partition.

    This function can be used to compare values with a predefined Q matrix.

    Parameters
    ----------
    y_arc : np.ndarray
        The observed (real) arc values (should be in a flattened array).
    y_est : np.ndarray
        The estimated arc values (should be in a flattened array).
    Q_dd : np.ndarray
        The variance-covariance matrix (Q_dd) of the residuals.
    bkps : list or np.ndarray
        Breakpoints indicating where the arc is divided into different partitions.

    Returns
    -------
    tuple
        - rmse_partition (list): List of RMSE values computed per partition.
        - std_est_partition (list): List of standard deviations of the residuals per partition.
        - q_per_partition (list): List of standard deviations from the Q_dd matrix for each partition.
    """
    start = 0
    rmse_partition = []
    std_est_partition = []
    q_per_partition = []

    y_arc = y_arc.flatten()
    y_est = y_est.flatten()

    # Loop over the partitions to comput residues per partition
    # and add the values to lists
    for i in range(len(bkps)):
        n = bkps[i] - start

        rmse_s = np.sqrt((np.sum((y_arc[start : bkps[i]] - y_est[start : bkps[i]]) ** 2)) / n)
        std_s = np.sqrt(np.var(y_arc[start : bkps[i]] - y_est[start : bkps[i]]))
        q_per_s = np.sqrt(Q_dd[start, start])  # the apriori defined quality per partition

        rmse_partition.append(rmse_s)
        std_est_partition.append(std_s)
        q_per_partition.append(q_per_s)

        start = bkps[i]

    return rmse_partition, std_est_partition, q_per_partition


def _flatten_arrays_in_dict(dictionary):
    """Flatten arrays in a dictionary.

    Function is required and used in arc_estimation_functions

    Args:
        dictionary (dict): dictionary with arrays

    Returns:
        dict: dictionary with flattend arrays
    """
    dictionary = {key: np.array(value) for key, value in dictionary.items()}

    for key, value in dictionary.items():
        if isinstance(value, np.ndarray):
            if value.ndim > 1:  # Only  flatten the array if the dimension is larger than 1
                dictionary[key] = value.ravel()  # get 1D
    return dictionary


def arc_estimation_xarray_input(
    stm_pnt_i,
    stm_pnt_j,
    bounds,
    m2ph,
    n_max_iter,
    partition_quality_label: str,
    x_crd_label: str = "rd_x",
    y_crd_label: str = "rd_y",
    coordinate_type: Literal["euclidean", "geographic"] = "euclidean",
    filter_length_complex=30,
    jump_percentage_2pi=0.85,
    vcm_complex_method="mad_median",
    test_stochastics=False,
    print_output=False,
):
    """Estimate parameters for the arc defined between the connection point j and control point i.

    This function performs a series of computations for and arc,
    including variance-covariance matrix computation, double-difference phase estimation, and parameter fitting
    in both the phase and complex domains. The results of these calculations are stored in structured dictionaries.

    Parameters
    ----------
    stm_pnt_i : Xarray.Dataset
        Input space time matrix for the reference point i
    stm_pnt_j : Xarray.Dataset
        Input space time matrix for connection point j
    bounds : tuple of lists
        Bounds for parameter estimation in the format (lower_bounds, upper_bounds).
    m2ph : float
        Conversion factor from meters to phase.
    n_max_iter : np.ndarray
        The maximum nr of iterations for non-linear lsq per arc
    partition_quality_label: str
        Layer name in the STM of the SLC quality
    x_crd_label: str, default "rd_x"
        Label of the x-coordinate in the STMs (for geographic, this is 'lon')
    y_crd_label: str, default "rd_y"
        Label of the y-coordinate in the STMs (for geographic, this is 'lat')
    coordinate_type: Literal["euclidean", "geographic"], default "euclidean"
        Whether to compute distances in Euclidean space (for RD) or geographic distance (for lat/lon)
    filter_length_complex : int, optional
        Length of the filter for phase unwrapping (default: 30).
    jump_percentage_2pi : float, optional
        Threshold for unwrapping phase jumps in terms of 2π (default: 0.85).
    vcm_complex_method : str, optional
        Method for variance-covariance matrix estimation in the complex domain.
        Options are "sigma_mean" or "mad_median" (default: "mad_median").
    test_stochastics : bool, optional
        Flag for performing stochastic testing (default: False).
    print_output : bool, optional
        Flag for enabling or disabling print statements (default: False).

    Returns
    -------
    results : dict
        Dictionary containing results for each arc, including:
            - 'unwrap_phases_arc': Unwrapped phases for each arc.
            - 'sigma_phases_arc': Phase variances for each arc.
            - 'estimated_phase': Estimated phases for each arc.
            - 'estimated_displ_phase': Displacement-related phases for each arc.
            - 'estimated_thermal': Estimated thermal expansion coefficients.
            - 'estimated_cross_range': Estimated cross-range components.
            - 'estimated_cross_range_sigma': Uncertainties of cross-range estimates.
            - 'estimated_thermal_sigma': Uncertainties of thermal estimates.
            - 'estimated_thermal_phase': Thermal-related phases for each arc.
            - 'estimated_cross_range_phase': Cross-range related phases for each arc.
            - 'cr2ph_arcs': Cross-range-to-phase conversion factors for each arc.
            - 'succeeded_arcs': List of arcs where parameter estimation succeeded.

    stochastic_results : dict, optional
        Dictionary containing stochastic testing results (if `test_stochastics=True`), including:
            - 'q_per_partition': Quality metrics for each partition.
            - 'std_residuals_partition': Standard deviations of residuals for each partition.
            - 'rmse_residuals_partition': RMSE of residuals for each partition.
            - 'mean_sigma_post_arc': Mean post-fit sigma values for each arc.
            - 'mean_a_priori_sigma_arc': Mean a priori sigma values for each arc.
            - 'arc_length': Lengths of the arcs.
            - 'mean_sigma_p_i': Mean quality metrics for the first point in each arc.
            - 'mean_sigma_p_j': Mean quality metrics for the second point in each arc.

    Notes
    -----
    1. The function uses deterministic assignment for the CR component, setting it to zero for one of the points.
    2. The estimation process includes fallback mechanisms to handle cases where optimal parameters cannot be found.
    3. Requires external utility functions for phase unwrapping, functional model construction, and lsq estimation.

    Raises
    ------
    ValueError
        If an unknown `vcm_complex_method` is specified.
    RuntimeError, ValueError
        If parameter estimation fails for an arc during optimization.
    """
    # If we want to do some tests on the stochastics
    if test_stochastics:
        stochastic_results = {
            "q_per_partition": [],
            "std_residuals_partition": [],
            "rmse_residuals_partition": [],
            "mean_sigma_post_arc": [],
            "mean_a_priori_sigma_arc": [],
            "arc_length": [],
            "mean_sigma_p_i": [],
            "mean_sigma_p_j": [],
        }

    # Dictionary to store results for one arc
    results = {
        "unwrap_phases_arc": [],
        "sigma_phases_arc": [],
        "estimated_phase": [],
        "estimated_displ_phase": [],
        "estimated_thermal": [],
        "estimated_cross_range": [],
        "estimated_cross_range_sigma": [],
        "estimated_thermal_sigma": [],
        "estimated_thermal_phase": [],
        "estimated_cross_range_phase": [],
        "cr2ph_arcs": [],
        "succeeded_arcs": [],
    }

    print(f"idx pnt i: {int(stm_pnt_i['space'].values)}")
    print(f"idx pnt j: {int(stm_pnt_j['space'].values)}")

    dates = stm_pnt_i["time"].values
    Btemporal = stm_pnt_i["years_since_first_img"].values
    temp = stm_pnt_i["temperature"].values

    # Extract information of the two points of the arc
    pnt_i_idx = int(stm_pnt_i["space"].values)
    sd_complex_i = stm_pnt_i["sd_complex"].values
    slc_quality_i = stm_pnt_i[partition_quality_label].values
    bkps_stm_i = stm_pnt_i["breakpoints"].values
    sigma_ampl_sd_i = stm_pnt_i["partition_sd_amplitude_sigma"].values
    mean_ampl_sd_i = stm_pnt_i["partition_sd_amplitude_mean"].values
    mad_ampl_sd_i = stm_pnt_i["partition_sd_mad"].values
    median_ampl_sd_i = stm_pnt_i["partition_sd_amplitude_median"].values

    pnt_j_idx = int(stm_pnt_j["space"].values)
    sd_complex_j = stm_pnt_j["sd_complex"].values
    slc_quality_j = stm_pnt_j[partition_quality_label].values
    cr2ph_j = stm_pnt_j["sd_cr2ph"].values
    bkps_stm_j = stm_pnt_j["breakpoints"].values
    sigma_ampl_sd_j = stm_pnt_j["partition_sd_amplitude_sigma"].values
    mean_ampl_sd_j = stm_pnt_j["partition_sd_amplitude_mean"].values
    mad_ampl_sd_j = stm_pnt_j["partition_sd_mad"].values
    median_ampl_sd_j = stm_pnt_j["partition_sd_amplitude_median"].values

    # Step 1: extract information on the arc
    # Compute the arc length
    arc_length = get_distance(
        [stm_pnt_i[x_crd_label], stm_pnt_i[y_crd_label]],
        [stm_pnt_j[x_crd_label], stm_pnt_j[y_crd_label]],
        mode=coordinate_type,
    )

    # Combine breakpoints to have breakpoints per arc
    bkps_arc_stm = bkps_stm_i + bkps_stm_j
    # Define the indexes of the breakpoints for the arc
    bkps = [index for index, value in enumerate(bkps_arc_stm) if value > 0]
    bkps.append(len(dates) - 1)

    # Compute the cr2ph for the arc, equals to point j since we determinsitcally set the value for point i to zero
    cr2ph_arc = cr2ph_j

    # Step 2. Compute the DD phases for the arc
    # point i is the reference point and is subtracted from point j:
    # Note that the output of Qyy_diagonal are actually sigmas and NO variances. Therefore we need to square the values
    dd_arc, Qyy_diagonal = _compute_dd(sd_complex_i, sd_complex_j, slc_quality_i, slc_quality_j)

    # Compute the variance covariance matrix of the DD based on the NMAD for the arc
    Qyy = np.identity(len(dates)) * Qyy_diagonal**2
    Qyy_inv = np.linalg.inv(Qyy)

    # Step 3. Estimate parameters in the phase domain
    # This step is required to get proper intial estimates for the parameter estimation in the complex domain
    # Unwrap the phases based on the filtered real and imaginary part
    phase_arc_unwrap, _, _, _, _ = _unwrap_phases_filter(
        filter_length_complex, dd_arc, np.angle(dd_arc), jump_percentage_2pi
    )

    # Compute different columns for the A matrices and construct to one A matrix
    A_cr = md.a_cross_range(cr2ph_arc)
    A_off = md.a_offset(len(Btemporal))
    A_lin = md.a_velocity(Btemporal, m2ph)
    A_temp = md.a_temperature(temp)
    A_arc = np.column_stack((A_cr, A_temp, A_off, A_lin))

    # Define the observation vector for the arc, which is based on the 'unwrapped' phase based on the filter
    y_arc = np.reshape(phase_arc_unwrap, (len(phase_arc_unwrap), 1))

    # Estimate parameters in the phase domain
    x_hat_arc_ph, Qx_hat_arc_ph = est.blue_q_yy_inv(A_arc, y_arc, Qyy_inv)

    # Step 4. Create VCM in the complex domain
    # Here we will compute the VCM for the complex domain.
    # It is possible to choose between the mean and sigma or mad and median amplitude per partition.

    # Estimate sigma of the DD phases
    if vcm_complex_method == "sigma_mean":
        sigma_dd = np.abs(mean_ampl_sd_i * mean_ampl_sd_j) * np.sqrt(
            (sigma_ampl_sd_i / mean_ampl_sd_i) ** 2 + (sigma_ampl_sd_j / mean_ampl_sd_j) ** 2
        )

    elif vcm_complex_method == "mad_median":
        sigma_dd = np.abs(median_ampl_sd_i * median_ampl_sd_j) * np.sqrt(
            (mad_ampl_sd_i * 1.4826 / median_ampl_sd_i) ** 2 + (mad_ampl_sd_j * 1.4826 / median_ampl_sd_j) ** 2
        )
    else:
        raise ValueError(
            f"You specified an unknown vcm complex method. The method -- {vcm_complex_method} -- does not exist"
        )

    # Compute VCM in the complex domain
    sigma_complex = np.append(
        sigma_dd, sigma_dd
    )  # Real and Imag are stacked together since we use both of the observations
    Q_dd_cmplx = np.identity(len(sigma_complex))
    np.fill_diagonal(Q_dd_cmplx, sigma_complex**2)

    # Step 5. Parameter estimation in the complex domain
    # Complex data preparation for the arc
    re_arc = dd_arc.real
    im_arc = dd_arc.imag
    arc_obs = np.append(re_arc, im_arc)

    # Combine all the independent variables in one independent variable
    x_data = (bkps, Btemporal, temp, cr2ph_arc)  # used in phase estimation
    X_data = (Btemporal, temp, cr2ph_arc)  # used in curve fit

    # Create arrays with initial values
    x0_2_p = np.zeros(3 * len(bkps) + 3)  # Create empty array for the bounds
    x0_2_p[0 : len(bkps)] = np.ones(len(bkps)) * np.max(re_arc)  # The amplitude to be estimated
    x0_2_p[len(bkps) + 1] = x_hat_arc_ph[2, 0]  # Interception of the dispalcement polynomial
    x0_2_p[len(bkps) + 1 : 2 * len(bkps) + 1] = (
        np.ones(len(bkps)) * x_hat_arc_ph[3, 0]
    )  # Value related displacement velocity
    x0_2_p[2 * len(bkps) + 1 : 3 * len(bkps) + 1] = np.zeros(
        len(bkps)
    )  # Value related to the second compont of polynomial
    x0_2_p[-2] = x_hat_arc_ph[0, 0]  # Cross range
    x0_2_p[-1] = x_hat_arc_ph[1, 0]  # Thermal expansion

    # Define the bounds
    (
        A_lower,
        a_lower,
        b_lower,
        c_lower,
        CR_lower,
        exp_lower,
        A_upper,
        a_upper,
        b_upper,
        c_upper,
        CR_upper,
        exp_upper,
    ) = bounds

    # define bounds for second order polynomial with partitions
    bounds_upper_2_p = np.ones(len(bkps) * 3 + 3)
    bounds_upper_2_p[0 : len(bkps)] = A_upper
    bounds_upper_2_p[len(bkps)] = a_upper
    bounds_upper_2_p[len(bkps) + 1 : 2 * len(bkps) + 1] = np.ones(len(bkps)) * b_upper
    bounds_upper_2_p[2 * len(bkps) + 1 : 3 * len(bkps) + 1] = np.ones(len(bkps)) * c_upper
    bounds_upper_2_p[-2] = CR_upper
    bounds_upper_2_p[-1] = exp_upper

    bounds_lower_2_p = np.ones(len(bkps) * 3 + 3)
    bounds_lower_2_p[0 : len(bkps)] = A_lower
    bounds_lower_2_p[len(bkps)] = a_lower
    bounds_lower_2_p[len(bkps) + 1 : 2 * len(bkps) + 1] = np.ones(len(bkps)) * b_lower
    bounds_lower_2_p[2 * len(bkps) + 1 : 3 * len(bkps) + 1] = np.ones(len(bkps)) * c_lower
    bounds_lower_2_p[-2] = CR_lower
    bounds_lower_2_p[-1] = exp_lower

    bounds_2_p = (list(bounds_lower_2_p), list(bounds_upper_2_p))

    # Curvefit with 2nd order displacement polynomial, partitions, and bounds
    # Estimate the unknown parameters:
    try:
        x_hat_2_p_b, pcov_2_p_b = _scipy_fit_partition_2nd_order_bounds(
            bkps,
            X_data,
            arc_obs,
            x0_2_p,
            bounds_2_p,
            Q_dd_cmplx,
            n_max_iter,
        )

    except (RuntimeError, ValueError) as e:
        print(f"Optimal parameters not found. Skipping arc {(pnt_i_idx, pnt_j_idx)}")
        print(f"Encountered error: {e}")

        # Fill everything with nans
        ts_length = len(Btemporal)

        for key in [
            "unwrap_phases_arc",
            "sigma_phases_arc",
            "estimated_phase",
            "estimated_displ_phase",
            "estimated_thermal_phase",
            "estimated_cross_range_phase",
            "cr2ph_arcs",
        ]:
            results[key].append(np.full([ts_length], np.nan))
        for key in [
            "estimated_thermal",
            "estimated_cross_range",
            "estimated_cross_range_sigma",
            "estimated_thermal_sigma",
        ]:
            results[key].append(np.nan)

        results["succeeded_arcs"].append((np.nan, np.nan))

        if test_stochastics:
            for key in [
                "q_per_partition",
                "std_residuals_partition",
                "rmse_residuals_partition",
                "mean_sigma_post_arc",
            ]:
                stochastic_results[key].append(np.nan)
            stochastic_results["mean_a_priori_sigma_arc"].append(np.mean(Qyy_diagonal))
            stochastic_results["arc_length"].append(arc_length)
            stochastic_results["mean_sigma_p_i"].append(np.mean(slc_quality_i))
            stochastic_results["mean_sigma_p_j"].append(np.mean(slc_quality_j))

    else:
        # Estimate the phases:
        phase_est_2_p_b, phase_th_2_p_b, phase_cross_range_2_p_b, phase_disp_2_p_b, _, _ = (
            _model_arc_partitions_2nd_order_phases(x_data, x_hat_2_p_b)
        )
        # Unwrap the observed phases:
        phase_unwrap_2_p_b = _unwrap_phases(np.angle(dd_arc), phase_est_2_p_b)
        # Estimate 'residual' phase:
        phase_res_2_p_b = phase_unwrap_2_p_b - phase_est_2_p_b

        # Add the results for the arc to the dictionary
        results["unwrap_phases_arc"].append(phase_unwrap_2_p_b)
        results["sigma_phases_arc"].append(Qyy_diagonal)
        results["estimated_phase"].append(phase_est_2_p_b)
        results["estimated_displ_phase"].append(phase_disp_2_p_b)
        results["estimated_thermal"].append(x_hat_2_p_b[-1] * 1000 / m2ph)
        results["estimated_cross_range"].append(x_hat_2_p_b[-2])
        results["estimated_cross_range_sigma"].append(np.sqrt(pcov_2_p_b[-2, -2]))
        results["estimated_thermal_sigma"].append(np.sqrt(pcov_2_p_b[-1, -1]) * 1000 / m2ph)
        results["estimated_thermal_phase"].append(phase_th_2_p_b)
        results["estimated_cross_range_phase"].append(phase_cross_range_2_p_b)
        results["cr2ph_arcs"].append(cr2ph_arc)
        results["succeeded_arcs"].append((pnt_i_idx, pnt_j_idx))
        # Get the dictionaries in the right shape and format
        results = _flatten_arrays_in_dict(results)

        if test_stochastics:
            rmse_res_partition, std_res_partition, q_per_part = _compute_residuals_per_partition_stm(
                phase_unwrap_2_p_b, phase_est_2_p_b, Qyy, bkps
            )
            stochastic_results["q_per_partition"].append(q_per_part)
            stochastic_results["std_residuals_partition"].append(std_res_partition)
            stochastic_results["rmse_residuals_partition"].append(rmse_res_partition)
            stochastic_results["mean_sigma_post_arc"].append(np.std(phase_res_2_p_b))
            stochastic_results["mean_a_priori_sigma_arc"].append(np.mean(Qyy_diagonal))
            stochastic_results["arc_length"].append(arc_length)
            stochastic_results["mean_sigma_p_i"].append(np.mean(slc_quality_i))
            stochastic_results["mean_sigma_p_j"].append(np.mean(slc_quality_j))

            # Get the dictionaries in the right shape and format
            stochastic_results = _flatten_arrays_in_dict(stochastic_results)

        # Step 6. Printing

        if print_output:
            print(
                "Estimated cross_range (phase domain NMAD):",
                np.around(x_hat_arc_ph[0, 0], 2),
                "+/-",
                np.around((np.sqrt(Qx_hat_arc_ph[0, 0])) / (-1 * m2ph), 2),
            )
            print(
                "Estimated cross_range (2nd order + partitions and bounds):",
                np.around(x_hat_2_p_b[-2], 2),
                "+/-",
                np.around((np.sqrt(pcov_2_p_b[-2, -2])) / (-1 * m2ph), 2),
            )
            print(
                "Estimated thermal expansion (phase domain NMAD):",
                np.around(x_hat_arc_ph[1, 0] * 1000 / m2ph, 4),
                "+/-",
                np.around(np.sqrt(Qx_hat_arc_ph[1, 1]) * 1000 / m2ph, 2),
            )
            print(
                "Estimated thermal expansion (2nd order + partitions and bounds):",
                np.around(x_hat_2_p_b[-1] * 1000 / m2ph, 4),
                np.around(np.sqrt(pcov_2_p_b[-1, -1]) * 1000 / m2ph, 2),
            )
            print("")
            print("")

    if test_stochastics:
        return results, stochastic_results

    else:
        return results


def arc_estimation_control_network(
    arcs_to_analyse,
    bounds,
    m2ph,
    n_max_iter,
    Btemporal,
    dates,
    temp,
    sd_complex,
    slc_quality,
    cr2ph,
    ampl_ts,
    bkps_stm,
    mean_ampl_sd,
    sigma_ampl_sd,
    mad_ampl_sd,
    median_ampl_sd,
    x_coordinates,
    y_coordinates,
    coordinate_type: Literal["euclidean", "geographic"] = "euclidean",
    filter_length_complex=30,
    jump_percentage_2pi=0.85,
    vcm_complex_method="mad_median",
    test_stochastics=False,
    print_output=False,
):
    """Estimate parameters for arcs in a control network based on input time series and geodetic measurements.

    This function performs a series of computations for each arc in the control network,
    including variance-covariance matrix computation, double-difference phase estimation, and parameter fitting
    in both the phase and complex domains. The results of these calculations are stored in structured dictionaries.

    Parameters
    ----------
    arcs_to_analyse : list of tuples
        List of arcs, where each arc is defined as a tuple (i, j) representing indices of two points.
    bounds : tuple of lists
        Bounds for parameter estimation in the format (lower_bounds, upper_bounds).
    m2ph : float
        Conversion factor from meters to phase.
    n_max_iter : np.ndarray
        The maximum nr of iterations for non-linear lsq per arc
    Btemporal : numpy.ndarray
        Array of decimal years corresponding to the time series epochs.
    dates : numpy.ndarray
        Array of date indices or timestamps corresponding to the time series.
    temp : numpy.ndarray
        Array of temperature values for thermal expansion modeling.
    sd_complex : numpy.ndarray
        Complex-valued standard deviations of the signal for all points.
    slc_quality : numpy.ndarray
        Quality metric for single-look complex (SLC) data.
    cr2ph : numpy.ndarray
        Cross-range to phase conversion factors for the points.
    ampl_ts : numpy.ndarray
        Amplitude time series for each point.
    bkps_stm : numpy.ndarray
        Breakpoints for state transition modeling.
    mean_ampl_sd : numpy.ndarray
        Mean amplitudes for each point, used in variance modeling.
    sigma_ampl_sd : numpy.ndarray
        Standard deviations of amplitudes for each point.
    mad_ampl_sd : numpy.ndarray
        Median absolute deviations (MAD) of amplitudes.
    median_ampl_sd : numpy.ndarray
        Median amplitudes for each point.
    x_coordinates : numpy.ndarray
        X-coordinates of the points in the control network.
    y_coordinates : numpy.ndarray
        Y-coordinates of the points in the control network.
    coordinate_type: Literal["euclidean", "geographic"], default "euclidean"
        Whether the given coordinates are in Euclidean space (such as RD) or in geographic space (such as lon/lat)
    filter_length_complex : int, optional
        Length of the filter for phase unwrapping (default: 30).
    jump_percentage_2pi : float, optional
        Threshold for unwrapping phase jumps in terms of 2π (default: 0.85).
    vcm_complex_method : str, optional
        Method for variance-covariance matrix estimation in the complex domain.
        Options are "sigma_mean" or "mad_median" (default: "mad_median").
    test_stochastics : bool, optional
        Flag for performing stochastic testing (default: False).
    print_output : bool, optional
        Flag for enabling or disabling print statements (default: False).

    Returns
    -------
    results : dict
        Dictionary containing results for each arc, including:
            - 'unwrap_phases_arc': Unwrapped phases for each arc.
            - 'sigma_phases_arc': Phase variances for each arc.
            - 'estimated_phase': Estimated phases for each arc.
            - 'estimated_displ_phase': Displacement-related phases for each arc.
            - 'estimated_thermal': Estimated thermal expansion coefficients.
            - 'estimated_cross_range': Estimated cross-range components.
            - 'estimated_cross_range_sigma': Uncertainties of cross-range estimates.
            - 'estimated_thermal_sigma': Uncertainties of thermal estimates.
            - 'estimated_thermal_phase': Thermal-related phases for each arc.
            - 'estimated_cross_range_phase': Cross-range related phases for each arc.
            - 'cr2ph_arcs': Cross-range-to-phase conversion factors for each arc.
            - 'succeeded_arcs': List of arcs where parameter estimation succeeded.

    stochastic_results : dict, optional
        Dictionary containing stochastic testing results (if `test_stochastics=True`), including:
            - 'q_per_partition': Quality metrics for each partition.
            - 'std_residuals_partition': Standard deviations of residuals for each partition.
            - 'rmse_residuals_partition': RMSE of residuals for each partition.
            - 'mean_sigma_post_arc': Mean post-fit sigma values for each arc.
            - 'mean_a_priori_sigma_arc': Mean a priori sigma values for each arc.
            - 'arc_length': Lengths of the arcs.
            - 'mean_sigma_p_i': Mean quality metrics for the first point in each arc.
            - 'mean_sigma_p_j': Mean quality metrics for the second point in each arc.

    Notes
    -----
    1. The function uses deterministic assignment for the CR component, setting it to zero for one of the points.
    2. The estimation process includes fallback mechanisms to handle cases where optimal parameters cannot be found.
    3. Requires external utility functions for phase unwrapping, functional model construction, and lsq estimation.

    Raises
    ------
    ValueError
        If an unknown `vcm_complex_method` is specified.
    RuntimeError, ValueError
        If parameter estimation fails for an arc during optimization.
    """
    # If we want to do some tests on the stochastics
    if test_stochastics:
        stochastic_results = {
            "q_per_partition": [],
            "std_residuals_partition": [],
            "rmse_residuals_partition": [],
            "mean_sigma_post_arc": [],
            "mean_a_priori_sigma_arc": [],
            "arc_length": [],
            "mean_sigma_p_i": [],
            "mean_sigma_p_j": [],
        }

    # Dictionary to store results
    results = {
        "unwrap_phases_arc": [],
        "sigma_phases_arc": [],
        "estimated_phase": [],
        "estimated_displ_phase": [],
        "estimated_thermal": [],
        "estimated_cross_range": [],
        "estimated_cross_range_sigma": [],
        "estimated_thermal_sigma": [],
        "estimated_thermal_phase": [],
        "estimated_cross_range_phase": [],
        "cr2ph_arcs": [],
        "succeeded_arcs": [],
    }

    # Counting needed for saving data
    p = 0

    for a in arcs_to_analyse:
        pnt_i_idx, pnt_j_idx = a

        print(f"idx pnt i: {pnt_i_idx}")
        print(f"idx pnt j: {pnt_j_idx}")

        # Extract information of the two points of the arc
        sd_complex_i = sd_complex[pnt_i_idx, :]
        slc_quality_i = slc_quality[pnt_i_idx, :]
        ampl_i = ampl_ts[pnt_i_idx, :]
        bkps_stm_i = bkps_stm[pnt_i_idx, :]
        sigma_ampl_sd_i = sigma_ampl_sd[pnt_i_idx, :]
        mean_ampl_sd_i = mean_ampl_sd[pnt_i_idx, :]
        mad_ampl_sd_i = mad_ampl_sd[pnt_i_idx, :]
        median_ampl_sd_i = median_ampl_sd[pnt_i_idx, :]

        sd_complex_j = sd_complex[pnt_j_idx, :]
        slc_quality_j = slc_quality[pnt_j_idx, :]
        cr2ph_j = cr2ph[pnt_j_idx]
        bkps_stm_j = bkps_stm[pnt_j_idx, :]
        sigma_ampl_sd_j = sigma_ampl_sd[pnt_j_idx, :]
        mean_ampl_sd_j = mean_ampl_sd[pnt_j_idx, :]
        mad_ampl_sd_j = mad_ampl_sd[pnt_j_idx, :]
        median_ampl_sd_j = median_ampl_sd[pnt_j_idx, :]

        # Step  1. Extract information for the ARC
        # Compute the arc length

        arc_length = get_distance(
            [x_coordinates[pnt_i_idx], y_coordinates[pnt_i_idx]],
            [x_coordinates[pnt_j_idx], y_coordinates[pnt_j_idx]],
            mode=coordinate_type,
        )

        # Extract the breakpoints for the arc
        bkps_arc_stm = bkps_stm_i + bkps_stm_j

        # Define the indexes of the breakpoints for the arc
        bkps = [index for index, value in enumerate(bkps_arc_stm) if value > 0]
        bkps.append(len(dates) - 1)

        # The value for the cross range component is cr2ph of point j, since we deterministically set the
        # cross-range component of point i to zero
        cr2ph_arc = cr2ph_j

        # Step 2. Compute the DD phases for the arc
        # point i is the reference point and is subtracted from point j:
        dd_arc, Q_yy_diagonal = _compute_dd(sd_complex_i, sd_complex_j, slc_quality_i, slc_quality_j)

        # Compute the variance covariance matrix of the DD based on the NMAD for the arc
        Q_yy = np.identity(len(dates)) * Q_yy_diagonal**2
        Q_yy_inv = np.linalg.inv(Q_yy)

        # Step 3. Estimate parameters in the phase domain
        # Step is required to get proper intial estimates for the parameter estimation in the complex domain
        # Unwrap the phases based on the filtered real and imaginary part
        phase_arc_unwrap, _, _, _, _ = _unwrap_phases_filter(
            filter_length_complex, dd_arc, np.angle(dd_arc), jump_percentage_2pi
        )

        # Contruct the A matrices for functional model
        # Compute different columns for the A matrices and construct to one A matrix
        A_cr = md.a_cross_range(cr2ph_arc)
        A_off = md.a_offset(len(Btemporal))
        A_lin = md.a_velocity(Btemporal, m2ph)
        A_temp = md.a_temperature(temp)
        A_arc = np.column_stack((A_cr, A_temp, A_off, A_lin))

        # Define the observation vector for the arc, which is based on the 'unwrapped' phase based on the filter
        y_arc = np.reshape(phase_arc_unwrap, (len(phase_arc_unwrap), 1))

        # Estimate parameters in the phase domain
        x_hat_arc_ph, Q_x_hat_arc_ph = est.blue_q_yy_inv(A_arc, y_arc, Q_yy_inv)

        # Step 4. VCM in the complex domain
        # Here we will compute the VCM for the complex domain.
        # It is possible to choose between the mean and sigma or mad and median amplitude per partition.

        # Estimate the DD sigma
        if vcm_complex_method == "sigma_mean":
            sigma_dd = np.abs(mean_ampl_sd_i * mean_ampl_sd_j) * np.sqrt(
                (sigma_ampl_sd_i / mean_ampl_sd_i) ** 2 + (sigma_ampl_sd_j / mean_ampl_sd_j) ** 2
            )

        elif vcm_complex_method == "mad_median":
            sigma_dd = np.abs(median_ampl_sd_i * median_ampl_sd_j) * np.sqrt(
                (mad_ampl_sd_i * 1.4826 / median_ampl_sd_i) ** 2 + (mad_ampl_sd_j * 1.4826 / median_ampl_sd_j) ** 2
            )
        else:
            raise ValueError(
                f"You specified an unknown vcm complex method. The method -- {vcm_complex_method} -- does not exist"
            )

        # Compute VCM in the complex domain
        sigma_complex = np.append(
            sigma_dd, sigma_dd
        )  # Real and Imag are stacked together since we use both of the observations
        Q_dd_cmplx = np.identity(len(sigma_complex))
        np.fill_diagonal(Q_dd_cmplx, sigma_complex**2)

        # Step 5. Parameter estimation in complex domain
        # Complex data preparation for the arc
        re_arc = dd_arc.real
        im_arc = dd_arc.imag
        arc_obs = np.append(re_arc, im_arc)

        # Combine all the independent variables in one independent variable
        x_data = bkps, Btemporal, temp, cr2ph_arc
        xx_data = Btemporal, temp, cr2ph_arc

        # Create initial value arrays
        x0_2_p = np.zeros(3 * len(bkps) + 3)  # Create empty array for the bounds
        x0_2_p[0 : len(bkps)] = np.ones(len(bkps)) * np.max(re_arc)  # The amplitude to be estimated
        x0_2_p[len(bkps) + 1] = x_hat_arc_ph[
            2, 0
        ]  # Interception of the dispalcement polynomial. We use estimated values in the phase domain
        x0_2_p[len(bkps) + 1 : 2 * len(bkps) + 1] = (
            np.ones(len(bkps)) * x_hat_arc_ph[3, 0]
        )  # Value related displacement velocity in the displacement polynomial
        x0_2_p[2 * len(bkps) + 1 : 3 * len(bkps) + 1] = np.zeros(
            len(bkps)
        )  # Value related to the second compont of the dispalcement polynomial
        x0_2_p[-2] = x_hat_arc_ph[0, 0]  # Cross range
        x0_2_p[-1] = x_hat_arc_ph[1, 0]  # Thermal expansion

        # Define bounds
        (
            amp_lower,
            a_lower,
            b_lower,
            c_lower,
            cr_lower,
            exp_lower,
            amp_upper,
            a_upper,
            b_upper,
            c_upper,
            cr_upper,
            exp_upper,
        ) = bounds

        # define bounds for second order polynomial with partitions
        bounds_upper_2_p = np.ones(len(bkps) * 3 + 3)
        bounds_upper_2_p[0 : len(bkps)] = amp_upper
        bounds_upper_2_p[len(bkps)] = a_upper
        bounds_upper_2_p[len(bkps) + 1 : 2 * len(bkps) + 1] = np.ones(len(bkps)) * b_upper
        bounds_upper_2_p[2 * len(bkps) + 1 : 3 * len(bkps) + 1] = np.ones(len(bkps)) * c_upper
        bounds_upper_2_p[-2] = cr_upper
        bounds_upper_2_p[-1] = exp_upper

        bounds_lower_2_p = np.ones(len(bkps) * 3 + 3)
        bounds_lower_2_p[0 : len(bkps)] = amp_lower
        bounds_lower_2_p[len(bkps)] = a_lower
        bounds_lower_2_p[len(bkps) + 1 : 2 * len(bkps) + 1] = np.ones(len(bkps)) * b_lower
        bounds_lower_2_p[2 * len(bkps) + 1 : 3 * len(bkps) + 1] = np.ones(len(bkps)) * c_lower
        bounds_lower_2_p[-2] = cr_lower
        bounds_lower_2_p[-1] = exp_lower

        bounds_2_p = (list(bounds_lower_2_p), list(bounds_upper_2_p))

        # Estiamte parameters in the complex domain with a 2nd order displacement polynomial
        # bounds and partitions
        try:
            x_hat_2_p_b, pcov_2_p_b = _scipy_fit_partition_2nd_order_bounds(
                bkps, xx_data, arc_obs, x0_2_p, bounds_2_p, Q_dd_cmplx, n_max_iter
            )
        except (RuntimeError, ValueError):
            print(f"Optimal parameters not found. Skipping arc {(pnt_i_idx, pnt_j_idx)}")

            # Fill everything with nans
            ts_length = len(ampl_i)

            for key in [
                "unwrap_phases_arc",
                "sigma_phases_arc",
                "estimated_phase",
                "estimated_displ_phase",
                "estimated_thermal_phase",
                "estimated_cross_range_phase",
                "cr2ph_arcs",
            ]:
                results[key].append(np.full([ts_length], np.nan))
            for key in [
                "estimated_thermal",
                "estimated_cross_range",
                "estimated_cross_range_sigma",
                "estimated_thermal_sigma",
            ]:
                results[key].append(np.nan)

            results["succeeded_arcs"].append((np.nan, np.nan))

            if test_stochastics:
                for key in [
                    "q_per_partition",
                    "std_residuals_partition",
                    "rmse_residuals_partition",
                    "mean_sigma_post_arc",
                ]:
                    stochastic_results[key].append(np.nan)
                stochastic_results["mean_a_priori_sigma_arc"].append(np.mean(Q_yy_diagonal))
                stochastic_results["arc_length"].append(arc_length)
                stochastic_results["mean_sigma_p_i"].append(np.mean(slc_quality_i))
                stochastic_results["mean_sigma_p_j"].append(np.mean(slc_quality_j))

        else:
            # Estimate the phases:
            (
                phase_est_2_p_b,
                phase_th_2_p_b,
                phase_cross_range_2_p_b,
                phase_disp_2_p_b,
                _,
                _,
            ) = _model_arc_partitions_2nd_order_phases(x_data, x_hat_2_p_b)
            # Unwrap the observed phases:
            phase_unwrap_2_p_b = _unwrap_phases(np.angle(dd_arc), phase_est_2_p_b)
            # Estimate 'residual' phase:
            phase_res_2_p_b = phase_unwrap_2_p_b - phase_est_2_p_b

            # Add the results for the arc to the dictionary
            results["unwrap_phases_arc"].append(phase_unwrap_2_p_b)
            results["sigma_phases_arc"].append(Q_yy_diagonal)
            results["estimated_phase"].append(phase_est_2_p_b)
            results["estimated_displ_phase"].append(phase_disp_2_p_b)
            results["estimated_thermal"].append(x_hat_2_p_b[-1] * 1000 / m2ph)
            results["estimated_cross_range"].append(x_hat_2_p_b[-2])
            results["estimated_cross_range_sigma"].append(np.sqrt(pcov_2_p_b[-2, -2]))
            results["estimated_thermal_sigma"].append(np.sqrt(pcov_2_p_b[-1, -1]) * 1000 / m2ph)
            results["estimated_thermal_phase"].append(phase_th_2_p_b)
            results["estimated_cross_range_phase"].append(phase_cross_range_2_p_b)
            results["cr2ph_arcs"].append(cr2ph_arc)
            results["succeeded_arcs"].append((pnt_i_idx, pnt_j_idx))
            # Get the dictionaries in the right shape and format
            # results = flatten_arrays_in_dict(results)

            if test_stochastics:
                rmse_res_partition, std_res_partition, q_per_part = _compute_residuals_per_partition_stm(
                    phase_unwrap_2_p_b, phase_est_2_p_b, Q_yy, bkps
                )
                stochastic_results["q_per_partition"].append(q_per_part)
                stochastic_results["std_residuals_partition"].append(std_res_partition)
                stochastic_results["rmse_residuals_partition"].append(rmse_res_partition)
                stochastic_results["mean_sigma_post_arc"].append(np.std(phase_res_2_p_b))
                stochastic_results["mean_a_priori_sigma_arc"].append(np.mean(Q_yy_diagonal))
                stochastic_results["arc_length"].append(arc_length)
                stochastic_results["mean_sigma_p_i"].append(np.mean(slc_quality_i))
                stochastic_results["mean_sigma_p_j"].append(np.mean(slc_quality_j))

                # Get the dictionaries in the right shape and format
                # stochastic_results = flatten_arrays_in_dict(stochastic_results)

            if print_output:
                print(
                    "Estimated cross_range (phase domain NMAD):",
                    np.around(x_hat_arc_ph[0, 0], 2),
                    "+/-",
                    np.around((np.sqrt(Q_x_hat_arc_ph[0, 0])) / (-1 * m2ph), 2),
                )
                print(
                    "Estimated cross_range (2nd order + partitions and bounds):",
                    np.around(x_hat_2_p_b[-2], 2),
                    "+/-",
                    np.around((np.sqrt(pcov_2_p_b[-2, -2])) / (-1 * m2ph), 2),
                )
                print(
                    "Estimated thermal expansion (phase domain NMAD):",
                    np.around(x_hat_arc_ph[1, 0] * 1000 / m2ph, 4),
                    "+/-",
                    np.around(np.sqrt(Q_x_hat_arc_ph[1, 1]) * 1000 / m2ph, 2),
                )
                print(
                    "Estimated thermal expansion (2nd order + partitions and bounds):",
                    np.around(x_hat_2_p_b[-1] * 1000 / m2ph, 4),
                    np.around(np.sqrt(pcov_2_p_b[-1, -1]) * 1000 / m2ph, 2),
                )
                print("")
                print("")

        p = p + 1

    # tranform from list to np array
    results = {key: np.array(value) for key, value in results.items()}

    if test_stochastics:
        return results, stochastic_results

    return results


def periodogram(
    stm: xr.Dataset,
    key_dphase: str,
    key_h2ph: str,
    key_Btemporal: str,
    std_obs: float = 1.0,
    std_height: float = 50.0,
    std_vel: float = 0.02,
    init_height: float = 0.0,
    init_vel: float = 0.0,
    init_step_height: float = 1.0,
    init_step_vel: float = 1e-3,
    min_steps: int = 10,
):
    """Periodogram algorithm.

    This function performs periodogram unwrapping on arcs.

    It uses a deformation model with two parameters: height and velocity to estimate the unwrapped phase.

    For computation efficiency, the design matrix is constructed only once for all arcs, utilizing the average
    height-to-phase conversion factor (h2ph) across all arcs. The effect of using this average is corrected later.

    Parameters
    ----------
    stm : xr.Dataset
        Input Space-Time Matrix (STM) containing the wrapped phase, height-to-phase conversion factor, and year-time.
    key_dphase : str
        Key for the wrapped differential phase data variable in the STM.
    key_h2ph : str
        Key for the height-to-phase conversion factor in the STM.
    key_Btemporal : str
        Key for the temporal baseline in the STM.
        The value should be in decimal years.
    std_obs : float, optional
        A-poriori standard deviation of the observations in rads, by default 1.0.
        This value is used to construct the stochastic model (Qyy) of the observations.
    std_height : float, optional
        A-priori standard deviation of the height in meters, by default 50.0.
        This value is used to construct the boundaries of the initial search space for the height parameter.
    std_vel : float, optional
        A-priori standard deviation of the velocity in meters per year, by default 0.02.
        This value is used to construct the boundaries of the initial search space for the velocity parameter.
    init_height : float, optional
        Initial value for the height parameter in meters, by default 0.0.
    init_vel : float, optional
        Initial value for the velocity parameter in meters per year, by default 0.0.
    init_step_height : float, optional
        Initial step size for the height parameter in meters, by default 1.0.
        This value sets the resolution of the initial search space for the height parameter.
        After every search, the step size will be reduced by a factor of 10.
    init_step_vel : float, optional
        Initial step size for the velocity parameter in meters per year, by default 1e-3.
        This value sets the resolution of the initial search space for the velocity parameter.
        After every search, the step size will be reduced by a factor of 10.
    min_steps : int, optional
        Minimum number of steps in the search space for the height and velocity parameters, by default 10.
        If the number of steps in the initial search space is smaller than this value, it will be set to this value.
        After the first search, the number of steps will be set to this value.

    Returns
    -------
    Tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]
        Returns the unwrapped phase, ambiguities, estimated height, estimated velocity, and temporal coherence.
        - Unwrapped phase: in rads, shape (n_arcs, n_obs), dtype np.float64.
        - Ambiguities: unitless, shape (n_arcs, n_obs), dtype np.float64.
        - Estimated height: in meters, shape (n_arcs,), dtype np.float64.
        - Estimated velocity: in meters per year, shape (n_arcs,), dtype np.float64.
        - Temporal coherence: unitless float number, norm of the complex coherence, scalar, dtype np.float64.
    """
    # Compute m2ph (meters to phase) conversion factor from wavelength
    if "wavelength" not in stm.attrs:
        raise ValueError(
            "Wavelength is not provided and not found in attributes of STM."
            "Please make sure it is provided."
            "For example: stm = stm.assign_attrs({'wavelength': wavelength})"
        )
    wavelength = stm.attrs["wavelength"]
    m2ph = -4 * np.pi / wavelength

    # Make sure year time only contains the time dimension
    assert (len(stm[key_Btemporal].dims) == 1) and (
        "time" in stm[key_Btemporal].dims
    ), "year time should and only should contain the 'time' dimension."

    # Load year time in memory
    Btemporal = stm[key_Btemporal].values

    # Set up functional and stochastic model for all arcs
    # Here we use the same h2ph (average over all arcs) for all arcs and correct the effect later
    # Doing this avoids perform matrix inversion for each arc
    h2ph_approx = stm[key_h2ph].mean(dim="space").values  # Mean h2ph of all arcs

    # Design matrix B, size n_obs x n_params
    # In B, h2ph should also be multiplied by m2ph since it did not when it was created
    B = np.stack([h2ph_approx * m2ph, Btemporal * m2ph]).T

    # Stochastic model Qyy, size n_obs x n_obs
    # This is the covariance matrix of the observations
    n_obs = stm[key_dphase].sizes["time"]  # number of observations
    Qyy = np.diag(np.repeat(std_obs**2, n_obs))

    # Normal matrix N and rhs for the least squares solution
    N = B.T @ np.linalg.inv(Qyy) @ B  # B.T * Qyy^-1 * B , size n_params x n_params
    # Solve N * x = B.T * Qyy^-1, then rhs =  N^-1 * B.T * Qyy^-1
    rhs = np.linalg.inv(N) @ B.T @ np.linalg.inv(Qyy)

    # check if the time dimension is not chunked, and unchunk it if necessary
    chunk_sizes = dict(zip(list(stm.sizes), stm.chunks, strict=False))
    if "time" in chunk_sizes.keys():
        if chunk_sizes["time"] != 1:
            stm = stm.chunk({"time": -1})

    # Build initial search space for height and velocity
    n_steps_height = max(round(2 * std_height / init_step_height), min_steps)
    n_steps_vel = max(round(2 * std_vel / init_step_vel), min_steps)
    init_search_space = _build_periodogram_search_space(
        init_height, init_vel, init_step_height, init_step_vel, n_steps_height, n_steps_vel
    )

    # Perform one search for all arcs, and get the best initial estimates for height and velocity per arc
    # This is motivated by the fact that the initial search space is the largest, and can be vectorized for all arcs
    # First iteration, candidate modeled phases are identical for all arcs
    # The redisuals phase_residual_all_arcs is a large array with n_arcs x n_obs x n_search
    # so use .data to avoid loading it into memory if it is a dask array
    dphase_obs = stm[key_dphase].data[:, :, None]  # n_arcs x n_obs x 1
    phs_model = B @ init_search_space.T  # n_obs x n_search
    if isinstance(dphase_obs, da.Array):
        # If dphase_obs is a dask array
        # Also chunk the phs_model in the search dimension, making each chunk about 100 MB
        # No chunk in observation dimension since we need to do sum in that dimension
        chunksize_search = max(1, 100 * 1024**2 // (dphase_obs.chunks[0][0] * n_obs * dphase_obs.dtype.itemsize))
        phs_model = da.from_array(phs_model, chunks=(n_obs, chunksize_search))
    phs_model = phs_model[None, :, :]  # 1 x n_obs x n_search
    phase_residual_all_arcs = dphase_obs - phs_model
    coh_search_space_all_arcs = (
        np.cos(phase_residual_all_arcs).sum(axis=1) + 1j * np.sin(phase_residual_all_arcs).sum(axis=1)
    ) / stm[key_dphase].sizes["time"]  # n_arcs x n_search
    coh_idx_all_arcs = np.argmax(np.abs(coh_search_space_all_arcs), axis=1)  # n_arcs

    # Build xr.DataArray for the initial height and velocity of all arcs
    da_init_height_all_arcs = xr.DataArray(
        init_search_space[coh_idx_all_arcs, 0],
        dims=["space"],
    )
    da_init_vel_all_arcs = xr.DataArray(
        init_search_space[coh_idx_all_arcs, 1],
        dims=["space"],
    )

    # Apply the _periodogram_arc on stm[key_dphase] along "space" dimension
    # Set up input core dimensions, which are the dimensions _periodogram_arc will be applied to
    # We are broadcasting _periodogram_arc on stm[key_dphase] and stm[key_h2ph] along the space dimension
    # The height and velocity are scalars
    # Therefore, we are only calling it on the "time" dimension for the first two parameters
    # So we have the input_core_dims as  [["time"], ["time", [], []]
    input_core_dims = [["time"], ["time"], [], []]

    # There are 5 outputs from _periodogram_arc
    # The first two are np arrays with time dimension
    # The other three are scalars, so they have no dimensions
    output_core_dims = [["time"], ["time"], [], [], []]

    results = xr.apply_ufunc(
        _periodogram_arc,
        stm[key_dphase],
        stm[key_h2ph],
        da_init_height_all_arcs,
        da_init_vel_all_arcs,
        input_core_dims=input_core_dims,
        output_core_dims=output_core_dims,
        kwargs={
            "h2ph_approx": h2ph_approx,
            "B": B,
            "Qyy": Qyy,
            "N": N,
            "rhs": rhs,
            "init_step_height": init_step_height,
            "init_step_vel": init_step_vel,
            "min_steps": min_steps,
        },
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64, np.float64, np.float64, np.float64, np.float64],
    )

    return results


def _periodogram_arc(
    phs_obs_wrapped: np.ndarray,
    h2ph: np.ndarray,
    init_height: float,
    init_vel: float,
    h2ph_approx: np.ndarray,
    B: np.ndarray,
    Qyy: np.ndarray,
    N: np.ndarray,
    rhs: np.ndarray,
    init_step_height: float,
    init_step_vel: float,
    min_steps: float,
):
    """Periodogram unwrapping for a single arc.

    Parameters
    ----------
    phs_obs_wrapped : np.ndarray
        Wrapped phase observations in radians, shape (n_obs,).
    h2ph : np.ndarray:
        Height-to-phase factor of the arc, shape (n_obs,).
    init_height : float
        Initial value for the height parameter in meters.
    init_vel : float
        Initial value for the velocity parameter in meters per year.
    h2ph_approx : np.ndarray
        Approximate height-to-phase factor calculated by spatial average of all h2ph, shape (n_obs,).
    B : np.ndarray
        Design matrix, size n_obs x n_params, where n_params = 2 (height and velocity).
    Qyy : np.ndarray
        Stochastic model of the observations, size n_obs x n_obs.
    N : np.ndarray
        Normal matrix, size n_params x n_params.
    rhs : np.ndarray
        Right-hand side matrix for the least squares solution, size n_params x n_obs.
    init_step_height : float
        Initial step size for the height parameter in meters.
    init_step_vel : float
        Initial step size for the velocity parameter in meters per year.
    min_steps : float
        Minimum number of steps in the search space for the height and velocity parameters.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, float, float, float]
        Returns the unwrapped phase, ambiguities, estimated height, estimated velocity, and temporal coherence.
        - Unwrapped phase: in rads, shape (n_obs,), dtype np.float64.
        - Ambiguities: unitless, shape (n_obs,), dtype np.float64.
        - Estimated height: in meters, scalar, dtype np.float64.
        - Estimated velocity: in meters per year, scalar, dtype np.float64.
        - Temporal coherence: unitless float number, norm of the complex coherence, scalar, dtype np.float64.
    """
    step_height = init_step_height
    step_vel = init_step_vel
    param_height = init_height
    param_vel = init_vel
    count = 0

    # Calculate the initial temporal coherence for the initial height and velocity,
    # in case the search loop is not entered
    phs_model = B @ np.array([param_height, param_vel])  # size n_obs
    phase_residual = phs_obs_wrapped[:, None] - phs_model
    coh_best = (np.cos(phase_residual).sum() + 1j * np.sin(phase_residual).sum()) / phs_obs_wrapped.shape[0]

    # Search loop
    while step_height > STOP_HEIGHT and step_vel > STOP_VEL and count < MAX_COUNT:
        # Build search space
        search_space = _build_periodogram_search_space(
            param_height, param_vel, step_height, step_vel, min_steps, min_steps
        )

        # Calculate the wrapped model phase for all candidates
        phs_model = wrap_phase(B @ search_space.T)  # size n_obs x n_search

        # Calculate the temporal coherence for all search candidates
        # Expand dimension of phs_obs_wrapped to facilitate broadcasting
        # No need to repeat phs_obs_wrapped since the minus operation will broadcast to the shape of phs_model
        # Sum along axis=0 which is the observation axis
        # Reference: van Leijen 2014, Eq. 4.55
        # The following implementation equivalent to:
        # np.exp(1j * (np.expand_dims(phs_obs_wrapped, axis=1) - phs_model)).sum(axis=0) / phs_obs_wrapped.shape[0]
        phase_residual = phs_obs_wrapped[:, None] - phs_model
        coh_search_space = (
            np.cos(phase_residual).sum(axis=0) + 1j * np.sin(phase_residual).sum(axis=0)
        ) / phs_obs_wrapped.shape[0]

        # Get the best temporal coherence value and its index
        coh_idx = np.argmax(np.abs(coh_search_space))
        coh_best = coh_search_space[coh_idx]

        # Update values needed for search space
        # Reduce step size to 1/10
        param_height = search_space[coh_idx, 0]
        param_vel = search_space[coh_idx, 1]
        step_height /= 10
        step_vel /= 10

        count += 1

    # Correct the height parameter for using h2ph_approx
    # Method copied from MATLAB DePSI code
    factor = np.median(h2ph / h2ph_approx)  # correct factor
    param_height = param_height / factor

    # Calculate the modelled phase and unwrapped phase
    model_est = B @ np.array([param_height, param_vel]) + np.angle(coh_best)  # Absolute modelled phase
    dphase_new = wrap_phase(phs_obs_wrapped - model_est)  # Wrapped modelled phase
    ambiguities = np.round((model_est + dphase_new - phs_obs_wrapped) / (2 * np.pi))  # Ambiguities
    phs_obs_unwrapped = 2 * np.pi * ambiguities + phs_obs_wrapped  # Unwrapped phase
    param = rhs @ phs_obs_unwrapped  # [height_est, velocity_est]

    return phs_obs_unwrapped, ambiguities, param[0], param[1], np.abs(coh_best)


def _build_periodogram_search_space(init_height, init_vel, step_height, step_vel, n_steps_height, n_steps_vel):
    """Construct the periodogram search space for height and velocity parameters.

    For both height and velocity, the candidates are generated around the initial values according to the step size
    and the number of steps. On each side of the initial value, N candidates are generated with a step size, where
    N is specified by `n_steps_height` and `n_steps_vel`, and the step size is specified by `step_height` and
    `step_vel`.

    Then all possible combinations of height and velocity candidates are created to form
    the search space.

    Parameters
    ----------
    init_height : float
        Initial height parameter in meters.
    init_vel : float
        Initial velocity parameter in meters per year.
    step_height : int
        Step size of search for height parameter, in meters.
    step_vel : int
        Step size of search for velocity parameter, in meters per year.
    n_steps_height : int
        Number of steps for height parameter on each side of the initial value.
    n_steps_vel : int
        Number of steps for velocity parameter on each side of the initial value.

    Returns
    -------
    np.ndarray
        Search space for height and velocity parameters, shape (n_candidates_vel * n_candidates_height, 2)
    """
    height_candidates = np.arange(
        init_height - n_steps_height * step_height,
        init_height + n_steps_height * step_height + step_height,
        step_height,
    )

    vel_candidates = np.arange(
        init_vel - n_steps_vel * step_vel, init_vel + n_steps_vel * step_vel + step_vel, step_vel
    )

    # All possible combinations of height and velocity
    search_space = np.array(np.meshgrid(height_candidates, vel_candidates)).T.reshape(-1, 2)

    return search_space
