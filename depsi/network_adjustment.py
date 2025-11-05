"""Functions for the network adjustment in the CORG network."""

from itertools import combinations, product

import numpy as np
import xarray as xr
from scipy.sparse import csr_matrix

import depsi.arc_estimation as arc_est
import depsi.deformation_models as dm
import depsi.estimation as est
import depsi.network as dn


def network_adjustment_control_network(A_adjustment, y_obs, sigma_obs, threshold_Qyy_change=0.005, epsilon=0):
    """Perform a network adjustment to estimate the parameters of points using the estimated parameters for the arc.

    The function:
    1. Reshapes the observation vector.
    2. Constructs the variance-covariance matrix (VCM) for the observations.
    3. Inverts the VCM to obtain the precision matrix.
    4. Uses the BLUE method to compute the adjusted parameters (x_hat) and their associated covariance matrix (Qx_hat).

    Parameters
    ----------
    A_adjustment : numpy.ndarray
        The design matrix that links the observations (arcs) to the unknown parameters (points).
    y_obs : numpy.ndarray
        The observation vector (or array) containing the estimated parameters for the arcs.
    sigma_obs : numpy.ndarray
        The standard deviations (uncertainties) of the observations.
    threshold_Qyy_change : float, optional
        Threshold for acceptable changes to the corrected VCM, by default 0.005.
    epsilon : float, optional
        Smallest allowable eigenvalue for ensuring the VCM is positive-definite, by default 0.

    Returns
    -------
    x_hat : numpy.ndarray
        The adjusted estimates of the unknown parameters (points).
    Qx_hat : numpy.ndarray
        The variance-covariance matrix of the adjusted parameters.
    Qyy : numpy.ndarray
        The variance-covariance matrix of the observations.
    Qyy_inv : numpy.ndarray
        The inverse of the variance-covariance matrix (precision matrix) of the observations.
    y : numpy.ndarray
        The reshaped observation vector.
    y_hat : numpy.ndarray
        The adjusted observation vector.
    e_hat : numpy.ndarray
        The residuals between the observations and the adjusted observations.
    """
    # Construct the observation vector
    y = y_obs.reshape(len(y_obs), 1)

    # Construct the VCM of the observations
    q_yy = np.identity(len(y_obs))
    np.fill_diagonal(q_yy, sigma_obs**2)

    # Compute eigenvectors and eigenvalues to check whether matrix is positive-definite, if the matrix is
    # positive-definite we need to correct it
    eigvals, eigvecs = np.linalg.eigh(q_yy)
    if np.any(eigvals < 0):
        # print(f"There are non negative eigenvalues: Eigenvalues: {eigvals}")
        print("There are non negative eigenvalues: Qyy will be corrected...")

        eigvals[eigvals < epsilon] = epsilon
        Qyy_corrected = eigvecs @ np.diag(eigvals) @ eigvecs.T
        # Check the change
        change = np.linalg.norm(Qyy_corrected - q_yy, ord="fro")
        if change < threshold_Qyy_change:
            q_yy = Qyy_corrected
            condition_number = np.linalg.cond(q_yy)
            Qyy_inv = np.linalg.pinv(q_yy) if condition_number > 1e12 else np.linalg.inv(q_yy)
        else:
            raise ValueError(f"Qyy correction too large: {change} > {threshold_Qyy_change}")

    else:
        Qyy_inv = np.linalg.inv(q_yy)

    # Integrate the network
    x_hat, Qx_hat = est.blue_q_yy_inv(A_adjustment, y, Qyy_inv)

    # Estimate the y_hat and e_hat
    y_hat = A_adjustment @ x_hat
    e_hat = y - y_hat

    return x_hat, Qx_hat, q_yy, Qyy_inv, y, y_hat, e_hat


def apply_w_test_control_network(m, A, Qx_hat, Qyy, Qyy_inv, y, e_hat):
    """Apply the w-test to detect and re-wrap arc outliers during a network adjustment.

    This function identifies the observation with the largest deviation (outlier) among a set of observations
    and adjusts it by re-wrapping the phase. The w-test calculates the standardized residuals of the observations,
    enabling the identification of outliers based on their significance.

    The function:
    1. Computes the estimated variance-covariance matrix of the observations.
    2. Applies the w-test to each observation to identify the one with the largest standardized deviation.
    3. Re-wraps the identified outlier's observation by adjusting its value to minimize its deviation.

    Parameters
    ----------
    m : int
        The number of observations (outliers) to be tested.
    A : numpy.ndarray
        The design matrix of the adjustment model, linking observations to the unknown parameters.
    Qx_hat : numpy.ndarray
        The variance-covariance matrix of the estimated parameters from the network adjustment.
    Qyy : numpy.ndarray
        The variance-covariance matrix of the observations.
    Qyy_inv : numpy.ndarray
        The inverse of the variance-covariance matrix of the observations.
    y : numpy.ndarray
        The vector of observations (arc measurements).
    e_hat : numpy.ndarray
        The residuals (differences between observed and estimated values).

    Returns
    -------
    corrected_y : numpy.ndarray
        The corrected vector of observations after re-wrapping the biggest outlier.
    idx_biggest_w : int
        The index of the observation with the largest w-test value, i.e., the biggest detected outlier.
    correct : int
        The correction applied to the outlier observation:
        -1 for decreasing the value,
        1 for increasing the value,
        0 for no change.

    Example
    -------
    corrected_y, idx_outlier, correction_type = apply_w_test_control_network(m, A, Qx_hat, Qyy, Qyy_inv, y, e_hat)
    """
    # Construct vector where we save the w-test results per detected outlier
    w_result = np.zeros(m)

    # Compute Qyy_hat and Qee
    Qyy_hat = A @ Qx_hat @ A.transpose()
    Qee = Qyy - Qyy_hat

    # Apply the w-test for every observation
    # the observation with the highest result, will be flaged as a potential outlier
    # Remark that we apply a simplification below since Qyy is often a diagonal matrix
    for w in range(m):
        # c_i = np.zeros((m, 1))
        # c_i[w, 0] = 1
        # A = (c_i.transpose() @ Qyy_inv @ e_hat)[0, 0]
        # B = np.sqrt((c_i.transpose() @ Qyy_inv @ Qee @ Qyy_inv @ c_i)[0, 0])
        # w_result[w] = A / B

        w_result[w] = e_hat[w] / np.sqrt(
            Qee[w, w]
        )  # Simplification since Qyy is a diagonal matrix, see above for full equation

    w_result = np.abs(w_result)

    # Re-wrap the observation with the highest w-test result
    idx_biggest_w = np.argmax(w_result)
    corrected_y = np.copy(y)

    # Compute whether we need to shift up or down
    mean_excl_outlier = np.mean(np.delete(y, idx_biggest_w))

    # Compute the difference
    mean_difference = y[idx_biggest_w] - mean_excl_outlier

    if mean_difference > 0:
        # if e_hat[idx_biggest_w] > 0:
        corrected_y[idx_biggest_w] = corrected_y[idx_biggest_w] - 2 * np.pi
        correct = -1
    elif mean_difference < 0:
        # if e_hat[idx_biggest_w] < 0:
        corrected_y[idx_biggest_w] = corrected_y[idx_biggest_w] + 2 * np.pi
        correct = 1
    else:
        e_hat[idx_biggest_w] = corrected_y[idx_biggest_w]
        correct = 0

    return corrected_y, idx_biggest_w, correct


def update_estimated_parameters(
    estimated_values,
    estimated_time_series,
    adjustment_points,
    arc_add,
    nr_epochs,
    point_cross_range,
    Qx_cross_range,
    point_thermal,
    Qx_thermal,
    point_time_series,
    point_time_series_corrected,
    point_time_series_vcm,
):
    """Update the estimated values, time series, and VCM for the given points at a specific iteration.

    This function iterates over the control points, updating the estimated heights, thermal components, and VCM.
    We do multiple network adjustments by adding additional arcs. For every adjustment we get new estimates for
    the unknown parameters. Here we update the dictionaries per iteration.

    The function:
    1. Checks if the point is new or already exists in the dictionaries.
    2. For new points, initializes the dictionaries with `None` values for previous iterations.
    3. For each iteration, appends the estimated values, time series data, and their variances to the dictionaries.

    Parameters
    ----------
    estimated_values : dict
        Dictionary storing the estimated CR, thermal components, and their variance for all points across iterations.
    estimated_time_series : dict
        Dictionary storing the time series data, corrected time series, and their vars for all points across iterations.
    adjustment_points : numpy.ndarray
        Indices of the points being adjusted at the current iteration.
    arc_add : int
        The current iteration index (arc being added).
    nr_epochs : int
        The number of epochs for the time series.
    point_cross_range : numpy.ndarray
        The new estimated cross range values for the adjustment points at the current iteration.
    Qx_cross_range : numpy.ndarray
        The variance-covariance matrix corresponding to the cross range estimates.
    point_thermal : numpy.ndarray
        The new estimated thermal components for the adjustment points at the current iteration.
    Qx_thermal : numpy.ndarray
        The variance-covariance matrix corresponding to the thermal estimates.
    point_time_series : numpy.ndarray
        The time series data for the adjustment points at the current iteration.
    point_time_series_corrected : numpy.array
        The time series data for the adjustment points if it was re-wrapped
    point_time_series_vcm : numpy.ndarray
        The variance-covariance matrix for the time series data at the current iteration.

    Returns
    -------
    estimated_values : dict
        Updated dictionary storing the estimated heights, thermal components, and their variances across iterations.
    estimated_time_series : dict
        Updated dictionary storing the time series data, corrected time series, and their variances across iterations.

    Example
    -------
    estimated_values, estimated_time_series = update_estimated_parameters(
        estimated_values,
        estimated_time_series,
        adjustment_points,
        arc_add,
        nr_epochs,
        point_heights,
        Qx_height,
        point_thermal,
        Qx_thermal,
        point_time_series,
        point_time_series_vcm
    )
    """
    # Loop over each point in the control network and update estimated values and time series
    for i, point in enumerate(adjustment_points):
        if point not in estimated_values:
            # Initialize estimated values and time series for new points
            estimated_values[point] = {
                "cross_range": [None] * arc_add,
                "cross_range_variances": [None] * arc_add,
                "thermal_comp": [None] * arc_add,
                "thermal_comp_variances": [None] * arc_add,
            }
            estimated_time_series[point] = {
                "time_series": [None] * arc_add,
                "time_series_corrected": [None] * arc_add,
                "time_series_variances": [None] * arc_add,
            }
            # Initialize arrays for previous iterations
            for j in range(arc_add):
                estimated_time_series[point]["time_series"][j] = [None] * nr_epochs
                estimated_time_series[point]["time_series_corrected"][j] = [None] * nr_epochs
                estimated_time_series[point]["time_series_variances"][j] = [None] * nr_epochs

        # Append the current iteration's values to the estimated values
        estimated_values[point]["cross_range"].append(point_cross_range[i, 0])
        estimated_values[point]["cross_range_variances"].append(Qx_cross_range[i, i])
        estimated_values[point]["thermal_comp"].append(point_thermal[i, 0])
        estimated_values[point]["thermal_comp_variances"].append(Qx_thermal[i, i])

        # Append the time series data for the current iteration
        estimated_time_series[point]["time_series"].append(point_time_series[i, :])
        estimated_time_series[point]["time_series_corrected"].append(point_time_series_corrected[i, :])
        estimated_time_series[point]["time_series_variances"].append(point_time_series_vcm[i, i, :])

    return estimated_values, estimated_time_series


def adjustment_matrix_control_network(arcs, nr_pnts, points_network):
    """Construct a sparse design matrix (A matrix) for network adjustment in the 'control' network.

    This function creates the A matrix, which describes the relationships between arcs (observations) and points
    in the network.

    The function:
    1. Defines the number of equations based on the number of arcs.
    2. Creates a sparse matrix (in CSR format) that links the arcs to the corresponding points.
    3. Converts the sparse matrix to a dense matrix for further processing.
    4. Extracts a submatrix (`A_small`) that represents only the specified points in the network.

    Parameters
    ----------
    arcs : numpy.ndarray
        A 2D array where each row represents an arc. Each arc is described by two point indices (start and end point).
        The first point is denoted as point i, and the second point as point j.
    nr_pnts : int
        The total number of points in the network.
    points_network : numpy.ndarray
        A 1D array of indices corresponding to the points involved in the current adjustment. The reference point
        should not be included in this array.

    Returns
    -------
    A_small : numpy.ndarray
        A dense submatrix of the design matrix that represents the part of the matrix corresponding to the selected
        points in the network.
    """
    neq = len(arcs)  # The number of arcs define the number of equations in the A matrix
    npt = nr_pnts
    eqs = np.concatenate((np.ones(neq), -1 * np.ones(neq)))
    rows = np.concatenate((np.arange(neq), np.arange(neq)))  # Rows in A matrix equal to nr of eqs.
    cols = np.concatenate((arcs[:, 1], arcs[:, 0]))  # the columnds are defined by the points in adjustment

    A = csr_matrix((eqs, (rows, cols)), shape=(neq, npt))
    A_dense = A.toarray()

    # A_dense will contain all points from the stm matrix (e.g., can be up to 10000)
    # For a lot of points we have no arcs, therefore we can remove these points from the matrix
    # resulting in A_small
    A_small = A_dense[:, points_network]

    return A_small


def fill_covariances_from_shared_points(q_yy, arcs, sigma_points, reference_point=None):
    """Fill the covariances based on shared points between observations.

    Parameters
    ----------
    q_yy : np.ndarray
        The already initialized variance-covariance matrix (VCM) with variances on the diagonal.
    arcs : list of tuples
        List of (start, end) points for each observation.
    sigma_points : dict
        Standard deviations (sigmas) of the individual points.
    reference_point : int or None, optional
        The reference point with a fixed value (e.g., 0), default is None.
        If a shared point equals the reference point, no covariance is added.

    Returns
    -------
    q_yy : np.ndarray
        The updated variance-covariance matrix.
    """
    n = len(arcs)

    for i in range(n):
        a_i, b_i = arcs[i]
        for j in range(i + 1, n):
            a_j, b_j = arcs[j]

            cov = 0.0
            if a_i == a_j and a_i != reference_point:
                cov += sigma_points.get(a_i, 0.0) ** 2
            if a_i == b_j and a_i != reference_point:
                # cov -= sigma_points.get(a_i, 0.0) ** 2
                cov += sigma_points.get(a_i, 0.0) ** 2
            if b_i == a_j and b_i != reference_point:
                # cov -= sigma_points.get(b_i, 0.0) ** 2
                cov += sigma_points.get(b_i, 0.0) ** 2
            if b_i == b_j and b_i != reference_point:
                cov += sigma_points.get(b_i, 0.0) ** 2

            # if a_i == a_j and a_i != reference_point:
            #     cov += sigma_points.get(a_i, 0.0)
            # if a_i == b_j and a_i != reference_point:
            #     cov -= sigma_points.get(a_i, 0.0)
            # if b_i == a_j and b_i != reference_point:
            #     cov -= sigma_points.get(b_i, 0.0)
            # if b_i == b_j and b_i != reference_point:
            #     cov += sigma_points.get(b_i, 0.0)

            q_yy[i, j] = cov
            q_yy[j, i] = cov  # Symmetric

    return q_yy


def network_adjustment_control_network_displ(
    A_adjustment, y_obs, sigma_obs, sigma_points, adjustment_arcs, ref_pnt, threshold_Qyy_change=0.005, epsilon=0
):
    """Perform a network adjustment to estimate the parameters of points using the estimated parameters for the arc.

    The function:
    1. Reshapes the observation vector.
    2. Constructs the variance-covariance matrix (VCM) for the observations.
    3. Inverts the VCM to obtain the precision matrix.
    4. Uses the BLUE method to compute the adjusted parameters (x_hat) and their associated covariance matrix (Qx_hat).

    Parameters
    ----------
    A_adjustment : numpy.ndarray
        The design matrix that links the observations (arcs) to the unknown parameters (points).
    y_obs : numpy.ndarray
        The observation vector (or array) containing the estimated parameters for the arcs.
    sigma_obs : numpy.ndarray
        The standard deviations (uncertainties) of the observations.
    threshold_Qyy_change : float, optional
        Threshold for acceptable changes to the corrected VCM, by default 0.005.
    epsilon : float, optional
        Smallest allowable eigenvalue for ensuring the VCM is positive-definite, by default 0.
    adjustment_arcs:
        ?
    ref_pnt:
        Reference point
    sigma_points:
        ?

    Returns
    -------
    x_hat : numpy.ndarray
        The adjusted estimates of the unknown parameters (points).
    Qx_hat : numpy.ndarray
        The variance-covariance matrix of the adjusted parameters.
    Qyy : numpy.ndarray
        The variance-covariance matrix of the observations.
    Qyy_inv : numpy.ndarray
        The inverse of the variance-covariance matrix (precision matrix) of the observations.
    y : numpy.ndarray
        The reshaped observation vector.
    y_hat : numpy.ndarray
        The adjusted observation vector.
    e_hat : numpy.ndarray
        The residuals between the observations and the adjusted observations.
    """
    # Construct the observation vector
    y = y_obs.reshape(len(y_obs), 1)

    # fill the diagonal with the arc sigma values
    q_yy_1 = np.identity(len(y_obs))
    np.fill_diagonal(q_yy_1, sigma_obs**2)

    # fill the off diagonal with correlations between the arcs
    q_yy = fill_covariances_from_shared_points(q_yy_1, adjustment_arcs, sigma_points, ref_pnt)

    # Compute eigenvectors and eigenvalues to check whether matrix is positive-definite, if the matrix is
    # positive-definite we need to correct it
    eigvals, eigvecs = np.linalg.eigh(q_yy)
    if np.any(eigvals < 0):
        # print(f"There are non negative eigenvalues: Eigenvalues: {eigvals}")
        print("There are non negative eigenvalues: Eigenvalues: Qyy will be corrected...")

        eigvals[eigvals < epsilon] = epsilon
        Qyy_corrected = eigvecs @ np.diag(eigvals) @ eigvecs.T
        # Check the change
        change = np.linalg.norm(Qyy_corrected - q_yy, ord="fro")
        if change < threshold_Qyy_change:
            q_yy = Qyy_corrected
            condition_number = np.linalg.cond(q_yy)
            Qyy_inv = np.linalg.pinv(q_yy) if condition_number > 1e12 else np.linalg.inv(q_yy)
        else:
            raise ValueError(f"Qyy correction too large: {change} > {threshold_Qyy_change}")

    else:
        Qyy_inv = np.linalg.inv(q_yy)

    # Integrate the network
    x_hat, Qx_hat = est.blue_q_yy_inv(A_adjustment, y, Qyy_inv)

    # Estimate the y_hat and e_hat
    y_hat = A_adjustment @ x_hat
    e_hat = y - y_hat

    return x_hat, Qx_hat, q_yy, Qyy_inv, y, y_hat, e_hat


def _estimate_connection_point_stm(variable, stm_control, arc_variable, arc_sigma, nr_conn, control_conn_points):
    """Estimate a variable (e.g., cross-range or displacement) for a conneciton point relative to the control network.

    This function uses BLUE to calculate the value of a variable for a connection point.
    It combines information from arc observations and previously estimated values for the control points (grondslag).
    Variance-covariance matrices (VCMs) are constructed to ensure proper error propagation

    Parameters
    ----------
    variable : str
        The name of the variable in the STM to estimate (e.g., 'cross_range', 'displacement').
    stm_control : xarray.DataArray
        State Transition Matrix (STM) of the control points, including the reference point.
    arc_variable : numpy.ndarray
        Observations related to the arcs in the network (e.g., arc displacements).
    arc_sigma : numpy.ndarray
        Standard deviations (uncertainties) of the arc observations.
    nr_conn : int
        Number of connection points used for the estimation.
    control_conn_points : list
        List of indices representing the control points connected to the connection point.

    Returns
    -------
    x_hat : numpy.ndarray
        Estimated value of the variable for the connection point.
    Qx_hat : numpy.ndarray
        Variance-covariance matrix of the estimated variable.
    Qyy : numpy.ndarray
        Variance-covariance matrix of the observations.
    Qyy_inv : numpy.ndarray
        Inverse of the variance-covariance matrix (precision matrix).
    A : numpy.ndarray
        Design matrix used in the adjustment process.
    y_obs : numpy.ndarray
        Observation vector combining arc and point variables.
    y_hat : numpy.ndarray
        Adjusted observation vector after the estimation process.
    e_hat : numpy.ndarray
        Residuals between the observed and adjusted observations.

    Notes
    -----
    - The observations vector is constructed by adding the estimated values of the control points
      to the arc observations.
    - The VCM for the observations accounts for uncertainties in both the arc observations and the
      previously estimated control values.
    """
    # variable_variance = variable + "_variance"

    # Select the right values from the stm_control_solved stm
    variable_estimate_control = stm_control.sel(space=control_conn_points)[variable]
    # variable_variance_control = stm_control.sel(space=control_conn_points)[variable_variance]

    # Construct the observations vector
    y_obs = arc_variable + variable_estimate_control.values
    y_obs = np.reshape(y_obs, (len(y_obs), 1))

    # Construct the VCM for the observations that need to be solved.
    # This is the sum of a matrix with the 'estimated' variances for the connection point and a matrix with the
    # variances of the control points on the diagonal
    # The is the same as filling the diagonal of the matrix with the arc variances

    # Since the variance for the cross_range and thermal component is not known, it is calculated
    # variance_connection_point = arc_sigma**2 - variable_variance_control.values

    # Get the variances of the 'connection point' into the full matrix,
    # since the variances of the connection point are causing the covariance terms
    # Qyy = np.ones((nr_conn, nr_conn)) * np.mean(variance_connection_point)

    # Method where we use the identity matrix
    Qyy = np.identity(nr_conn)

    # Get the arc variances on the diagonal
    np.fill_diagonal(Qyy, arc_sigma**2)

    # Compute the inverse
    Qyy_inv = np.linalg.inv(Qyy)

    A = np.ones((nr_conn, 1))

    # Estimate the unknown parameter for the connection point
    x_hat, Qx_hat = est.blue_q_yy_inv(A, y_obs, Qyy_inv)

    # Estimate residue
    y_hat = A @ x_hat
    e_hat = y_obs - y_hat

    return x_hat, Qx_hat, Qyy, Qyy_inv, A, y_obs, y_hat, e_hat


def moving_average(data, window):
    """Calculate a moving average over the data with a specified window size.

    Parameters
    ----------
    data: np.ndarray
        1D array with the data
    window: int
        Windowsize

    Returns
    -------
    np.ndarray
        Moving average of `data` (only the valid part)

    """
    return np.convolve(data, np.ones(window) / window, mode="valid")


def _detect_ambiguous_series2(time_series, window_length=30):
    """Detect and correct time series that differ by a multiple of pi.

    Parameters
    ----------
    time_series : np.ndarray
        Array of shape (N, T), with N time seris of length T.
    window_length : int
        Moving average time window.

    Returns
    -------
    ts_aligned : np.ndarray
        Aligned time series after correction.
    """
    ts_all = time_series
    N, T = ts_all.shape

    # Toegestane shifts in pi
    shift_options = np.array([-2, 0, 2]) * np.pi

    # Smoothe tijdseries
    ts_smoothed = np.array([moving_average(ts, window_length) for ts in ts_all])
    means = np.mean(ts_smoothed, axis=1)

    max_shiftable = N // 2
    best_score = np.inf
    best_shifts = np.zeros(N)

    for k in range(1, max_shiftable + 1):
        for shift_indices in combinations(range(N), k):
            for shift_combo in product(shift_options, repeat=k):
                current_shifts = np.zeros(N)
                for idx, shift in zip(shift_indices, shift_combo, strict=True):
                    current_shifts[idx] = shift

                shifted_means = means + current_shifts
                score = np.std(shifted_means)

                if score < best_score:
                    best_score = score
                    best_shifts = current_shifts.copy()

    # Pas de shifts toe
    ts_aligned = np.array([ts_all[i] + best_shifts[i] for i in range(N)])
    return ts_aligned


def _detect_ambiguous_series(timeseries, threshold=0.85 * 2 * np.pi):
    """Detect ambiguous time series that may require correction due to phase wrapping and applys corrections.

    Sometimes when the time series for an arc is estimated, it is accidentally shifted by +pi or -pi wrt
    the other estimated time series. We need to test this, and potentially shift a time series up or down

    Parameters
    ----------
    timeseries : numpy.ndarray
        A 2D array of shape (n, m), where `n` is the number of time series and `m` is the length of each time series.
    threshold : float, optional
        Tolerance for detecting deviations, default is `0.85 * 2 * np.pi`.

    Returns
    -------
    timeseries : numpy.ndarray
        Corrected time series with adjustments applied to resolve ambiguity.
    mean_values : list of tuple
        A list of tuples representing the mean differences between pairs of time series.
        Each tuple is of the form `(i, j, mean_diff)` where `i` and `j` are indices of the time series compared,
        and `mean_diff` is their mean difference.
    corrections : list of tuple
        A list of tuples specifying the index of the corrected time series and the applied correction value.
        Each tuple is of the form `(index, correction_value)`.
    """
    n_series = timeseries.shape[0]
    mean_values = []
    correction_counts = {i: 0 for i in range(n_series)}  # Counter for correction suggestions

    # Define mean differences for a pair of time series
    for i in range(n_series):
        for j in range(i + 1, n_series):
            verschil = timeseries[i] - timeseries[j]
            mean_diff = np.mean(verschil)
            mean_values.append((i, j, mean_diff))

    for i, j, mean_diff in mean_values:
        if np.abs(mean_diff) > threshold:
            correction_counts[j] += 1
            correction_counts[i] += 1

    # Identify arc that has most correction suggestions
    # That arc time series need to be shifted up or down
    max_correction_point = max(correction_counts, key=correction_counts.get)

    # Define the direction of the correction (up or down)
    corrections = []
    for i, j, mean_diff in mean_values:
        if max_correction_point in (i, j) and np.abs(mean_diff) > threshold:
            if max_correction_point == j and mean_diff < 0:
                if max_correction_point not in corrections:
                    corrections.append((max_correction_point, np.pi))  # j needs to go down
            if max_correction_point == j and mean_diff > 0:
                if max_correction_point not in corrections:
                    corrections.append((max_correction_point, np.pi))  # j needs to go up
            if max_correction_point == i and mean_diff < 0:
                if max_correction_point not in corrections:
                    corrections.append((max_correction_point, np.pi))  # i needs to go down
            elif max_correction_point == i and mean_diff > 0:
                if max_correction_point not in corrections:
                    corrections.append((max_correction_point, np.pi))  # j needs to go up

    if corrections:
        timeseries[max_correction_point, :] += corrections[0][1]

    return timeseries, mean_values, corrections


def _estimate_connection_point_displ_stm_input(
    stm_control, arc_displ, arc_sigma_displ, nr_conn, control_conn_points, t, threshold_Qyy_change=0.005, epsilon=0
):
    """Estimate the displacement of a connection point using the space-time matrix (STM) and arc displacement estimates.

    This method assumes that the arc noise is predominantly attributed to the connection point and considers
    covariances between arc observations.
    The function performs the following steps:
    1. Constructs the displacement vector for the `control network` points using STM.
    2. Forms the observation vector by combining displ. of the arcs and the `control network` points.
    3. Constructs the VCM for the observations using the provided arc variances and the variances of control pnts.
    4. Ensures the VCM is positive-definite by correcting its eigenvalues if necessary.
    5. Uses BLUE to compute the displacement estimate (`x_hat`) and its associated covariance matrix (`Qx_hat`).
    6. Computes the adjusted observation vector (`y_hat`) and the residuals (`e_hat`).

    Parameters
    ----------
    stm_control : xarray.DataArray
        Space-time matrix containing the `control network` points and their estimated time series, variances, and
        orrected displacement values.
    arc_displ : numpy.ndarray
        The estimated displacements for the arcs.
    arc_sigma_displ : numpy.ndarray
        Standard deviations (uncertainties) of the arc displacements.
    nr_conn : int
        Number of connections to the new point.
    control_conn_points : numpy.ndarray
        Indices of the `control network` points to which the new point is connected.
    t : int
        Time epoch for which the displacement is being estimated.
    threshold_Qyy_change : float, optional
        Threshold for acceptable changes to the corrected VCM, by default 0.005.
    epsilon : float, optional
        Smallest allowable eigenvalue for ensuring the VCM is positive-definite, by default 0.

    Returns
    -------
    x_hat : numpy.ndarray
        Estimated displacement of the new connection point (1D array).
    Qx_hat : numpy.ndarray
        Variance-covariance matrix of the estimated displacement (1x1 matrix).
    Qyy : numpy.ndarray
        Variance-covariance matrix of the observations (nxn matrix, where `n=nr_conn`).
    Qyy_inv : numpy.ndarray
        Inverse (or precision matrix) of the variance-covariance matrix of the observations (nxn matrix).
    A : numpy.ndarray
        Design matrix linking the observations to the unknown displacement (nx1 matrix).
    y_obs : numpy.ndarray
        Observation vector for the displacements (nx1 array).
    y_hat : numpy.ndarray
        Adjusted observation vector (nx1 array).
    e_hat : numpy.ndarray
        Residuals between the observed and adjusted displacements (nx1 array).

    Raises
    ------
    ValueError
        If the correction applied to the VCM (`Qyy`) exceeds the specified `threshold_Qyy_change`.

    Notes
    -----
    - The function ensures numerical stability by correcting VCM to be positive-definite if eigenvalues are negative.
    - The corrected VCM is validated to ensure that the changes are within the specified threshold.
    """
    # Extract the values from the stm of the control points
    time_series_control = stm_control.sel(space=control_conn_points)["displ_time_series_corrected"]
    time_series_variance_control = stm_control.sel(space=control_conn_points)["displ_time_series_variances"]

    # Extract the right epoch for the displacement adjustment
    displ_control = time_series_control.values[:, t]
    displ_variance_control = time_series_variance_control.values[:, t]

    # Construct the observations vector which is the sum of the arc estimations and the control points
    y_obs = arc_displ + displ_control
    y_obs = np.reshape(y_obs, (len(y_obs), 1))

    # Construct the VCM for the observations that need to be solved.
    # This is the sum of a matrix with the 'estimated' variances for the new connection point and an identity matrix
    # with the variances of the control points on the diagonal
    # The is the same as filling the diagonal of the matrix with the arc variances

    # We estimate the variance of the new point
    variance_connection_point = arc_sigma_displ**2 - displ_variance_control

    # Get the variances of the 'new point' into the full matrix, since the variances of the new point are
    # causing the covariance terms
    Qyy = np.ones((nr_conn, nr_conn)) * np.mean(variance_connection_point)
    # Get the arc variances on the diagonal
    np.fill_diagonal(Qyy, arc_sigma_displ**2)

    # Compute eigenvectors and eigenvalues to check whether matrix is positive-definite, if the matrix is
    # positive-definite we need to correct it
    eigvals, eigvecs = np.linalg.eigh(Qyy)
    if np.any(eigvals < 0):
        # print(f"There are non negative eigenvalues at epoch {t}: Eigenvalues: {eigvals}")
        print(f"There are non negative eigenvalues at epoch {t}: Qyy will be corrected...")

        eigvals[eigvals < epsilon] = epsilon
        Qyy_corrected = eigvecs @ np.diag(eigvals) @ eigvecs.T
        # Check the change
        change = np.linalg.norm(Qyy_corrected - Qyy, ord="fro")
        if change < threshold_Qyy_change:
            Qyy = Qyy_corrected
            condition_number = np.linalg.cond(Qyy)
            Qyy_inv = np.linalg.pinv(Qyy) if condition_number > 1e12 else np.linalg.inv(Qyy)
        else:
            raise ValueError(f"Qyy correction too large: {change} > {threshold_Qyy_change}, happens at epoch {t}")

    else:
        Qyy_inv = np.linalg.inv(Qyy)

    # Design matrix
    A = np.ones((nr_conn, 1))

    # Estimate the unknown parameter for the new point
    x_hat, Qx_hat = est.blue_q_yy_inv(A, y_obs, Qyy_inv)

    # Estimate residue
    y_hat = A @ x_hat
    e_hat = y_obs - y_hat

    return x_hat, Qx_hat, Qyy, Qyy_inv, A, y_obs, y_hat, e_hat


def _estimate_connection_point_full_phase_stm_input(
    stm_control, arc_displ, arc_sigma_displ, nr_conn, control_conn_points, t, threshold_Qyy_change=0.005, epsilon=0
):
    """Estimate the displacement of a connection point using the space-time matrix (STM) and arc phases.

    This method assumes that the arc noise is predominantly attributed to the connection point and considers
    covariances between arc observations.
    The function performs the following steps:
    1. Constructs the displacement vector for the `control network` points using STM.
    2. Forms the observation vector by combining displ. of the arcs and the `control network` points.
    3. Constructs the VCM for the observations using the provided arc variances and the variances of control pnts.
    4. Ensures the VCM is positive-definite by correcting its eigenvalues if necessary.
    5. Uses BLUE to compute the displacement estimate (`x_hat`) and its associated covariance matrix (`Qx_hat`).
    6. Computes the adjusted observation vector (`y_hat`) and the residuals (`e_hat`).

    Parameters
    ----------
    stm_control : xarray.DataArray
        Space-time matrix containing the `control network` points and their estimated time series, variances, and
        orrected displacement values.
    arc_displ : numpy.ndarray
        The estimated displacements for the arcs.
    arc_sigma_displ : numpy.ndarray
        Standard deviations (uncertainties) of the arc displacements.
    nr_conn : int
        Number of connections to the new point.
    control_conn_points : numpy.ndarray
        Indices of the `control network` points to which the new point is connected.
    t : int
        Time epoch for which the displacement is being estimated.
    threshold_Qyy_change : float, optional
        Threshold for acceptable changes to the corrected VCM, by default 0.005.
    epsilon : float, optional
        Smallest allowable eigenvalue for ensuring the VCM is positive-definite, by default 0.

    Returns
    -------
    x_hat : numpy.ndarray
        Estimated displacement of the new connection point (1D array).
    Qx_hat : numpy.ndarray
        Variance-covariance matrix of the estimated displacement (1x1 matrix).
    Qyy : numpy.ndarray
        Variance-covariance matrix of the observations (nxn matrix, where `n=nr_conn`).
    Qyy_inv : numpy.ndarray
        Inverse (or precision matrix) of the variance-covariance matrix of the observations (nxn matrix).
    A : numpy.ndarray
        Design matrix linking the observations to the unknown displacement (nx1 matrix).
    y_obs : numpy.ndarray
        Observation vector for the displacements (nx1 array).
    y_hat : numpy.ndarray
        Adjusted observation vector (nx1 array).
    e_hat : numpy.ndarray
        Residuals between the observed and adjusted displacements (nx1 array).

    Raises
    ------
    ValueError
        If the correction applied to the VCM (`Qyy`) exceeds the specified `threshold_Qyy_change`.

    Notes
    -----
    - The function ensures numerical stability by correcting VCM to be positive-definite if eigenvalues are negative.
    - The corrected VCM is validated to ensure that the changes are within the specified threshold.
    """
    # Extract the values from the stm of the control points
    time_series_control = stm_control.sel(space=control_conn_points)["displ_time_series_corrected"]
    cr_phase_timeseries = stm_control.sel(space=control_conn_points)["cross_range_phase_timeseries"]
    thermal_phase_timeseries = stm_control.sel(space=control_conn_points)["thermal_phase_timeseries"]
    time_series_variance_control = stm_control.sel(space=control_conn_points)["displ_time_series_variances"]

    # Extract the right epoch for the displacement adjustment
    # displ_control = time_series_control.values[:, t] # Old version
    phases_control = (
        time_series_control.values[:, t] + thermal_phase_timeseries.values[:, t] + cr_phase_timeseries.values[:, t]
    )
    displ_variance_control = time_series_variance_control.values[:, t]

    # Construct the observations vector which is the sum of the arc estimations and the control points
    y_obs = arc_displ + phases_control
    y_obs = np.reshape(y_obs, (len(y_obs), 1))

    # Construct the VCM for the observations that need to be solved.
    # This is the sum of a matrix with the 'estimated' variances for the new connection point and an identity matrix
    # with the variances of the control points on the diagonal
    # The is the same as filling the diagonal of the matrix with the arc variances

    # We estimate the variance of the new point
    variance_connection_point = arc_sigma_displ**2 - displ_variance_control

    # Get the variances of the 'new point' into the full matrix, since the variances of the new point are
    # causing the covariance terms
    Qyy = np.ones((nr_conn, nr_conn)) * np.mean(variance_connection_point)
    # Get the arc variances on the diagonal
    np.fill_diagonal(Qyy, arc_sigma_displ**2)

    # Compute eigenvectors and eigenvalues to check whether matrix is positive-definite, if the matrix is
    # positive-definite we need to correct it
    eigvals, eigvecs = np.linalg.eigh(Qyy)
    if np.any(eigvals < 0):
        # print(f"There are non negative eigenvalues at epoch {t}: Eigenvalues: {eigvals}")
        print("There are non negative eigenvalues at epoch {t}: Qyy will be corrected...")

        eigvals[eigvals < epsilon] = epsilon
        Qyy_corrected = eigvecs @ np.diag(eigvals) @ eigvecs.T
        # Check the change
        change = np.linalg.norm(Qyy_corrected - Qyy, ord="fro")
        if change < threshold_Qyy_change:
            Qyy = Qyy_corrected
            condition_number = np.linalg.cond(Qyy)
            Qyy_inv = np.linalg.pinv(Qyy) if condition_number > 1e12 else np.linalg.inv(Qyy)
        else:
            raise ValueError(f"Qyy correction too large: {change} > {threshold_Qyy_change}, happens at epoch {t}")

    else:
        Qyy_inv = np.linalg.inv(Qyy)

    # Design matrix
    A = np.ones((nr_conn, 1))

    # Estimate the unknown parameter for the new point
    x_hat, Qx_hat = est.blue_q_yy_inv(A, y_obs, Qyy_inv)

    # Estimate residue
    y_hat = A @ x_hat
    e_hat = y_obs - y_hat

    return x_hat, Qx_hat, Qyy, Qyy_inv, A, y_obs, y_hat, e_hat


def _estimate_connection_point_outlier(y_corrected, A, Qyy_inv):
    """Estimate the displacement for a connection point using a corrected observation vector.

    This function applies BLUE method to estimate the
    unknown displacement parameter for a connection point, along with its variance-covariance
    matrix, the adjusted observation vector, and the residuals.

    Parameters
    ----------
    y_corrected : numpy.ndarray
        observation vector
    A : numpy.ndarray
        The design matrix linking the observations to the unknown parameter(s).
    Qyy_inv : numpy.ndarray
        The inverse of the variance-covariance matrix of the observations.

    Returns
    -------
    x_hat : numpy.ndarray
        Estimated value of the unknown parameter (e.g., displacement) for the connection point.
    Qx_hat : numpy.ndarray
        Variance-covariance matrix of the estimated parameter(s).
    y_hat : numpy.ndarray
        Adjusted observation vector based on the estimated parameter(s).
    e_hat : numpy.ndarray
        Residuals, calculated as the difference between the corrected observation vector
        (`y_corrected`) and the adjusted observation vector (`y_hat`).
    """
    # Estimate the unknown parameter for the new point
    x_hat, Qx_hat = est.blue_q_yy_inv(A, y_corrected, Qyy_inv)

    # Estimate residue
    y_hat = A @ x_hat
    e_hat = y_corrected - y_hat

    return x_hat, Qx_hat, y_hat, e_hat


def connect_point_to_control_network(
    stm_1_point,
    stm_control,
    stm_ref_control_solved,
    ref_pnt_idx,
    dist_to_quality,
    nr_conn,
    sigma_post_over_sigma_prior,
    alpha,
    bounds,
    m2ph,
    n_max_iter,
    correct_epochs_arc,
):
    """Estimate parameters for a new point relative to the control network points.

    This function performs the following steps:
    1. Identifies the most suitable arcs between the new point (`point_add`) and the control network.
    2. Estimates parameters for these arcs. A time limit can be set for this step, as the estimation involves solving
       a non-linear problem that can be computationally expensive.
    3. Adjusts the new point based on the estimated parameters for the arcs and estimates its parameters relative to
       the control points:
       - Handles cases where the Qyy matrix is not positive-definite by modifying the matrix to ensure solvability.
       - Adjustments are performed for the cross-range, thermal component, and displacement values per epoch.
       - If the T-OMT test statistic exceeds the critical value, a w-test is applied to detect incorrectly unwrapped
         phases for specific arcs and epochs. These arcs are rewrapped as needed.
    4. Stores the results in a dictionary and STM format.

    Parameters
    ----------
    stm_1_point : xarray.DataArray
        STM of the point to be estimated, containint SLC phase values and more.
    stm_control : xarray.DataArray
        STM of the control network points (including the reference point), containing SLC phase values.
    stm_ref_control_solved : xarray.DataArray
        STM of the control network points (including the reference point) with estimated values, such as cross-range,
        thermal component, and displacement time series.
    ref_pnt_idx : integer
        index of the reference point
    dist_to_quality : numpy.ndarray
        Value representing additional variability in the stochastic model for longer arcs.
    nr_conn : numpy.ndarray
        Number of arcs that the new point needs to connect to within the control network.
    sigma_post_over_sigma_prior:
        Upper limit on the allowed change from apriori sigma to aposteriori sigma before an arc is rejected
    alpha : float
        Level of significance for overall model test
    bounds : tuple
        Bounds for parameter estimation in the complex domain. Example format:
        (A_lower, a_lower, b_lower, c_lower, H_lower, therm_lower,
        A_upper, a_upper, b_upper, c_upper, H_upper, therm_upper),
        where:
        - `A` is the amplitude,
        - `a`, `b`, `c` are polynomial parameters,
        - `CR` is the cross-range,
        - `therm` is the thermal component.
    m2ph : float
        Factor relating the phase to meters, which varies per mission.
    n_max_iter : np.ndarray
        The maximum nr of iterations for non-linear lsq per arc
    correct_epochs_arc : int
        Specify whether to re-wrap time series if the T-OMT exceeds critical values (1 = yes, 0 = no).

    Returns
    -------
    stm_point_solved : xarray.DataArray
        STM for the point with estimated values.
    stm_arc_results : xarray.DataArray
        STM containing results for each arc.
    estimated_values_pnt_add : dict
        Dictionary with estimation results for the new point.
    arc_results_pnt_add : dict
        Dictionary with results for the arcs.
    solved_the_point : int
        Indicates whether the point could be solved (1 = yes, 0 = no).
    arcs_closing_variable : numpy.ndarray
        Closing variable for arcs during estimation.
    sigma_phases_arc : numpy.ndarray
        Variance of phase values for the arcs.
    """
    # Make an empty dictionary where the estimation results for the point will be stored
    point_add = int(stm_1_point["pnt_idx"].values)
    estimated_values_pnt_add = {}
    arc_results_pnt_add = {}

    nr_epochs = len(stm_1_point["time"])

    # Compute the ordered-arcs between the connection point and control points
    arcs, sorted_quality_values = dn.ordered_arcs_connection_point_and_control_network(
        stm_1_point["rd_x"],
        stm_1_point["rd_y"],
        stm_1_point["slc_quality"],
        stm_1_point["pnt_idx"].astype(int),
        stm_control["rd_x"],
        stm_control["rd_y"],
        stm_control["slc_quality"],
        stm_control["pnt_idx"].astype(int),
        dist_to_quality,
    )

    # To count how many 'succesfull' estimations of the connection points we have
    succes_arcs = np.full(len(arcs), False)
    arc_results = {}

    for i, arc_add in enumerate(arcs):
        # STM of the control points still consists of all points, so only take out the control point that is in the arc
        stm_control_1_point = stm_control.where(stm_control["pnt_idx"].isin(arc_add[0]), drop=True)
        stm_control_1_point = stm_control_1_point.squeeze()

        # Estimate the unknown parameters for the arc.
        # This occurs within a wrapper that makes sure it does not take too much time
        # arc_results_1_arc = _run_with_timeout(
        #     #arc_est.arc_estimation_xarray_input,
        #     arc_estimation_xarray_input_v2,
        #     max_time_arc_estimation,  # The maximum allowed time in seconds for this function
        #     stm_control_1_point,
        #     stm_1_point,
        #     bounds,
        #     m2ph,
        #     test_stochastics=0,
        #     print_output=0,
        # )
        arc_results_1_arc = arc_est.arc_estimation_xarray_input(
            stm_control_1_point,
            stm_1_point,
            bounds,
            m2ph,
            n_max_iter,
            test_stochastics=0,
            print_output=0,
        )

        if arc_results_1_arc is None:
            print(f"Computation for {arc_add} failed")
        else:
            if "succeeded_arcs" in arc_results_1_arc and not np.isnan(arc_results_1_arc["succeeded_arcs"]).any():
                # We need to test whether the arc solution that we found matches the a priori quality
                # Only if the posterior sigma is not to much bigger than the prior sigma, the arc is considered
                est_displ_phase = (
                    arc_results_1_arc["unwrap_phases_arc"]
                    - arc_results_1_arc["estimated_cross_range_phase"]
                    - arc_results_1_arc["estimated_thermal_phase"]
                )
                displ_phase = arc_results_1_arc["estimated_displ_phase"]
                residues_per_arc = est_displ_phase - displ_phase
                sigma_post_arc = np.std(residues_per_arc)
                mean_sigma_prior_arc = np.mean(arc_results_1_arc["sigma_phases_arc"])

                if sigma_post_arc < sigma_post_over_sigma_prior * mean_sigma_prior_arc:
                    succes_arcs[i] = True
                    print(f"Computation succesfull for {arc_add}")
                    # Add results for this arc to dictionary
                    for key, value in arc_results_1_arc.items():
                        if key in arc_results:
                            arc_results[key] = (
                                np.vstack([arc_results[key], value])
                                if isinstance(value, np.ndarray)
                                else arc_results[key] + [value]
                            )
                        else:
                            arc_results[key] = [value] if not isinstance(value, np.ndarray) else value
                else:
                    print(f"Solution for arc {arc_add} is too noisey")
            else:
                print(f"No solution was found for {arc_add}")

        if np.sum(succes_arcs) >= nr_conn:
            print("There are enough successful arcs")
            break

    if np.sum(succes_arcs) < nr_conn:
        print(f"There were not enough arcs that could be solved for point {point_add}")
        solved_the_point = 0
        stm_point_solved = []
        stm_arc_results = []
        estimated_values_pnt_add = []
        arc_results_pnt_add = []

    else:
        # Get information from the dictionaries
        unwrap_phases_arc = arc_results["unwrap_phases_arc"]
        sigma_phases_arc = arc_results["sigma_phases_arc"]
        estimated_cross_range_arc = arc_results["estimated_cross_range"].flatten()
        estimated_thermal_arc = arc_results["estimated_thermal"].flatten()
        estimated_cross_range_sigma_arc = arc_results["estimated_cross_range_sigma"].flatten()
        estimated_thermal_sigma_arc = arc_results["estimated_thermal_sigma"].flatten()
        arcs_closing_variable = (
            unwrap_phases_arc - arc_results["estimated_cross_range_phase"] - arc_results["estimated_thermal_phase"]
        )  # Version where we test on the displacement phase
        # arcs_closing_variable = (
        #     unwrap_phases_arc
        # ) # Version where we close on the full phase

        # Get the succeeded arcs
        arcs_connection_point = arc_results["succeeded_arcs"]
        control_conn_points = arcs_connection_point[:, 0]

        print("arcs")
        print(arcs_connection_point)
        print("Control network connection points")
        print(control_conn_points)

        # Apply OMT for cross_range and thermal component
        m_omt = nr_conn  # The nr of observations is defined by the nr of arcs
        n_omt = 1  # value for OMT is always equal to 1

        # Estimate time series for the point
        point_time_series = np.zeros(nr_epochs)
        point_time_series_corrected = np.zeros(nr_epochs)
        point_time_series_sigma = np.zeros(nr_epochs)
        T_omt_displ = np.zeros(nr_epochs)
        T_omt_displ_correction = np.zeros(nr_epochs)

        y_arcs = np.zeros((nr_conn, nr_epochs))
        y_arcs_corrected = np.zeros((nr_conn, nr_epochs))
        y_hat_arcs = np.zeros((nr_conn, nr_epochs))
        y_hat_arcs_corrected = np.zeros((nr_conn, nr_epochs))

        count_omt_reject = 0
        arc_idx_reject = np.full(nr_epochs, np.nan)

        # Test whether we need to correct a time series by 2pi or not
        # Sometimes one of the arc will have the timeseries one full cycle above the other two
        arcs_closing_variable = _detect_ambiguous_series2(arcs_closing_variable)
        print("Check whether we need to shift a time series ")

        for t in range(nr_epochs):
            # Estimate displacement for the new point
            point_epoch, Qx_epoch, Qyy, Qyy_inv, A, y_obs, y_hat_obs, e_hat_epoch = (
                _estimate_connection_point_displ_stm_input(
                    # _estimate_connection_point_full_phase_stm_input(
                    stm_ref_control_solved,
                    arcs_closing_variable[:, t],
                    sigma_phases_arc[:, t],
                    nr_conn,
                    control_conn_points,
                    t,
                )
            )

            point_time_series[t] = point_epoch[0, 0]
            point_time_series_sigma[t] = np.sqrt(np.abs(Qx_epoch[0, 0]))
            point_time_series_corrected[t] = np.copy(point_time_series[t])
            y_arcs[:, t] = y_obs.flatten()
            y_hat_arcs[:, t] = y_hat_obs.flatten()
            y_arcs_corrected[:, t] = np.copy(y_arcs[:, t])
            y_hat_arcs_corrected[:, t] = np.copy(y_hat_arcs[:, t])

            # Apply OMT
            k, T_omt_displ[t] = est.overall_model_test(alpha, e_hat_epoch, Qyy_inv, m_omt - n_omt, 0)
            T_omt_displ_correction[t] = np.copy(T_omt_displ[t])

            # Apply w-test if the OMT is rejected
            if correct_epochs_arc == 1 and T_omt_displ[t] > k:
                count_omt_reject += 1

                corrected_y, idx_biggest_w, _ = apply_w_test_control_network(
                    m_omt, A, Qx_epoch, Qyy, Qyy_inv, y_obs, e_hat_epoch
                )

                # Re-estimate the displacement with the corrected y
                point_epoch, _, y_hat_corr, e_hat_epoch_corr = _estimate_connection_point_outlier(
                    corrected_y, A, Qyy_inv
                )
                # # Compute the OMT for this particular epoch
                k, T_omt_temp = est.overall_model_test(alpha, e_hat_epoch_corr, Qyy_inv, m_omt - n_omt, 0)

                # if T_omt_temp < T_omt_displ[t] and T_omt_temp < k:
                if T_omt_temp < T_omt_displ[t]:  # correct if the new T is better than the old one
                    T_omt_displ_correction[t] = T_omt_temp

                    y_hat_arcs_corrected[:, t] = y_hat_corr.flatten()
                    point_time_series_corrected[t] = point_epoch[0, 0]
                    y_arcs_corrected[:, t] = corrected_y.flatten()
                    arc_idx_reject[t] = idx_biggest_w

                # TODO possible to be added that we consider other arcs if the omt is rejected over and over again

        # if there are more thatn 20% omt displ reject count which arc is causing the issue
        # counts_bad_arc = np.array([(arc_idx_reject == i).sum() for i in range(nr_conn + 1)])
        count_bad_arcs = np.array([(arc_idx_reject == i).sum() for i in range(nr_conn)])
        print(f"Counts of the bad arcs {count_bad_arcs}")
        print(f"The worst arc is therefore {np.argmax(count_bad_arcs)}")

        # Estimate the cross range and thermal component
        # Compute the thermal component and cross range for the point
        (
            cross_range_pnt,
            Qx_cross_range,
            _,
            Qyy_cross_range_inv,
            A,
            _,
            _,
            e_cross_range,
        ) = _estimate_connection_point_stm(
            "cross_range",
            stm_ref_control_solved,
            estimated_cross_range_arc,
            estimated_cross_range_sigma_arc,
            nr_conn,
            control_conn_points,
        )
        thermal_pnt, Qx_thermal, _, Qyy_thermal_inv, A, _, _, e_thermal = _estimate_connection_point_stm(
            "thermal_comp",
            stm_ref_control_solved,
            estimated_thermal_arc,
            estimated_thermal_sigma_arc,
            nr_conn,
            control_conn_points,
        )

        k, T_omt_cross_range = est.overall_model_test(alpha, e_cross_range, Qyy_cross_range_inv, m_omt - n_omt, 0)
        k, T_omt_thermal = est.overall_model_test(alpha, e_thermal, Qyy_thermal_inv, m_omt - n_omt, 0)

        estimated_values_pnt_add["T_omt_cross_range"] = T_omt_cross_range
        estimated_values_pnt_add["T_omt_thermal"] = T_omt_thermal
        estimated_values_pnt_add["k_omt"] = k

        estimated_values_pnt_add["cross_range"] = cross_range_pnt.flatten()[0]
        estimated_values_pnt_add["cross_range_sigma"] = np.sqrt(Qx_cross_range.flatten()[0])

        estimated_values_pnt_add["thermal_comp"] = thermal_pnt.flatten()[0]
        estimated_values_pnt_add["thermal_comp_sigma"] = np.sqrt(Qx_thermal.flatten()[0])

        # Estimate the mean a priori an posteriori sigma values per arc
        mean_priori_quality = np.mean(sigma_phases_arc, axis=1)
        arc_residues_displacement = y_arcs_corrected - y_hat_arcs_corrected
        sigma_residues = np.std(arc_residues_displacement, axis=1)

        # Add the estimations to the dictionaries for the point (estimated result) and arc estimated
        estimated_values_pnt_add["time_series"] = point_time_series
        estimated_values_pnt_add["time_series_corrected"] = point_time_series_corrected
        estimated_values_pnt_add["point_time_series_sigma"] = point_time_series_sigma
        estimated_values_pnt_add["T_omt_displ"] = T_omt_displ
        estimated_values_pnt_add["T_omt_displ_corrected"] = T_omt_displ_correction
        estimated_values_pnt_add["connection points"] = control_conn_points
        estimated_values_pnt_add["omt_reject"] = count_omt_reject

        arc_results_pnt_add["arcs"] = arcs_connection_point
        arc_results_pnt_add["y_arcs"] = y_arcs
        arc_results_pnt_add["y_arcs_corrected"] = y_arcs_corrected
        arc_results_pnt_add["y_hat_arcs"] = y_hat_arcs
        arc_results_pnt_add["y_hat_arcs_corrected"] = y_hat_arcs_corrected
        arc_results_pnt_add["sigma arc"] = sigma_phases_arc
        arc_results_pnt_add["sigma_residues"] = sigma_residues
        arc_results_pnt_add["mean_priori_quality"] = mean_priori_quality

        # Create stm for the arcs
        # Create new stm where we solve the estimated values
        # The stm_solved will have the same dimensions as the stm that we already have
        stm_arc_results = xr.Dataset(coords={"space": np.ones(nr_conn) * point_add, "time": stm_1_point.time.values})

        stm_arc_results["pnt_idx"] = (["space"], control_conn_points)
        stm_arc_results["y_arcs"] = (["space", "time"], y_arcs)
        stm_arc_results["y_arcs_corrected"] = (["space", "time"], y_arcs_corrected)
        stm_arc_results["y_hat_arcs"] = (["space", "time"], y_hat_arcs)
        stm_arc_results["y_hat_arcs_corrected"] = (["space", "time"], y_hat_arcs_corrected)
        stm_arc_results["sigma arc"] = (["space", "time"], sigma_phases_arc)
        stm_arc_results["sigma_residues"] = (["space"], sigma_residues)
        stm_arc_results["mean_priori_quality"] = (["space"], mean_priori_quality)

        # Create new stm where we solve the estimated values
        # The stm_solved will have the same dimensions as the stm that we already have
        stm_point_solved = xr.Dataset(coords={"space": [point_add], "time": stm_1_point.time.values})

        # Add the values to the STM
        stm_point_solved["ref_pnt"] = (["space"], [ref_pnt_idx])
        stm_point_solved["control_or_not"] = (["space"], [0])
        stm_point_solved["conn_points"] = (["space"], [0])
        stm_point_solved["pnt_idx"] = (["space"], [int(point_add)])

        estimated_values_pnt_add["T_omt_cross_range"] = T_omt_cross_range

        stm_point_solved["cross_range"] = (["space"], [cross_range_pnt.flatten()[0]])
        stm_point_solved["cross_range_variance"] = (["space"], [Qx_cross_range.flatten()[0]])
        stm_point_solved["thermal_comp"] = (["space"], [thermal_pnt.flatten()[0]])
        stm_point_solved["thermal_comp_variance"] = (["space"], [Qx_thermal.flatten()[0]])
        stm_point_solved["omt_cross_range"] = (["space"], [T_omt_cross_range])
        stm_point_solved["omt_thermal"] = (["space"], [T_omt_thermal])
        stm_point_solved["k_omt"] = (["space"], [k])

        stm_point_solved["omt_reject"] = (["space"], [count_omt_reject])
        stm_point_solved["omt_displ"] = (["space", "time"], T_omt_displ.reshape(1, nr_epochs))
        stm_point_solved["omt_displ_corrected"] = (["space", "time"], T_omt_displ_correction.reshape(1, nr_epochs))

        stm_point_solved["displ_time_series"] = (["space", "time"], point_time_series.reshape(1, nr_epochs))
        stm_point_solved["displ_time_series_corrected"] = (
            ["space", "time"],
            point_time_series_corrected.reshape(1, nr_epochs),
        )
        stm_point_solved["displ_time_series_variances"] = (
            ["space", "time"],
            (point_time_series_sigma**2).reshape(1, nr_epochs),
        )

        solved_the_point = 1

    return stm_point_solved, stm_arc_results, estimated_values_pnt_add, arc_results_pnt_add, solved_the_point


def connect_point_to_control_network_full_phase(
    stm_1_point,
    stm_control,
    stm_ref_control_solved,
    ref_pnt_idx,
    dist_to_quality,
    nr_conn,
    sigma_post_over_sigma_prior,
    alpha,
    bounds,
    m2ph,
    n_max_iter,
    correct_epochs_arc,
):
    """Estimate parameters for a new point relative to the control network points.

    This function performs the following steps:
    1. Identifies the most suitable arcs between the new point (`point_add`) and the control network.
    2. Estimates parameters for these arcs. A time limit can be set for this step, as the estimation involves solving
       a non-linear problem that can be computationally expensive.
    3. Adjusts the new point based on the estimated parameters for the arcs and estimates its parameters relative to
       the control points:
       - Handles cases where the Qyy matrix is not positive-definite by modifying the matrix to ensure solvability.
       - Adjustments are performed for the cross-range, thermal component, and displacement values per epoch.
       - If the T-OMT test statistic exceeds the critical value, a w-test is applied to detect incorrectly unwrapped
         phases for specific arcs and epochs. These arcs are rewrapped as needed.
    4. Stores the results in a dictionary and STM format.

    Parameters
    ----------
    stm_1_point : xarray.DataArray
        STM of the point to be estimated, containint SLC phase values and more.
    stm_control : xarray.DataArray
        STM of the control network points (including the reference point), containing SLC phase values.
    stm_ref_control_solved : xarray.DataArray
        STM of the control network points (including the reference point) with estimated values, such as cross-range,
        thermal component, and displacement time series.
    ref_pnt_idx : integer
        index of the reference point
    dist_to_quality : numpy.ndarray
        Value representing additional variability in the stochastic model for longer arcs.
    nr_conn : numpy.ndarray
        Number of arcs that the new point needs to connect to within the control network.
    sigma_post_over_sigma_prior:
        Upper limit on the allowed change from apriori sigma to aposteriori sigma before an arc is rejected
    alpha : float
        Level of significance for overall model test
    bounds : tuple
        Bounds for parameter estimation in the complex domain. Example format:
        (A_lower, a_lower, b_lower, c_lower, H_lower, therm_lower,
        A_upper, a_upper, b_upper, c_upper, H_upper, therm_upper),
        where:
        - `A` is the amplitude,
        - `a`, `b`, `c` are polynomial parameters,
        - `CR` is the cross-range,
        - `therm` is the thermal component.
    m2ph : float
        Factor relating the phase to meters, which varies per mission.
    n_max_iter : np.ndarray
        The maximum nr of iterations for non-linear lsq per arc
    correct_epochs_arc : int
        Specify whether to re-wrap time series if the T-OMT exceeds critical values (1 = yes, 0 = no).

    Returns
    -------
    stm_point_solved : xarray.DataArray
        STM for the point with estimated values.
    stm_arc_results : xarray.DataArray
        STM containing results for each arc.
    estimated_values_pnt_add : dict
        Dictionary with estimation results for the new point.
    arc_results_pnt_add : dict
        Dictionary with results for the arcs.
    solved_the_point : int
        Indicates whether the point could be solved (1 = yes, 0 = no).
    arcs_closing_variable : numpy.ndarray
        Closing variable for arcs during estimation.
    sigma_phases_arc : numpy.ndarray
        Variance of phase values for the arcs.
    """
    # Make an empty dictionary where the estimation results for the point will be stored
    point_add = int(stm_1_point["pnt_idx"].values)
    solved_the_point = 0
    estimated_values_pnt_add = {}
    arc_results_pnt_add = {}
    cr2ph_arc = stm_1_point["cr2ph"].values
    temperature = stm_1_point["temperature"].values
    years = stm_1_point["years"].values

    nr_epochs = len(stm_1_point["time"])

    # Compute the ordered-arcs between the connection point and control points
    arcs, sorted_quality_values = dn.ordered_arcs_connection_point_and_control_network(
        stm_1_point["rd_x"],
        stm_1_point["rd_y"],
        stm_1_point["slc_quality"],
        stm_1_point["pnt_idx"].astype(int),
        stm_control["rd_x"],
        stm_control["rd_y"],
        stm_control["slc_quality"],
        stm_control["pnt_idx"].astype(int),
        dist_to_quality,
    )

    # To count how many 'succesfull' estimations of the connection points we have
    succes_arcs = np.full(len(arcs), False)
    arc_results = {}

    for i, arc_add in enumerate(arcs):
        # STM of the control points still consists of all points, so only take out the control point that is in the arc
        stm_control_1_point = stm_control.where(stm_control["pnt_idx"].isin(arc_add[0]), drop=True)
        stm_control_1_point = stm_control_1_point.squeeze()

        # Estimate the unknown parameters for the arc.
        # This occurs within a wrapper that makes sure it does not take too much time
        # arc_results_1_arc = _run_with_timeout(
        #     #arc_est.arc_estimation_xarray_input,
        #     arc_estimation_xarray_input_v2,
        #     max_time_arc_estimation,  # The maximum allowed time in seconds for this function
        #     stm_control_1_point,
        #     stm_1_point,
        #     bounds,
        #     m2ph,
        #     test_stochastics=0,
        #     print_output=0,
        # )
        arc_results_1_arc = arc_est.arc_estimation_xarray_input(
            stm_control_1_point,
            stm_1_point,
            bounds,
            m2ph,
            n_max_iter,
            test_stochastics=0,
            print_output=0,
        )

        if arc_results_1_arc is None:
            print(f"Computation for {arc_add} failed")
        else:
            if "succeeded_arcs" in arc_results_1_arc and not np.isnan(arc_results_1_arc["succeeded_arcs"]).any():
                # We need to test whether the arc solution that we found matches the a priori quality
                # Only if the posterior sigma is not to much bigger than the prior sigma, the arc is considered
                est_displ_phase = (
                    arc_results_1_arc["unwrap_phases_arc"]
                    - arc_results_1_arc["estimated_cross_range_phase"]
                    - arc_results_1_arc["estimated_thermal_phase"]
                )
                displ_phase = arc_results_1_arc["estimated_displ_phase"]
                residues_per_arc = est_displ_phase - displ_phase
                sigma_post_arc = np.std(residues_per_arc)
                mean_sigma_prior_arc = np.mean(arc_results_1_arc["sigma_phases_arc"])

                if sigma_post_arc < sigma_post_over_sigma_prior * mean_sigma_prior_arc:
                    succes_arcs[i] = True
                    print(f"Computation succesfull for {arc_add}")
                    # Add results for this arc to dictionary
                    for key, value in arc_results_1_arc.items():
                        if key in arc_results:
                            arc_results[key] = (
                                np.vstack([arc_results[key], value])
                                if isinstance(value, np.ndarray)
                                else arc_results[key] + [value]
                            )
                        else:
                            arc_results[key] = [value] if not isinstance(value, np.ndarray) else value
                else:
                    print(f"Solution for arc {arc_add} is too noisey")
            else:
                print(f"No solution was found for {arc_add}")

        if np.sum(succes_arcs) >= nr_conn:
            print("There are enough successful arcs")
            break

    if np.sum(succes_arcs) < nr_conn:
        print(f"There were not enough arcs that could be solved for point {point_add}")
        solved_the_point = 0
        stm_point_solved = []
        stm_arc_results = []
        estimated_values_pnt_add = []
        arc_results_pnt_add = []

    else:
        # Get information from the dictionaries
        unwrap_phases_arc = arc_results["unwrap_phases_arc"]
        sigma_phases_arc = arc_results["sigma_phases_arc"]
        estimated_cross_range_arc = arc_results["estimated_cross_range"].flatten()
        estimated_thermal_arc = arc_results["estimated_thermal"].flatten()
        # estimated_cross_range_sigma_arc = arc_results["estimated_cross_range_sigma"].flatten()
        # estimated_thermal_sigma_arc = arc_results["estimated_thermal_sigma"].flatten()
        # arcs_closing_variable = (
        #     unwrap_phases_arc - arc_results["estimated_cross_range_phase"] - arc_results["estimated_thermal_phase"]
        # )   # Version where we test on the displacement phase
        arcs_closing_variable = unwrap_phases_arc  # Version where we close on the full phase

        # Get the succeeded arcs
        arcs_connection_point = arc_results["succeeded_arcs"]
        control_conn_points = arcs_connection_point[:, 0]

        print("arcs")
        print(arcs_connection_point)
        print("Grondslag connection points")
        print(control_conn_points)

        # Define settings for the OMT
        m_omt = nr_conn  # The nr of observations is defined by the nr of arcs
        n_omt = 1  # value for OMT is always equal to 1

        # # Compute the thermal component and cross range for the point
        # (
        #     cross_range_pnt,
        #     Qx_cross_range,
        #     _,
        #     Qyy_cross_range_inv,
        #     A,
        #     _,
        #     _,
        #     e_cross_range,
        # ) = _estimate_connection_point_stm(
        #     "cross_range",
        #     stm_ref_control_solved,
        #     estimated_cross_range_arc,
        #     estimated_cross_range_sigma_arc,
        #     nr_conn,
        #     control_conn_points,
        # )
        # thermal_pnt, Qx_thermal, _, Qyy_thermal_inv, A, _, _, e_thermal = _estimate_connection_point_stm(
        #     "thermal_comp",
        #     stm_ref_control_solved,
        #     estimated_thermal_arc,
        #     estimated_thermal_sigma_arc,
        #     nr_conn,
        #     control_conn_points,
        # )

        # # Apply OMT for cross_range and thermal component
        # m_omt = nr_conn  # The nr of observations is defined by the nr of arcs
        # n_omt = 1  # value for OMT is always equal to 1

        # k, T_omt_cross_range = est.overall_model_test(alpha, e_cross_range, Qyy_cross_range_inv, m_omt - n_omt, 0)
        # estimated_values_pnt_add["T_omt_cross_range"] = T_omt_cross_range

        # k, T_omt_thermal = est.overall_model_test(alpha, e_thermal, Qyy_thermal_inv, m_omt - n_omt, 0)

        # print ('Test results thermal component and cross range:')
        # print (f'k value: {k}. T_omt_cross_range: {T_omt_cross_range}. and T_omt_thermal: {T_omt_thermal}')

        # estimated_values_pnt_add["T_omt_thermal"] = T_omt_thermal
        # estimated_values_pnt_add["k_omt"] = k

        # estimated_values_pnt_add["cross_range"] = cross_range_pnt.flatten()[0]
        # estimated_values_pnt_add["cross_range_sigma"] = np.sqrt(Qx_cross_range.flatten()[0])

        # estimated_values_pnt_add["thermal_comp"] = thermal_pnt.flatten()[0]
        # estimated_values_pnt_add["thermal_comp_sigma"] = np.sqrt(Qx_thermal.flatten()[0])

        # Estimate time series for the point
        point_time_series = np.zeros(nr_epochs)
        point_time_series_corrected = np.zeros(nr_epochs)
        point_time_series_sigma = np.zeros(nr_epochs)
        T_omt_displ = np.zeros(nr_epochs)
        T_omt_displ_correction = np.zeros(nr_epochs)

        y_arcs = np.zeros((nr_conn, nr_epochs))
        y_arcs_corrected = np.zeros((nr_conn, nr_epochs))
        y_hat_arcs = np.zeros((nr_conn, nr_epochs))
        y_hat_arcs_corrected = np.zeros((nr_conn, nr_epochs))

        count_omt_reject = 0

        # Test whether we need to correct a time series by 2pi or not
        # Sometimes one of the arc will have the timeseries one full cycle above the other two
        arcs_closing_variable = _detect_ambiguous_series2(arcs_closing_variable)
        print("Check whether we need to shift a time series ")

        for t in range(nr_epochs):
            # Estimate displacement for the new point
            point_epoch, Qx_epoch, Qyy, Qyy_inv, A, y_obs, y_hat_obs, e_hat_epoch = (
                #  _estimate_connection_point_displ_stm_input(  # old version
                _estimate_connection_point_full_phase_stm_input(
                    stm_ref_control_solved,
                    arcs_closing_variable[:, t],
                    sigma_phases_arc[:, t],
                    nr_conn,
                    control_conn_points,
                    t,
                )
            )

            point_time_series[t] = point_epoch[0, 0]
            point_time_series_sigma[t] = np.sqrt(np.abs(Qx_epoch[0, 0]))
            point_time_series_corrected[t] = np.copy(point_time_series[t])
            y_arcs[:, t] = y_obs.flatten()
            y_hat_arcs[:, t] = y_hat_obs.flatten()
            y_arcs_corrected[:, t] = np.copy(y_arcs[:, t])
            y_hat_arcs_corrected[:, t] = np.copy(y_hat_arcs[:, t])

            # Apply OMT
            k, T_omt_displ[t] = est.overall_model_test(alpha, e_hat_epoch, Qyy_inv, m_omt - n_omt, 0)
            T_omt_displ_correction[t] = np.copy(T_omt_displ[t])

            # Apply w-test if the OMT is rejected
            if correct_epochs_arc == 1 and T_omt_displ[t] > k:
                count_omt_reject += 1

                corrected_y, _, _ = apply_w_test_control_network(m_omt, A, Qx_epoch, Qyy, Qyy_inv, y_obs, e_hat_epoch)

                # Re-estimate the displacement with the corrected y
                point_epoch, _, y_hat_corr, e_hat_epoch_corr = _estimate_connection_point_outlier(
                    corrected_y, A, Qyy_inv
                )
                # # Compute the OMT for this particular epoch
                k, T_omt_temp = est.overall_model_test(alpha, e_hat_epoch_corr, Qyy_inv, m_omt - n_omt, 0)

                # if T_omt_temp < T_omt_displ[t] and T_omt_temp < k:
                if T_omt_temp < T_omt_displ[t]:  # correct if the new T is better than the old one
                    T_omt_displ_correction[t] = T_omt_temp

                    y_hat_arcs_corrected[:, t] = y_hat_corr.flatten()
                    point_time_series_corrected[t] = point_epoch[0, 0]
                    y_arcs_corrected[:, t] = corrected_y.flatten()

                # TODO possible to be added that we consider other arcs if the omt is rejected over and over again

        thermal_comp_arcs = np.zeros(nr_conn)
        cross_range_values_arcs = np.zeros(nr_conn)

        thermal_comp_arcs_sigma = np.zeros(nr_conn)
        cross_range_arcs_sigma = np.zeros(nr_conn)

        for arc_estimate in range(nr_conn):
            # estimate the cross range and thermal component based on the correctly unwrapped phases
            # Construct the A matrices for functional model
            # Compute different columns for the A matrices and construct to one A matrix
            A_cr = dm.a_cross_range(cr2ph_arc)
            A_lin = dm.a_linear(years)
            A_temp = dm.a_temperature(temperature)
            A_arc = np.column_stack((A_cr, A_temp, A_lin))

            # Define the observation vector for the arc, which is based on the 'unwrapped' phase based on the filter
            y_arc = np.reshape(y_arcs[arc_estimate, :], (len(y_arcs[arc_estimate, :]), 1))
            Q_yy_arc_inv = np.identity(nr_epochs) * 1 / (sigma_phases_arc[arc_estimate, :] ** 2)

            # Estimate parameters in the phase domain
            x_hat_arc_ph, Q_x_hat_arc_ph = est.blue_q_yy_inv(A_arc, y_arc, Q_yy_arc_inv)

            thermal_comp_arcs[arc_estimate] = x_hat_arc_ph[1, 0] * 1000 / m2ph
            cross_range_values_arcs[arc_estimate] = x_hat_arc_ph[0, 0]
            thermal_comp_arcs_sigma[arc_estimate] = np.sqrt(Q_x_hat_arc_ph[1, 1]) * 1000 / m2ph
            cross_range_arcs_sigma[arc_estimate] = np.sqrt(Q_x_hat_arc_ph[0, 0])

        print("The comparison for the cross range and thermal comp")
        print(f"Estimated cross range values based on the phase values are {np.around(cross_range_values_arcs,1)}")
        print(f"Estimated cross range values based on the complex values are {np.around(estimated_cross_range_arc,1)}")

        print(f"Estimated thermal values based on the phase values are {np.around(thermal_comp_arcs,3)}")
        print(f"Estimated thermal values based on the complex values are {np.around(estimated_thermal_arc,3)}")

        # Compute the thermal component and cross range for the point
        (
            cross_range_pnt,
            Qx_cross_range,
            _,
            Qyy_cross_range_inv,
            A,
            _,
            _,
            e_cross_range,
        ) = _estimate_connection_point_stm(
            "cross_range",
            stm_ref_control_solved,
            cross_range_values_arcs,
            cross_range_arcs_sigma,
            nr_conn,
            control_conn_points,
        )
        thermal_pnt, Qx_thermal, _, Qyy_thermal_inv, A, _, _, e_thermal = _estimate_connection_point_stm(
            "thermal_comp",
            stm_ref_control_solved,
            thermal_comp_arcs,
            thermal_comp_arcs_sigma,
            nr_conn,
            control_conn_points,
        )

        # Apply OMT for cross_range and thermal component
        m_omt = nr_conn  # The nr of observations is defined by the nr of arcs
        n_omt = 1  # value for OMT is always equal to 1

        k, T_omt_cross_range = est.overall_model_test(alpha, e_cross_range, Qyy_cross_range_inv, m_omt - n_omt, 0)
        estimated_values_pnt_add["T_omt_cross_range"] = T_omt_cross_range

        k, T_omt_thermal = est.overall_model_test(alpha, e_thermal, Qyy_thermal_inv, m_omt - n_omt, 0)

        print("Test results thermal component and cross range:")
        print(f"k value: {k}. T_omt_cross_range: {T_omt_cross_range}. and T_omt_thermal: {T_omt_thermal}")

        estimated_values_pnt_add["T_omt_thermal"] = T_omt_thermal
        estimated_values_pnt_add["k_omt"] = k

        estimated_values_pnt_add["cross_range"] = cross_range_pnt.flatten()[0]
        estimated_values_pnt_add["cross_range_sigma"] = np.sqrt(Qx_cross_range.flatten()[0])

        estimated_values_pnt_add["thermal_comp"] = thermal_pnt.flatten()[0]
        estimated_values_pnt_add["thermal_comp_sigma"] = np.sqrt(Qx_thermal.flatten()[0])

        # Estimate the mean a priori an posteriori sigma values per arc
        mean_priori_quality = np.mean(sigma_phases_arc, axis=1)
        arc_residues_displacement = y_arcs_corrected - y_hat_arcs_corrected
        sigma_residues = np.std(arc_residues_displacement, axis=1)

        # Add the estimations to the dictionaries for the point (estimated result) and arc estimated
        estimated_values_pnt_add["time_series"] = point_time_series
        estimated_values_pnt_add["time_series_corrected"] = point_time_series_corrected
        estimated_values_pnt_add["point_time_series_sigma"] = point_time_series_sigma
        estimated_values_pnt_add["T_omt_displ"] = T_omt_displ
        estimated_values_pnt_add["T_omt_displ_corrected"] = T_omt_displ_correction
        estimated_values_pnt_add["connection points"] = control_conn_points
        estimated_values_pnt_add["omt_reject"] = count_omt_reject

        arc_results_pnt_add["arcs"] = arcs_connection_point
        arc_results_pnt_add["y_arcs"] = y_arcs
        arc_results_pnt_add["y_arcs_corrected"] = y_arcs_corrected
        arc_results_pnt_add["y_hat_arcs"] = y_hat_arcs
        arc_results_pnt_add["y_hat_arcs_corrected"] = y_hat_arcs_corrected
        arc_results_pnt_add["sigma arc"] = sigma_phases_arc
        arc_results_pnt_add["sigma_residues"] = sigma_residues
        arc_results_pnt_add["mean_priori_quality"] = mean_priori_quality

        # Create stm for the arcs
        # Create new stm where we solve the estimated values
        # The stm_solved will have the same dimensions as the stm that we already have
        stm_arc_results = xr.Dataset(coords={"space": np.ones(nr_conn) * point_add, "time": stm_1_point.time.values})

        stm_arc_results["pnt_idx"] = (["space"], control_conn_points)
        stm_arc_results["y_arcs"] = (["space", "time"], y_arcs)
        stm_arc_results["y_arcs_corrected"] = (["space", "time"], y_arcs_corrected)
        stm_arc_results["y_hat_arcs"] = (["space", "time"], y_hat_arcs)
        stm_arc_results["y_hat_arcs_corrected"] = (["space", "time"], y_hat_arcs_corrected)
        stm_arc_results["sigma arc"] = (["space", "time"], sigma_phases_arc)
        stm_arc_results["sigma_residues"] = (["space"], sigma_residues)
        stm_arc_results["mean_priori_quality"] = (["space"], mean_priori_quality)

        # Create new stm where we solve the estimated values
        # The stm_solved will have the same dimensions as the stm that we already have
        stm_point_solved = xr.Dataset(coords={"space": [point_add], "time": stm_1_point.time.values})

        # Add the values to the STM
        stm_point_solved["ref_pnt"] = (["space"], [ref_pnt_idx])
        stm_point_solved["control_or_not"] = (["space"], [0])
        stm_point_solved["conn_points"] = (["space"], [0])
        stm_point_solved["pnt_idx"] = (["space"], [int(point_add)])

        estimated_values_pnt_add["T_omt_cross_range"] = T_omt_cross_range

        stm_point_solved["cross_range"] = (["space"], [cross_range_pnt.flatten()[0]])
        stm_point_solved["cross_range_variance"] = (["space"], [Qx_cross_range.flatten()[0]])
        stm_point_solved["thermal_comp"] = (["space"], [thermal_pnt.flatten()[0]])
        stm_point_solved["thermal_comp_variance"] = (["space"], [Qx_thermal.flatten()[0]])
        stm_point_solved["omt_cross_range"] = (["space"], [T_omt_cross_range])
        stm_point_solved["omt_thermal"] = (["space"], [T_omt_thermal])
        stm_point_solved["k_omt"] = (["space"], [k])

        stm_point_solved["omt_reject"] = (["space"], [count_omt_reject])
        stm_point_solved["omt_displ"] = (["space", "time"], T_omt_displ.reshape(1, nr_epochs))
        stm_point_solved["omt_displ_corrected"] = (["space", "time"], T_omt_displ_correction.reshape(1, nr_epochs))

        stm_point_solved["time_series"] = (["space", "time"], point_time_series.reshape(1, nr_epochs))
        stm_point_solved["time_series_corrected"] = (
            ["space", "time"],
            point_time_series_corrected.reshape(1, nr_epochs),
        )
        stm_point_solved["time_series_variances"] = (
            ["space", "time"],
            (point_time_series_sigma**2).reshape(1, nr_epochs),
        )

        solved_the_point = 1

    return stm_point_solved, stm_arc_results, estimated_values_pnt_add, arc_results_pnt_add, solved_the_point
