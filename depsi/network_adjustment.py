import numpy as np
from scipy.sparse import csr_matrix

import depsi.estimation as est


def network_adjustment_control_network(A_adjustment, y_obs, sigma_obs):
    """Perform a network adjustment to estimate the parameters of points using the estimated parameters for the arc.

    This function applies BLUE to adjust the parameters of the points in the
    network, minimizing the residuals between the observed and adjusted values of the arcs. The adjustment process also
    provides estimates of the variance-covariance matrix of the adjusted parameters and the observations.

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

    # Compute the matrix inverse
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

    # Apply the w-test per detected outlier
    for w in range(m):
        c_i = np.zeros((m, 1))
        c_i[w, 0] = 1

        A = (c_i.transpose() @ Qyy_inv @ e_hat)[0, 0]  # Haal de enkele waarde op
        B = np.sqrt((c_i.transpose() @ Qyy_inv @ Qee @ Qyy_inv @ c_i)[0, 0])

        w_result[w] = A / B

    w_result = np.abs(w_result)

    # Re-wrap the biggest outlier
    idx_biggest_w = np.argmax(w_result)
    corrected_y = np.copy(y)

    if e_hat[idx_biggest_w] > 0:
        corrected_y[idx_biggest_w] = corrected_y[idx_biggest_w] - 2 * np.pi
        correct = -1
    elif e_hat[idx_biggest_w] < 0:
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

    This function iterates over the adjustment points, updating the estimated heights, thermal components, and VCM.
    Additionally, it updates the time series data, corrected time series, and their variance-covariance
    matrices for each point.

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
    # Loop over each point in the adjustment and update estimated values and time series
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
    in the network. The matrix is used for adjusting the network to minimize the residuals between the observed
    arcs and the estimated positions of the points.

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
    neq = len(arcs)
    npt = nr_pnts
    eqs = np.concatenate((np.ones(neq), -1 * np.ones(neq)))
    rows = np.concatenate((np.arange(neq), np.arange(neq)))
    cols = np.concatenate((arcs[:, 1], arcs[:, 0]))

    A = csr_matrix((eqs, (rows, cols)), shape=(neq, npt))
    A_dense = A.toarray()

    A_small = A_dense[:, points_network]

    return A_small
