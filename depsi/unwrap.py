"""modules for phase unwrapping."""

import numpy as np
import xarray as xr

from depsi.utils import wrap_phase

# Constants
STOP_HEIGHT = 1e-4  # Stop search step for height
STOP_VEL = 1e-7  # Stop search step for velocity
MAX_COUNT = 10  # Maximum number of search iterations
WAVELENGTH_S1 = 0.055465763  # m, sentinel-1 wavelength


def periodogram(
    stm: xr.Dataset,
    key_phs: str,
    key_h2ph: str,
    key_yeartime: str,
    std_obs: float = 1.0,
    std_height: float = 50.0,
    std_vel: float = 0.02,
    init_height: float = 0.0,
    init_vel: float = 0.0,
    init_step_height: float = 2.0,
    init_step_vel: float = 1e-3,
    min_searches: int = 11,
    wavelength: float = None,
):
    """Periodogram unwrapping algorithm.

    This function performs periodogram unwrapping on arcs as Space-Time Matrices (STMs).

    It uses a deformation model with two parameters: height and velocity to estimate the unwrapped phase.

    For computation efficiency, the design matrix is constructed only once for all arcs, utilizing the average
    height-to-phase conversion factor (h2ph) across all arcs. The effect of using this average is corrected later.

    Parameters
    ----------
    stm : xr.Dataset
        Input Space-Time Matrix (STM) containing the wrapped phase, height-to-phase conversion factor, and year-time.
    key_phs : str
        Key for the wrapped phase data variable in the STM.
    key_h2ph : str
        Key for the height-to-phase conversion factor in the STM.
    key_yeartime : str
        key for the year-time data variable in the STM.
    std_obs : float, optional
        A-poriori standard deviation of the observations in rads, by default 1.0.
        This value is used to construct the stochastic model (Qyy) of the observations.
    std_height : float, optional
        A-poriori standard deviation of the height in meters, by default 50.0.
        This value is used to construct the boundaries initial search space for the height parameter.
    std_vel : float, optional
        A-poriori standard deviation of the velocity in meters per year, by default 0.02.
        This value is used to construct the boundaries initial search space for the velocity parameter.
    init_height : float, optional
        Initial guess for the height parameter in meters, by default 0.0.
    init_vel : float, optional
        Initial guess for the velocity parameter in meters per year, by default 0.0.
    init_step_height : float, optional
        Initial step size for the height parameter in meters, by default 2.0.
        This value is used to construct the resolution of the initial search space for the height parameter.
        After every search, the step size will be reduced by a factor of 10.
    init_step_vel : float, optional
        Initial step size for the velocity parameter in meters per year, by default 1e-3.
        This value is used to construct the resolution of the initial search space for the velocity parameter.
        After every search, the step size will be reduced by a factor of 10.
    min_searches : int, optional
        Minimum number of searches for the height and velocity parameters, by default 11.
        If the number of the initial search space is smaller than this value, it will be increased to this value.
        After the first search, the number of searches will be set to this value.
    wavelength : float, optional
        Wavelength of the sensor in meters, by default None.
        If not provided, the default wavelength for Sentinel-1 will be used.
        This value is used to calculate the height-to-phase conversion factor (m2ph).

    Returns
    -------
    Tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]
        Returns the unwrapped phase, ambiguities, estimated height, estimated velocity, and temporal coherence.
        - Unwrapped phase: in rads, shape (n_arcs, n_obs), dtype np.float64.
        - Ambiguities: unitless, shape (n_arcs, n_obs), dtype np.float64.
        - Estimated height: in meters, shape (n_arcs,), dtype np.float64.
        - Estimated velocity: in meters per year, shape (n_arcs,), dtype np.float64.
        - Temporal coherence: unitless complex number, shape (n_arcs,), dtype np.complex128.
    """
    # If wavelength is not provided, use the default sentinel-1 wavelength
    # TODO: get wavelength from metadata stm.attrs
    if wavelength is None:
        m2ph = -4 * np.pi / WAVELENGTH_S1

    # Set up functional and stochastic model for all arcs
    # Here we use the same h2ph (average over all arcs) for all arcs and correct the effect later
    # Doing this avoids perform matrix inversion for each arc
    h2ph_mean = stm[key_h2ph].mean(dim="space").data  # Mean h2ph of all arcs

    # Design matrix B, size n_obs x n_params
    # In B, h2ph should also be multiplied by m2ph since it did not when it was created
    B = np.stack([h2ph_mean * m2ph, stm[key_yeartime] * m2ph]).T.compute()

    # Stochastic model Qyy, size n_obs x n_obs
    # This is the covariance matrix of the observations
    Qyy = np.diag(np.repeat(std_obs**2, stm[key_h2ph].sizes["time"]))

    # Covenience matrix R and rhs for the least squares solution
    R = B.T @ np.linalg.inv(Qyy) @ B  # B.T * Qyy^-1 * B , size n_params x n_params
    rhs = np.linalg.inv(R) @ B.T @ np.linalg.inv(Qyy)  # (B.T * Qyy^-1 * B)^-1 * B.T * Qyy^-1, size n_params x n_obs

    # Set up core dimensions, which are the dimensions _periodogram_single will be applied to
    # We are broadcasting _periodogram_single on stm[key_phs] along the space dimension
    # Threfore, we are calling it on the "time" dimension of every space entry
    # So we have the input_core_dims as ["time"]
    input_core_dims = [
        ["time"],
    ]
    # There are 5 outputs from _periodogram_single
    # The first two are np arrays with time dimension
    # The other three are scalars, so they have no dimensions
    output_core_dims = [["time"], ["time"], [], [], []]

    # Apply the _periodogram_single on stm[key_phs] along "space" dimension
    # Other parameters are duplicated for each space entry
    # Therefore they can be passed as kwargs
    results = xr.apply_ufunc(
        _periodogram_single,
        stm[key_phs],
        input_core_dims=input_core_dims,
        output_core_dims=output_core_dims,
        kwargs={
            "B": B,
            "Qyy": Qyy,
            "R": R,
            "rhs": rhs,
            "std_height": std_height,
            "std_vel": std_vel,
            "init_height": init_height,
            "init_vel": init_vel,
            "init_step_height": init_step_height,
            "init_step_vel": init_step_vel,
            "min_searches": min_searches,
        },
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64, np.float64, np.float64, np.float64, np.complex128],
    )

    return results


def _periodogram_single(
    phs_obs_wrapped: np.ndarray,
    B: np.ndarray,
    Qyy: np.ndarray,
    R: np.ndarray,
    rhs: np.ndarray,
    std_height: float,
    std_vel: float,
    init_height: float,
    init_vel: float,
    init_step_height: float,
    init_step_vel: float,
    min_searches: float,
):
    """Periodogram unwrapping for a single arc.

    Parameters
    ----------
    phs_obs_wrapped : np.ndarray
        Wrapped phase observations in radians, shape (n_obs,).
    B : np.ndarray
        Design matrix, size n_obs x n_params, where n_params = 2 (height and velocity).
    Qyy : np.ndarray
        Stochastic model of the observations, size n_obs x n_obs.
    R : np.ndarray
        Covariance matrix of the parameters, size n_params x n_params.
    rhs : np.ndarray
        Right-hand side matrix for the least squares solution, size n_params x n_obs.
    std_height : float
        A-poriori standard deviation of the height in meters.
    std_vel : float
        A-poriori standard deviation of the velocity in meters per year.
    init_height : float
        Initial guess for the height parameter in meters.
    init_vel : float
        Initial guess for the velocity parameter in meters per year.
    init_step_height : float
        Initial step size for the height parameter in meters.
    init_step_vel : float
        Initial step size for the velocity parameter in meters per year.
    min_searches : float
        Minimum number of searches for the height and velocity parameters.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, float, float, complex]
        Returns the unwrapped phase, ambiguities, estimated height, estimated velocity, and temporal coherence.
        - Unwrapped phase: in rads, shape (n_obs,), dtype np.float64.
        - Ambiguities: unitless, shape (n_obs,), dtype np.float64.
        - Estimated height: in meters, scalar, dtype np.float64.
        - Estimated velocity: in meters per year, scalar, dtype np.float64.
        - Temporal coherence: unitless complex number, scalar, dtype np.complex128.
    """
    # Build initial search space for height and velocity
    param_height = init_height
    param_vel = init_vel
    step_height = init_step_height
    step_vel = init_step_vel
    n_search_height = max(round(2 * std_height / step_height), min_searches)
    n_search_vel = max(round(2 * std_vel / step_vel), min_searches)

    # Search loop
    count = 0
    while step_height > STOP_HEIGHT and step_vel > STOP_VEL and count < MAX_COUNT:
        # Build search space
        search_space = _build_search_space(
            param_height, param_vel, step_height, step_vel, n_search_height, n_search_vel
        )

        # Calculate the wrpped model phase for all candidates
        phs_model = wrap_phase(B @ search_space.T)  # size n_obs x n_search

        # Calculate the temporal coherence for all search candidates
        # Expand dimension of phs_obs_wrapped to facilitate broadcasting
        # No need to repeat phs_obs_wrapped since the minus operation will broadcast to the shape of phs_model
        # Sum along axis=0 which is the observation axis
        # Reference: van Leijen 2014, Eq. 4.55
        coh_search_space = (
            np.exp(1j * (np.expand_dims(phs_obs_wrapped, axis=1) - phs_model)).sum(axis=0) / phs_obs_wrapped.shape[0]
        )

        # Get the best temporal coherence value and its index
        coh_idx = np.argmax(np.abs(coh_search_space))
        coh_best = coh_search_space[coh_idx]

        # Update values needed for search space
        # Reduce step size to 1/10
        param_height = search_space[coh_idx, 0]
        param_vel = search_space[coh_idx, 1]
        step_height /= 10
        step_vel /= 10
        n_search_height = min_searches
        n_search_vel = min_searches

        count += 1

    phs_model_abs = B @ np.array([param_height, param_vel])  # Absolute modelled phase
    phs_model_wrapped = wrap_phase(phs_model_abs)  # Wrapped modelled phase
    ambigs = np.round((phs_model_abs + phs_model_wrapped - phs_obs_wrapped) / (2 * np.pi))  # Ambiguities
    phs_obs_unwrapped = 2 * np.pi * ambigs + phs_obs_wrapped  # Unwrapped phase
    param = rhs @ phs_obs_unwrapped  # [height_est, velocity_est]

    return phs_obs_unwrapped, ambigs, param[0], param[1], coh_best


def _build_search_space(param_height, param_vel, step_height, step_vel, n_search_height, n_search_vel):
    """Construct the search space for height and velocity parameters.

    For both height and velocity, the candidates are generated around the initial values according to the step size
    and the number of searches. On each side of the initial value, N candidates are generated with a step size, where
    N is specified by `n_search_height` and `n_search_vel`, and the step size is specified by `step_height` and
    `step_vel`.

    Then all possible combinations of height and velocity candidates are created to form
    the search space.

    Parameters
    ----------
    param_height : float
        Initial height parameter in meters.
    param_vel : _type_
        Initial velocity parameter in meters per year.
    step_height : _type_
        Step size of search for height parameter, in meters.
    step_vel : _type_
        Step size of search for velocity parameter, in meters per year.
    n_search_height : _type_
        Number of searches for height parameter on each side of the initial value.
    n_search_vel : _type_
        Number of searches for velocity parameter on each side of the initial value.

    Returns
    -------
    _type_
        _description_
    """
    height_candidates = np.arange(
        param_height - n_search_height * step_height,
        param_height + n_search_height * step_height + step_height,
        step_height,
    )

    vel_candidates = np.arange(
        param_vel - n_search_vel * step_vel, param_vel + n_search_vel * step_vel + step_vel, step_vel
    )

    # All possible combinations of height and velocity
    search_space = np.array(np.meshgrid(height_candidates, vel_candidates)).T.reshape(-1, 2)

    return search_space
