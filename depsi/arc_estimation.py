"""arc estimation algorithms."""

import dask.array as da
import numpy as np
import xarray as xr

from depsi.utils import wrap_phase

# Constants
STOP_HEIGHT = 1e-4  # Stop search step for height [m]
STOP_VEL = 1e-7  # Stop search step for velocity [m/y]
MAX_COUNT = 10  # Maximum number of search iterations
THRES_TEMP_COH_MEMORY = 500  # unit MB. Threshold to use dask array for temporal coherence calculation in periodogram.


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
    init_step_height: float = 3.0,
    init_step_vel: float = 2e-3,
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
        Initial step size for the height parameter in meters, by default 3.0.
        This value sets the resolution of the initial search space for the height parameter.
        After every search, the step size will be reduced by a factor of 10.
    init_step_vel : float, optional
        Initial step size for the velocity parameter in meters per year, by default 2e-3.
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
    assert (len(stm[key_Btemporal].dims) == 1) and ("time" in stm[key_Btemporal].dims), (
        "year time should and only should contain the 'time' dimension."
    )

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
    if "time" in stm.chunks.keys():
        if len(stm.chunks["time"]) != 1:
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
    # The residuals phase_residual_all_arcs is a large array with n_arcs x n_obs x n_search
    # so use .data to avoid loading it into memory if it is a dask array
    dphase_obs = stm[key_dphase].data  # n_arcs x n_obs x 1

    # If the memory size of phase_residual_all_arcs will exceed the threshold
    # chunk dphase_obs and init_search_space to enable computation
    mem_size_estimation = (
        dphase_obs.shape[0] * dphase_obs.shape[1] * init_search_space.shape[0] * dphase_obs.dtype.itemsize
    ) / (1024**2)  # in MB
    if mem_size_estimation > THRES_TEMP_COH_MEMORY:
        dphase_obs, init_search_space = _chunk_for_temp_coh_compute(dphase_obs, init_search_space)

    # Compute modelled phase for all arcs and all search candidates
    phs_model = B @ init_search_space.T  # n_obs x n_search

    # Expand dimensions and compute the phase residuals for all arcs and all search candidates
    dphase_obs = dphase_obs[:, :, None]  # n_arcs x n_obs x 1
    phs_model = phs_model[None, :, :]  # 1 x n_obs x n_search
    phase_residual_all_arcs = dphase_obs - phs_model

    # Find the best initial height and velocity based on the temporal coherence
    coh_search_space_all_arcs = (
        np.cos(phase_residual_all_arcs).sum(axis=1) + 1j * np.sin(phase_residual_all_arcs).sum(axis=1)
    ) / stm[key_dphase].sizes["time"]  # n_arcs x n_search
    coh_idx_all_arcs = np.argmax(np.abs(coh_search_space_all_arcs), axis=1)  # n_arcs

    # Implicitly compute best coh index if dask array
    coh_idx_all_arcs = coh_idx_all_arcs.compute() if isinstance(coh_idx_all_arcs, da.Array) else coh_idx_all_arcs

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
    # So we have the input_core_dims as [["time"], ["time"], [], []]
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
    # Assign initial values for the search
    step_height = init_step_height
    step_vel = init_step_vel
    param_height = init_height
    param_vel = init_vel

    # Calculate the initial temporal coherence for the initial height and velocity,
    # in case the search loop is not entered
    phs_model = B @ np.array([param_height, param_vel])  # size n_obs
    phase_residual = phs_obs_wrapped[:, None] - phs_model
    coh_best = (np.cos(phase_residual).sum() + 1j * np.sin(phase_residual).sum()) / phs_obs_wrapped.shape[0]

    # Search loop
    count = 0
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


def _chunk_for_temp_coh_compute(phs_obs_wrapped, search_space):
    """Chunk observations and search space for temporal coherence computation."""
    if isinstance(phs_obs_wrapped, da.Array):  # Existing chunk size for the arc dimension
        chunk_arcs = phs_obs_wrapped.chunks[0][0]
    else:
        # If phs_obs_wrapped is not a dask array, chunk it in the arc dimension, making each chunk about 10 MB
        chunk_arcs = max(1, 10 * 1024**2 // (phs_obs_wrapped.shape[1] * phs_obs_wrapped.dtype.itemsize))
        phs_obs_wrapped = da.from_array(phs_obs_wrapped, chunks=(chunk_arcs, phs_obs_wrapped.shape[1]))
    # Decide the chunk size for the search space dimension
    # making each n_arcs x n_obs x n_search chunk about 100 MB
    chunk_searches = max(1, 100 * 1024**2 // (chunk_arcs * phs_obs_wrapped.shape[1] * phs_obs_wrapped.dtype.itemsize))
    search_space = da.from_array(search_space, chunks=(chunk_searches, 2))

    return phs_obs_wrapped, search_space


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
