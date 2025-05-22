"""modules for phase unwrapping."""

import numpy as np

from depsi.utils import wrap_phase

# Constants
STOP_HEIGHT = 1e-4  # Stop search step for height
STOP_VEL = 1e-7  # Stop search step for velocity
MAX_COUNT = 10  # Maximum number of search iterations
WAVELENGTH = 0.055465763  # m, sentinel-1 wavelength


def _periodogram_single(
    phs_wrapped: np.ndarray,
    Btemp: np.ndarray,
    h2ph: np.ndarray,
    std_obs: float,
    std_height: float,
    std_vel: float,
    init_height: float = 0.0,
    init_vel: float = 0.0,
    init_step_height: float = 1.0,
    init_step_vel: float = 1e-4,
    min_searches=11,
):
    # This function estimates DEM error and linear deformation for a single persistent scatterer point
    # using periodogram search.

    n_obs = phs_wrapped.shape[0]  # Number of observations
    m2ph = -4 * np.pi / WAVELENGTH
    B = np.stack([h2ph * m2ph, Btemp * m2ph]).T  # Design matrix, TODO: does h2ph need to be multiplied by m2ph?
    Qyy = np.diag(np.repeat(std_obs**2, n_obs))  # Construct covariance matrix for observations
    R = B.T @ np.linalg.inv(Qyy) @ B  # B.T * Qyy^-1 * B
    rhs = np.linalg.inv(R) @ B.T @ np.linalg.inv(Qyy)  # (B.T * Qyy^-1 * B)^-1 * B.T * Qyy^-1

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
        # Expand dimension to make later calculations easier
        # No need to repeat the phase model since the minus operation will handle it
        # Sum along axis=0 which is the observation axis
        # Reference: van Leijen 2014, Eq. 4.55
        coh_search_space = np.exp(1j * (np.expand_dims(phs_wrapped, axis=1) - phs_model)).sum(axis=0) / n_obs

        # Get the best temporal coherence value and its index
        coh_idx = np.argmax(np.abs(coh_search_space))
        coh_best = coh_search_space[coh_idx]

        # Update values needed for search space
        # Reduce step size to 1/10
        # TODO: check if we need to update the search space according to the post-priori Qxx
        param_height = search_space[coh_idx, 0]
        param_vel = search_space[coh_idx, 1]
        step_height /= 10
        step_vel /= 10

        count += 1

    phs_model_abs = B @ np.array([param_height, param_vel])  # Absolute modelled phase
    phs_model_wrapped = wrap_phase(phs_model_abs)  # Wrapped modelled phase
    ambigs = np.round((phs_model_abs + phs_model_wrapped - phs_wrapped) / (2 * np.pi))  # Ambiguities
    phs_wrapped_unw = 2 * np.pi * ambigs + phs_wrapped  # Unwrapped phase
    param = rhs @ phs_wrapped_unw  # [height_est, velocity_est]

    return phs_wrapped_unw, ambigs, param[0], param[1], coh_best


def _build_search_space(param_height, param_vel, step_height, step_vel, n_search_height, n_search_vel):
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
