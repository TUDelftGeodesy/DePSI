"""stochastic model related functions."""

import numpy as np
import xarray as xr

from depsi.arc_estimation import periodogram
from depsi.network import _independent_arcs, form_network
from depsi.utils import get_m2ph

SIGMA_MOTHER_ATMO = 15.0  # std of mother atmosphere when its included in the stochastic model
SIGMA_OTHERS = 20.0  # std of other interferograms when master
SIGMA_OVERALL = 30  # std when mother atmosphere is in functional model instead of stochastic model
SIGMA_MINIMUM = 10  # minimum value for variance components, used to avoid negative values, in degrees


def vce_temporal(
    stm: xr.Dataset,
    key_phase: str,
    key_Btemporal: str,
    key_h2ph: str,
    key_x="lon",
    key_y="lat",
    max_length: float = 0.01,  # TODO: change distance from degree to meters
    include_mother_atmo: bool = False,
) -> np.ndarray:
    """Estimate variance components per epoch.

    This estimation is performed on an Space-Time Matrix (STM) of points.
    Independent arcs are formed from the STM and unwrapped as redundancies of the estimation.

    Parameters
    ----------
    stm : xr.Dataset
        Input space-time matrix (STM).
    key_phase : str
        key for the phase data variable in the STM.
    key_Btemporal : str
        key for the Btemp data variable in the STM.
    key_h2ph : str
        key for the h2ph data variable in the STM.
    key_x : str, optional
        x coordinate for network formation, by default "lon"
    key_y : str, optional
        y coordinate for network formation, by default "lat"
    max_length : float, optional
        maximum length of the arcs, in degrees, by default 0.01
    include_mother_atmo : bool, optional
        whether to include mother atmosphere in the stochastic model.
        when False, it is assumed that the mother atmosphere is included in the functional model,
        by default False

    Returns
    -------
    np.ndarray
        Estimated variance components per epoch, in radians squared.
        This is an array of shape (Nifgs + 1,) where Nifgs is the number of interferograms.
        The extra element is for the mother atmosphere.
        The first element is the variance component for the mother atmosphere.
        If `include_mother_atmo` is False, the first element is 0.0.
    """
    # Generate a delaunay network
    arcs = form_network(
        stm, key_phase=key_phase, key_h2ph=key_h2ph, key_Btemporal=key_Btemporal, network_method="delaunay"
    )

    # Select independent arcs
    arcs_source_target = _independent_arcs(np.stack([arcs["source"].data, arcs["target"].data], axis=1))
    arcs_source_target_set = set(map(tuple, arcs_source_target))
    pairs = np.stack([arcs["source"].data, arcs["target"].data], axis=1)
    mask = np.array([tuple(x) in arcs_source_target_set for x in pairs])
    arcs = arcs.isel(space=mask)

    # Unwrap the arcs, arcs stm has standard data vars 'd_phase', 'h2ph', 'Btemp'
    phase_unwrapped, _, _, _, _ = periodogram(arcs, key_dphase="d_phase", key_h2ph="h2ph", key_Btemporal="Btemp")

    # Intiate variance components
    Nifgs = stm.sizes["time"]
    if include_mother_atmo:
        # include mother atmosphere in the stochastic model
        Qy1, Qy = _q_with_mother_atmo(Nifgs)
    else:
        # estimate variance components in the functional model
        Qy1, Qy = _q_no_mother_atmo(Nifgs)

    Qyinv = np.linalg.inv(Qy)

    # Compute Pao and QP
    Btemp = stm[key_Btemporal].values
    h2ph_approx = stm[key_h2ph].mean(dim="space").values  # Mean h2ph of all arcs
    m2ph = get_m2ph()  # Convert meters to phase
    B = np.stack([h2ph_approx * m2ph, Btemp * m2ph]).T  # Design matrix of the functional model
    Pao = np.eye(Nifgs) - B @ np.linalg.inv(B.T @ Qyinv @ B) @ B.T @ Qyinv
    QP = Qyinv @ Pao

    # Compute QPQy1QP and N, optimized using einsum
    # The following code is equivalent to this nested loop:
    # Nsig = Qy1.shape[2]  # Number of sigmas, i.e. the number of components to estimate
    # Narcs_vce = phase_unwrapped.shape[0]  # Number of arcs for VCE
    # QPQy1QP = np.full((Nifgs, Nifgs, Nsig), np.nan)
    # N = np.full((Nsig, Nsig), np.nan)
    # for k in range(Nsig):
    #     QPQy1QP[:, :, k] = QP @ Qy1[:, :, k] @ QP
    #     for j in range(Nsig):
    #         N[k, j] = np.trace(QPQy1QP[:, :, k] @ Qy1[:, :, j])
    QPQy1QP = np.einsum("ij,jlk,lm ->imk", QP, Qy1, QP, optimize=True)
    N = np.einsum("abk,baj->kj", QPQy1QP, Qy1, optimize=True)

    Ninv = np.linalg.inv(N)

    # Estimate variance components from all independent arcs
    # The following code is equivalent to this nested loop:
    # sig2 = np.full((Nsig, Narcs_vce), np.nan)
    # l = np.full((Nsig, 1), np.nan)
    # for v in range(Narcs_vce):
    #     y = phase_unwrapped[v, :].reshape(-1, 1)
    #     for k in range(Nsig):
    #         l[k, 0] = (y.T @ QPQy1QP[:, :, k] @ y).squeeze()
    #     sig2[:, v] = (Ninv @ l).flatten()
    l_vec = np.einsum("ij,jmk,mi ->ki", phase_unwrapped, QPQy1QP, phase_unwrapped.T, optimize=True)  # l vector
    sig2_all_arcs = Ninv @ l_vec
    sig2_est = np.mean(sig2_all_arcs, axis=1)

    # Apply threshold to avoid small and negative values
    threshold = (np.pi * SIGMA_MINIMUM / 180) ** 2
    sig2_est[sig2_est < threshold] = threshold

    # if atmosphere is not included in the stochastic model, add a zero at the beginning
    if not include_mother_atmo:
        sig2_est = np.insert(sig2_est, 0, 0.0)

    return sig2_est


def _q_with_mother_atmo(Nifgs: int) -> tuple:
    """Build Qy1 and Qy matrices with mother atmosphere."""
    sig0_mother = (np.pi * SIGMA_MOTHER_ATMO / 180) ** 2
    sig0_others = (np.pi * SIGMA_OTHERS / 180) ** 2

    # Build "Design matrix" for variance components estimation
    Qy1 = np.zeros((Nifgs, Nifgs, Nifgs + 1), dtype=np.int16)  # Build Qy1 matrix, has 0 or 2
    Qy1[:, :, 0] = 2  # Mother epoch, all 2
    for v in range(0, Nifgs):
        Qy1[v, v, v + 1] = 2  # Other epochs, on location set 2
    Qy = np.full((Nifgs, Nifgs), 2 * sig0_mother)  # Build Qy matrix has 2 * sig0_mother in background
    for v in range(Nifgs):
        Qy[v, v] += 2 * sig0_others  # Per epoch add 2 * sig0_others

    return Qy1, Qy


def _q_no_mother_atmo(Nifgs: int) -> tuple:
    """Build Qy1 and Qy matrices without mother atmosphere."""
    sig0 = (np.pi * SIGMA_OVERALL / 180) ** 2

    # Build "Design matrix" for variance components estimation
    Qy1 = np.zeros((Nifgs, Nifgs, Nifgs), dtype=np.int16)  # Build Qy1 matrix, has 0 or 2
    for v in range(Nifgs):
        Qy1[v, v, v] = 2  # Other epochs, on location set 2
    Qy = np.diag(np.full(Nifgs, 2 * sig0))  # Build Qy matrix has 2 * sig0 in background
    return Qy1, Qy
