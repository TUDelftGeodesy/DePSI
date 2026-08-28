import numpy as np
import pytest

from depsi.model_definition import construct_design_matrix


def test_construct_design_matrix_basic_case():
    """Test construct_design_matrix against hand-computed offset + velocity columns."""
    # Arrange
    m2ph = 2.0
    time = np.array([0.0, 1.0, 2.0])

    # Act
    A = construct_design_matrix(["offset", "velocity"], m2ph, n_epochs=3, time=time)

    # Assert
    expected = np.array(
        [
            [1.0, 0.0],
            [1.0, 2.0],
            [1.0, 4.0],
        ]
    )
    np.testing.assert_allclose(A, expected, atol=1e-8)


def test_construct_design_matrix_respects_models_order():
    """Test that the output columns follow the order of the `models` argument, not insertion order."""
    # Arrange
    m2ph = 2.0
    time = np.array([0.0, 1.0, 2.0])

    # Act
    A_offset_first = construct_design_matrix(["offset", "velocity"], m2ph, n_epochs=3, time=time)
    A_velocity_first = construct_design_matrix(["velocity", "offset"], m2ph, n_epochs=3, time=time)

    # Assert
    np.testing.assert_allclose(A_offset_first, A_velocity_first[:, ::-1], atol=1e-8)


def test_construct_design_matrix_rejects_unsupported_model():
    """Test that an unsupported model name raises a NotImplementedError."""
    with pytest.raises(NotImplementedError):
        construct_design_matrix(["not_a_real_model"], m2ph=1.0)


def test_construct_design_matrix_empty_models_raises():
    """Test the current boundary behavior for an empty `models` list."""
    with pytest.raises(ValueError):
        construct_design_matrix([], m2ph=1.0)


def test_construct_design_matrix_rejects_missing_model_inputs():
    """Test that a ValueError is raised when a requested model's required input is not provided."""
    with pytest.raises(ValueError):
        construct_design_matrix(["velocity"], m2ph=1.0, n_epochs=3)
