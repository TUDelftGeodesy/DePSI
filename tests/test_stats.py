import pytest

from depsi.stats import lambda0, pretest


@pytest.mark.parametrize(
    "b, alpha0, gamma0, expected",
    [
        (2, 0.005, 0.80, (13.31268332, 2.807033768, 8.7698881, 0.01246359)),
        (5, 0.01, 0.90, (14.87938717, 2.575829303, 9.96265261, 0.07630041)),
        (10, 0.05, 0.95, (12.99470909, 1.95996398, 10.76181669, 0.37637336)),
    ],
)
def test_pretest(b, alpha0, gamma0, expected):
    lam0, k1, kb, alphab = pretest(b, alpha0, gamma0)
    assert pytest.approx(lam0, rel=1e-8) == expected[0]
    assert pytest.approx(k1, rel=1e-6) == expected[1]
    assert pytest.approx(kb, rel=1e-6) == expected[2]
    assert pytest.approx(alphab, rel=1e-4) == expected[3]


def test_pretest_zero_redundancy(caplog):
    # Expect a warning in logger output
    with caplog.at_level("WARNING"):
        _ = pretest(0, 0.01, 0.95)
        assert "zero redundancy" in caplog.text


@pytest.mark.parametrize(
    "gam0, df, cv",
    [
        (0.80, 0, 3.84),  # df must be >= 1
        (0.90, 5, -1.1),  # cv must be >= 0
        (1.2, 10, 18.31),  # gam0 must be in [0, 1)
        (0.001, 2, 10),  # negative lambda
    ],
)
def test_lambda0_error(gam0, df, cv):
    with pytest.raises(ValueError):
        lambda0(gam0, df, cv)
