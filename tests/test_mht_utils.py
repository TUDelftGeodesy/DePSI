import pytest

from depsi.mht_utils import pretest


@pytest.mark.parametrize(
    "b, alpha0, gamma0, expected",
    [
        (1, 0.001, 0.80, (17.07464681, 3.29052673, 10.8275661, 0.001)),
        (5, 0.01, 0.90, (14.87938717, 2.575829303, 9.96265261, 0.0763)),
        (10, 0.05, 0.95, (12.99470909, 1.95996398, 10.76181669, 0.37637)),
    ],
)
def test_pretest(b, alpha0, gamma0, expected):
    lam0, k1, kb, alphab = pretest(b, alpha0, gamma0)
    assert pytest.approx(lam0, rel=1e-5) == expected[0]
    assert pytest.approx(k1, rel=1e-5) == expected[1]
    assert pytest.approx(kb, rel=1e-5) == expected[2]
    assert pytest.approx(alphab, rel=1e-3) == expected[3]


def test_pretest_zero_redundancy(caplog):
    # Expect a warning in logger output
    with caplog.at_level("WARNING"):
        lam0, k1, kb, alphab = pretest(0, 0.01, 0.95)
        assert "zero redundancy" in caplog.text
