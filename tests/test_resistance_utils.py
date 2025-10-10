import math
import pytest

from FetoFlow.resistance_utils import (
    calculate_viscosity_factor_from_radius,
    calculate_convolute_resistance,
    calculate_capillary_equivalent_resistance,
)


def test_viscosity_factor_monotonic():
    # Larger radius should give a viscosity factor <= smaller radius (empirical expectation)
    small = calculate_viscosity_factor_from_radius(1e-6)
    large = calculate_viscosity_factor_from_radius(5e-6)
    assert isinstance(small, float) and isinstance(large, float)
    assert large <= small or math.isclose(large, small)


def test_convolute_resistance_positive():
    R = calculate_convolute_resistance(n_series=3, n_parallel=6, mu=0.33600e-02, visc_factor=1, capillary_length=0.0005, capillary_convolute_radius=7.2e-6)
    assert R > 0


def test_capillary_equivalent_resistance_invalid_model():
    with pytest.raises(ValueError):
        calculate_capillary_equivalent_resistance(radius_in_artery=1e-5, radius_in_vein=2e-5, capillary_model="not_a_model")


def test_capillary_equivalent_resistance_consistent():
    R1 = calculate_capillary_equivalent_resistance(radius_in_artery=1e-5, radius_in_vein=2e-5, viscosity_model="constant")
    R2 = calculate_capillary_equivalent_resistance(radius_in_artery=1e-5, radius_in_vein=2e-5, viscosity_model="pries_vessel")
    assert R1 > 0 and R2 > 0


def test_viscosity_factor_extremes():
    # Very small radius should produce a viscosity factor > 0
    small = calculate_viscosity_factor_from_radius(1e-7)
    large = calculate_viscosity_factor_from_radius(1e-3)
    assert small > 0
    assert large > 0


def test_capillary_resistance_varied_parameters():
    for n_series in [1, 3]:
        for n_parallel in [1, 6]:
            R = calculate_convolute_resistance(n_series=n_series, n_parallel=n_parallel, mu=0.33600e-02, visc_factor=1, capillary_length=0.0005, capillary_convolute_radius=7.2e-6)
            assert R > 0
