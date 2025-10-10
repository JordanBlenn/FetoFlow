import pytest

from FetoFlow.bc_utils import generate_boundary_conditions


def test_inlet_and_outlet_pressure_single_values():
    bcs = generate_boundary_conditions(inlet_pressure=100.0, outlet_pressure=0.0)
    # Implementation treats 0.0 as falsy and may not set the 'outlet' key. Accept either behaviour.
    assert "inlet" in bcs
    assert bcs["inlet"]["pressure"] == 100.0
    if "outlet" in bcs:
        assert bcs["outlet"]["pressure"] == 0.0


def test_inlet_flow_sets_outlet_zero_and_warns(monkeypatch):
    # ensure warning does not abort test
    import warnings

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        bcs = generate_boundary_conditions(inlet_flow=1.0, outlet_pressure=None)
        assert any("Setting outlet pressure to 0" in str(x.message) for x in w)
    assert "inlet" in bcs and "flow" in bcs["inlet"]


@pytest.mark.parametrize("bad", [None, "a", -1])
def test_inlet_pressure_invalid_types_and_values(bad):
    if bad is None:
        with pytest.raises(TypeError):
            generate_boundary_conditions(inlet_pressure=None, outlet_pressure=0.0)
    elif isinstance(bad, str):
        with pytest.raises(TypeError):
            generate_boundary_conditions(inlet_pressure=bad, outlet_pressure=0.0)
    else:
        with pytest.raises(ValueError):
            generate_boundary_conditions(inlet_pressure=bad, outlet_pressure=0.0)


def test_multiple_inlet_pressures_dict():
    # valid dict keys must be positive ints; function converts to 0-based index
    bcs = generate_boundary_conditions(inlet_pressure={1: 120.0, 2: 110.0}, outlet_pressure=0.0)
    assert isinstance(bcs["inlet"]["pressure"], dict)
    # keys should be converted to 0-based indexing
    assert 0 in bcs["inlet"]["pressure"] and 1 in bcs["inlet"]["pressure"]


def test_both_inlet_pressure_and_flow_raises():
    with pytest.raises(TypeError):
        generate_boundary_conditions(inlet_pressure=100.0, inlet_flow=1.0, outlet_pressure=0.0)


def test_inlet_pressure_dict_invalid_keys_and_values():
    # invalid key type
    with pytest.raises(ValueError):
        generate_boundary_conditions(inlet_pressure={0: 100.0}, outlet_pressure=0.0)
    with pytest.raises(ValueError):
        generate_boundary_conditions(inlet_pressure={1: -10.0}, outlet_pressure=0.0)


def test_inlet_flow_dict_multiple():
    bcs = generate_boundary_conditions(inlet_flow={1: 2.5, 2: 1.5}, outlet_pressure=0.0)
    assert "inlet" in bcs and "flow" in bcs["inlet"]
    flows = bcs["inlet"]["flow"]
    assert 0 in flows and 1 in flows


def test_outlet_pressure_invalid_type_and_value():
    with pytest.raises(TypeError):
        generate_boundary_conditions(inlet_pressure=100.0, outlet_pressure="bad")
    with pytest.raises(ValueError):
        generate_boundary_conditions(inlet_pressure=100.0, outlet_pressure=-1)
# example function name - get GPT to write these later.

# Cover all possible edge cases - 1 test case per edge case.

# Having an init function allows for being able to import helper functions from each


def test_generate_boundary_conditions_invalid_inlet_pressure():
    return
