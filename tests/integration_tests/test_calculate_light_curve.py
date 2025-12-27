import numpy as np

from binara import internal_calculate_light_curve


def test_calculate_light_curve():
    # Parameters and sample light curve taken from a run of Gioula Kalapotharakos for TIC ID 110602878 sector 34.
    model_parameters = np.array([
        0.385886490907, 0.0664402461459, 0.728154, 0.398693367346, 0.0250731943869, -2.58470776396,
        0.986066718147, -0.0495279552116, -1.02057860777, 0.368375172834, 0.275439743308,
        0.0288241638398, 0.374759311771, 0.922036025043, 1.61396897845, 0.113311357046, 0.242166052741,
        -1.42801281267, 1.40666171312, 0.650303472667, 0.999811553496, 0.0909320263842
    ], dtype=np.float64)
    model_phases = np.array([
        0.0534753500, 1.0019700296, 1.2655832833, 4.7593061500, 5.0445080167], dtype=np.float64)
    expected_model_fluxes = np.array([
        1.0000618467, 1.0873321167, 0.9824356398, 0.9997781792, 0.9364980438], dtype=np.float64)
    model_fluxes = internal_calculate_light_curve(
        model_phases,
        model_parameters,
    )
    assert np.allclose(model_fluxes, expected_model_fluxes)
