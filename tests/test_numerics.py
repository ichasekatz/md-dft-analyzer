"""Tests for the shared polynomial fitting utilities."""

from __future__ import annotations

from unittest import TestCase

import numpy as np
import pytest

from md_dft_analyzer.numerics import gauss_pivot, poly_eval, poly_fit, r_squared


class TestPolyFit(TestCase):
    """Tests for poly_fit and gauss_pivot."""

    def test_linear_fit_recovers_slope(self):
        """Degree-1 fit recovers slope and intercept of a linear signal."""
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = 3.0 * x + 7.0
        c = poly_fit(x, y, degree=1)
        assert c[0] == pytest.approx(7.0, abs=1e-8)
        assert c[1] == pytest.approx(3.0, abs=1e-8)

    def test_quadratic_fit_recovers_coefficients(self):
        """Degree-2 fit recovers all three coefficients of a quadratic."""
        x = np.linspace(0, 5, 20)
        y = 2.0 * x**2 - x + 4.0
        c = poly_fit(x, y, degree=2)
        assert c[0] == pytest.approx(4.0, abs=1e-6)
        assert c[1] == pytest.approx(-1.0, abs=1e-6)
        assert c[2] == pytest.approx(2.0, abs=1e-6)

    def test_poly_eval_at_known_points(self):
        """poly_eval matches manual evaluation."""
        c = np.array([1.0, 2.0, 3.0])  # 1 + 2x + 3x²
        assert poly_eval(c, 0.0) == pytest.approx(1.0)
        assert poly_eval(c, 1.0) == pytest.approx(6.0)
        assert poly_eval(c, 2.0) == pytest.approx(17.0)

    def test_r_squared_perfect_fit(self):
        """R² = 1.0 when y_pred == y_true."""
        y = np.array([1.0, 2.0, 3.0, 4.0])
        assert r_squared(y, y) == pytest.approx(1.0)

    def test_r_squared_mean_prediction(self):
        """R² = 0.0 when predicting the mean for all points."""
        y = np.array([1.0, 2.0, 3.0, 4.0])
        y_mean = np.full_like(y, np.mean(y))
        assert r_squared(y, y_mean) == pytest.approx(0.0)
