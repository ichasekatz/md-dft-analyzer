"""Shared numerical utilities: polynomial least-squares fitting and error statistics.

All functions here are used by the MD and DFT analysis modules to fit
thermodynamic properties versus temperature (Cp, Cv, thermal expansion, etc.)
or DFT convergence curves (ENCUT, KPOINTS).
"""

from __future__ import annotations

import numpy as np


def gauss_pivot(a: np.ndarray, b: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """Solve a linear system via Gaussian elimination with partial pivoting.

    Args:
        a (np.ndarray): Coefficient matrix of shape ``(n, n)``. Modified in place.
        b (np.ndarray): Right-hand side vector of length ``n``. Modified in place.
        tol (float): Pivot tolerance below which the matrix is treated as singular.
            Defaults to ``1e-12``.

    Returns:
        Solution vector of length ``n``.

    Raises:
        ValueError: If the pivot element falls below ``tol``.
    """
    n = len(b)
    s = np.array([np.max(np.abs(a[i, :])) for i in range(n)])

    for k in range(n - 1):
        p = int(np.argmax(np.abs(a[k:n, k]) / s[k:n])) + k
        if p != k:
            b[[k, p]] = b[[p, k]]
            s[[k, p]] = s[[p, k]]
            a[[k, p], :] = a[[p, k], :]
        for i in range(k + 1, n):
            if a[i, k] != 0.0:
                lam = a[i, k] / a[k, k]
                a[i, k + 1 : n] -= lam * a[k, k + 1 : n]
                b[i] -= lam * b[k]

    b[n - 1] /= a[n - 1, n - 1]
    for k in range(n - 2, -1, -1):
        b[k] = (b[k] - np.dot(a[k, k + 1 : n], b[k + 1 : n])) / a[k, k]
    return b


def poly_fit(x_data: np.ndarray, y_data: np.ndarray, degree: int) -> np.ndarray:
    """Fit a polynomial of given degree to (x, y) data via least squares.

    Uses the normal-equation formulation solved by :func:`gauss_pivot`.

    Args:
        x_data (np.ndarray): Independent variable values of length ``n``.
        y_data (np.ndarray): Dependent variable values of length ``n``.
        degree (int): Polynomial degree (number of coefficients = degree + 1).

    Returns:
        Coefficient array ``[c0, c1, ..., c_degree]`` such that the
        polynomial is ``c0 + c1*x + ... + c_degree*x**degree``.
    """
    m = degree
    a = np.zeros((m + 1, m + 1))
    b = np.zeros(m + 1)
    s = np.zeros(2 * m + 1)

    for xi, yi in zip(x_data, y_data, strict=True):
        term = yi
        for j in range(m + 1):
            b[j] += term
            term *= xi
        term = 1.0
        for j in range(2 * m + 1):
            s[j] += term
            term *= xi

    for i in range(m + 1):
        for j in range(m + 1):
            a[i, j] = s[i + j]

    return gauss_pivot(a, b)


def poly_eval(coeffs: np.ndarray, x: float | np.ndarray) -> float | np.ndarray:
    """Evaluate a polynomial given its coefficients.

    Args:
        coeffs (np.ndarray): Coefficients ``[c0, c1, ..., cm]`` from
            :func:`poly_fit` (ascending power order).
        x (float | np.ndarray): Evaluation point(s).

    Returns:
        Polynomial value(s) at ``x``.
    """
    result = np.zeros_like(np.asarray(x, dtype=float))
    for i, c in enumerate(coeffs):
        result = result + c * np.asarray(x, dtype=float) ** i
    return float(result) if result.ndim == 0 else result


def r_squared(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Compute the coefficient of determination R².

    Args:
        y_true (np.ndarray): Observed values.
        y_pred (np.ndarray): Predicted values (same length as ``y_true``).

    Returns:
        R² value in ``[0, 1]``; 1.0 indicates a perfect fit.
    """
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return float(1.0 - ss_res / ss_tot)
