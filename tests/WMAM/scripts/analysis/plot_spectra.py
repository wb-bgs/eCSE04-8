#!/usr/bin/env python
"""
Functions (top) and script (bottom) to plot power spectral of model results,
and differences to BGS "true" model.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def coef2spec(coef, *, n_max=None):
    """
    Compute power spectrum of vector of spherical harmonic coefficients,
    assumed to be for internal source at Earth's surface.

    Parameters
    ----------
    coef : ndarray, shape (N*(N+2),)
        Spherical harmonic cofficients to degree `N`.
    n_max : int, optional
        Maximum spherical harmonic degree to calculate spectrum to
        (defaults to `N`). Spectra nan padded to n_max if needed.

    Returns
    -------
    P_n : ndarray, shape (n_max,)
        Power spectrum for to degree `n_max`.
    """

    # Assume model and plot at Earth's surface
    rho = 1

    N = int(np.sqrt(coef.shape[-1] + 1) - 1)
    if n_max is None:
        n_max = N

    if N > n_max:
        N = n_max

    def factor(n, rho):
        # Internal source
        return (n+1)*rho**(2*n+4)

    P_n = np.empty(coef.shape[:-1] + (n_max,))
    P_n[:] = np.nan

    for n in range(1, N+1):
        idx_m_min = n**2 - 1
        idx_m_max = idx_m_min + (2*n + 1)
        P_n[..., n-1] = factor(n, rho)*np.sum(coef[...,
                         idx_m_min:idx_m_max]**2, axis=-1)

    return P_n


def plot_power_spectrum(spectrum, labels, title):
    """
    Plot spherical harmonic spectrum.

    Parameters
    ----------
    spectrum : ndarray, shape (N,)
        Spherical harmonics spectrum of degree `N`.
    labels : str
        Strings to label each series plotted.
    title : str
        String to title plot with.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        Matplotlib figure.
    axes : :class:`matplotlib.axes.Axes`
        A single axes instance.
    """

    degrees = np.arange(1, spectrum.shape[0] + 1, step=1.0)

    # Make sure spectrum is positive if differences given
    #spectrum[spectrum == 0] = np.nan  # remove non-positive values for log
    spectrum = np.abs(spectrum)

    # create axis handle
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(8,6))

    h = ax.semilogy(degrees, spectrum)
    ax.set_title(title)
    ax.grid(visible=True, which='minor', linestyle=':')
    ax.grid(visible=True, which='major', linestyle='-', axis='both')
    ax.set(ylabel="nT^2", xlabel='degree')

    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.xlim((0, degrees[-1]))
    plt.tight_layout()
    plt.legend(h, labels)

    return fig, ax


def load_starred_csv(fname):
    """
    Load spherical harmonic model file e.g. model_No_P.out when file is a csv.
    Specifically a csv that is missing the comma after the 1st column and
    occasionally has a merged 2nd and 3rd column entry separated by a * for
    some reason... Fortran write statement wrong width?

    Expects 4 line header then rows of "type deg, ord, coef, error".

    Parameters
    ----------
    fname : str
        File to load.

    Returns
    -------
    coef : ndarray, shape (num coeffs,)
        Array of coefficients.

    """
    with open(fname, 'rb') as f:
        for _ in range(4):
            next(f)
        clean_lines = (line.replace(b'*',b',') for line in f)
        coef = np.genfromtxt(clean_lines, delimiter=',', usecols=2)

    return coef


def load_tabbed(fname):
    """
    Load spherical harmonic model file e.g. model_No_P.out when file is tab
    delimited.

    Expects 4 line header then rows of "type deg ord coef error".

    Parameters
    ----------
    fname : str
        File to load.

    Returns
    -------
    coef : ndarray, shape (num coeffs,)
        Array of coefficients.
    """
    with open(fname, 'rb') as f:
        for _ in range(4):
            next(f)
        coef = np.genfromtxt(f, dtype=None, usecols=3)

    return coef


# Load files
# BGS degree 1440 model to use as benchmark
truth = load_tabbed("BGS1440/model_No_P.out")
# Output files located in "/home/ecsead08/ecsead08/shared/arc/"
coef1 = load_starred_csv("eCSE_short/model_No_P.out")
coef2 = load_starred_csv("eCSE_medium/model_No_P.out")
coef3 = load_starred_csv("eCSE_long/model_No_P.out")
coef4 = load_tabbed("720/model_No_P.out")
coef5 = load_starred_csv("eCSE_large_problem/model_No_P.out")
coef6 = load_starred_csv("mod_wdmam_4Nick_1024/model_No_P.out")
coef7 = load_starred_csv("mod_wdmam_4Nick_2048/model_No_P.out")
coef8 = load_starred_csv("mod_wdmam_4Nick_4096/model_No_P.out")
coef9 = load_starred_csv("mod_wdmam_4Nick_8192/model_No_P.out")

# Calculate power spectra, padded to deg 2000
N = 2000
T = coef2spec(truth, n_max=N)
P1 = coef2spec(coef1, n_max=N)
P2 = coef2spec(coef2, n_max=N)
P3 = coef2spec(coef3, n_max=N)
P4 = coef2spec(coef4, n_max=N)
P5 = coef2spec(coef5, n_max=N)
P6 = coef2spec(coef6, n_max=N)
P7 = coef2spec(coef7, n_max=N)
P8 = coef2spec(coef8, n_max=N)
P9 = coef2spec(coef9, n_max=N)

# Plot spectra
fig, _ = plot_power_spectrum(np.column_stack((P1, P2, P3, P4, P5, P6, P7, P8, P9, T)),
                    ["eCSE_short", "eCSE_medium", "eCSE_long", "720",
                     "eCSE_large_problem",
                     "1024", "2048", "4096", "8192", "BGS1440"],
                    "Power spectra")
fig.savefig("power_spectra.png", format="png", dpi=300)

# Plot absolute spectra differences
fig, _ = plot_power_spectrum(np.column_stack((P1-T, P2-T, P3-T, P4-T, P5-T, P6-T, P7-T, P8-T, P9-T)),
                    ["eCSE_short", "eCSE_medium", "eCSE_long", "720",
                     "eCSE_large_problem",
                     "1024", "2048", "4096", "8192"],
                    "Power spectra realtive to BGS1440")
fig.savefig("power_spectra_diff.png", format="png", dpi=300)
