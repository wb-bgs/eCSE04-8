#!/usr/bin/env python
"""
Functions (top) and script (bottom) to plot power spectral of model results,
and differences to BGS "true" model.

Credit to ChaosMagPy for basis of power spectrum and correlation functions,
stripped down here to lighten environment dependencies.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import glob
import os
import sys


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


def degree_correlation(coef1, coef2):
    """
    Correlation per spherical harmonic degree between model coef1 and model
    coef2.

    Parameters
    ----------
    coef1, coef2 : ndarray, shape (N,)
        Two sets of coefficients of equal length `N`.

    Returns
    -------
    C_n : ndarray, shape (nmax,)
        Degree correlation of the two models. There are `N = nmax(nmax+2)`
        coefficients.
    """

    if coef1.ndim != 1:
        raise ValueError(f'Only 1-D input allowed {coef1.ndim} != 1')

    if coef2.ndim != 1:
        raise ValueError(f'Only 1-D input allowed {coef2.ndim} != 1')

    # Use shorter of two models to set degree limit
    if coef1.size <= coef2.size:
        ncoef = coef1.size
    else:
        ncoef = coef2.size

    nmax = int(np.sqrt(ncoef + 1) - 1)

    C_n = np.zeros((nmax,))
    R_n = np.zeros((nmax,))  # elements are prop. to power spectrum of coef1
    S_n = np.zeros((nmax,))  # elements are prop. to power spectrum of coef2

    coef12 = coef1[:ncoef]*coef2[:ncoef]

    for n in range(1, nmax+1):
        m_min = n**2 - 1
        m_max = m_min + (2*n + 1)
        R_n[n-1] = np.sum(coef1[m_min:m_max]**2)
        S_n[n-1] = np.sum(coef2[m_min:m_max]**2)
        C_n[n-1] = (np.sum(coef12[m_min:m_max]) / np.sqrt(R_n[n-1]*S_n[n-1]))

    return C_n


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
    ax : :class:`matplotlib.axes.Axes`
        A single axes instance.
    """

    """Expecting list of different length arrays now, but keep compatibility if
    given array"""
    if type(spectrum) is not list:
        spectrum = [spectrum]
        labels = [labels]

    # create axis handle
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(11,6))

    h = []
    max_deg = 0
    for i, s in enumerate(spectrum):
        degrees = np.arange(1, s.shape[0] + 1, step=1.0)
        max_deg = max([max_deg, max(degrees)])
        # Make sure spectrum is positive if differences given
        s = np.abs(s)
        h.append(ax.semilogy(degrees, s, label=labels[i]))

    ax.set_title(title)
    ax.grid(visible=True, which='minor', linestyle=':')
    ax.grid(visible=True, which='major', linestyle='-', axis='both')
    ax.set(ylabel="nT^2", xlabel='degree')

    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.xlim((0, max_deg))
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, l, loc='center left', ncol=2,
               bbox_to_anchor=(1, 0.5))
    plt.tight_layout()

    return fig, ax


def plot_degree_correlation(corr, labels, title):
    """
    Plot correlation per spherical harmonic degree between models.

    Parameters
    ----------
    corr : ndarray, shape (N,)
        Spherical harmonic correlation of degree `N`.
    labels : str
        Strings to label each series plotted.
    title : str
        String to title plot with.

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        Matplotlib figure.
    ax : :class:`matplotlib.axes.Axes`
        A single axes instance.
    """

    """Expecting list of different length arrays now, but keep compatibility if
    given array"""
    if type(corr) is not list:
        corr = [corr]
        labels = [labels]

    # create axis handle
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(11,6))

    h = []
    max_deg = 0
    for i, c in enumerate(corr):
        degrees = np.arange(1, c.shape[0] + 1, step=1.0)
        max_deg = max([max_deg, max(degrees)])
        h.append(ax.plot(degrees, c, label=labels[i]))

    ax.set_title(title)
    ax.grid(visible=True, which='minor', linestyle=':')
    ax.grid(visible=True, which='major', linestyle='-', axis='both')
    ax.set(ylabel="C_n", xlabel='degree')

    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.xlim((0, max_deg))
    h, l = ax.get_legend_handles_labels()
    plt.legend(h, l)
    plt.tight_layout()

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


###############################################################################
# v3.1 degree 2000 results analysis

#dir_root = '/users/globmod/eCSE_2021/'
dir_root = '/work/ecsead08/ecsead08/mrbad08/tests/WMAM/'

# Find all model output files below a given dir
dir_results = dir_root+'outputs/3.1'  # don't use "~" for home!

dir_figs = dir_root+'outputs/figures_3.1'
if not os.path.exists(dir_figs):
    os.makedirs(dir_figs)

paths = glob.glob(os.path.join(dir_results ,"**", "model_No_P.out"),
                  recursive=True)

# BGS degree 1440 model to use as benchmark
truth = load_tabbed(dir_root+"outputs/BGS1440/model_No_P.out")
T = coef2spec(truth)

P = []
labels = []
P.append(T)
labels.append("BGS1440")
for path in paths:
    coef = load_starred_csv(path)
    P.append(coef2spec(coef))
    """Will only work for current path setup! Should use regex if this
    directory naming convention is consistent."""
    labels.append('_'.join([path.split('/')[i] for i in [5,6,9,10]]))

fig, _ = plot_power_spectrum(P, labels, "Power spectra")
fig.savefig(os.path.join(dir_figs, "power_spectra.png"),
            format="png", dpi=300)

D = P[1:]
for i, d in enumerate(D):
    if len(d) <= len(T):
        D[i] = d - T[:len(d)]
    else:
        D[i] = d[:len(T)] - T
fig, _ = plot_power_spectrum(D, labels[1:], "Power spectra relative to BGS1440")
fig.savefig(os.path.join(dir_figs, "power_spectra_diff.png"),
            format="png", dpi=300)

C = degree_correlation(truth, load_starred_csv(paths[0]))
labels = "BGS1440 : "+labels[1]
fig, _ = plot_degree_correlation(C, labels, "Degree correlation")
fig.savefig(os.path.join(dir_figs, "degree_correlation.png"),
            format="png", dpi=300)

sys.exit()

###############################################################################
# v3.0 results analysis

# Find all model output files below a given dir
dir_results = dir_root+'outputs/3.0'  # don't use "~" for home!

dir_figs = dir_root+'outputs/figures_3.0'
if not os.path.exists(dir_figs):
    os.makedirs(dir_figs)

paths = glob.glob(os.path.join(dir_results ,"**", "model_No_P.out"),
                  recursive=True)

# BGS degree 1440 model to use as benchmark
truth = load_tabbed(dir_root+"outputs/BGS1440/model_No_P.out")
T = coef2spec(truth)

P = []
labels = []
P.append(T)
labels.append("BGS1440")
for path in paths:
    coef = load_starred_csv(path)
    P.append(coef2spec(coef))
    """Will only work for current path setup! Should use regex if this
    directory naming convention is consistent."""
    labels.append('_'.join([path.split('/')[i] for i in [5,6,9,10]]))

fig, _ = plot_power_spectrum(P, labels, "Power spectra")
fig.savefig(os.path.join(dir_figs, "power_spectra.png"),
            format="png", dpi=300)

D = P[1:]
for i, d in enumerate(D):
    if len(d) <= len(T):
        D[i] = d - T[:len(d)]
    else:
        D[i] = d[:len(T)] - T
fig, _ = plot_power_spectrum(D, labels[1:], "Power spectra relative to BGS1440")
fig.savefig(os.path.join(dir_figs, "power_spectra_diff.png"),
            format="png", dpi=300)

sys.exit()

###############################################################################
# Initial model files analysis plots:

dir_figs = dir_root+"outputs/figures_initial_runs"
if not os.path.exists(dir_figs):
    os.makedirs(dir_figs)

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
fig.savefig(os.path.join(dir_figs, "power_spectra.png"),
            format="png", dpi=300)

# Plot absolute spectra differences
fig, _ = plot_power_spectrum(np.column_stack((P1-T, P2-T, P3-T, P4-T, P5-T, P6-T, P7-T, P8-T, P9-T)),
                    ["eCSE_short", "eCSE_medium", "eCSE_long", "720",
                     "eCSE_large_problem",
                     "1024", "2048", "4096", "8192"],
                    "Power spectra realtive to BGS1440")
fig.savefig(os.path.join(dir_figs, "power_spectra_diff.png"),
            format="png", dpi=300)

sys.exit()
