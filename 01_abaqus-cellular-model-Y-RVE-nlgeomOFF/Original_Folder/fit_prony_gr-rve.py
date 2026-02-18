from math import *
import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
from pathlib import Path
import re
import sys

# --------------------------------------
# Set working directory
# --------------------------------------
"""script_path = Path(__file__).resolve() # get path of current script
cwd = script_path.parent # get directory"""
cwd = Path(sys.argv[1])

# --------------------------------------
# Define functions
# --------------------------------------
def prony_response(comp_i, tdata, tau_0):
    # Compute the Prony series response.  
    # Args:
    # - comp_i : 1D array of prony coefficients
    # - tdata  : 1D array of time data
    # - tau_0  : 1D array of fixed retardation times
    # Returns:
    # - 1D numpy array with the response for each time point.
    comp_i = np.array(comp_i)
    return np.sum(comp_i) - np.sum(comp_i[:, None] * np.exp(-tdata / tau_0[:, None]), axis=0)

def fit_column(ydata, time, tau, maxiter=500):
    # Fit a compliance data series to the Prony series model.
    # Args:
    # - ydata  : 1D numpy array of compliance data
    # - time   : 1D numpy array of time values
    # - tau    : 1D numpy array of fixed retardation times
    # - maxiter: maximum number of iterations for the optimizer     
    # Returns:
    # - sol   : the optimized prony coefficients (accounting for sign flip)
    # - R2    : coefficient of determination for the fit
    
    y = np.array(ydata, float)
    t = np.array(time, float)

    # Handle sign‐flip
    if y.mean() < 0:
        y_fit, flip = -y, -1
    else:
        y_fit, flip = y.copy(), 1

    N = len(tau)
    p0 = np.full(N, 0.5)
    lower = np.full(N, 1e-6)
    upper = np.full(N, np.inf)

    popt, _ = curve_fit(
        lambda t, *ci: prony_response(ci, t, tau),
        time_data, y_fit, p0=p0,
        bounds=(lower, upper), max_nfev=maxiter)
    sol = popt * flip

    # Calculate R-square coefficients
    fitted = prony_response(sol, t, tau)
    SS_res = np.sum((y - fitted) ** 2)
    SS_tot = np.sum((y - y.mean()) ** 2)
    R2 = 1 - SS_res/SS_tot

    return sol, R2

# --------------------------------------
# Main code execution
# --------------------------------------
tissue = "GR"
tissue_folder = cwd / tissue

# --------------------------------------
# Load compliances
# --------------------------------------
elastic_path = tissue_folder / (tissue + "_elastic_compliance_coeffs.csv")
C0inv = pd.read_csv(elastic_path)
creep_path = tissue_folder / (tissue + "_creep_compliance_coeffs.csv")
Jc = pd.read_csv(creep_path)
# Subtract the elastic C0 values from total compliance
D_cols = [c for c in Jc.columns if c.startswith("D")]
time_cols = [c for c in Jc.columns if c.startswith("time")]
Jc[D_cols] = Jc[D_cols].subtract(C0inv.loc[0, D_cols], axis='columns')
# Prepend a zero‐row
zero_row = {**{c: 0 for c in time_cols}, **{c: 0 for c in D_cols}}
Jc = pd.concat([pd.DataFrame([zero_row]), Jc], ignore_index=True)

# --------------------------------------
# Fit creep compliances
# --------------------------------------
tau = np.array([0.1, 1.0, 10.0, 100.0])
prony_coeffs = {}
gamma_values = {}
R2_values    = {}
pairs = ['11','22','12','13','23','33','44','55','66']
for pair in pairs:
    time_data = Jc[f"time{pair}"].values
    ydata = Jc[f"D{pair}"].values
    # Identify NaN values in either time or compliance data
    nan_mask = np.isnan(time_data) | np.isnan(ydata)
    if nan_mask.any():
        # Indices where NaNs occur
        nan_indices = np.where(nan_mask)[0]
        # If NaNs are not all consecutive at the end, something went wrong
        if not np.all(nan_indices == np.arange(nan_indices[0], len(ydata))):
            raise ValueError(f"NaN found in the middle of data for {tissue} {pair}")
        # Remove trailing NaNs
        first_nan = nan_indices[0]
        time_data = time_data[:first_nan]
        ydata = ydata[:first_nan]
    sol, R2 = fit_column(ydata, time_data, tau)
    prony_coeffs[f"D{pair}"] = sol
    R2_values[f"D{pair}"] = R2
    baseline = C0inv[f"D{pair}"].iloc[0]
    gamma_values[f"D{pair}"] = baseline / np.array(sol)

# --------------------------------------
# Save gamma values
# --------------------------------------
gamma_df = pd.DataFrame(columns=['gamma1','gamma2','gamma3','gamma4'])
for pair in pairs:
    gamma_df.loc[f"D{pair}", :] = gamma_values[f"D{pair}"]
savepath = tissue_folder / (tissue + "_gamma.csv")
gamma_df.to_csv(savepath, index=False)