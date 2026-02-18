from math import *
import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
from pathlib import Path
import sys

# --------------------------------------
# Set working directory
# --------------------------------------
#script_path = Path(__file__).resolve() # get path of current script
#cwd = script_path.parent # get directory
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
tau = np.array([0.1, 1.0, 10.0, 100.0]) # Fixed retardation times [h]
pairs = ['11', '22', '33', '44', '55', '66', '12', '13', '23']

# Repeat for each tissue
for tissue in ["EW", "TW", "LW"]:
    # --------------------------------------
    # Load data
    # --------------------------------------
    tissue_folder = cwd / tissue
    elastic_path = tissue_folder / (tissue + "_elastic_compliance_coeffs.csv")
    C0inv = pd.read_csv(elastic_path)
    creep_path = tissue_folder / (tissue + "_creep_compliance_coeffs.csv")
    Jc = pd.read_csv(creep_path)
    # From that first row of D-columns, subtract the elastic C0 values
    D_cols = [c for c in Jc.columns if c.startswith("D")]
    time_cols = [c for c in Jc.columns if c.startswith("time")]
    Jc[D_cols] = Jc[D_cols].subtract(C0inv.loc[0, D_cols], axis='columns')
    # Prepend a zero‐row
    zero_row = {**{c: 0 for c in time_cols}, **{c: 0 for c in D_cols}}
    Jc = pd.concat([pd.DataFrame([zero_row]), Jc], ignore_index=True)

    # --------------------------------------
    # Fit creep compliances
    # --------------------------------------
    # Dictionaries to store results
    prony_coeffs = {}
    gamma_values = {}
    R2_values = {}
    # Loop over each pair to process its corresponding time and D data.
    for pair in pairs:
        # Construct the column names for time and compliance
        time_col = "time" + pair
        D_col = "D" + pair
        # Extract time data and compliance (D) data for the current pair.
        time_data = Jc[time_col].values
        ydata = Jc[D_col].values
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
        # Fit the data using your fit_column function, passing the paired time_data.
        sol, R2 = fit_column(ydata, time_data, tau)
        # Store the fitted coefficients and the R2 value.
        prony_coeffs[D_col] = sol
        R2_values[D_col] = R2
        # Retrieve the baseline elastic compliance from C0inv for the corresponding D coefficient.
        baseline = C0inv[D_col].iloc[0]
        # Calculate gamma values as baseline divided by the fitted prony coefficients.
        gamma_values[D_col] = baseline / np.array(sol)

    # --------------------------------------
    # Assemble results
    # --------------------------------------
    prony_cols = ["C1", "C2", "C3", "C4"]
    gamma_cols = ["gamma1", "gamma2", "gamma3", "gamma4"]
    columns = prony_cols + gamma_cols + ["R2"]
    row_labels, rows = [], []
    for col in prony_coeffs.keys():
        row_labels.append(col)
        row_data = list(prony_coeffs[col]) + list(gamma_values[col]) + [R2_values[col]]
        rows.append(row_data)

    df_results = pd.DataFrame(rows, index=row_labels, columns=columns)

    # --------------------------------------
    # Calculate and save mean gamma values
    # --------------------------------------
    gamma_means = df_results[gamma_cols].mean()
    df_gamma_new = gamma_means.to_frame().T
    savepath = tissue_folder / (tissue + "_gamma.csv")
    if not os.path.exists(savepath):
        df_gamma_new.to_csv(savepath, index=False)
    else:
        df_gamma_new.to_csv(savepath, mode='a', header=False, index=False)