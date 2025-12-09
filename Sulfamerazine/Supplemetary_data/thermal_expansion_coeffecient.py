# import required libraries
import pandas as pd
import numpy as np
import sys
from scipy.stats import t

# User parameters
INPUT_CSV = sys.argv[1]  # Change as needed
TEMP_COL = 'Temperature/ K'
VOLUME_COL = 'volume / Angstrom-3'
NOTES_COL = 'notes'
INCREMENT = 10  # Kelvin
if "Form" in INPUT_CSV:
    name = f"Form_{INPUT_CSV.split('_')[1]}_"
else:
    name = ""
# Output file names
OUTPUT_SCXRD = f'thermal_expansion_coefficient_{name}SCXRD.csv'
OUTPUT_PXRD = f'thermal_expansion_coefficient_{name}PXRD.csv'
OUTPUT_ALL = f'thermal_expansion_coefficient_{name}ALL.csv'
# Read the data
try:
    df = pd.read_csv(INPUT_CSV)
except Exception as e:
    print(f"Error reading '{INPUT_CSV}': {e}")
    exit(1)

df.columns = df.columns.str.strip()

required_cols = [TEMP_COL, VOLUME_COL, NOTES_COL]
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    print(f"Error: The following required columns are missing from '{INPUT_CSV}': {missing_cols}")
    print("Please ensure your CSV contains columns named exactly as follows:")
    for col in required_cols:
        print(f"  - {col}")
    print("Or, modify this script to match your column names.")
    exit(1)
# Helper function to process a subset
def process_and_save(sub_df, output_csv):
    # Drop rows with missing temperature or volume
    sub_df = sub_df.dropna(subset=[TEMP_COL, VOLUME_COL])
    if len(sub_df) < 3:
        print(f"Not enough data for {output_csv}")
        return
    temps = sub_df[TEMP_COL].astype(float).values
    vols = sub_df[VOLUME_COL].astype(float).values
    coeffs = np.polyfit(temps, vols, 2)
    p = np.poly1d(coeffs)

    def expansion_coefficient(T):
        return 2*coeffs[0]*T + coeffs[1]

    T_min = int(np.ceil(temps.min() / INCREMENT) * INCREMENT)
    T_max = int(np.floor(temps.max() / INCREMENT) * INCREMENT)
    T_points = np.arange(T_min, T_max + 1, INCREMENT)
    exp_coeffs = expansion_coefficient(T_points)
    volumes_at_points = np.polyval(coeffs, T_points)
    norm_exp_coeffs = exp_coeffs / volumes_at_points

    # Confidence band calculation
    # Fit quadratic: vols = a*T^2 + b*T + c
    # Get predicted values
    y_pred = p(temps)
    residuals = vols - y_pred
    dof = len(temps) - 3  # 3 parameters for quadratic
    s_err = np.sqrt(np.sum(residuals**2) / dof)

    # For each T_point, calculate confidence interval
    # Design matrix for quadratic
    X = np.vstack([T_points**2, T_points, np.ones_like(T_points)]).T
    X_data = np.vstack([temps**2, temps, np.ones_like(temps)]).T
    # (X^T X)^-1
    XT_X_inv = np.linalg.inv(X_data.T @ X_data)
    # t-value for 95% confidence
    t_val = t.ppf(0.975, dof)

    conf_band = []
    for i, T in enumerate(T_points):
        x0 = np.array([T**2, T, 1])
        # Standard error of prediction
        se_pred = s_err * np.sqrt(x0 @ XT_X_inv @ x0.T)
        conf_band.append(t_val * se_pred)
    conf_band = np.array(conf_band)
    norm_conf_band = conf_band / volumes_at_points
    output_df = pd.DataFrame({
        'Identifier': np.arange(1, len(T_points) + 1),
        'Property': [f"Volumetric thermal expansion coefficient @{int(T)}K" for T in T_points],
#        'Thermal Expansion Coefficient (dV/dT)': exp_coeffs,
        'Value (1/K)': norm_exp_coeffs,
#        'Regression Fit (Volume)': volumes_at_points,
#        'Confidence Band (+/-)': conf_band,
        'Std': norm_conf_band,
        'N': "",
        'Name': "",
        'Reference': "",
        'Comment': ""
    })
    output_df.to_csv(output_csv, index=False)

# Separate SC-XRD and PXRD data
df[NOTES_COL] = df[NOTES_COL].astype(str)
scxrd_df = df[df[NOTES_COL].str.contains('SC-XRD', case=False, na=False)]
pxrd_df = df[df[NOTES_COL].str.contains('PXRD', case=False, na=False)]

has_scxrd = not scxrd_df.empty
has_pxrd = not pxrd_df.empty

if has_scxrd or has_pxrd:
    if has_scxrd:
        process_and_save(scxrd_df, OUTPUT_SCXRD)
    else:
        print("No SC-XRD data found. Explicitly labelled 'SC-XRD' in the Notes column")
    if has_pxrd:
        process_and_save(pxrd_df, OUTPUT_PXRD)
    else:
        print("No PXRD data found. Explicitly labelled 'PXRD' in the Notes column")
else:
    print("No Explicitly labelled SC-XRD or PXRD data found. Processing all data together.")
    process_and_save(df, OUTPUT_ALL)
