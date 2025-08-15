"""
Code for precipitation data simulation preparation
Author: Vaclav Steinbach
Date: 15.08.2025
Dissertation Work
"""
import pandas as pd
import numpy as np
import os

data_dir = 'data/'
out_dir = "out/"
os.makedirs(out_dir, exist_ok=True)
out_file = "rain_tree.in"

# Read the data
filename = "canopy_rain_25072025_30072025_raw"
# filename = "free_rain_25072025_30072025_raw"
df = pd.read_csv(data_dir+filename+".csv")

# Rename columns
df.columns = ["Time", "precip"]  

# Reformat time column to datetime
df["Time"] = pd.to_datetime(df["Time"])
# Convert precipitation from mm -> m
df["precip"] = df["precip"] / 1000.0

# Calculate seconds since start
start_time = df["Time"].iloc[0]
df["seconds"] = (df["Time"] - start_time).dt.total_seconds().astype(int)

# Reorder columns to seconds + precip
rain_arr = df[["seconds", "precip"]].to_numpy()

# Save to .in file
np.savetxt(
    data_dir+out_file,
    rain_arr,
    fmt="%d %.10f",
    header="seconds precipitation [m]",
    comments="# "
)

print("Data converted!")

