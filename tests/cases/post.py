import pandas as pd
import pdb
import matplotlib.pyplot as plt
import numpy as np

out_name = "out.csv"
res = pd.read_csv(out_name)

variables = ['flow:INFLOW:branch0_seg0',
             'flow:branch0_seg0:OUT',
             'pressure:INFLOW:branch0_seg0',
             'pressure:branch0_seg0:OUT',
             'r:branch0_seg0',
             'v:branch0_seg0',
             'S:branch0_seg0',
             'tau:branch0_seg0',
             'V:branch0_seg0']
fig, axs = plt.subplots(5, 2, figsize=(12, 15))
axs = axs.ravel()

for idx, var in enumerate(variables):
    name = f"{var}"
    ids = res.name == name
    axs[idx].plot(np.array(res[ids].time), np.array(res[ids].y), label=var)
    axs[idx].set_xlabel('Time')
    axs[idx].set_ylabel(var)
    axs[idx].legend()

plt.tight_layout()
plt.savefig('output_plot.pdf')

# # # Load out.csv for time values
# out_name = "out.csv"
# res_out = pd.read_csv(out_name)

# # Load screen_out.csv for y values (assume it has no titles)
# screen_out_name = "screen_out.csv"
# res_screen = pd.read_csv(screen_out_name, header=None)  # Load without headers

# # Extract time values from out.csv
# time_name = "tau:branch0_seg0"
# time_ids = res_out.name == time_name
# time_values = np.array(res_out[time_ids].time)

# # Extract y values directly from the first column of screen_out.csv
# y_values = np.array(res_screen.iloc[0::2,0])  # Use the first column and every other row

# # Plot time (x-axis) vs y (y-axis)
# plt.plot(time_values, y_values, label="Screen Output")
# plt.xlabel("Time")
# plt.ylabel("f")
# plt.legend()
# plt.show()