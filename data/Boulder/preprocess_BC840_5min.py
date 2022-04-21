#%% SETUP
import numpy as np
from datetime import datetime
import time as time
import matplotlib.pyplot as plt
import scipy.interpolate as spi

import sys
sys.path.insert(0, '../')
import data_utils as du

plot_figs = True
save_python = True
save_matlab = True

#%% LOAD RAW DATA
# Read the whole file
data_name = 'BC840_2018-02-01_2018-02-13'
fname = data_name + '.txt'

day = 288
h_min, h_max = 100, 500
# Keep 13 days (the last 6 days are missing too many values for this file)
num_keep = 13

f = open(fname, 'r')
Lines = f.readlines()
f.close()

# Create profile lists
num_entries = len(Lines)//4    # There are 4 lines per entry (including 1 blank line in between)
timestamps = [None] * num_entries
heights = [None] * num_entries
freqs = [None] * num_entries

# Populate profiles
line_idx = 0
prof_idx = 0
missing = 0
while prof_idx<num_entries:
    timestamps[prof_idx] = '_'.join(Lines[line_idx].strip().split()[0::2])
    if Lines[line_idx+1].strip().split()[0] == 'NO':
        heights[prof_idx] = [np.nan]
        freqs[prof_idx] = [np.nan]
        missing += 1
    else:
        heights[prof_idx] = np.asarray(Lines[line_idx+1].strip().split(), dtype=float)
        freqs[prof_idx] = np.asarray(Lines[line_idx+2].strip().split(), dtype=float)
    line_idx += 4
    prof_idx += 1

#%% INTERPOLATE TO REGULAR GRID
# Create nan-filled matrix from the raw data
# Interpolate to a regular height grid
num_height = 1000
num_time = len(freqs)
ne_mat = np.ones(shape=(num_height, num_time))*np.nan
prof_heights = np.linspace(1, num_height, num_height)
for ii in range(num_time):
    if not np.isnan(heights[ii][0]):
        h_new = np.arange(np.ceil(np.min(heights[ii])), np.floor(np.max(heights[ii])), 1, dtype=int)
        interp = spi.interp1d(heights[ii], freqs[ii])
        f_new = interp(h_new)
        for jj in range(len(h_new)):
            ne_mat[h_new[jj], ii] = f_new[jj]

ne_mat = ne_mat[h_min:h_max, :num_keep*day]
timestamps = timestamps[:num_keep*day]
freqs = freqs[:num_keep*day]
heights = heights[:num_keep*day]

if plot_figs:
    # Verify the nan-filled data matrix
    plt.figure(1, figsize=(20,10))
    plt.contourf(ne_mat[:, :])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('freq (MHz)')
    plt.title('profiles interpolated to regular grid with nans above/below')

#%% IMPUTE MISSING VALUES
# Create daily sinusoid patterns:
phase = np.zeros((2,len(freqs)))
_times = np.array([time.mktime(time.strptime(ts.split('_')[-1],'%H:%M:%S')) for ts in timestamps])/60
_times -= _times[0]

phase[0,:] = np.sin(_times/24/60*2*np.pi)
phase[1,:] = np.cos(_times/24/60*2*np.pi)

if plot_figs:
    fig = plt.figure(2, figsize=(20,10))
    ax1 = plt.subplot(1, 2, 1)
    ax1.plot(_times, '.')
    ax1.set_xlabel('index'); ax1.set_ylabel('minute of day')
    ax2 = plt.subplot(1, 2, 2)
    ax2.plot(phase[0,:day], phase[1,:day], 'o')
    ax2.set_title('time of day phase')
    plt.show()

# Find the phase of the diurnal pattern, for imputation
avg = np.nanmean(ne_mat,0); avg[avg!=avg] = np.nanmean(avg);
ph = np.array([avg.dot(phase[0,:]),avg.dot(phase[1,:])])
ph /= np.sqrt((ph**2).sum())
print(ph)

# Impute
X = ne_mat.T
XX = np.hstack((ph[0]*phase[0:1,:].T + ph[1]*phase[1:2,:].T, X))
Xim, params = du.imputeMissing(XX, method='gausscopula')
ne_impute = Xim[:, 1:].T

# Plot the imputed data
xticks = np.arange(0, ne_mat.shape[1], day)
xticklabels = np.arange(0, num_keep, 1)

if plot_figs:
    fig = plt.figure(3, figsize=(20,10))
    ax1 = plt.subplot(2, 1, 1)
    con1 = ax1.contourf(X.T)
    cbar = plt.colorbar(con1)
    cbar.ax.set_ylabel('freq (MHz)')
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticklabels)
    ax2 = plt.subplot(2, 1, 2)
    con2 = ax2.contourf(ne_impute)
    cbar = plt.colorbar(con2)
    cbar.ax.set_ylabel('freq (MHz)')
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels)
    plt.show()


#%% SAVE 
ne_dict = {"ne": ne_impute,
    "heights" : [x for x in range(h_min, h_max)],
    "times" : timestamps}

if save_python:
    import pickle
    pickle.dump(ne_dict, open(data_name + '.pkl', 'wb'))

if save_matlab:
    from scipy.io import savemat
    savemat(data_name + '.mat', ne_dict)

