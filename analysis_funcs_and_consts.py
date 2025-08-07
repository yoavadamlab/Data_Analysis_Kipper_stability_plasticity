
#################### IMPORTS & CONSTANTS ####################
import sys
import os
import pandas as pd
import numpy as np
import scipy as sp
from scipy.signal import hilbert
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append(r'Z:\Adam-Lab-Shared\FromExperimentToAnalysis\Rotem')
from utils import files_paths as paths
from utils import pipeline_constants as consts
from utils import data_utils as data_utils
from utils import pipeline_utils as pipe_utils
import ast
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st
import json
import glob 
import warnings
warnings.filterwarnings('ignore')
from scipy.signal import butter, lfilter,filtfilt,firwin
from scipy.fft import fft
import pickle
import math
import pyppt as ppt
from statannotations.Annotator import Annotator
from scipy.optimize import curve_fit
from datetime import datetime
from scipy import signal


### columns names ###
EXPERIMENT_DATE = "experiment_date"
CAGE = "cage"
MOUSE_NAME = "mouse_name"
SEQ = "pipeline_seq"
GOOD_CELLS = "good_cells"
CELL_TYPE = "cell_type"
STRAIN = "strain"
FOV = "FOV"
FRAME_RATE = "frame_rate"
REMAPPING = "remapping"
REMOVED_LAPS = "removed_laps"
COMMENTS = "comments"
VIDEO_PATH = "video_path"
SESSIONS_COUNTER = "sessions_counter"
CELLS_NUM = "number_of_cells"
### analysis options ###
CELLS_ACTIVITY = "Cells' activity"
FR_AND_SUB = "Firing rate & subthreshold"
LAP_FR = "Firing per lap"
REMAPP_ANALYSIS = "Remapping analysis"
ACTIVITY_PER_LAP = "Activity per lap"
LONGITUDIAL_ANALYSIS = "Longitudinal analysis"
FR_POPULATION = "Population firing rate"
### analysis columns  #####
BIN = "binned_position"
WORLD="current_World"
TIME_IN_BIN_PER_LAP = 'time_in_bin_per_lap'
SPEED_IN_BIN_PER_LAP = 'speed_in_bin_per_lap'
FR_PREFIX = "fr_cell_"
SPIKES_PREFIX = "spikes_cell_"
MEAN_FR_PREFIX = "mean_fr_over_laps_cell_"
SEM_FR_PREFIX = "sem_fr_over_laps_cell_"
NON_SPIKED_TRACE_PREFIX = "non_spiked_trace_cell_"
SUB_ACTIVITY_PREFIX = "sub_activity_cell_"
SEM_SUB_ACTIVITY_PREFIX = "sem_sub_activity_cell_"
TRACE_PREFIX="spatial_component_on_mc_cell_"
OTHER_TRACE_PREFIX="sf_on_mc_cell_"


### analysis params ###
BINS_NUM = 48

###plotting params ###
PYR_COLOR = "tab:red"
IN_COLOR = "royalblue"
PYR_COLOR_weak='indianred'
IN_COLOR_weak='cornflowerblue'
SVG_TICKS = 6
SVG_LABELS = 7






#################### FUNCTIONS ####################
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = convolve1d(y, box)
    return y_smooth
from scipy.ndimage.filters import convolve1d

def normalize(vector):
    norm_vec=(vector-vector.min())/(vector.max()-vector.min())
    return norm_vec

def save_fr_matrices(result_dict):
    result_path = r'Z:\Adam-Lab-Shared\Code\data_analysis\Analysis_Results\FR_mats'
    for exp_name, fr_matrices in result_dict.items():
        if len(fr_matrices) == 1:
            for cell_name, fr_mat in fr_matrices[0].items():
                file_name = os.path.join(result_path, exp_name + '_' + cell_name + '.npy')
                np.save(file_name, fr_mat)
        else:
            for cell_name, fr_mat in fr_matrices[0].items():
                file_name = os.path.join(result_path, exp_name + '_' + cell_name + '_fam.npy')
                np.save(file_name, fr_mat)
            for cell_name, fr_mat in fr_matrices[1].items():
                file_name = os.path.join(result_path, exp_name + '_' + cell_name + '_nov.npy')
                np.save(file_name, fr_mat)


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def smoothdata(y,kernel_size=10):
    y=np.real(y)
    kernel = np.ones(kernel_size) / kernel_size
    data_convolved = np.convolve(y, kernel, mode='same')
    return data_convolved

def psd(sig,fr,nuMax=100,nuMin=1):
    sig=sig-np.mean(sig)
    freq=np.arange(0,sig.size)*fr/sig.size
    sMax = np.where(freq > nuMax)[0][0]
    sMin = np.where(freq > nuMin)[0][0]
    psd=np.multiply(fft(sig),np.conj(fft(sig)))
    return psd[sMin:sMax],freq[sMin:sMax]


def sort_rows_by_max_index(array):
    """
    Sort the rows of a NumPy array based on the index of the maximal value in each row.

    Parameters:
        array (numpy.ndarray): Input array.

    Returns:
        numpy.ndarray: Sorted array.
    """
    max_indices = np.argmax(array, axis=1)  # Find the indices of the maximum values in each row
    sorted_indices = np.argsort(max_indices)  # Sort the indices

    sorted_array = array[sorted_indices]  # Sort the rows based on the sorted indices

    return sorted_array,sorted_indices

def sort_rows_by_min_index(array):
    """
    Sort the rows of a NumPy array based on the index of the minimal value in each row.

    Parameters:
        array (numpy.ndarray): Input array.

    Returns:
        numpy.ndarray: Sorted array.
    """
    min_indices = np.argmin(array, axis=1)  # Find the indices of the maximum values in each row
    sorted_indices = np.argsort(min_indices)  # Sort the indices

    sorted_array = array[sorted_indices]  # Sort the rows based on the sorted indices

    return sorted_array,sorted_indices

def tolerant_mean(lists):
    # Determine the maximum length of the lists
    max_len = max(len(lst) for lst in lists)
    
    # Initialize a list to store the sums and counts of each position
    sums = [0] * max_len
    counts = [0] * max_len
    
    # Iterate through each list and accumulate sums and counts
    for lst in lists:
        for i, value in enumerate(lst):
            if value is None:
                sums[i] += 0
            else:
                sums[i] += value
            counts[i] += 1
    
    # Calculate the mean for each position
    means = [sums[i] / counts[i] if counts[i] != 0 else 0 for i in range(max_len)]
    
    return means

def tolerant_sem(lists):
    # Determine the maximum length of the lists
    max_len = max(len(lst) for lst in lists)
    
    # Initialize lists to store sums, sums of squares, and counts
    sums = [0] * max_len
    sums_of_squares = [0] * max_len
    counts = [0] * max_len
    
    # Iterate through each list to accumulate sums, sums of squares, and counts
    for lst in lists:
        for i, value in enumerate(lst):
            if value is not None:
                sums[i] += value
                sums_of_squares[i] += value ** 2
                counts[i] += 1
    
    # Calculate the SEM for each position
    sems = []
    for i in range(max_len):
        if counts[i] > 1:
            mean = sums[i] / counts[i]
            variance = (sums_of_squares[i] - counts[i] * mean ** 2) / (counts[i] - 1)
            sem = math.sqrt(variance) / math.sqrt(counts[i])
            sems.append(sem)
        else:
            sems.append(0)
    
    return sems

def preprocess_behavioral_data(df, speed_threshold=50, bin_size=4):
    behav_df = df.copy(deep=True)  # Creates a deep copy
    # Shift the position by 30 virmen units to make it positive
    behav_df['position'] = behav_df['position'] + 30
    # Remove the black screen
    behav_df = behav_df[behav_df[WORLD] != 5]
    # Remove frames with speed below the threshold
    behav_df = behav_df[behav_df['speed'] >= speed_threshold]
    # Bin the position data
    bins = np.arange(0, 190 + bin_size, bin_size)
    digitized = np.digitize(behav_df['position'], bins)
    behav_df['binned_position'] = digitized
    # behav_df = behav_df.reset_index(drop=True)
    return behav_df, bins

def get_avg_spike_height_per_time(cell, df=None, time_window=15):
    """
    Calculate the average spike height every `time_window` seconds for a given cell.

    Parameters:
    - cell: Object containing cell data and metadata (e.g., frame rate).
    - df: DataFrame containing experiment data. If None, uses `cell.exp.data`.
    - time_window: Time window in seconds for averaging spike heights.

    Returns:
    - DataFrame with a new column 'mean_spike_height' containing averaged spike heights per time window.
    """
    if df is None:
        df = cell.exp.data
    df = df.copy(deep=True)
    # Prepare the trace and spikes
    trace_col_prefix = data_utils.get_trace_col_prefix(df)
    trace = df[trace_col_prefix + str(cell.cell_num)]
    # Detrend the trace
    trace = trace - butter_lowpass_filter(trace, cutoff=0.5, fs=cell.metadata['frame_rate'], order=5)
    df[trace_col_prefix + str(cell.cell_num)] = trace
    spikes_col = SPIKES_PREFIX + str(cell.cell_num)
    spikes = df[spikes_col]

    # Identify single spikes with ISI > 10ms
    ISI_ms = 10
    ISI_frames = ISI_ms * cell.metadata['frame_rate'] / 1000
    single_spikes = np.where(spikes)[0]
    single_spikes = single_spikes[np.where(np.diff(single_spikes)>ISI_frames)]

    # Time window in frames
    frames_per_window = int(time_window * cell.metadata['frame_rate'])

    # Initialize column for spike heights
    df['mean_spike_height'] = np.nan

    # Calculate average spike height per time window
    total_frames = len(df)
    for start_frame in range(0, total_frames, frames_per_window):
        end_frame = start_frame + frames_per_window

        # Handle edge case for the final window
        if end_frame > total_frames:
            end_frame = total_frames

        window_indices = range(start_frame, end_frame)

        spikes_in_window = [spike for spike in single_spikes if start_frame <= spike < end_frame]
        spike_heights = []

        for spike in spikes_in_window:
            spike_height = trace.iloc[spike] - np.min(trace.iloc[max(0, spike - 3):spike])
            spike_heights.append(spike_height)

        if spike_heights:
            mean_spike_height = np.mean(spike_heights)
        else:
            mean_spike_height = np.nan

        df.loc[window_indices, 'mean_spike_height'] = mean_spike_height
    # check if there is a time window with no spikes, and fill it with the mean of the previous and next time windows
    df['mean_spike_height'] = df['mean_spike_height'].interpolate(method='linear')
    return df

def get_hilbert_trace_on_df(cell,df=None,with_norm=False,freqs=[6,10]): #returns a list of in and out theta power, instead of a ratio
    if df is None:
        df=cell.exp.data
    sub_trace=delete_spikes_bittner_2017(cell,df)
    theta_trace=butter_bandpass_filter(sub_trace, freqs[0], freqs[1], cell.metadata[FRAME_RATE], order=3)
    if with_norm:
        theta_trace = theta_trace/df['mean_spike_height'] #make sure i ran df=get_avg_spike_height_per_time(df) before somewhere in the script
    hilbert_trace=np.abs(hilbert(theta_trace))
    #make into a series
    hilbert_trace=pd.Series(hilbert_trace)
    return hilbert_trace


def low_pass_filter_bittner2017(trace, cutoff_freq=3.0, window_size=0.2, sampling_rate=1000):
    """
    Apply a low-pass FIR filter with a given cutoff frequency and Hamming window size.

    Parameters:
    trace (array-like): Input signal trace to be filtered.
    cutoff_freq (float): Cutoff frequency of the low-pass filter in Hz.
    window_size (float): Size of the Hamming window in seconds.
    sampling_rate (int): Sampling rate of the input signal in Hz.

    Returns:
    filtered_trace (array-like): Filtered signal trace with delay compensated.
    """
    # Number of taps (filter coefficients)
    num_taps = int(window_size * sampling_rate)

    # Design the FIR filter
    fir_coeff = firwin(num_taps, cutoff=cutoff_freq, window='hamming', fs=sampling_rate)

    # Apply the FIR filter to the trace
    filtered_trace = lfilter(fir_coeff, 1.0, trace)
    
    # Compensate for the filter delay
    delay = (num_taps - 1) // 2
    filtered_trace = np.roll(filtered_trace, -delay)
    
    # Set the first 'delay' elements to zero (or use the original trace if you prefer)
    filtered_trace[:delay] = 0

    return filtered_trace

def butter_lowpass_filter(data, cutoff, fs, order):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return filtfilt(b, a, data)

def delete_spikes_bittner_2017(cell,df=None):
    """
    Delete spikes from the cell's trace. until now I used a symetric window, but they deleted 0.26ms before and 3.4ms after the spike
    I don't have this resolution, so I will delete 2ms before and 4ms after the spike
    """
    if df is None:
        df = cell.exp.data
    # Define the window size in seconds
    window_size_before = 0.002  # 2ms
    window_size_after = 0.004   # 4ms
    window_size_before_frames = int(window_size_before * cell.metadata[FRAME_RATE]) #should be 1 frame
    window_size_after_frames = int(window_size_after * cell.metadata[FRAME_RATE])   #should be 2 frames
    trace_col_prefix = data_utils.get_trace_col_prefix(df)
    spikes_col = SPIKES_PREFIX + str(cell.cell_num)
    trace = df[trace_col_prefix + str(cell.cell_num)]
    # Find the indices of the spikes
    spike_inds = np.where(df[spikes_col])[0]
    # Create a copy of the trace
    trace_clean = trace.copy()
    # Delete the spikes and interpulate
    for spike_ind in spike_inds:
        # Delete the spike and the window around it
        trace_clean[spike_ind - window_size_before_frames:spike_ind + window_size_after_frames] = np.nan
        # Interpolate the deleted part
        trace_clean = np.interp(np.arange(len(trace_clean)), np.arange(len(trace_clean))[~np.isnan(trace_clean)], trace_clean[~np.isnan(trace_clean)])
    return trace_clean


def calc_spatial_info_Pval_on_FR_mat_new(FR_mat, time_mat):    #instead of calculating Zscore, I will calculate the p value by counting the number of times 
    SI=calc_spatial_info_on_mat_new(FR_mat, time_mat,smooth_size=3)         #the shuffled SI values are greater than the original SI value
    shuffled_si_values=create_n_shuffled_SI_values_new(FR_mat, time_mat, num_of_permutations=1000)  # to keep simulated FR vectors of every perm
    SI_95 = np.percentile(shuffled_si_values, 95)
    return SI,SI_95

def calc_spatial_info_Zscore_on_FR_mat_new(FR_mat, time_mat):    #instead of calculating Zscore, I will calculate the p value by counting the number of times 
    SI=calc_spatial_info_on_mat_new(FR_mat, time_mat,smooth_size=3)         #the shuffled SI values are greater than the original SI value
    shuffled_si_values=create_n_shuffled_SI_values_new(FR_mat, time_mat, num_of_permutations=1000)  # to keep simulated FR vectors of every perm
    SI_Z = (SI-np.mean(shuffled_si_values))/np.std(shuffled_si_values)
    return SI_Z

def calc_spatial_info_on_mat_new(FR_mat, time_mat,smooth_size=3):
    FR_vec = np.nanmean(FR_mat, 0)
    FR_vec = np.convolve(FR_vec, np.ones(smooth_size), 'same') / smooth_size
    mean_time_vec = np.nanmean(time_mat, 0)
    pi_vec = mean_time_vec/mean_time_vec.sum() 
    r = np.dot(FR_vec,pi_vec)
    si = 0
    for ri,pi in zip(FR_vec, pi_vec):
            if ri !=0:
                    si += pi*ri*np.log2(ri/r) # running this formula without ri/r outside the log
    return si

def create_n_shuffled_SI_values_new(original_FR_matrix,original_time_matrix, num_of_permutations):
    lap_number=original_FR_matrix.shape[0]
    #circshift each lap and calc SI
    shuffled_SI_values=[]
    for i in range(num_of_permutations):
        shuffled_FR_matrix=np.zeros(original_FR_matrix.shape)
        for lap in range(lap_number):
            shuffled_FR_matrix[lap,:]=np.roll(original_FR_matrix[lap,:], np.random.randint(original_FR_matrix.shape[1]))
        shuffled_SI_values.append(calc_spatial_info_on_mat_new(shuffled_FR_matrix, original_time_matrix))
    return shuffled_SI_values


def get_FR_and_time_mats(cell,df=None, speed_threshold=50, bin_size=4, smooth_size=3):
    FR_mat = np.array([]) ; Time_mat = np.array([])
    if df is None:
        df = cell.exp.data
    # data preprocessing
    behav_df = df.copy()
    behav_df = behav_df[behav_df[WORLD]!=5] # remove black screen
    behav_df.position = behav_df.position+30 # shift the position by 30 virmen units to make it positive, don't think this is necessary
    behav_df = behav_df[behav_df['speed']>=speed_threshold] # remove frames with speed<speed_threshold
    bins = np.arange(0, 190+bin_size, bin_size) # bin the position from 0 to 190 into bins of size bin_size
    digitized = np.digitize(behav_df.position, bins) # digitize the position
    behav_df['binned_position'] = digitized # add binned_position column to behav_df
    bins_num = len(bins)-1
    spikes_col = SPIKES_PREFIX+str(cell.cell_num)

    for lap in behav_df.lap_counter.unique():
        lap_df = behav_df[behav_df['lap_counter']==lap] # get the data for each lap
        num_spikes_per_position = lap_df.groupby('binned_position')[spikes_col].sum() # sum the spikes in each bin
        num_frames_per_position = lap_df.groupby('binned_position')[spikes_col].count() # count the number of frames in each bin

        if len(num_frames_per_position)!=bins_num:
            num_frames_per_position = num_frames_per_position.reindex(range(bins_num), fill_value=0) # fill the missing bins with 0
            num_spikes_per_position = num_spikes_per_position.reindex(range(bins_num), fill_value=0) # fill the missing bins with 0
        
        fr_by_position_per_lap = np.array(num_spikes_per_position)/np.array(num_frames_per_position)*cell.metadata[FRAME_RATE] # calculate the firing rate in each bin
        #replace nan with 0
        fr_by_position_per_lap = np.nan_to_num(fr_by_position_per_lap) 

        time_spent_by_position = np.array(num_frames_per_position)/cell.metadata[FRAME_RATE] # calculate the time spent in each bin
        
        if lap==np.min(behav_df.lap_counter.unique()): # if it's the first lap
            FR_mat = fr_by_position_per_lap            # initialize FR_mat
            Time_mat = time_spent_by_position        # initialize Time_mat
        else:
            FR_mat = np.vstack((FR_mat,fr_by_position_per_lap)) # stack the firing rate of each lap
            Time_mat = np.vstack((Time_mat,time_spent_by_position)) # stack the time spent in each bin of each lap
    return FR_mat,Time_mat

def get_FR_and_time_mats_after_preprocessing(cell,bins,df=None):
    FR_mat = np.array([]) ; Time_mat = np.array([])
    if df is None:
        df = cell.exp.data
    
    bins_num = len(bins)-1
    spikes_col = SPIKES_PREFIX+str(cell.cell_num)

    for lap in df.lap_counter.unique():
        lap_df = df[df['lap_counter']==lap] # get the data for each lap
        num_spikes_per_position = lap_df.groupby('binned_position')[spikes_col].sum() # sum the spikes in each bin
        num_frames_per_position = lap_df.groupby('binned_position')[spikes_col].count() # count the number of frames in each bin

        if len(num_frames_per_position)!=bins_num:
            num_frames_per_position = num_frames_per_position.reindex(range(bins_num), fill_value=0) # fill the missing bins with 0
            num_spikes_per_position = num_spikes_per_position.reindex(range(bins_num), fill_value=0) # fill the missing bins with 0
        
        fr_by_position_per_lap = np.array(num_spikes_per_position)/np.array(num_frames_per_position)*cell.metadata[FRAME_RATE] # calculate the firing rate in each bin
        #replace nan with 0
        fr_by_position_per_lap = np.nan_to_num(fr_by_position_per_lap) 

        time_spent_by_position = np.array(num_frames_per_position)/cell.metadata[FRAME_RATE] # calculate the time spent in each bin
        
        if lap==np.min(df.lap_counter.unique()): # if it's the first lap
            FR_mat = fr_by_position_per_lap            # initialize FR_mat
            Time_mat = time_spent_by_position        # initialize Time_mat
        else:
            FR_mat = np.vstack((FR_mat,fr_by_position_per_lap)) # stack the firing rate of each lap
            Time_mat = np.vstack((Time_mat,time_spent_by_position)) # stack the time spent in each bin of each lap
    return FR_mat,Time_mat

def get_PF_bins_on_FR_mat(FR_mat,smooth_size=3):
    FR_vec = np.nanmean(FR_mat, 0)
    FR_vec = smooth(FR_vec, smooth_size)
    threshold=FR_vec.mean() + 0.10 * (FR_vec.max() - FR_vec.mean())
    a,b=get_neighboring_indexes(FR_vec,threshold)
    if len(b)==0:
        in_field_indices=a
    else:
        in_field_indices=np.concatenate((a,b),axis=0)
    return in_field_indices


def calc_in_out_ratio_on_FR_mat(FR_mat,smooth_size=3):
    FR_vec = np.nanmean(FR_mat, 0)
    FR_vec = np.convolve(FR_vec, np.ones(smooth_size), 'same') / smooth_size
    threshold=FR_vec.mean() + 0.10 * (FR_vec.max() - FR_vec.mean())
    a,b=get_neighboring_indexes(FR_vec,threshold)
    norm_vec=(FR_vec-FR_vec.min())/(FR_vec.max()-FR_vec.min())#adding normalization
    if len(b)==0:
        in_field_norm_vec=norm_vec[a]
        in_field_indices=a
    else:
        in_field_norm_vec = np.concatenate((norm_vec[a],norm_vec[b]),axis=0)
        in_field_indices=np.concatenate((a,b),axis=0)
    out_field_norm_vec = [ elem for elem in norm_vec if elem not in in_field_norm_vec]
    in_fr=np.mean(in_field_norm_vec)
    out_fr=np.mean(out_field_norm_vec)
    ratio=in_fr/out_fr
    return ratio,in_field_indices,in_fr,out_fr

def calc_in_out_theta_ratio_on_df(cell,df,freqs=[6,10]):
    df = get_avg_spike_height_per_time(cell,df)
    hlibert_trace = get_hilbert_trace_on_df(cell,df,with_norm=True,freqs=freqs)
    df['hilbert_trace'] = hlibert_trace
    df,bins = preprocess_behavioral_data(df)
    FR_mat,Time_mat=get_FR_and_time_mats_after_preprocessing(cell,bins,df)
    in_field_bins=get_PF_bins_on_FR_mat(FR_mat)
    #find the indices of elements that are in the field using numpy, not for loop
    in_field_indices=np.where(np.isin(df.binned_position,in_field_bins))[0]
    out_field_indices=np.where(~np.isin(df.binned_position,in_field_bins))[0]
    in_field_trace = df.hilbert_trace.iloc[in_field_indices]
    out_field_trace = df.hilbert_trace.iloc[out_field_indices]
    # hilbert_trace=get_hilbert_trace_on_df(cell,with_norm=True,freqs=freqs)
    #calculate the theta power for in and out field indices
    mean_power_in_field=np.mean(in_field_trace)
    mean_power_out_field=np.mean(out_field_trace)
    return mean_power_in_field,mean_power_out_field


def calc_reliability_on_FR_mat(FR_mat,smooth_size=3):
    ratio,in_field_indices,in_fr,out_fr=calc_in_out_ratio_on_FR_mat(FR_mat,smooth_size)
    cell_pf=in_field_indices
    mean_fr=np.mean(FR_mat)
    std_fr=np.std(FR_mat)
    significant_bins=np.where(FR_mat>mean_fr+2*std_fr)
    bool_mask=np.isin(significant_bins[1], cell_pf)
    laps_with_sig_bins=significant_bins[0][bool_mask]
    laps_with_sig_bins=np.unique(laps_with_sig_bins)
    reliability=len(laps_with_sig_bins)/len(FR_mat)*100
    return reliability

def in_out_ratio_on_FR_vec(FR_vec):  
    threshold=FR_vec.mean() + 0.10 * (FR_vec.max() - FR_vec.mean())
    a,b=get_neighboring_indexes(FR_vec,threshold)
    norm_vec=(FR_vec-FR_vec.min())/(FR_vec.max()-FR_vec.min())#adding normalization
    if len(b)==0:
        in_field_norm_vec=norm_vec[a]
        in_field_indices=a
    else:
        in_field_norm_vec = np.concatenate((norm_vec[a],norm_vec[b]),axis=0)
        in_field_indices=np.concatenate((a,b),axis=0)
    out_field_norm_vec = [ elem for elem in norm_vec if elem not in in_field_norm_vec]
    in_fr=np.mean(in_field_norm_vec)
    out_fr=np.mean(out_field_norm_vec)
    ratio=in_fr/out_fr
    return ratio,in_field_indices,in_fr,out_fr

def get_neighboring_indexes(FR_vec,threshold):
            max_index = np.argmax(FR_vec)
            neighbor_indexes = [max_index]
            start_indexes=[] #i saparated the end and start of PFs that cross from end to start. this is temporary for plotting reasons
            for i in range(max_index + 1, len(FR_vec)):
                if FR_vec[i] < threshold:
                    break
                neighbor_indexes.append(i)
            if len(FR_vec)-1 in neighbor_indexes:
                for i in range(0 , max_index):
                    if FR_vec[i] < threshold:
                        break
                    start_indexes.append(i)
                            
            for i in range(max_index - 1, -1, -1):
                if FR_vec[i] < threshold:
                    break
                neighbor_indexes.append(i)
            return np.sort(neighbor_indexes),np.sort(start_indexes)



def calc_speed_corr_on_df(cell,df=None):
    #using 1s time bin (frame rate is cell.metadata['frame_rate']), calculate the firing rate in each bin (adopted the second approach in this notebook)
    if df is None:
        df=cell.exp.data
    # df=df[(~df['position'].between(107,128))] #removing RZ bins
    df=df[df['speed']>=20]
    spikes_col="spikes_cell_"+str(cell.cell_num)
    #bin the data into 100ms bins
    df=df.reset_index(drop=True)
    df['time']=df['TS_time']-df['TS_time'][0]
    df['binned_time']=pd.cut(df['time'],bins=round(len(df)/cell.metadata['frame_rate']),labels=False)
    #count the number of spikes in each bin
    df['spikes_per_bin']=df.groupby('binned_time')[spikes_col].transform('sum')
    df['speed_per_bin']=df.groupby('binned_time')['speed'].transform('mean')
    #bin the speed per bin into 5 bins
    df['binned_speed']=pd.cut(df['speed_per_bin'],bins=20,labels=False)
    #count the number of spikes in each bin
    spikes_per_bin=df.groupby('binned_speed')['spikes_per_bin'].sum()
    #count the number of frames in each bin
    frames_per_bin=df.groupby('binned_speed')['spikes_per_bin'].count()
    #take into account only bins with more than 1% of the frames
    spikes_per_bin=spikes_per_bin[frames_per_bin>0.01*len(df)]
    frames_per_bin=frames_per_bin[frames_per_bin>0.01*len(df)]
    #divide to get the firing rate in each bin
    FR_per_bin=(spikes_per_bin/frames_per_bin)
    #turn ndarray to pandas series
    FR_per_bin=pd.Series(FR_per_bin,index=spikes_per_bin.index)
    #preform linear regression on the data
    slope, intercept, r_value, p_value, std_err = stats.linregress(FR_per_bin.index,FR_per_bin)
    return r_value


def is_positively_speed_tuned_on_df(cell,p_val_threshold=0.001,df=None):  #need to adjust this function to work with "calc_FR_per_speed_bin_on_df" function
    #using 1s time bin (frame rate is cell.metadata['frame_rate']), calculate the firing rate in each bin (adopted the second approach in this notebook)
    if df is None:
        df=cell.exp.data
    # df=df[(~df['position'].between(107,128))] #removing RZ bins
    df=df[df['speed']>=20]
    spikes_col="spikes_cell_"+str(cell.cell_num)
    spike_times=cell.spikes
    #bin the data into 100ms bins
    df=df.reset_index(drop=True)
    df['time']=df['TS_time']-df['TS_time'][0]
    df['binned_time']=pd.cut(df['time'],bins=round(len(df)/cell.metadata['frame_rate']),labels=False)
    #count the number of spikes in each bin
    df['spikes_per_bin']=df.groupby('binned_time')[spikes_col].transform('sum')
    df['speed_per_bin']=df.groupby('binned_time')['speed'].transform('mean')
    #bin the speed per bin into 5 bins
    df['binned_speed']=pd.cut(df['speed_per_bin'],bins=20,labels=False)
    #count the number of spikes in each bin
    spikes_per_bin=df.groupby('binned_speed')['spikes_per_bin'].sum()
    #count the number of frames in each bin
    frames_per_bin=df.groupby('binned_speed')['spikes_per_bin'].count()
    #take into account only bins with more than 1% of the frames
    spikes_per_bin=spikes_per_bin[frames_per_bin>0.01*len(df)]
    frames_per_bin=frames_per_bin[frames_per_bin>0.01*len(df)]
    #divide to get the firing rate in each bin
    FR_per_bin=(spikes_per_bin/frames_per_bin)
    # FR_per_bin=smooth(FR_per_bin,3)
    #turn ndarray to pandas series
    FR_per_bin=pd.Series(FR_per_bin,index=spikes_per_bin.index)

    #preform linear regression on the data
    slope, intercept, r_value, p_value, std_err = stats.linregress(FR_per_bin.index,FR_per_bin)
    if p_value<p_val_threshold and r_value>0.3:
        return True
    else:
        return False

def is_positively_speed_tuned_on_cell(cell,p_val_threshold=0.001,laps_to_consider=10):
    is_tuned_all = is_positively_speed_tuned_on_df(cell,p_val_threshold)
    f_df,n_df = cell.get_fam_and_novel_df_partial(laps_to_consider=laps_to_consider)
    is_tuned_F=is_positively_speed_tuned_on_df(cell,p_val_threshold,df=f_df)
    is_tuned_N=is_positively_speed_tuned_on_df(cell,p_val_threshold,df=n_df)
    return [is_tuned_all,is_tuned_F,is_tuned_N]    



    
def get_shuffled_ci_of_FR_mat(FR_mat,num_shuffles=1000):
    #in each permutation, roll each row of the FR_mat by a random number of bins
    for i in range(num_shuffles):
        shuffled_FR_mat = np.array([np.roll(row, np.random.randint(100)) for row in FR_mat])
        shuffled_FR_vec = shuffled_FR_mat.mean(axis=0)
        # shuffled_FR_vec = smooth(shuffled_FR_vec,5)
        if i==0:
            shuffled_FR_vecs=shuffled_FR_vec
        else:
            shuffled_FR_vecs=np.vstack((shuffled_FR_vecs,shuffled_FR_vec))

    #calculate the 95% confidence interval for each bin
    ci = np.percentile(shuffled_FR_vecs, [10, 90], axis=0)
    return ci

def is_spatially_modulated_on_FR_mat(FR_mat,smooth_size=3):
    is_modulated = False
    ci = get_shuffled_ci_of_FR_mat(FR_mat,num_shuffles=1000)
    FR_vec = FR_mat.mean(axis=0)
    FR_vec = smooth(FR_vec,smooth_size)
    #mark with asteriks the bins that are outside the confidence interval, only if they are outside for more than 3 consecutive bins
    for i in range(len(FR_vec)):
        if FR_vec[i] < ci[0][i] or FR_vec[i] > ci[1][i]:
            if i==0:
                consecutive=1
            else:
                if FR_vec[i-1] < ci[0][i-1] or FR_vec[i-1] > ci[1][i-1]:
                    consecutive+=1
                else:
                    consecutive=1
            # if consecutive>=3:
            if consecutive>=5:  #changed from 3 to 4 to match the new number of bins (4.2.25)
                is_modulated = True
    return is_modulated

def bin_binary_data_by_time(cell,col_name,df=None,time_window=1):  #time window is in seconds
    if df is None:
        df=cell.exp.data
    df=df.reset_index(drop=True)
    df['time']=df['TS_time']-df['TS_time'][0]
    df['binned_time']=pd.cut(df['time'],bins=round(len(df)/cell.metadata['frame_rate']/time_window),labels=False)
    col_name=col_name+str(cell.cell_num)
    binned_data=df.groupby('binned_time')[col_name].sum()
    return binned_data



        









