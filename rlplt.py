"""Radial Profile Langmuir Probe Data Plotting Module

This module contains the functions used to plot electron number density vs
radial position for the electric propulsion systems in the Advanced Propulsion
Laboratory at the University of Washington.
"""

__author__ = 'Kaito Durkee'

import sys
import os
import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import butter, lfilter
from scipy.stats import maxwell
from scipy.interpolate import CubicSpline, splev, splrep
from scipy import interpolate as inter
import re
import scipy.constants as const

from PyQt5.QtWidgets import QCheckBox

from ErrorClasses import DataError, FileError
import warnings


def get_radial_position(filename):

    # Pull out radial position from filename as a float.
    position_match = re.search(
    	        '[-]?[0-9]?[0-9][ ]?cm',filename)

    if position_match == None:
        raise ValueError("Filename format is incorrect: %r" % filename)

    position_string = position_match.group()
    position_string = "".join(position_string.split())
    position_string = position_string[:-2]
    # if position_string[-1] == '-':
    #     position_string = position_string[:-1]

    position = float(position_string)

    return position


def check_folders_in_directory(id_list):
    flag = False
    for folder in os.listdir():
        if folder[0] in id_list:
            flag = True
    if flag == False:
        raise FileError("No folders in directory are named correctly: " +
                         "CHECK DIRECTORY")


def get_data(name):

    data = {}
    os.chdir(name)

    # id is L, R, D, or T (Left, Right, Double, or Triple)
    id_list = ['L','R','D','T']
    check_folders_in_directory(id_list)

    for id in id_list:
        for id_folder in os.listdir():
            if id_folder[0] == id:
                data[id] = {}
                for position_folder in os.listdir(id_folder):
                    for shot_file in os.listdir(
                    	        id_folder + '/' + position_folder):
                        position = get_radial_position(shot_file)
                        data[id][position] = {shot_file: np.ndfromtxt(
                                id_folder + '/' + position_folder +'/' +
                                shot_file, delimiter='\t')}
    if data == {}:
    	raise DataError("Data could not be read.")
    return data


def butter_filter(data, order, cutoff):
    buttered = {}
    sos = signal.butter(order, cutoff, btype='low', 
            analog=False, output='sos')
    for id in data:
        buttered[id] = {}
        for position in data[id]:
            buttered[id][position] = {}
            for shot_file in data[id][position]:
                corrected = np.array(
                    signal.sosfiltfilt(sos,
                            data[id][position][shot_file][:, 1]))

                buttered[id][position][shot_file] = corrected
    return buttered


def butter_avg(buttered):

    avg = {}

    for id in buttered.keys():
        avg[id] = {}
        for position in buttered[id]:
            avg[id][position] = {}
            for shot_file in buttered[id][position]:
                length = len(buttered[id][position][shot_file])
            average_vals = (np.sum(
                    buttered[id][position][shot_file], axis=0) / length)

            avg[id][position] = average_vals

    return avg

def get_max_vals(avg):
	max_vals = {}

	for id in avg.keys():
		max_vals[id] = {}
		for position in avg[id]:
				max_val_at_position = np.amax(
					np.absolute(avg[id][position]))
				max_vals[id][position] = max_val_at_position
	return max_vals

def density(max_vals):
    # Cross-sectional area of the Langmuir probe used
    area_of_probe = 1.749E-5 # m^2

    # assuming argon propellant
    ion_mass = 6.6335209E-26 # kg (per ion)

    k_B = const.Boltzmann # J/K; Boltzmann constant
    q_e = const.e # C; electron charge

    temp_estimate = 10  # eV
    temp_K = temp_estimate * 1.160451812E4 # K

    proportionality_constant = m.exp(-0.5) * q_e * area_of_probe \
            * np.sqrt(k_B * temp_K / ion_mass)

    n_e = {}

    # Create a dictionary of radial position mapped to
    # max electron number density values.
    for id in max_vals.keys():
        n_e[id] = {}
        for position in max_vals[id]:
            n_e[id][position] = max_vals[id][position] \
                    / proportionality_constant
    return n_e


if __name__ == 'main':
    print('Running rlplt')
