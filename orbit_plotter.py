#! /usr/bin/env python
# @file orbit_plotter.py
# @author John Keeling
# @date 6th December 2021
# Generates plots from odl20lite (UCL Orbit Prediction Software) output files 
# ...and converted CHAMP GPS files. Requires data files to be saved in 
# ...plot_input_files folder, two levels below orbit_plotter.py.

import matplotlib.pyplot as plt
from collections import OrderedDict
import math
import os
import numpy as np


def obtain_plot_title():

    # Returns plot title from first of either TLE or GPS filenames
    plot_title = ''
    title_count = 0
    filecheck = ((x for x in os.listdir("../../plot_input_files/") 
            if any(y in x for y in ('GPS', 'TLE'))
            and title_count == 0)
            )
    for x in filecheck:
        plot_title = (
            ' Simulated Orbits vs Actual, Starting Day ' 
            + x.split('_')[3] 
            + ' of ' + x.split('_')[2]
            )
        title_count += 1

    return plot_title


def load_orbit_data(y_axis):

    # Returns full contents of each file in plot_input_files as item in list
    file_data = OrderedDict()
    plot_data = OrderedDict()
    filecheck = ((x for x in os.listdir("../../plot_input_files/") 
            if not 'DS_Store' in x)
            )
    for x in filecheck:
        file_path = "../../plot_input_files/" + x
        with open(file_path) as f:
            lines = f.readlines()
        f.close()

        # Populates dict with x and y-axis plot data from each file
        if not 'GPS' in x:
            if not 'TLE' in x:
                plot_data[x] = filter_OPS_data(
                    lines, x, y_axis)
            else:
                plot_data[x] = filter_TLE_data(
                    lines, x, y_axis)
        else:
            plot_data[x] = filter_GPS_data(
                lines, x, y_axis)
    return plot_data


def apogee_perigee_check(extrema_check):

    # Checks central item of three neighbouring altitude or semi-major axis 
    # ...data points to see if maximum or minimum
    if (
            max(extrema_check) == extrema_check[1]
            or min(extrema_check) == extrema_check[1]
            ):
        apcheck = True
    else:
        apcheck = False
    
    return apcheck

def filter_OPS_data(data_lines, filename, y_axis):

    # Returns OPS x-axis (timestamp) and y-axis (alt / sma) data
    if y_axis == 'altitude':
        split_index = 16
    elif y_axis == 'semi-major':
        split_index = 8
    else:
        print('Exception: invalid y_axis value: ', y_axis)
        exit()
    base_time = data_lines[1].split()[1].rstrip(",")
    file_data = OrderedDict()
    neighbours_offset = [-1,0,1]
    for count in enumerate(data_lines[2:-1]):
        extrema_check = []
        for i in neighbours_offset:
            y_val = (
                data_lines[count[0] + i + 2].split()[split_index].rstrip(",")
                )
            extrema_check.append(float(y_val))
        if apogee_perigee_check(extrema_check) == True:
            str_mjd = data_lines[count[0] + i + 2].split()[1].rstrip(",")
            key = round((float(str_mjd) - float(base_time)) * 24, 6)
            file_data[key] = round(extrema_check[1],6)

    return file_data
    

def filter_TLE_data(data_lines, filename, y_axis):

    # Returns TLE x-axis (timestamp) and y-axis (alt / sma) data
    file_data = OrderedDict()
    base_time = data_lines[1].split()[1].rstrip(",")
    if y_axis == 'altitude':
        split_index = 14
    elif y_axis == 'semi-major':
        split_index = 8
    else:
        print('Exception: invalid y_axis value: ', y_axis)
        exit()
    for count in enumerate(data_lines[1:]):
        y_val = data_lines[count[0] + 1].split()[split_index].rstrip(",")
        str_mjd = data_lines[count[0] + 1].split()[1].rstrip(",")
        key = round((float(str_mjd) - float(base_time)) * 24, 6)
        file_data[key] = round(float(y_val),6)

    return file_data


def filter_GPS_data(data_lines, filename, y_axis):

    # Returns GPS x-axis (timestamp) and y-axis (alt / sma) data

    # Time at start
    base_mjd = float(data_lines[0][:5])
    base_deci = float(data_lines[0][6:15].strip())/86400
    base_time = base_mjd + base_deci

    # Extract y-axis data and check if apogee or perigee
    file_data = OrderedDict()
    neighbours_offset = [-1,0,1]
    for count in enumerate(data_lines[1:-1]):
        extrema_check = []
        for i in neighbours_offset:
            if 'ECI' in data_lines[count[0] + i + 1].split()[0]:
                split_index = 1
            else:
                split_index = 2
            if y_axis == 'semi-major':
                str_x, str_y, str_z, str_u, str_v, str_w = (
                    [data_lines[count[0] + i + 1].split()[split_index 
                    + j] for j in range(6)]
                    )
                y_val = compute_semi_maj_axis(
                    float(str_x), float(str_y), float(str_z), 
                    float(str_u), float(str_v), float(str_w)
                    )
            elif y_axis == 'altitude':
                y_val = data_lines[count[0] + i + 1].split()[split_index + 6]
            else:
                print('Exception: invalid y_axis value: ', y_axis)
                exit()
            extrema_check.append(y_val)

        # Pair timestamp with y-value for apogee and perigee data:
        if apogee_perigee_check(extrema_check) == True:
            mjd = float(data_lines[count[0] + 1][:5])
            deci = float(data_lines[count[0] + 1][6:15].strip()) / 86400
            time = mjd + deci
            key = round((time - base_time) * 24, 6)
            file_data[key] = round(extrema_check[1],6)

    return file_data


def compute_semi_maj_axis(x, y, z, u, v, w):

    # Calculates semi-major axis from position and velocity (route a)
    r = math.sqrt(x**2 + y**2 + z**2)
    v2 = u**2 + v**2 + w**2
    # Grav. scale parameter, from UCL OPS: Keplerian_elements.cpp:
    GM = 398600.4415  
    sma =  r * GM / (2.0 * GM - v2 * r)

    return sma


def compute_propoagation_of_error(init_sma, step_size, no_of_steps, 
        a, b, c, d, e, f, g):

    x = 6

    # Calculates stepwise error propagation for initial sma 
    # ...from uncertainties a, b, c,...


def plot_ols_reg(x_lab, y_lab, title = '', scatter = False, 
        start_aligned = False): 
    
    # Ordinary least squares regression plot with optional scatter dots

    # Fetch plot data and y-intercept values for predicted orbits
    if 'ajor' in y_lab:
        y_axis = 'semi-major'
    elif 'ltitude' in y_lab:
        y_axis = 'altitude'
    if title == '':
        title = obtain_plot_title()
    input_data = load_orbit_data(y_axis)
    y_intercept = 0.0
    filecheck = False
    gen = (j for j in input_data if 'GPS' not in j)
    for x in gen:
        y_intercept += input_data[x].values()[0]
        filecheck = True

    # For error bars, see 
    # https://jakevdp.github.io/PythonDataScienceHandbook/04.03-errorbars.html

    # Colour ordering
    colours = [
        'aqua',
        'crimson', 
        'blue',
        'blueviolet',
        'brown',
        'cadetblue',
        'chartreuse',
        'black',
        'darkgoldenrod',
        'darkgreen',
        'darkslateblue',
        'darkslategray'
        ]

    # Plot scatter points and reg line
    plt.figure(figsize=(18, 6), dpi=80).patch.set_facecolor('white')
    for file_key, colour in zip(input_data, enumerate(
            colours[:len(input_data)*2])
            ):
        filelabel = file_key            
        x = np.array([i for i in input_data[file_key].keys()])
        y = np.array([j for j in input_data[file_key].values()])
        m, b = np.polyfit(x, y, 1)
        if start_aligned == True and filecheck == True:
                b =  y_intercept / (len(input_data) -  1)
        if scatter == True:
            plt.scatter(x,y,s=12, color = colours[colour[0]])
        plt.plot(x, m*x + b, label = filelabel, color = colours[colour[0]])

    # Format plot axes, gridlines, title and legend
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.grid(axis = 'y', color = 'grey', linestyle = '--', linewidth = 0.5)
    plt.title(title + ' - OLS reg')
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position(
        [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9]
        )
    ax.legend(
        loc='upper center', bbox_to_anchor=(0.5, -0.11), 
        fancybox=True, shadow=True, ncol=5, prop={'size': 10}
        )

    #plt.xlim(741.5,742) # ZOOM-IN ON SECTION (MSIS 02-02-01)
    #plt.ylim(413.7265,413.7275)

    #plt.xlim(730,750)# ZOOM-IN ON SECTION (JB 02-02-01)
    #plt.ylim(408.65,408.85)

    plt.show()


plot_ols_reg(
    x_lab = 'Elapsed Time / hours', 
    y_lab = 'Semi-Major Axis / km', 
    start_aligned = True
    )
#xxxxxxxxxxxxx MAX LINE LENGTH = 79 CHARACTERS xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
