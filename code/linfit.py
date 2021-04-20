#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 12:07:33 2018

@author: GreenMonster
"""

from scipy import optimize
from scipy import stats
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import data_manipulation as dm
import panc_presets as ps
#from matplotlib.colors import LogNorm

#print("in main...")

# x and y coordinates of the plot
#plf = pd.read_csv(ps.pembro_PDL1p_PFS_filename)
#plf_OS = pd.read_csv(ps.pembro_Fuchs_OS_filename)
#time = plf.iloc[:, 0].values
#PFS =  plf.iloc[:, 1].values
#time_OS = plf_OS.iloc[:, 0].values
#OS = plf_OS.iloc[:, 1].values
#plt.plot(time, PFS, 'go')
#plt.plot(time_OS, OS, 'bo')

def slope_list(rise, run):
    # Takes a list of x and y distances from their 
    # respective averages and returns a slope
#    print(run)
    A = sum(np.multiply(rise, run))
    slope = A/sum(np.multiply(run, run))
    return slope
    
def intercept(x, y, m):
    # b = y - mx
    b = y - m * x
    return b

def mylinregress(x, y):
    # Takes x and y to produce a slope and intercept
    if (len(x) != len(y)):
        print("ERROR-- x_val and y_val have inconsistent lengths (in pw_linregress)")
        return
    
    x_avg = np.mean(x)
    y_avg = np.mean(y)
    A_x = [i-x_avg for i in x]
    A_y = [j-y_avg for j in y]
#    print(A_y)
    slope = slope_list(A_y, A_x)
    inter = intercept(x_avg, y_avg, slope)
    return slope, inter


#m, b = mylinregress(x_val, y_val)
        

#def pw_mylinregress(*args):
#    pw_slopes = []
#    pw_intercept = []
#    count = 1
#    for i in args:
#        if type(i) == int:
#            if count == 1:
#                slope, intercept = mylinregress(x_val[:i], y_val[:i])
#                pw_slopes.append(slope)
#                pw_intercept.append(intercept)
#                j = i
#            elif count == len(args):
#                slope, intercept = mylinregress(x_val[i:], y_val[:i])
#                pw_slopes.append(slope)
#                pw_intercept.append(intercept)
#            else:
#                slope, intercept = mylinregress(x_val[i:j], y_val[i:j])
#                pw_slopes.append(slope)
#                pw_intercept.append(intercept)
#                i = j
#        count += 1
#    pw = pd.DataFrame([pw_slopes, pw_intercept])
#    pw.index +=1
#    return pw_slopes

def pw_mylinregress(x_segments, y_segments):
    
    if (len(x_segments) != len(y_segments)):
        print("ERROR-- x_val and y_val have inconsistent lengths (in pw_linregress)")
        return
    

    slope = 0
    intercept = 0
    pw_slopes = list()
    pw_intercepts = list()
#    
##    for x_seg in x_segments, y_seg in y_segments:
    for x_seg, y_seg in zip(x_segments, y_segments):
        slope, intercept = mylinregress(x_seg, y_seg)
        pw_slopes.append(slope)
        pw_intercepts.append(intercept)
#    
    pw = pd.DataFrame([pw_slopes, pw_intercepts])
    return pw

def ensure_KM_bounds_consistency(x,y):
    
    if (x[0] < 0):
        x[0] = 0
        
    # check for y<0
        
    return x,y


def linpredict(x, m, b):
    p = m*x + b
    return p


def generate_segments(x,y, nodes):
    #exception
    if ( not (len(nodes) < len(x))): 
        print("Node length exceeds length of x")
    if ( not (len(x) == len(y))): #exception
        print("length of x does not match length of y")
    ind = 0
#    c = 0
    x_segments = list()
    y_segments = list()
    
    for i in range(len(nodes)-1):        
        x_segment = list()
        y_segment = list()
        for ind in range(len(x)):        
            if(x[ind] >= nodes[i] and x[ind] < nodes[i+1]):
                x_segment.append(x[ind])
                y_segment.append(y[ind])
#                print(ind)
        x_segments.append(x_segment)
        y_segments.append(y_segment)
        
    return x_segments, y_segments

def add_nodes(x_seg, y_seg, nodes, df): 
    for node in nodes:
        for i in range(1, len(x_seg)):
            if node not in x_seg[i] and node > x_seg[i-1][len(x_seg)-1] and node < x_seg[i][0]:
                x_seg[i].append(node)
                y_seg[i].append(linpredict(node, df.loc[0, i], df.loc[1, i]))
#                print(node)
            sorted(x_seg[i])
    return x_seg, y_seg


#pw = pw_mylinregress(11, 25)
# Takes a dataframe of slope and intercepts and plots them
# Assumes first row is slopes and second row is intercepts
def pw_plot(x, df, color):
    for i in df.columns:
        x[i] = np.array(x[i])
        pw_plot = plt.plot(x[i], df.loc[0, i]*x[i] + df.loc[1, i], color)
    return pw_plot


def pw_values(x, nodes, df):
    pw_values = []
    for j in range(len(nodes)-1):
        for ind in range(x):
            if ind >= nodes[j] and ind < nodes[j+1]:
                pw_value = df.loc[0, j]*ind + df.loc[1, j]
                pw_values.append(pw_value) 
    pw_values = np.array(pw_values)
    pw_values[0], pw_values[1] = 1, 1
    return pw_values
        


def r_squared(y, p):
    y_avg = sum(y)/len(y)
    A_y = []
    A_p = []
    for i in range(len(y)):
        d_y = y[i] - y_avg
        A_y.append(d_y)
        d_p = y[i] - p[i]
        A_p.append(d_p)
    ss_tot = sum(np.multiply(A_y, A_y))
    ss_res = sum(np.multiply(A_p, A_p))
    r_s = 1 - ss_res/ss_tot
    return r_s



def pw_r_squared(x, y, seg):
    df = pw_mylinregress(x, y)
    r_squared_list = []
    for i in range(len(x)):
        y_p = r_squared(y[i], linpredict(x[i], df.loc[0, i], df.loc[1, i]))
        r_squared_list.append(y_p)
    return
        

def extract_KM_rates(filename, node_1, node_2):
    KM_t, KM_value = dm.csv_to_lists(filename) 
    KM_t = np.asarray(KM_t)
    KM_value = np.asarray(KM_value)/100
#    x, y = [int(x) for x in input('Enter two nodes: ').split()]
    nodes = [min(KM_t), node_1, node_2, max(KM_t)+1]
#    print(nodes)
    KM_x, KM_y = generate_segments(KM_t, KM_value, nodes)
    KM_linfit_rates = pw_mylinregress(KM_x, KM_y)
    return KM_linfit_rates, nodes


def plot_linfit_KM(filename):
    KM_t, KM_value = dm.csv_to_lists(filename) 
    KM_t = np.asarray(KM_t) 
    KM_value = np.asarray(KM_value)/100
    x, y = [int(x) for x in input('Enter two nodes: ').split()]
    nodes = [min(KM_t), x, y, max(KM_t)+1]
    KM_x, KM_y = generate_segments(KM_t, KM_value, nodes)
    KM_linfit_rates = pw_mylinregress(KM_x, KM_y)
    plt.axis([0, 75, 0, 1])
    plt.plot(KM_t, KM_value, 'go')
    pw_plot(KM_x, KM_linfit_rates, 'r')
    return


def plot_KM(filename):
    KM_t, KM_value = dm.csv_to_lists(filename)
    KM_t = np.asarray(KM_t)/12/ps.CYCLE_LENGTH
    KM_value = np.asarray(KM_value)/100
    if "Fuchs" in str(filename):
        label = "Fuchs Pembro"
    if "Pembro" in str(filename):
        if "Shitara" in str(filename):
            label = "Shitara Pembro PDL1-"
        else: 
            label = "Shitara Pembro PDL1+"
    if "Wilkes" in str(filename):
        label = "Wilkes Pac"
    plt.plot(KM_t, KM_value, label=label)
    return


def plot_KMs(title, *filename):
    for name in filename:
        plot_KM(name)
    x, y = dm.combine_KM(range(27), ps.PDL1_PFS_path, ps.non_PDL1_PFS_path, .6622)
    y = y/100
    plt.plot(x, y, label="Shitara Combined")
    plt.legend()
    plt.title(title)
    plt.xlabel("time (Months)")
    plt.ylabel("Survival (%)")
    return


def make_line(slope, intercept, x):
    values = [slope*i + intercept for i in x]
    return values


def cancer_death_BSC(filename):
    KM_t, KM_value = dm.csv_to_lists(filename)
    KM_t = np.asarray(KM_t)/12/ps.CYCLE_LENGTH
    KM_value = np.asarray(KM_value)/100
    death_rate, death_intercept = mylinregress(
            [KM_t[0], KM_t[-1]], [KM_value[0], KM_value[-1]]
            )
    plt.axis([0, 20, 0, 1])
    plt.plot(KM_t, KM_value, 'go')
    plt.plot(KM_t, make_line(death_rate, death_intercept, KM_t), 'r')
    return death_rate, death_intercept



# normal distribution center at x=0 and y=5
#x = np.random.randn(100000)
#y = np.random.randn(100000) + 5
#
#plt.hist2d(x, y, bins=20, norm=LogNorm())
#plt.colorbar()
#plt.show()


# sanity checks:
#time,PFS = ensure_KM_bounds_consistency(time,PFS)

#print("time=", time)
#print("PFS=", PFS)
    
#nodes = [0, 2, 7, max(time)+1]   # why +1? (generate_segments is buggy..)
#nodes_OS = [0, 2, 7, max(time_OS)+1]
#x_val, y_val = generate_segments(time, PFS, nodes)
#x_val_OS, y_val_OS = generate_segments(time_OS, OS, nodes)
##PFS = y_val
#pwlin_KM_fit_OS = pw_mylinregress(x_val_OS, y_val_OS)
#pwlin_KM_fit = pw_mylinregress(x_val, y_val)
#print(pwlin_KM_fit_OS)
#print(pwlin_KM_fit)
#x_val, y_val = add_nodes(x_val, y_val, nodes, pwlin_KM_fit)
#x_val_OS, y_val_OS = add_nodes(x_val_OS, y_val_OS, nodes_OS, pwlin_KM_fit_OS)
#
##pw = pw_values(pwlin_KM_fit)
##print("x_val = " , x_val)
##print("y_val = " , y_val)
##plt.axis([0,29, 0, 100])
#pw_plot(x_val, pwlin_KM_fit, 'r')
#pw_plot(x_val_OS, pwlin_KM_fit_OS, 'm')


# =============================================================================
# m1, b1, m2, b2, m3, b3 = pw_mylinregress(11, 25)
# r_s = r_squared(x_val[11:25], y_val[11:25], (b2 + m2*x_val[11:25]))
# plt.plot(x, y)
# plt.plot(x_val[:11], b1 + m1*x_val[:11], 'r')
# plt.plot(x_val[11:25], b2 + m2*x_val[11:25], 'r')
# plt.plot(x_val[25:len(x)], b3 + m3*x_val[25:len(x)], 'r')
# plt.axis([0,29, 0, 100])
# print(r_s)
# 
# =============================================================================


# ====================================================
# # finding inflection point of the graph
# tck = interpolate.splrep(x, y, k=2, s=0)
# xnew = np.linspace(0, 30)
# dev_2 = interpolate.splev(x, tck, der=2)
# turning_point_mask = dev_2 == np.amax(dev_2)
# dev_2_dict = {}
# for i  in range(len(dev_2)):
#     dev_2_dict[dev_2[i]] = x[i]
# 
# # the slopes of the parts of the piecewise graph
# slope1, intercept1, r_value, p_value, std_err = stats.linregress(x[0:17], y[0:17])
# slope2, intercept2, r_value, p_value, std_err = stats.linregress(x[17:73], y[17:73])
# 
# # creates a piecewise linear function for a given graph found at https://stackoverflow.com/questions/29382903/how-to-apply-piecewise-linear-fit-in-python
# def piecewise_linear(x, x0, y0, k1, k2):
#     return np.piecewise(x, [x < x0 ], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
# 
# # fits piecewise linear function to a curve
# p , e = optimize.curve_fit(piecewise_linear, x, y)
# xd = np.linspace(0, 10, 100)
# plt.plot(x, y, "o")
# plt.plot(xd, piecewise_linear(xd, *p))
# 
# # Interpolation
# yinterp = np.interp(x, x, y)
# 
# plt.plot(x, y, 'go')
# plt.plot(x, intercept + slope2*x, 'r')
# plt.plot(x, yinterp, 'b-')
# plt.show()
# 
# 
# =============================================================================



#def main():
#    
#    print("in main...")
#    
#    # x and y coordinates of the plot
#    plf = pd.read_csv(ps.pembro_PDL1p_PFS_filename)
#    time = plf.iloc[:, 0].values
#    PFS =  plf.iloc[:, 1].values
#    plt.plot(time, PFS, 'go')
#
#    
#    print("time=", time)
#    print("PFS=", PFS)
#    
##    pwlin_KM_fit = pw_linregress(x_val, y_val)
##    
##    print(pwlin_KM_fit)
#
#
#
#if __name__ == "__main___":
#    main()
    
    
