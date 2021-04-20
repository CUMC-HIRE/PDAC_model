# Markov Plotting Functions

import matplotlib.pyplot as plt
import pandas as pd
import probability_functions as pf
import numpy as np
import panc_presets as ps
import data_manipulation as dm

def median_survival(D_matrix, state_subset):
    '''
    Gives the median overall survival or progression free survival of a model's
    distribution matrix
    
    Input: Distribution matrix of a Markov model, states of interest
    
    Output: Median survival
    '''
    # itereates over the states of interest in the distribution matrix to get
    # sum of population in those states
    Sum = [0] * len(D_matrix)
    for state in state_subset:
        Sum = Sum + D_matrix[:, state]
    # turns into percentage of people not in those states
    Sum = 1-Sum
    # if there is a .5 in the sum array then finds it. If not calculates the 
    # median survival point
    if Sum[Sum==.5].size > 0:
        med_srvl = np.where(Sum==.5)
    else:
        x = len(Sum[Sum>.5])
        slope = Sum[x-1] - Sum[x]
        intercept = Sum[x] - slope * x
        med_srvl = (.5 - intercept)/slope
    # rounds the final answer
    med_srvl = round(med_srvl, 2)
    return med_srvl#, Sum


def all_median_survivals(D_matrix_dict, state_subset):
    '''
    Gives all the median overall survivals or progression free survivals of 
    the arms in the model
    
    Input: Dictionary of all arms and corresponding distribution matrices, 
    states of interest
    
    Output: Median survivals for each arm
    '''
    # initialize med survival dictionary
    med_srvl_dict = {}
    
    for key in D_matrix_dict.keys():
        med_srvl = median_survival(D_matrix_dict[key], state_subset)
        med_srvl_dict.update({key: med_srvl})
    
    return med_srvl_dict
    
    
def survival_sum(D_matrix, states):
    '''
    Calculates the percentage surviving in all the states listed
    
    Input: states to be accounted for in the summation, distribution matrix
    
    Output: Percentage surviving in those states
    '''
    srvl_sum = [0] * len(D_matrix)
    for state in states:
        srvl_sum = srvl_sum + D_matrix[:, state]
    return srvl_sum
    
    


def OS_plot(D_matrix, death_states, arm, overwrite_file = False):
    '''
    Creates a plot of the overall survival of a model given the distribution
    matrix
    
    Input: distribution matrix of a Markov model
    
    Output: Plot of the overall survival of the Markov model
    '''
    OS_sum = [0] * len(D_matrix)
    for state in death_states:
        OS_sum = OS_sum + D_matrix[:, state]
    plt.title("PDAC Overall Survival " + str(arm.chemo_name))
    plt.plot(range(len(D_matrix)), 1-OS_sum, label="model")
    
    
    # create dictionary in presets for papers and arms
    if arm.chemo_name == "FOLF":
        exp = "Dhir 2018"
    elif arm.chemo_name == "AG":
        exp = "Dhir 2018"
    elif arm.chemo_name == "PEXG":
        exp = "Reni 2018"
    elif arm.chemo_name == "Natural History":
        exp = "Shapiro 2016"
    else:
        exp = None

    # plot expected
    if exp != None:
        KM_time_OS, KM_OS = dm.csv_to_lists(arm.OS_target)
        plt.plot(KM_time_OS, np.asarray(KM_OS)/100, label=exp)
    plt.xlabel("Months")
    plt.ylabel("Survival (%)")
   
    plt.legend()
    
    plot_name = "PDAC_OS_" + str(arm.chemo_name)
    if pf.check_valid_file(ps.graphs/plot_name) == False or overwrite_file == True:
        plt.savefig(ps.dump/plot_name)
    plt.show()
    return OS_sum
    

def PFS_plot(D_matrix, disease_states, arm, overwrite_file = False):
    '''
    Creates a plot of the progression-free survival of a model given 
    the distribution matrix
    
    Input: distribution matrix of a Markov model
    
    Output: Plot of the progression-free survival of the Markov model
    '''
    PFS_sum = [0] * len(D_matrix)
    for state in disease_states:
        PFS_sum = PFS_sum + D_matrix[:, state]
        
    plt.title("PDAC Progression-Free Survival " + str(arm.chemo_name))
    plt.plot(range(len(D_matrix)), 1-PFS_sum, label="model")
    
    if arm.chemo_name == "FOLF":
        exp = "Michelakos 2019"
    elif arm.chemo_name == "AG":
        exp = "Ielpo 2017"
    elif arm.chemo_name == "PEXG":
        exp = "Reni 2018"
    else:
        exp = None

    if exp != None:
        KM_time_PFS, KM_PFS = dm.csv_to_lists(arm.PFS_target)
        plt.plot(KM_time_PFS, np.asarray(KM_PFS)/100, label=exp)
    plt.xlabel("Months")
    plt.ylabel("Survival (%)")
 
    plt.legend()
    
    plot_name = "PDAC_PFS_" + str(arm.chemo_name)
    if pf.check_valid_file(ps.graphs/plot_name) == False or overwrite_file == True:
        plt.savefig(ps.dump/plot_name)
    plt.show()
    return PFS_sum


def comp_plot_OS(D_matrix_dict, death_states, overwrite_file = False):

    plt.title("PDAC Overall Survival")
    for key in D_matrix_dict.keys():
        OS_sum = [0] * len(D_matrix_dict[key])
        for state in death_states:
            OS_sum = OS_sum + D_matrix_dict[key][:, state]
        plt.plot(range(len(D_matrix_dict[key])), 1-OS_sum, label=key)
    plt.xlabel("Months")
    plt.ylabel("Survival (%)")
    plt.legend(labels=["FOLFIRINOX", "G-nP",  "G-nP-C", "Natural History",])
    plot_name = "PDAC_OS_all_strats"
    if pf.check_valid_file(ps.graphs/plot_name) == False or overwrite_file == True:
        plt.savefig(ps.dump/plot_name)
    plt.show()
    return


def comp_plot_PFS(D_matrix_dict, disease_states, overwrite_file = False):

    plt.title("PDAC Progression-Free Survival")
    for key in D_matrix_dict.keys():
        PFS_sum = [0] * len(D_matrix_dict[key])
        for state in disease_states:
            PFS_sum = PFS_sum + D_matrix_dict[key][:, state]
        plt.plot(range(len(D_matrix_dict[key])), 1-PFS_sum, label=key)
    plt.xlabel("Months")
    plt.ylabel("Survival (%)")
    plt.legend(labels=["FOLFIRINOX", "G-nP", "G-nP-C", "Natural History"])
    plot_name = "PDAC_PFS_all_strats"
    if pf.check_valid_file(ps.graphs/plot_name) == False or overwrite_file == True:
        plt.savefig(ps.dump/plot_name)
    plt.show()
    return
