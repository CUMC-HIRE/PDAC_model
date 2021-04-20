# simulator functions
'''
Functions used to create and run a Markov model including generating a 
transition matrix and making a distribution matrix form that transition matrix
'''

import pandas as pd
import numpy as np
import data_manipulation as dm
import panc_presets as ps
import probability_functions as pf
from scipy.stats import chisquare
from copy import copy
from copy import deepcopy
import markov_plot as mplt
from sklearn.preprocessing import MinMaxScaler 
import linfit as lf

# Turns dictionaries of states and connections into a dataframe
def dict_to_connect_matrix(states, connects):
    c_matrix = np.zeros((len(states), len(states)))
    for key in connects.keys():
        for connect in connects[key]:
            if connect in states.keys():
                c_matrix[key, connect] = 1
    return c_matrix


def switch_connect(c_matrix, row, column):
# =============================================================================
#     Switches a connection located at the given row and 
#     column on and all others in that row of connectivity matrix off 
# =============================================================================
    row_length = c_matrix.shape[-1]
    
    for connect in range(row_length):
        if connect == column:
            c_matrix[row, int(connect)] = 1
        else:
            c_matrix[row, int(connect)] = 0
        
    return c_matrix


def nonegative(array):
    new_array = []
    for element in array:
        if element <= 0:
            element = 0
        else:
            element = element
        new_array.append(element)
    return new_array


def preset_prob_matrix(states, connects, probs):
    '''
    Inserts the probabilities into the connectivity matrix
    
    Inputs: dict of states in a Markov model, 
    connectivity dict of an arm, and probability dict of an arm
    
    Outputs: A matrix with probabilites between states
    '''
    
    # creates connectiivty matrix
    prob_matrix = dict_to_connect_matrix(states, connects)
    # flips states_dictionary to interact with probs dictionary
    flip = dm.flip(states)
    
    # adds the probability in probs dictionary to its corresponding place in the
    # probability matrix
    for key in probs.keys():
        if prob_matrix[flip[key[0]], flip[key[1]]] == 1:
            prob_matrix[flip[key[0]], flip[key[1]]] = probs[key]
    
    return prob_matrix
    

def all_cause_prob(states, prob_matrix, morb_table, age):
    '''
    Adds the all_cause probabilites to the rows of the prob matrix that connect
    to all_cause death
    
    Inputs: dict of states in a Markov model, probability dict of an arm, age
    
    '''
    # increases readability of code 
    flip = dm.flip(states)
    
    # adds all cause probability based on age
    for key in states.keys():
        ac_prob = prob_matrix[key, flip["all_cause"]]
        if key == flip["all_cause"]:
            continue
        elif ac_prob != 0:
            prob_matrix[key, flip["all_cause"]] = pf.prob_conv(morb_table[age],1, ps.CYCLE_LENGTH)
            
    return prob_matrix


def normalize_select(prob_matrix):
    '''
    Normalizes select items from a matrix based on the criteria that they equal
    1. These elements are assumed to be unset from the connectivity matrix. 
    
    Inputs: Matrix of probabilities
    
    Output: Normalized probability matrix where every row adds up to 1
    '''
    # iterates over each row in the probability matrix
    count = 0
    for row in prob_matrix:    
        for j in range(len(row)):
            scaler = MinMaxScaler(feature_range=(0,1))
            scaler.fit(row.reshape(-1, 1))
            row = scaler.transform(row.reshape(-1, 1))
            if row[j] == 1:
    
                # replaces the values that still equal 1 after probability
                # insertation with a normalized number
                prob_matrix[count] = pf.normalize_target(row.flatten(), j)
        
        count += 1
            
    return prob_matrix
                

def generate_prob_matrix(states, connects, prob, morb_table, age):
    '''
    Creates one instance of a probability matrix for a certain age
    
    Inputs: states dictionary, connectivity dictionary, probability dictionary,
    age
    
    Outputs: Normalized probability matrix
    '''
    
    # adds the preset probabilities to the connectivity matrix
    prob_matrix = preset_prob_matrix(states, connects, prob)
    # adds all cause probabilities to the probability matrix
    prob_matrix = all_cause_prob(states, prob_matrix, morb_table, age)
    
    # normalizes the matrix so every probability row adds up to 1
    prob_matrix = normalize_select(prob_matrix)

    
    return prob_matrix


def matrix_check(prob_matrix):
    '''
    Checks if matrix adheres to the constraines of probabilites. 
    
    Input: Probability matrix
    
    Output: Probability matrix checked for errors
    '''
    # chekcs if probability x is < 0 or > 1
    for item in np.nditer(prob_matrix):
        if item <= 1 and item >=0:
            continue
        elif item > 1:
            raise ValueError("value of item is greater than 1")
        elif item < 0:
            print(prob_matrix)
            raise ValueError("value of item is negative")
    for row in prob_matrix:
        if int(round(sum(row), 8)) == 1 or int(round(sum(row), 8)) == 0.0:
            continue 
        else:
            print(round(sum(row), 8))
            raise ValueError("sum of row exceeds 1")
    return prob_matrix
 

# initializes the start state
def get_start_state(connects):
    '''
    Uses the connectivity dict of an arm to determine the starting state.
    Assumes that the first state of the keys of the dict is the start state of
    the model
    
    Input: Connectivity dictionary
    
    Outputs: vector corresponding to the starting state of the model
    '''
    start_state = np.zeros((1, len(connects.keys())))
    for key in connects.keys():
        if connects[key] == []:
            continue
        else:
            start_state[0][key] = 1
            break
        
    return start_state
        

def create_t_matrix(arm, run_time, morb_table, *args):
    '''
    Creates a transition martix for a Markov model for ages in the model
    
    Inputs: variable in the arm class and run_time of the model, and if available
    special connectivity dict (arg[0]) for a range of cycles (arg[1])
    
    Output: a 3-D matrix of a transition matrix for every age in the model
    '''
    
    for t in range(run_time):
        age = int(round(arm.START_AGE + t*ps.CYCLE_LENGTH, 1))
        
        if args:
            if t < args[1]:
                prob_matrix = generate_prob_matrix(arm.states, args[0], 
                                                   arm.params_dict, morb_table,
                                                   age)
#                prob_matrix = all_cause_prob(arm.states, 
#                                             dict_to_connect_matrix(arm.states, 
#                                                                    args[0]),
#                                                                    morb_table,
#                                                                    age)
#                prob_matrix = normalize_select(prob_matrix)
            else:
                prob_matrix = generate_prob_matrix(arm.states, 
                                           arm.CONNECTIVITY, 
                                           arm.params_dict, morb_table, age)
        else:
            prob_matrix = generate_prob_matrix(arm.states, 
                                           arm.CONNECTIVITY, 
                                           arm.params_dict, morb_table, age)
        prob_matrix = matrix_check(prob_matrix)
        
        if t == 0:
            t_matrix = prob_matrix
        else:
            t_matrix = np.vstack((t_matrix, prob_matrix))
    t_matrix = np.reshape(t_matrix, (run_time, len(arm.states), len(arm.states)))
    return t_matrix


def build_D_matrix(time, start_state, t_matrix):
    '''
    Recursively constructs a distribution matrix from a start state and
    transition matrix
    
    Input: Start time of the model, start state vector of the model, transition
    matrix
    
    Output: Distribution matrix of the model
    
    '''
    matrix_length = t_matrix.shape[-1]
    temp = np.transpose(t_matrix[time]) * start_state[-1]
    Distribution = [sum(temp[i, :]) for i in range(matrix_length)]
    D_matrix = np.vstack((start_state, Distribution))
    return D_matrix if time == ps.RUN_TIME-1 else build_D_matrix(time+1, D_matrix, t_matrix)


def run_markov(arm, run_time, morb_table, *args):
    '''
    Creates the distribution matrix for arm
    
    Input: variable of the arm class, run time of the model, relevant morbidity
    table, and if available special connectivity dict for a range of cycles
    
    Output: transition matrix and distribution matrix of the model
    '''
    # initializes distribution matrix
    start_state = get_start_state(arm.CONNECTIVITY)
    # creates the transition matrix
    t_matrix = create_t_matrix(arm, run_time, morb_table, *args)

    # creates the distribution matrix based on the transition matrix
    D_matrix = build_D_matrix(0, start_state, t_matrix)
    
    return t_matrix, D_matrix

# make all functions below into a different script

def randomize_params(params, params_list, bound):
    '''
    Randomizes select parameters for calibrations
    
    Input: parameter dictionary, list of select parameters and the bounds 
    (percentage in decimals) to randomly select from
    
    Output: new parameter dictionary with randomized select parameters
    '''
    P = params.copy()
    
    for prob in params_list:
        P[prob] = np.random.uniform(P[prob]*(1-bound), 
                                         P[prob]*(1+bound))
    
    return P

def get_survival(D_matrix, states):
    '''
    Gets the overall survival or progrssion free survival
    of the model from distribution matrix by summing
    over the death states
    
    Input: distribution matrix (D_matrix) and list of relevant state 
    
    Output: the survival per cycle of the relevant states of the model
    '''
    ssum = [0] * len(D_matrix)
    for state in states:
        ssum = ssum + D_matrix[:, state]
    return ssum


def gof(model, nodes, slopes):
    '''
    Gives goodness of fit measurement using chi-square test.
    
    Input: Model, nodes of expected data, and slope of expected data
    
    Output: goodness of fit measurement
    '''
    exp_len = len(lf.pw_values(len(model), nodes, slopes))
    p_value = chisquare(model[:exp_len], f_exp=lf.pw_values(len(model), nodes, slopes))[0]

    return p_value


#def gof(model, nodes, slopes):
#    '''
#    Gives the goodness of fit metric for the model results based on chi-square
#    test
#    
#    Input: observed model values, nodes and slopes of the expected values
#    
#    Output: GOF score
#    '''
#    
#    GOF = chisquare(model, lf.pw_values(len(model)+1, nodes, slopes))
#    
#    return GOF



def calibrate_markov(arm, run_time, morb_table, target_value, ans_num, *args):
    '''
    Loops over the markov model making process, randomizing select parameters
    until the p-value of the output of the markov goes below the target value and
    saves the number of solutions specified
    
    
    Input: variable of the arm class, run time of the model, relevant morbidity
    table, the target value that the p-value needs to be under/equal to, and 
    if available special connectivity dict for a range of cycles
    
    Output: possible t_matrices and D_matrices that fit the GOF parameters
    '''    
    # initialize answer list
    answers = []
    # gets pw slopes and nodes from KM curves (make into one function)
    if arm.PFS_target:
        slopes_PFS, nodes_PFS = lf.extract_KM_rates(arm.PFS_target, 
                                                arm.PFS_nodes[0], arm.PFS_nodes[1])
#    KM_time_PFS, KM_PFS = dm.csv_to_lists(arm.PFS_target)
    if arm.OS_target:
        slopes_OS, nodes_OS = lf.extract_KM_rates(arm.OS_target,
                                              arm.OS_nodes[0], arm.OS_nodes[1])
#    KM_times_OS, KM_OS = dm.csv_to_lists(arm.OS_target)
    # loops over markov creation until an adequete p value is found
    while len(answers) < ans_num:
        # initialize p values
        print(len(answers))
        if arm.OS_target:
            p_value_OS = target_value +1
        if arm.PFS_target:
            p_value_PFS = target_value +1
        # runs markov until a t_matrix is created that fits the expected values
        while (p_value_OS >= target_value) or (p_value_PFS >= target_value):
            A = deepcopy(arm)
            A.params_dict = randomize_params(A.params_dict, A.calib_params, .98)
            t_matrix, D_matrix = run_markov(A, run_time, morb_table, *args)
#            print(D_matrix)
            if arm.OS_target:
                OS_sum = get_survival(D_matrix, ps.death_states)
                p_value_OS = gof(1-OS_sum, nodes_OS, slopes_OS)
                
            if arm.PFS_target:
                PFS_sum = get_survival(D_matrix, ps.disease_states)
                p_value_PFS = gof(1-PFS_sum, nodes_PFS, slopes_PFS)
                
            print("p_value OS: " + str(p_value_OS))
            print("p_value PFS: " + str(p_value_PFS))
            # save params_dict
        answers.append([t_matrix, D_matrix, A.params_dict])
        
    return answers
        
    

# make this better
def calibrate_markov_single(arm, run_time, morb_table, target_value, ans_num, *args):
    '''
    Loops over the markov model making process, randomizing select parameters
    until the p-value of the output of the markov goes below the target value and
    saves the number of solutions specified. Uses only one graph instead of two
    
    
    Input: variable of the arm class, run time of the model, relevant morbidity
    table, the target value that the p-value needs to be under/equal to, and 
    if available special connectivity dict for a range of cycles
    
    Output: possible t_matrices and D_matrices that fit the GOF parameters
    '''    
    # initialize answer list
    answers = []
    # gets pw slopes and nodes from KM curves (make into one function)

    if arm.OS_target:
        slopes_OS, nodes_OS = lf.extract_KM_rates(arm.OS_target,
                                              arm.OS_nodes[0], arm.OS_nodes[1])
#    KM_times_OS, KM_OS = dm.csv_to_lists(arm.OS_target)
    # loops over markov creation until an adequete p value is found
    while len(answers) < ans_num:
        # initialize p values
#        print(len(answers))
        if arm.OS_target:
            p_value_OS = target_value +1
            
        p_value_med = target_value+1
        # runs markov until a t_matrix is created that fits the expected values
        while (p_value_OS >= target_value) or (p_value_med >= target_value):
            A = deepcopy(arm)
            A.params_dict = randomize_params(A.params_dict, A.calib_params, .98)
            t_matrix, D_matrix = run_markov(A, run_time, morb_table, *args)
#            print(D_matrix)
            if arm.OS_target:
                OS_sum = get_survival(D_matrix, ps.death_states)
                p_value_OS = gof(1-OS_sum, nodes_OS, slopes_OS)
                print(p_value_OS)
            p_value_med = 1-chisquare(mplt.median_survival(D_matrix, ps.death_states), f_exp=38.7)[1]
            # save params_dict
        answers.append([t_matrix, D_matrix, A.params_dict])
    return answers
        

def edit_t_matrix_morb(arm, t_matrix, morb_table):
    '''
    Edits the t_matrix in a t_matrix dictionary with new morbidity tables
    
    Inputs: arm of strategy, t_matrix of strategy, new morbidity tables
    
    Output: Updated t_matrix
    
    '''
    for t in range(t_matrix.shape[0]):
        age = int(round(arm.START_AGE + t*ps.CYCLE_LENGTH, 1))
        # adds all cause probabilities to the probability matrix
        prob_matrix = all_cause_prob(arm.states, t_matrix[t], morb_table, age)
        
        # normalizes the matrix so every probability row adds up to 1
        for row in prob_matrix:
            row = pf.normalize_static(row, 11)
#        print(prob_matrix)
        matrix_check(prob_matrix)
        if t == 0:
            t_matrix_new = prob_matrix
        else:
            t_matrix_new = np.vstack((t_matrix_new, prob_matrix))
    t_matrix_new = np.reshape(t_matrix_new, (t_matrix.shape[0], len(arm.states), len(arm.states)))
    return t_matrix_new


def update_matrix_dicts(strats, t_matrix_dict, morb):
    '''
    Updates the t_matrix dictionary based on the specified morbidity status
    
    Input: strats dfs, t_matrix dictionary, morbidity status
    
    Output: new t_matrix dict
    '''
    morb_mult = {"no": .7,
                 "low": .95,
                 "mod": 1.1,
                 "sev": 1.5}
    D_matrix_dict_new = {}
    for strat in strats:
        # introduces arm
        arm = ps.arm(strat)
        names = dm.flip(ps.ALL_STATES)
        t_matrix_old = t_matrix_dict[arm.chemo_name.lower()]
        # new morbidity values
        if morb in ps.morb_table_dict.keys():
            morb_table = ps.morb_table_dict[morb]
        # changes t_matrix with new all cause probs
        t_matrix_new = edit_t_matrix_morb(arm, t_matrix_old, morb_table)
#        t_matrix_dict_new.update({arm.chemo_name.lower():t_matrix_new})
        # initializes distribution matrix
        for pair in arm.calib_params:
            state_1 = names[pair[0]]
            state_2 = names[pair[1]]
            prob_matrix = t_matrix_new[:, state_1, state_2]*morb_mult[morb]
        
            while max(prob_matrix) > 1 or max(prob_matrix) < 0:
                mult = round((morb_mult[morb]-.014), 10)
                prob_matrix = t_matrix_new[:, state_1, state_2]*morb_mult[morb]
                print(mult)
    
            
            t_matrix_new[:, state_1, state_2] = prob_matrix
    #        print(prob_matrix)
            for i in range(t_matrix_new.shape[0]):
#                print(t_matrix_new[i, state_1])
    ##            temp[i] = pf.normalize_static(temp[i, state_1], state_2)
                t_matrix_new[i, state_1] = pf.normalize(t_matrix_new[i, state_1])
                matrix_check(t_matrix_new[i])
#        print(temp[:, state_1, state_2])
        start_state = get_start_state(arm.CONNECTIVITY)
        # updates D_matrix_dict
        D_matrix_new = build_D_matrix(0, start_state, t_matrix_new)
        D_matrix_dict_new.update({arm.chemo_name.lower():D_matrix_new})
        
        
    return D_matrix_dict_new

    


    

    


        
            
        
    
    
    
    
    
    
