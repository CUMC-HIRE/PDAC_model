# Probability functions

'''
Authors: Myles Ingran
Description: Functions to compute probabilities
Last update: 07.02.19
'''

import math
from scipy.stats import chisquare
import numpy as np
import pandas as pd
import data_manipulation as dm

#import lin_risk as lr
from pathlib import Path



# determines if x is between the bounds
def between(x, low_bound, up_bound):
    if x >= low_bound and x <= up_bound:
        return True
    else:
        return False

def to_df(arg1, arg2):
    df = pd.DataFrame(list(arg2), index=arg1)
    return df
        
        
# Returns flat probability independent of time
def prob_flat(p):
    flat = p
    return flat


# Turns rate into a probability
def rate_to_prob(rate, t):
    prob = 1-math.exp(-abs(rate)*t)
#    print(prob)
    return prob

def prob_to_rate(prob, t):
    rate = -(np.log(1-prob))/t
#    print(rate)
    return rate
    

def prob_conv(prob, old_t, new_t):
    temp = prob_to_rate(prob, old_t)
    new_prob = rate_to_prob(temp, new_t)
    return new_prob

# Creates a linear probability
def prob_line(m, b, x):
    y_list = []
    for i in x:
        y = m*x + b
        y_list.append(y)
    df = to_df(x, y)
    return df


# Creates a flat probability between the bounds
def prob_boxcar(p_min, p_max, low_bound, up_bound, time):
    p = list()
    for t in time:
        if between(t, low_bound, up_bound):
            p.append(p_max)
        else:
            p.append(p_min)
#    plt.plot(time, p)
    return p


# Gaussian probability distrubition
def prob_normal(x, mean, SD):
    prob_list = []
    for i in x:
        prob = (1/(np.sqrt(2*np.pi*SD**2)))*np.exp(-(i-mean)**2/(2*SD**2))
        prob_list.append(prob)
    df = to_df(x, prob_list)
    return df



# Gives probability based on pw function
def pw_prob(time, slope_array, nodes):
    for n in range(len(nodes)-1):
        if between(time, nodes[n], nodes[n+1]-1):
            prob = rate_to_prob(slope_array[n], 1)
#                                ((nodes[n+1]-nodes[n])*ps.CYCLE_LENGTH))
            print(prob)
            return prob
        elif time >= nodes[len(nodes)-1]:
            prob = rate_to_prob(slope_array[len(slope_array)-1], 1)
#                                ((nodes[len(nodes)-1])*ps.CYCLE_LENGTH))
            return prob

def pw_choose_prob(time, probs, nodes):
    for n in range(len(nodes)-1):
        if between(time, nodes[n], nodes[n+1]):
            prob = probs[n]
            return prob
        elif time >= nodes[len(nodes)-1]:
            prob = probs[len(probs)-1]
            return prob
        elif time <= nodes[0]:
            prob = probs[0]
            return prob


def weighted_avg(rate, weight):
    weighted_avg = np.average(rate, weights=weight, axis=0)
    return weighted_avg


     

#xl = pd.ExcelFile(ps.life_table)
#life_table = xl.parse('Table 2')
#life_table = life_table.iloc[2:, 1].reset_index(drop=True).dropna()

def prob_to_prob(prob, from_cycle_length, to_cycle_length, time): 
    # Converts prob per one cycle length to prob per another cycle length
    # Inputs: 
    # prob is list of probabilities at each time point of original cycle length
    # from_cycle_length is original cycle length
    # to_cycle_length is desired cycle length
    # tmin and tmax: minimum and maximum times
    # Output: list of probabilities at each time point at desired cycle length
    
    small_step = 0
    big_step = 0
    t = min(time)
    i = 0
    big_dt = max(from_cycle_length, to_cycle_length)
    small_dt = min(from_cycle_length, to_cycle_length)
    num_cycle_ratio = from_cycle_length / to_cycle_length
    to_prob = []
    
    while t < max(time):
        from_rate = prob_to_rate(prob[i], t)
        to_rate = from_rate / num_cycle_ratio
        to_prob.append(rate_to_prob(to_rate, t))
        if small_step < num_cycle_ratio:
            small_step += 1
        else:
            small_step = 1
            big_step += 1
            i += 1
        t = (big_step*big_dt) + (small_step*small_dt)
#        print(t)
    return to_prob


def discount(cost, t, input_type = 'cost'):
    
    # Discounts values over time
    # Inputs: values to discount, length of simulation, and discount rate
    # Output: new array with discounted values
    if input_type == 'cost':
        new_costs = cost/((1+ps.dRate)**t)
    else:
        new_costs = cost/((1+0.015)**t)
        
    return new_costs

def normalize_target(row, row_index):
    '''
    Normalizes value at positiion row_index using the sum of all the other 
    probabilities in that row. Assumes all positions are correct expect for the
    row_index position
    
    Input: np array and location of value in question (row_index)
    
    Output: Normalized array
    '''
    i = 0
    other_prob_sum = 0
    for i in range(0, len(row)):
        if i != row_index:
            other_prob_sum += row[i]
        
    row[row_index] = 1 - other_prob_sum
    return row
    
    
def normalize_static(row, static):
    '''
    keeps one value (static) static while normalizes all other values in a 
    given array
    '''
    if round(sum(row), 8) == 0:
        return row
    else:
        # sum without static
        new_sum = 1 - row[static]
        if new_sum == 0:
            return row
    #    print(row)
        # gives array excluding static
        new_row = [x for i, x in enumerate(row) if i!= static]
        exclude_sum = round(sum(new_row), 8)
        # the new_sum divided by the sum of all elements except static
        norm_ratio = new_sum/exclude_sum
        # normalizes all elements that are not statuc
        for i in range(len(row)):
            if i != static:
                row[i] = norm_ratio * row[i]
    return row

# normalizes to 1--for use with death states
def normalize(row):
    if max(row) == 0:
        return row
    else:
        arr = (row - min(row))/(max(row) - min(row))
        arr2 = arr/sum(arr)
    return arr2

# noramlizes to a specific number
def normalize_choose(array, number):
    if number == 0:
        print("number is equal to zero")
        return array
    if sum(array) > number:
        array = np.divide(array, sum(array)/number)
#        print(array)
        return array
    elif sum(array) < number:
        array = np.multiply(array, number/sum(array))
#        print(array)
        return array
    else:
        return array


def normalize_matrix(matrix):
    matrix_width = matrix.shape[-2]
    matrix_depth = matrix.shape[-3]
    
    for i in range(matrix_depth):
        for j in range(matrix_width):
            matrix[i, j] = normalize(matrix[i, j])
    return matrix
    

def chi_square(observed, expected):
    chi = float(((observed-expected)**2)/abs(expected))
    if type(chi) == float:
        chi_square = chi
    else:    
        chi_square = sum(chi)
    print(chi_square)
    return chi_square


def mean_probs(probs_dict, nodes):
    prob_segs = list()
    for key in probs_dict.keys():
        prob_seg = list()
        prob_seg_1 = np.mean(probs_dict[key][:nodes[0]])
        prob_seg_2 = np.mean(probs_dict[key][nodes[0]:nodes[1]])
        prob_seg_3 = np.mean(probs_dict[key][nodes[1]:])
        prob_seg.append(key)
        prob_seg.append(prob_seg_1)
        prob_seg.append(prob_seg_2)
        prob_seg.append(prob_seg_3)
        prob_segs.append(prob_seg)
    return prob_segs
    

def risk_increase(prob, risk):
    inc_prob = prob/(1-risk)
    return inc_prob


def minus(rate_1, rate_2):
# =============================================================================
#     Subtract rate 1 from rate 2
# =============================================================================
    new_rate = rate_2 - rate_1
    return new_rate


def annual_rate_to_prob(rates, age):
    annual_prob = [rate_to_prob(rates[i], 1/minus(age[j], age[j+1]))
                    for i, j in zip(range(len(rates)), 
                                    range(len(age)-1))]
    return annual_prob


def cumul_prob_to_annual(risks, gene, multiplier):
# =============================================================================
#     Converts cumulative probability to annual probability
# =============================================================================
    age, cumul_prob = dm.excel_to_lists(risks, gene)
#    print(cumul_prob)
#    print(age)
    cumul_rate = [prob_to_rate(prob/100, 1) for prob in cumul_prob]
    annual_rate = [minus(cumul_rate[i], cumul_rate[i+1]) 
                    for i in range(len(cumul_rate)-1)]
    new_annual_rate = [rate*multiplier for rate in annual_rate]
    annual_prob = [rate_to_prob(new_annual_rate[i], 1/minus(age[j], age[j+1]))
                    for i, j in zip(range(len(new_annual_rate)), 
                                    range(len(age)-1))]
    #print(annual_prob)
    #print(age)
    return age, annual_prob, new_annual_rate



# merge this function and function above later
def cumul_prob_to_annual_raw(cumul_prob, age):
    if any(prob > 1 for prob in cumul_prob):
        cumul_prob = np.array(cumul_prob)/100
    cumul_rate = [prob_to_rate(prob, 1) for prob in cumul_prob]
#    print(cumul_rate)
    annual_rate = [minus(cumul_rate[i], cumul_rate[i+1]) 
                    for i in range(len(cumul_rate)-1)]
#    print(annual_rate)
    annual_prob = [rate_to_prob(annual_rate[i], 1/minus(age[j], age[j+1]))
                    for i, j in zip(range(len(annual_rate)), 
                                    range(len(age)-1))]
    return age, annual_prob, annual_rate


def cumul_prob_KM(filename):
    KM_df, nodes = lr.extract_KM_rates(filename)
    values = lr.pw_node_value(KM_df, nodes)
    age, annual_prob, annual_rate = cumul_prob_to_annual_raw(values, nodes)
    return age, annual_prob, annual_rate

def get_cancer_death_probs(age, gender):
    AC_death_rate = prob_to_rate(gender.lynch_ac_mortality[age], 1)
    if ps.CRC_death_rate.loc[2, 'death_rate'] > AC_death_rate:
        stage_2 = ps.CRC_death_rate.loc[2, 'death_rate'] - AC_death_rate
    else:
        stage_2 = 0
    if ps.CRC_death_rate.loc[3, 'death_rate'] > AC_death_rate:
        stage_3 = ps.CRC_death_rate.loc[3, 'death_rate'] - AC_death_rate
    else:
        stage_3 = 0
    if ps.CRC_death_rate.loc[4, 'death_rate'] > AC_death_rate:
        stage_4 = ps.CRC_death_rate.loc[4, 'death_rate'] - AC_death_rate
    else:
        stage_4 = 0
    return stage_2, stage_3, stage_4
    
def check_valid_file(filename):
    config = Path(filename)
    if config.is_file():
        #print('file is valid')
        return True
    else:
        #print('file is not valid')
        return False


def gamma_dist(mean, sd):
#   Gives scale (theta) and shape (kappa) for a gamma distribution
    if mean > 0:
        theta = sd**2 / mean
        kappa = mean / theta
    elif mean == 1:
        return mean
    else:
        return mean
    output = np.random.gamma(kappa, theta)
    return output


def beta_dist(mean, sd):
#   Gives a beta distribution based on mean and sd
    if mean > 0 and mean < 1:
        alpha = ((1-mean)/sd**2-(1/mean))*mean**2
        beta = alpha * (1/mean - 1)
        output = np.random.beta(alpha, beta)
        return output
    else:
        return mean

 

def unif_dist(lower, upper):
    '''
    Gives a uniform distribution based on upper and lower bound
    
    Input: Upper and lower bound of the distribution
    
    Output: number from unif distribution
    '''
    output = np.random.uniform(lower, upper)
    return output


def exp_half_life_rate(t_half_life):
    '''
    Given the time of the half life and using .5 as the half life value, finds
    the rate of decay
    
    Inputs: half life time
    
    Outputs: exponential rate of decay
    '''
    
    decay_rate = np.log(.5)/(-t_half_life)
    
    return decay_rate

