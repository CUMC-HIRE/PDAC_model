# markov_icer

import numpy as np
import pandas as pd
import data_manipulation as dm
import panc_presets as ps
import markov_plot as mplt


'''
Calculates the incremental cost-effective ratio of a model including strategies
'''

def create_weight_vector(weights_dict, states):
    '''
    Creates vector with weights that correspond to a state
    
    Input: weights dictionary with weights as the values and states as the keys
    and states of the model
    
    Output: weight vector
    '''
    
    # finds the longest set of items in value in case value in dictionary as
    # more than one item
    val_lens = [len(val) for val in weights_dict.values()]
    
    vector = np.zeros((len(states), max(val_lens)))
    
    # flips states dictionary to help with readability
    flip = dm.flip(states)
    
    # populates the vector with the values of the weights dictionary
    for key in weights_dict.keys():
    
        if vector.shape[1] > 1:
            for i in range(vector.shape[1]):
                
                vector[flip[key], i] = weights_dict[key][i]
        else:
            vector[flip[key]] = weights_dict[key]
            
    return vector


def create_strat_ce_df(D_matrix, weight_vector, strat):
    '''
    Multipies the weight vector by each row of the distribution matrix to get 
    the total cost and utility of the model. Also adds the life states on each
    row of the distribution matrix to get the total life-years of the model
    
    Input: distribution matrix, weight_vector (1st column cost, 
    2nd column utilities), and name of strategy (string)
    
    Output: cost-effective dataframe with life-years, costs, and QALYs
    '''
    # initializes the life_year, costs, and utilities counts
    life_years = 0
    cost = 0
    if strat in ps.neoadj_list:
        utilities = 0.8 # progression-free when starting treatment for neoadj
    else:
        utilities = 0
        
    lf_list = []
    cost_list = []
    util_list = []
#    print(strat)
    for row in D_matrix:
        # accumulates the life years, costs, and utitilies per row of the
        # distribution matrix
        life_years += sum(row[min(ps.life_states):(max(ps.life_states)+1)])
        cost += sum(row * weight_vector[:, 0])
        utilities += sum(row * weight_vector[:, 1])
        lf_list.append(life_years)
        cost_list.append(cost)
        util_list.append(utilities)
    
    df = pd.DataFrame(data={"Life_Years": lf_list, 
                                "Costs": cost_list, "Util": util_list})
    df.to_csv(strat+ "_icer_values.csv")

        
    # consolidates the three values into one array for easy icer calculations
    strat_ce_df = pd.DataFrame(data={"life_years":life_years*ps.CYCLE_LENGTH, 
                                     "cost":cost, "QALYs":utilities*ps.CYCLE_LENGTH}, 
                             index=[strat])
#    print(strat_ce_df)
    return strat_ce_df

def create_all_ce_matrices(D_matrix_dict, states):
    '''
    Aggrates all strategies in the model into one dataframe to make the icer
    calulations easier
    
    Inputs: Dictionary of strategies as keys and their respective distrubution
    matrices, arm of the strategy, and states of the model
    '''
    # initializes the cost-effective dataframe for all genes
    all_ce_df = pd.DataFrame(columns=["life_years", "cost", "QALYs"])
    
    # stacks strat cost-effective dataframes on top of each other
    for strat in ps.strats():
        arm = ps.arm(strat)
        weight_vector = create_weight_vector(arm.weights_dict, states)
        strat_ce_df = create_strat_ce_df(D_matrix_dict[strat], weight_vector, strat)
        all_ce_df = pd.concat([all_ce_df, strat_ce_df])
    
    return all_ce_df
        

def calculate_icer(cost_1, cost_2, qaly_1, qaly_2):
    '''
    Calculates icers given the costs and qalys of two strategies
    
    Inputs: costs and qalys of two strategies
    
    Output: icer value
    '''
    delta_cost = cost_2 - cost_1
    delta_qaly = qaly_2 - qaly_1
    
    icer = round(delta_cost/delta_qaly, 5)
    
    return icer
    


def get_icer(all_ce_df):
    '''
    Calculates icers from the cost-effective dataframe containing all the
    strategies
    
    Inputs: cost-effective dataframe with all strategies
    
    Output: cost-effective dataframe with icers
    '''
    
    ce_icers = all_ce_df.sort_values(by=["cost"])
    
    strat_len = len(ce_icers)
    count = 0
    
    # Eliminates strongly dominated strats (higher cost, lower qalys)
    while count < strat_len-1:
        if (ce_icers["QALYs"][count+1] < ce_icers["QALYs"][count]):
            ce_icers = ce_icers.drop([ce_icers.index[count+1]])
            strat_len = len(ce_icers)
            count = 0
        else:
            count += 1
    
    # Initialize icers column
    ce_icers.loc[:, 'icers'] = 0
    
    # Eliminate weakly dominated strategies
    if len(ce_icers) > 1:
        strat_len = len(ce_icers)
        count = 1
        while count < strat_len:
            ce_icers.loc[ce_icers.index[count], "icers"] = (
            calculate_icer(ce_icers["cost"][count], ce_icers["cost"][count-1],
                           ce_icers["QALYs"][count], ce_icers["QALYs"][count-1]))
            # drops weakly dominated strat
            if(ce_icers.loc[ce_icers.index[count], "icers"] < 
               ce_icers.loc[ce_icers.index[count-1], "icers"]):
                ce_icers = ce_icers.drop(ce_icers.index.values[count-1])
                strat_len = len(ce_icers)
                count = count -1
            else:
                count += 1
                
    return ce_icers


def get_icer_w_dom(all_ce_df):
    '''
    Calculates icers from the cost-effective dataframe containing all the
    strategies
    
    Inputs: cost-effective dataframe with all strategies
    
    Output: cost-effective dataframe with icers
    '''
    
    ce_icers = all_ce_df.sort_values(by=["cost"])
    
    strat_len = len(ce_icers)
    count = 0
    
    # Initialize icers column
    ce_icers.loc[:, 'icers'] = 0
    
    # Eliminates strongly dominated strats (higher cost, lower qalys)
    while count < strat_len-1:
        if (ce_icers["QALYs"][count+1] < ce_icers["QALYs"][count]):
            ce_icers.loc[ce_icers.index[count+1], 'icers'] = "dominated"
        count += 1
    
#    print(ce_icers)
    
    # Eliminate weakly dominated strategies
    if len(ce_icers) > 1:
        strat_len = len(ce_icers)
        count = 1
        comp_count = 0
        while count < strat_len:
#            print(count)
            if ce_icers.loc[ce_icers.index[count], "icers"] != "dominated" and ce_icers.loc[ce_icers.index[comp_count], "icers"] != "dominated":
                ce_icers.loc[ce_icers.index[count], "icers"] = (
                        calculate_icer(ce_icers["cost"][count], ce_icers["cost"][comp_count],
                           ce_icers["QALYs"][count], ce_icers["QALYs"][comp_count]))
#                print(ce_icers.loc[ce_icers.index[count], "icers"])
            # drops weakly dominated strat
                if(ce_icers.loc[ce_icers.index[count], "icers"] < 
                    ce_icers.loc[ce_icers.index[comp_count], "icers"]):
                    ce_icers.loc[ce_icers.index[comp_count], 'icers'] = "dominated"
                    count = count - 1
                else:
                    count += 1
                    comp_count +=1
            else:
                count += 1
                
    return ce_icers



#def 


def run_icers(strats, states, D_matrix_dict):
    '''
    Calculates icers for the strategies
    
    Inputs: strats list, states dicitonary, dictionary of distribution matrices
    with the strats as the keys
    
    Output: Dataframe containing the icers 
    '''
    # initialize the dataframe for the costs and qalys of the strategies
    all_ce_df = pd.DataFrame(columns=["life_years", "cost", "QALYs"])
    for strat in strats:
        arm = ps.arm(strat)
#        print(arm.params_dict)
#        print(arm.cu_dict)
        # creates weight vector to multiply the rows of the distribution matrix by
        wv = create_weight_vector(arm.cu_dict, states)
        wv = wv * ps.discount
#        print(arm.chemo_name)
#        print(wv)
        # creates the strat ce series
#        print(arm.chemo_name.lower())
        strat_ce_df = create_strat_ce_df(D_matrix_dict[arm.chemo_name.lower()], 
                                                       wv, arm.chemo_name.lower())
        # adds the strat series to the overall dataframe
        all_ce_df = pd.concat([all_ce_df, strat_ce_df])
    # gets the icers from the complete ce dataframe
    ce_icers = get_icer(all_ce_df)
#    print(ce_icers)
    
    return ce_icers

def run_icers_w_dom(strats, states, D_matrix_dict):
    '''
    Calculates icers for the strategies
    
    Inputs: strats list, states dicitonary, dictionary of distribution matrices
    with the strats as the keys
    
    Output: Dataframe containing the icers 
    '''
    # initialize the dataframe for the costs and qalys of the strategies
    all_ce_df = pd.DataFrame(columns=["life_years", "cost", "QALYs"])
    for strat in strats:
        arm = ps.arm(strat)
#        print(arm.params_dict)
#        print(arm.cu_dict)
        # creates weight vector to multiply the rows of the distribution matrix by
        wv = create_weight_vector(arm.cu_dict, states)
        wv = wv * ps.discount
#        print(arm.chemo_name)
#        print(wv)
        # creates the strat ce series
#        print(arm.chemo_name.lower())
        strat_ce_df = create_strat_ce_df(D_matrix_dict[arm.chemo_name.lower()], 
                                                       wv, arm.chemo_name.lower())
        # adds the strat series to the overall dataframe
        all_ce_df = pd.concat([all_ce_df, strat_ce_df])
    # gets the icers from the complete ce dataframe
    ce_icers = get_icer_w_dom(all_ce_df)
#    print(ce_icers)
    
    return ce_icers

def run_icers_arms(arms, states, D_matrix_dict):
    '''
    Calculates icers for the strategies
    
    Inputs: arms list, states dicitonary, dictionary of distribution matrices
    with the strats as the keys
    
    Output: Dataframe containing the icers 
    '''
    # initialize the dataframe for the costs and qalys of the strategies
    all_ce_df = pd.DataFrame(columns=["life_years", "cost", "QALYs"])
    for arm in arms:
        # creates weight vector to multiply the rows of the distribution matrix by
        wv = create_weight_vector(arm.cu_dict, states)
        wv = wv * ps.discount
#        print(arm.chemo_name)
#        print(wv)
        # creates the strat ce series
#        print(arm.chemo_name.lower())
        strat_ce_df = create_strat_ce_df(D_matrix_dict[arm.chemo_name.lower()], 
                                                       wv, arm.chemo_name.lower())
        # adds the strat series to the overall dataframe
        all_ce_df = pd.concat([all_ce_df, strat_ce_df])
    # gets the icers from the complete ce dataframe
    ce_icers = get_icer_w_dom(all_ce_df)
#    print(ce_icers)
    
    return ce_icers


def get_x_year_srvl(x, D_matrix_dict=ps.D_matrix_dict, srvl="OS"):
    '''
    Gets the x year survival 
    
    Input: x number of years
    
    Output: dictionary of arm and percentage alive at x years
    '''
    if srvl == "OS":
        states = ps.death_states
    elif srvl == "PFS":
        states = ps.disease_states
    x_srvl_dict = {}
    for key in D_matrix_dict.keys():
        d_matrix = D_matrix_dict[key]
        months = int((x*12) + 1)
        x_d_matrix = d_matrix[months]
        temp = [x_d_matrix[state] for state in states]
        x_year_srvl = 1-sum(temp)
        x_srvl_dict.update({key:x_year_srvl})
    x_srvl_df = pd.DataFrame.from_dict(x_srvl_dict, orient="index", 
                                       columns=[str(x) + " Year " + srvl])
    return x_srvl_df


def get_x_year_death(x, D_matrix_dict=ps.D_matrix_dict, death_state = "cancer_death"):
    '''
    Gets the x year percentage of people dead in a certain death state
    
    Input: x number of years, name of specific death state
    
    Output: dictionary of arm and percentage alive at x years
    '''
    
    names_dict = dm.flip(ps.ALL_STATES)
    d_state = names_dict[death_state]
    x_srvl_dict = {}
    for key in D_matrix_dict.keys():
        d_matrix = D_matrix_dict[key]
        months = int((x*12) + 1)
        if months == 145:
            months = 144
        x_year_srvl = d_matrix[months][d_state]
        x_srvl_dict.update({key:x_year_srvl})
    x_srvl_df = pd.DataFrame.from_dict(x_srvl_dict, orient="index", 
                                       columns=[str(x) + " Year " + death_state])
    return x_srvl_df
    

def R0_status(D_matrix_dict=ps.D_matrix_dict, status = 3):
    '''
    Creates a df of the percentage of R0 status population after 
    surgery for each strat
    
    Input: R0 status 3 or R1 status 4
    
    Output: R1 df or R0 df per strat
    '''
    R0_dict = {}
    for key in D_matrix_dict.keys():
        d_matrix = D_matrix_dict[key]
        # 3 is R0 status state 
        R0 = d_matrix[ps.status_mnth][status]
        R0_dict.update({key:R0})
    R0_df = pd.DataFrame.from_dict(R0_dict, orient="index", columns=["R0 status"])
    return R0_df


def surg_status(D_matrix_dict=ps.D_matrix_dict):
    '''
    Creates a df of the percentage of population recieving surgeryfor each strat
    '''
    surg_dict = {}
    for key in D_matrix_dict.keys():
        d_matrix = D_matrix_dict[key]
        # 3 is R0 status state 
        surg = d_matrix[ps.surg_mnth][1]
        surg_dict.update({key:surg})
    surg_df = pd.DataFrame.from_dict(surg_dict, orient="index", columns=["% Recieving Surgery"])
    return surg_df
    
    
    


def get_surg_srvl(D_matrix, srvl='OS', time='post'):
    '''
    Gets the OS or PFS of pre or post surgery strategies
    
    Input: Distribution matrix, OS or PFS for srvl and pre or post for time
    
    Output: survival based on the preset conditions
    '''
    if time == 'post':
        surg_D_matrix = D_matrix[ps.post_surg_mnth:]
    elif time == 'pre':
        surg_D_matrix = D_matrix[:ps.post_surg_mnth]
        denom = 1
    if srvl == 'OS':
        surg_death = mplt.survival_sum(surg_D_matrix, ps.death_states)
        if time == 'post':
            temp = [surg_D_matrix[state] for state in ps.death_states]
            denom = 1-sum(temp)
    elif srvl == 'PFS':
        surg_death = mplt.survival_sum(surg_D_matrix, ps.disease_states)
        if time == 'post':
            temp = [surg_D_matrix[state] for state in ps.disease_states]
            denom = 1-sum(temp)
    surg_srvl = 1-surg_death
    srvl_perc = surg_srvl/denom
    return surg_srvl, srvl_perc
    

def get_all_surg_srvl(D_matrix_dict=ps.D_matrix_dict, srvls='OS', times='post'):
    '''
    Gets the surgical survival information and percentage of people 
    '''
    surg_dict = {}
    for key in D_matrix_dict.keys():
        surg_srvl, srvl_perc = get_surg_srvl(D_matrix_dict[key], srvl=srvls, time=times)
        surg_dict.update({key:(surg_srvl, srvl_perc)})
    return surg_dict
        
        
def all_x_year_srvl(D_matrix_dict=ps.D_matrix_dict):
    all_yr_srvl = pd.DataFrame()
    srvls = ["OS", "PFS"]
    yrs = [5, 10]
    for srvl in srvls:
        for yr in yrs:
            x_yr_srvl = get_x_year_srvl(yr, D_matrix_dict=D_matrix_dict, srvl=srvl)
            all_yr_srvl = pd.concat([all_yr_srvl, x_yr_srvl], axis=1)
    return all_yr_srvl


def all_x_year_death(D_matrix_dict=ps.D_matrix_dict):
    all_yr_death = pd.DataFrame()
    deaths = ["cancer_death", "all_cause"]
    yrs = [5, 10]
    for death in deaths:
        for yr in yrs:
            x_yr_death = get_x_year_death(yr, D_matrix_dict=D_matrix_dict, death_state=death)
            all_yr_death = pd.concat([all_yr_death, x_yr_death], axis=1)
    return all_yr_death
    
    
    

def add_clinical_factors(df, D_mat_dict=ps.D_matrix_dict_edit):
    '''
    Add clinical factors to the icers data 
    '''
    # median OS and PFS 
    med_OS = mplt.all_median_survivals(D_mat_dict, ps.death_states)
    med_OS_df = pd.DataFrame.from_dict(med_OS, orient="index", columns=["Median OS"])
    med_PFS = mplt.all_median_survivals(D_mat_dict, ps.disease_states)
    med_PFS_df = pd.DataFrame.from_dict(med_PFS, orient="index", columns=["Median PFS"])
    # 5 and 10 year OS and PFS
    all_srvl_df = all_x_year_srvl(D_matrix_dict=D_mat_dict)
    # 5 and 10 year cancer death and all cause death
    all_death_df = all_x_year_death(D_matrix_dict=D_mat_dict)
    R0_df = R0_status(D_matrix_dict=D_mat_dict)
    R1_df = R0_status(D_matrix_dict=D_mat_dict, status=4)
    surg_df = surg_status(D_matrix_dict=D_mat_dict)
    
    df_add_list = [med_OS_df, med_PFS_df, all_srvl_df, all_death_df, R0_df,
                   R1_df, surg_df, df]
    clinical_df = pd.concat(df_add_list, axis=1, sort=False)
    
    return clinical_df



def save_clinical_df():
    clinical_df = add_clinical_factors(ps.bc_table)
    clinical_df.to_csv(ps.icer/"bc_table_clinical_v2.csv")
    return
    
    
    


            
            
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        