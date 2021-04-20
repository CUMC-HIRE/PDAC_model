#

#markov_owsa
import pandas as pd
import numpy as np
import panc_presets as ps
import markov_icer as mic
import markov_plot as mplt
import simulator_functions as sim
import probability_functions as pf
from copy import deepcopy

def owsa_cu(strats_dfs, cu_name, strat, bound):
    '''
    Iterates over the strategy dataframe to 
    produce a table of ICERs based on the lower and upper bounds of each
    paramter.
 
    Input: Strategy dataframe to be changed, name of the parameter to be changed,
    the upper or lower bound of the paramter
    
    Output: Table of upper and lower ICERs 
    '''
    strats_copy_dict = deepcopy(strats_dfs)
    # make more generazible in the future
    if strat == "folf":
        strat_df = strats_copy_dict[0]
    elif strat == "ag":
        strat_df = strats_copy_dict[1]
    elif strat == "natural history":
        strat_df = strats_copy_dict[2]
#    base_value = strat_df.loc[cu_name, "Value"]
#    owsa_value = strat_df.loc[cu_name, bound]
#        print(base_value)
#        print(owsa_value)
    strat_df.loc[cu_name, "Value"] = strat_df.loc[cu_name, bound]
    owsa_icer = mic.run_icers_w_dom(strats_copy_dict, ps.ALL_STATES, ps.D_matrix_dict_edit)
#    for strat_df, real_strat_df in zip(strats_copy_dict, ps.strats_dfs):
##        print("Before " + str(strat_df.loc[cu_name, "Value"]))
#        strat_df.loc[cu_name, "Value"] = real_strat_df.loc[cu_name, "Value"]
##        print("After " + str(strat_df.loc[cu_name, "Value"]))
    
    return owsa_icer

#e = owsa_cu(ps.strats_dfs, "Chemo cost per cycle", "ag", "Low")
#print(e)

def create_owsa_cu_df(strats_dfs, param_list, strat):
    owsa_df = pd.DataFrame(index=param_list)
    # think of more succinct way of doing this
    base_values = []
    low_values = []
    high_values = []
    # contains, in order, ly, cost, qalys, icers
    low_endpoints = []
    high_endpoints = []
    # iterates over each parameter in the param list
    for cu_name in param_list:
        # gets bost low and high values
        for bound in ["Low", "High"]:
            if bound == "Low":
                base_values.append(ps.arm_to_df[strat].loc[cu_name, "Value"])
                owsa_icer = owsa_cu(strats_dfs, cu_name, strat, bound)
   
                # adds base value to list
                # adds low_value to list
                low_values.append(ps.arm_to_df[strat].loc[cu_name, bound])
                # adds all the endpoints from this variation to list
#                print(owsa_icer.loc["folf", :])
                
                low_endpoints.append(owsa_icer.loc[strat, :].values)
                
            elif bound == "High":
                owsa_icer = owsa_cu(strats_dfs, cu_name, strat, bound)
                # adds high value to list
                high_values.append(ps.arm_to_df[strat].loc[cu_name, bound])
                # adds all the endpoints from this variation to list
                high_endpoints.append(owsa_icer.loc[strat, :].values)
#                if cu_name == "Toxicity cost per month":
#                    print(owsa_icer) 

                
                    
    # create more succinct way of doing this especially
    owsa_df.loc[:, "Base Value"] = base_values
    owsa_df.loc[:, "Low Value"] = low_values
    owsa_df.loc[:, "High Value"] = high_values
#    print(low_endpoints)
    idx = 3
    for name in ["icers", "QALYs", "cost", "life_years"]:
        

        reformat_dict = {"icers": "ICER", "QALYs": "QALYs", "cost": "Costs",
                         "life_years": "Life Years"}
        for bound in ["Low", "Base", "High"]:
            if bound == "Low":
                array = np.array(low_endpoints)[:, idx]
            elif bound == "Base":
                array = ps.bc_table.loc[strat, name]
            elif bound == "High":
                array = np.array(high_endpoints)[:, idx]
            owsa_df.loc[:, (bound + " " + reformat_dict[name])] = array
        
        idx = idx - 1
#        print(idx)
    
    return owsa_df


def gen_owsa_cu(strats_dfs, gen_name, strat, bound):
    '''
    Iterates over the general parameters dataframe to 
    produce a table of ICERs based on the lower and upper bounds of each
    paramter.
   
    Input: Strategy dataframe to be changed, name of the parameter to be changed,
    the upper or lower bound of the paramter
    
    Output: Table of upper and lower ICERs 
    '''
    arms = []
    strat_dict = {0:"folf", 1:"ag", 2:"natural history"}
    owsa_params = deepcopy(ps.gen_params)
    owsa_value = owsa_params.loc[gen_name, bound]
    owsa_params.loc[gen_name, "Value"] = owsa_params.loc[gen_name, bound]
    ind = 0
    for df in strats_dfs:
        if strat_dict[ind] == strat:
            arm = ps.arm(df, gen_params=owsa_params)
        else: 
            arm = ps.arm(df)
        arms.append(arm)
        ind += 1
    owsa_icer = mic.run_icers_arms(arms, ps.ALL_STATES, ps.D_matrix_dict_edit)
    owsa_params.loc[gen_name, "Value"] = ps.gen_params.loc[gen_name, "Value"]
    
    return owsa_icer

#e = gen_owsa_cu(ps.strats_dfs, "Palliative care cost", "ag", "Low")
#print(e)

def create_owsa_gen_df(strats_dfs, param_list, strat):
    owsa_df = pd.DataFrame(index=param_list)
    # think of more succinct way of doing this
    base_values = []
    low_values = []
    high_values = []
    # contains, in order, ly, cost, qalys, icers
    low_endpoints = []
    high_endpoints = []
    # iterates over each parameter in the param list
    for gen_name in param_list:
        # gets bost low and high values
        for bound in ["Low", "High"]:
            if bound == "Low":
                base_values.append(ps.gen_params.loc[gen_name, "Value"])
                owsa_icer = gen_owsa_cu(strats_dfs, gen_name, strat, bound)
                print(owsa_icer)
                # adds base value to list
                # adds low_value to list
                low_values.append(ps.gen_params.loc[gen_name, bound])
                # adds all the endpoints from this variation to list
#                print(owsa_icer.loc["folf", :])
                
                low_endpoints.append(owsa_icer.loc[strat, :].values)
                
            elif bound == "High":
                owsa_icer = gen_owsa_cu(strats_dfs, gen_name, strat, bound)
                print(owsa_icer)
                # adds high value to list
                high_values.append(ps.gen_params.loc[gen_name, bound])
                # adds all the endpoints from this variation to list
                high_endpoints.append(owsa_icer.loc[strat, :].values)
                
                    
    # create more succinct way of doing this especially
    owsa_df.loc[:, "Base Value"] = base_values
    owsa_df.loc[:, "Low Value"] = low_values
    owsa_df.loc[:, "High Value"] = high_values
#    print(low_endpoints)
    idx = 3
    for name in ["icers", "QALYs", "cost", "life_years"]:
        

        reformat_dict = {"icers": "ICER", "QALYs": "QALYs", "cost": "Costs",
                         "life_years": "Life Years"}
        for bound in ["Low", "Base", "High"]:
            if bound == "Low":
                array = np.array(low_endpoints)[:, idx]
            elif bound == "Base":
                array = ps.bc_table.loc[strat, name]
            elif bound == "High":
                array = np.array(high_endpoints)[:, idx]
            owsa_df.loc[:, (bound + " " + reformat_dict[name])] = array
        
        idx = idx - 1
#        print(idx)
    
    return owsa_df

#g = create_owsa_cu_df(ps.strats_dfs, ps.cost_pname_list, "ag")
#d = create_owsa_gen_df(ps.strats_dfs, ps.cost_gname_list, "folf")
#d.to_csv(ps.dump/"Folf_gen_cost_owsa.csv")
        
#def all_owsa_cu(strats): 
    

def owsa_prob(strats, state_1, state_2, bound):
    '''
    Iterates over the strategy dataframe and general parameters dataframe to 
    produce a table of ICERs based on the lower and upper bounds of each
    paramter.
    
    Input: Strategy list (strs), name of from state the prob, 
    name of the to state of the prob to be changed,
    the upper or lower bound of the paramter
    
    Output: Table of upper or lower ICERs 
    '''
    D_matrix_dict = {}
    if bound == "High":
        mult = 2
    elif bound == "Low":
        mult = 0.5
    for strat in strats:

        arm = ps.arm(ps.arm_to_df[strat])
        start_state = sim.get_start_state(arm.CONNECTIVITY)
        temp = deepcopy(ps.t_matrix_dict_edit[strat])
        
                
        prob_matrix = temp[:, state_1, state_2]*mult
        
        while max(prob_matrix) > 1 or max(prob_matrix) < 0:
            mult = round((mult-.014), 10)
            prob_matrix = ps.t_matrix_dict[strat][:, state_1, state_2]*mult
            print(mult)

        
        temp[:, state_1, state_2] = prob_matrix
#        print(prob_matrix)
        for i in range(temp.shape[0]):
            temp[i, state_1] = pf.normalize(temp[i, state_1])
#            temp[i] = sim.normalize_select(temp[i])
            sim.matrix_check(temp[i])
#        print(temp[:, state_1, state_2])
        D_matrix = sim.build_D_matrix(0, start_state, temp)
        D_matrix_dict.update({strat:D_matrix})
#    owsa_icer = mic.run_icers(ps.strats_dfs, ps.ALL_STATES, D_matrix_dict)
    owsa_icer = mic.run_icers_w_dom(ps.strats_dfs, ps.ALL_STATES, D_matrix_dict)
    
    return owsa_icer, D_matrix_dict, mult


def vary_resection_rate(strats, t_matrix_dict, save=False):
    '''
    Updates the t_matrix dictionary based on the specified states status
    
    Input: strats dfs, t_matrix dictionary, states status
    
    Output: new t_matrix dict
    '''
    D_matrix_dict = {}
    states_list = [(0, 6), (0, 7), (0, 10), (6, 10), (7, 10), (8,9), (9,10)]
    for strat in strats:
        t_matrix = t_matrix_dict[strat]
        print(t_matrix)
        for states in states_list:
            state_1 = states[0]
            state_2 = states[1]

            arm = ps.arm(ps.arm_to_df[strat])
            start_state = sim.get_start_state(arm.CONNECTIVITY)
            temp = deepcopy(t_matrix)
            
            if arm.chemo_name.lower() == "folf":
                mult_dict = {(0, 6): 1.14, #1.34
                             (0, 7): 2.1, #1
                             (0, 10): 2.1, #1.5
                             (6, 10): .9, #.9
                             (7, 10): .9, #.9
                             (8, 9): .6, #1.3
                             (9, 10): .6 #1.3
                        }
            elif arm.chemo_name.lower() == "ag":
                mult_dict = {(0, 6): 1.2, #1.1
                             (0, 7): 2.8, #2.3
                             (0, 10): 2.3, #1.5
                             (6, 10): .8, #.8
                             (7, 10): .8, #.8
                             (8, 9): .65, #1.3
                             (9, 10): .6 #1.3
                        }
            elif arm.chemo_name.lower() == "agc":
                mult_dict = {(0, 6): 1.3,
                             (0, 7): 1.3,
                             (0, 10): 1.5,
                             (6, 10): 1,
                             (7, 10): 1,
                             (8, 9): 1.3,
                             (9, 10): 1.3
                        }
            else:
                mult_dict = {(0, 6): 1,
                             (0, 7): 1,
                             (0, 10): 1,
                             (6, 10): 1,
                             (7, 10): 1,
                             (8, 9): 1,
                             (9, 10): 1
                             
                             }
            
            mult_state = mult_dict[(state_1, state_2)]
            prob_matrix = temp[:, state_1, state_2]*mult_state
            
            while max(prob_matrix) > 1 or max(prob_matrix) < 0:
                mult_state = round((mult_state-.014), 10)
                prob_matrix = ps.t_matrix_dict[strat][:, state_1, state_2]*mult_state
            print(mult_state)
    
            
            temp[:, state_1, state_2] = prob_matrix
    #        print(prob_matrix)
            for i in range(temp.shape[0]):
                temp[i, state_1] = pf.normalize(temp[i, state_1])
    #            temp[i] = sim.normalize_select(temp[i])
                sim.matrix_check(temp[i])
    #        print(temp[:, state_1, state_2])
            t_matrix = temp
        D_matrix = sim.build_D_matrix(0, start_state, temp)
        if save == True:
            filename_D = str(arm.chemo_name) +"_D_matrix_edit.npy"
            filename_t = str(arm.chemo_name) +"_t_matrix_edit.npy"
            np.save(ps.matrices/filename_D, D_matrix)
            np.save(ps.matrices/filename_t, t_matrix)
        D_matrix_dict.update({strat:D_matrix})
#    owsa_icer = mic.run_icers(ps.strats_dfs, ps.ALL_STATES, D_matrix_dict)
    owsa_icer = mic.run_icers_w_dom(ps.strats_dfs, ps.ALL_STATES, D_matrix_dict)
    
    return owsa_icer, D_matrix_dict

#l = vary_resection_rate(ps.strats, ps.t_matrix_dict_agc)
#print(l[0])
#print(l[1]["folf"][8])
#print(l[1]["ag"][8])
#print(l[1])
#print(l[1]["agc"][7])
#mplt.OS_plot(l[1]["ag"], ps.death_states, ps.arm(ps.arm_to_df["ag"]))
#mplt.OS_plot(l[1]["folf"], ps.death_states, ps.arm(ps.arm_to_df["folf"]))
#mplt.PFS_plot(l[1]["ag"], ps.disease_states, ps.arm(ps.arm_to_df["ag"]))
#mplt.PFS_plot(l[1]["folf"], ps.disease_states, ps.arm(ps.arm_to_df["folf"]))

def create_owsa_prob_df(strats, strat):
    param_list = [ps.ALL_STATES[key] + '_to_' + ps.ALL_STATES[value] 
                  for key in ps.OWSA_CONNECT.keys() for value in ps.OWSA_CONNECT[key]]
    owsa_df = pd.DataFrame(index=param_list)
    # think of more succinct way of doing this
    base_values = []
    low_values = []
    high_values = []
    # contains, in order, ly, cost, qalys, icers
    low_endpoints = []
    high_endpoints = []
    # iterates over each parameter in the param list
    for state_1 in ps.OWSA_CONNECT.keys():
        # gets bost low and high values
        for state_2 in ps.OWSA_CONNECT[state_1]:
            print(state_1)
            print(state_2)
            for bound in ["Low", "High"]:
                if bound == "Low":
                    base_values.append(1)
                    owsa_icer,l, mult = owsa_prob(strats, state_1, state_2, bound)
                    # adds base value to list
                    # adds low_value to list
                    low_values.append(mult)
                    # adds all the endpoints from this variation to list
#                    print(owsa_icer.loc["folf", :])
                
                    low_endpoints.append(owsa_icer.loc[strat, :].values)
                
                elif bound == "High":
                    owsa_icer,l, mult = owsa_prob(strats, state_1, state_2, bound)
                    # adds high value to list
                    high_values.append(mult)
                    # adds all the endpoints from this variation to list
                    high_endpoints.append(owsa_icer.loc[strat, :].values)
                
                    
    # create more succinct way of doing this especially
    owsa_df.loc[:, "Base Value"] = base_values
    owsa_df.loc[:, "Low Value"] = low_values
    owsa_df.loc[:, "High Value"] = high_values
#    print(low_endpoints)
    idx = 3
    for name in ["icers", "QALYs", "cost", "life_years"]:
        

        reformat_dict = {"icers": "ICER", "QALYs": "QALYs", "cost": "Costs",
                         "life_years": "Life Years"}
        for bound in ["Low", "Base", "High"]:
            if bound == "Low":
                array = np.array(low_endpoints)[:, idx]
            elif bound == "Base":
                array = ps.bc_table.loc[strat, name]
            elif bound == "High":
                array = np.array(high_endpoints)[:, idx]
            owsa_df.loc[:, (bound + " " + reformat_dict[name])] = array
        
        idx = idx - 1
#        print(idx)
    
    return owsa_df

#h =  owsa_prob(ps.strats, 0, 6, "High")
#print(h[1]["folf"][7])
#print(h[0])
#f = create_owsa_prob_df(ps.strats, "ag")
#f.to_csv("test.csv")
#print(f)

def save_owsa_results():
    for strat in ps.strats:
        # all costs csv
        g_costs = create_owsa_gen_df(ps.strats_dfs, ps.cost_gname_list, strat)
        s_costs = create_owsa_cu_df(ps.strats_dfs, ps.cost_pname_list, strat)
        costs_owsa = pd.concat([s_costs,g_costs])
        c_owsa_name = strat + "_cost_owsa.csv"
        costs_owsa.to_csv(ps.owsa/c_owsa_name)
        # all utility csv
        g_util = create_owsa_gen_df(ps.strats_dfs, ps.util_gname_list, strat)
        s_util = create_owsa_cu_df(ps.strats_dfs, ps.util_pname_list, strat)
        util_owsa = pd.concat([s_util, g_util])
        u_owsa_name = strat + "_util_owsa.csv"
        util_owsa.to_csv(ps.owsa/u_owsa_name)
        # all prob csv
        prob_owsa = create_owsa_prob_df(ps.strats, strat)
        prob_owsa_name = strat + "_prob_owsa.csv"
        prob_owsa.to_csv(ps.owsa/prob_owsa_name)
    return
        