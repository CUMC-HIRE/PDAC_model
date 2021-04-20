#PSA

"""
Created on Thu Jun 11 13:47:00 2020

@author: mai2125
"""
import pandas as pd
import numpy as np
import panc_presets as ps
import simulator_functions as sim
import markov_icer as mic
import probability_functions as pf
import matplotlib.pyplot as plt
import csv as csv
import pickle
from scipy.ndimage.filters import gaussian_filter1d

def psa_matrix_check(prob_matrix):
    for item in np.nditer(prob_matrix):
        if item <= 1 and item >=0:
            return True
        elif item > 1:
            return False
        elif item < 0:
            return False
    for row in prob_matrix:
        if int(round(sum(row), 8)) == 1 or int(round(sum(row), 8)) == 0.0:
            return True
        else:
            return False


def params_psa(bound, dist="beta"):
    '''
    Uses a named distribution to randomize the probability parameters of the 
    t_matrix to generate a D_matrix for each strategy
    
    Input: bound is the standard deviation of the mean (base) value
    
    Output: D_matrix dictionary
    '''
    D_matrix_dict = {}
    for strat in ps.strats:

        arm = ps.arm(ps.arm_to_df[strat])
        start_state = sim.get_start_state(arm.CONNECTIVITY)
        temp = ps.t_matrix_dict[strat]
        prob_check = False
        
        while prob_check == False:
            # use lambda for gamma dist
            if dist == "beta":
                beta_dist = lambda x: pf.beta_dist(x, x*bound)
                psa_matrix = np.array([beta_dist(x_i) for x_i in temp.flatten()])
                temp = psa_matrix.reshape(temp.shape[0], temp.shape[1], temp.shape[2])
    
            for i in range(temp.shape[0]):
    
    #            temp[i] = pf.normalize_static(temp[i, state_1], state_2)
                temp[i] = sim.normalize_select(temp[i])

                prob_check = psa_matrix_check(temp[i])
    #        print(temp[:, state_1, state_2])
        D_matrix = sim.build_D_matrix(0, start_state, temp)
        D_matrix_dict.update({strat:D_matrix})
    
    return D_matrix_dict
    
def weight_psa(vec_col, bound, dist):
    '''
    Uses a named distribution to randomize the column of a weight vector
    
    Input: 1-D np.array weight vector column, bound is the standard deviation of 
    the mean (base) value, name of distribution used
    
    Output: randomized weight vector column
    '''
    if dist == "gamma":  
        psa_dist = lambda x: pf.gamma_dist(x, x*bound)
    elif dist == "uniform":
        psa_dist = lambda x: pf.unif_dist(x*(1-bound), x*(1+bound))
    new_vec_col = np.array([psa_dist(x_i) for x_i in vec_col])
    
    return new_vec_col

def wv_psa(arm, states, bound_u, bound_c, dist_u="uniform", dist_c="gamma"):
    '''
    Creates a weight vector (wv) based on the distribution of the utilities
    and costs.
    
    Input: strategy arm, model, states, bound (standard deviation in decimal) 
    used to calculate cost (bound_c) and utilities (bound_u), distribution for 
    utility (dist_u) is uniform and distribution for cost (dist_c) is gamma
    
    Output: randomized weight vector
    '''
    wv = mic.create_weight_vector(arm.cu_dict, states)
    # cost column
    wv[:, 0] = weight_psa(wv[:, 0], bound_c, dist_c)
    # utitl column
    wv[:, 1] = weight_psa(wv[:, 1], bound_u, dist_u)
    
    return wv
    


def run_psa_icers(dom_strats=False):
    '''
    Calculates icers for the strategies for each run of the PSA
    
    Inputs: 
    
    Output: Dataframe containing the icers 
    '''
    # randomizes probabilities with .1 standard deviation from
    D_matrix_dict = params_psa(.02)
    # initialize the dataframe for the costs and qalys of the strategies
    all_ce_df = pd.DataFrame(columns=["life_years", "cost", "QALYs"])
    for strat in ps.strats_dfs:
        arm = ps.arm(strat)
        # randomizes utilities and cost, .2 standard deviation for both
        wv = wv_psa(arm, ps.ALL_STATES, .05, .1)
        wv = wv * ps.discount
        # creates the strat ce series
        strat_ce_df = mic.create_strat_ce_df(D_matrix_dict[arm.chemo_name.lower()], 
                                                       wv, arm.chemo_name.lower())
        # adds the strat series to the overall dataframe
        all_ce_df = pd.concat([all_ce_df, strat_ce_df])
    # gets the icers from the complete ce dataframe
    if dom_strats == False:
        ce_icers = mic.get_icer(all_ce_df)
    else:
        ce_icers = mic.get_icer_w_dom(all_ce_df)
#    print(ce_icers)
    
    return ce_icers



def icer_plot(inc_cost_list, inc_eff_list, pt_color, WTP, save=False):
    

    plt.figure(figsize=(8, 8))
    plt.ylim(0, 200000)
    plt.xlim(0, 3)
    plt.scatter(inc_eff_list, inc_cost_list, c=pt_color, 
                s=5, alpha=0.6, facecolors = None)
    lower = plt.axis()[0]
    upper = plt.axis()[1]
    
    plt.plot([lower, upper], np.multiply([lower, upper], WTP), "r--", alpha=0.5, label='WTP threshold')
    plt.xlabel("Incremental Effectiveness")
    plt.ylabel("Incremental Cost (USD)")
    plt.title("FOLFIRINOX vs. G-nP PSA")
    if save == True:
        file_name = 'PSA_scatter_plot.png'
        plt.savefig(ps.psa/file_name, dpi=400)
    plt.show()
    return

def percent_ce(strat_dict, run_times):
# turns the count into a percentage of cost-effective for all runs
    strat_perc_dict = strat_dict.copy()
    for key in strat_dict.keys():
        percentage = strat_dict[key]/run_times
        print(str(percentage) + "% of iterations were greater than the WTP threshold for " + key)
        strat_perc_dict[key] = percentage
        
    return strat_perc_dict


def WTP_check(df, strat_dict, WTP):
    count = 0
#    print(df)
    while df.iloc[count, -1] > WTP:
        count += 1
        continue
    strat = df.index[count]
    strat_dict[strat] += 1
    return strat_dict

    
def count_best_strat(df, strat_dict, WTP):
    QALYs_sort = df.sort_values("QALYs", ascending=False)
    strat_dict = WTP_check(QALYs_sort, strat_dict, WTP)
    
    return strat_dict


def dom_to_fnl(dom_df, WTP):
    '''
    Converts a dataframe with dominated strategies into a dataframe with only
    non-dominated strategies below the WTP
    
    Input: dataframe with dominated strategies, willingness to pay threshold
    
    Output: final dataframe without dominated strategies
    '''
    # removes dominated strategies
    this_df = dom_df[dom_df.icers != "dominated"]
    # removes strategies over the WTP
    this_df = this_df[this_df.icers <= WTP]
    
    return this_df


def raw_ce_psa(this_df, strat_1, strat_2, WTP):
    '''
    Gives the raw incremental cost and effectiveness for the strategy with the 
    highest QALYs for a given run of the PSA
    
    Input: df with dominated strategies, strat_1 is optimal strat, 
    strat_2 is the next optimal strat
    
    Output: incremental cost and effectiveness from strat 1 to strat 2,
    optimal strat of the run
    '''
    
    inc_cost = this_df.loc[strat_1, "cost"] - this_df.loc[strat_2, "cost"]
    inc_eff = this_df.loc[strat_1, "QALYs"] - this_df.loc[strat_2, "QALYs"]
    # removes dominated strategies above WTP
    this_df = dom_to_fnl(this_df, WTP)
    dom_strat = this_df.index[-1]  
        
    return inc_cost, inc_eff, dom_strat






def run_best_psa(run_times, strat_list, WTP):
# inputs are run_times and strat_list is from the strategies of df of run_CEA
    result_dict = {"dom_strat":[], "inc_cost":[], "inc_eff":[]}
    zeros = [0] * len(strat_list)
    strat_count = dict(zip(strat_list, zeros))
#    print(strat_count)
    print("in run_psa...")
    for i in range(run_times):
        print(i)
        dom_df = run_psa_icers(dom_strats=True)
        df = dom_to_fnl(dom_df, WTP)
        strat_count = count_best_strat(df, strat_count, WTP)
        inc_cost, inc_eff, dom_strat = raw_ce_psa(dom_df, "folf", "ag", WTP)
        result_dict["dom_strat"].append(dom_strat)
        result_dict["inc_cost"].append(inc_cost)
        result_dict["inc_eff"].append(inc_eff)
#            print(strat_count)
    strat_perc_dict = percent_ce(strat_count, run_times)
    return strat_perc_dict, result_dict


# =============================================================================
# def comp_psa_gene(run_times, gene):
#     inc_cost_list = []
#     inc_eff_list = []
#     if gene == "MLH1":
#         comp_list = ps.MLH1_psa_strat
#     elif gene == "MSH2":
#         comp_list = ps.MSH2_psa_strat
#     elif gene == "MSH6":
#         comp_list = ps.MSH6_psa_strat
#     elif gene == "PMS2":
#         comp_list = ps.PMS2_psa_strat
#     else:
#         print("ERROR: Invalid gene entry")
#     
#     count = 0
#     for i in range(run_times):
#         print(i)
#         df = ic.run_CEA()
#         this_df = df[df['gene'] == gene]
#         this_df = dm.selection(this_df, comp_list, "Strategy").reset_index(drop=True)
# #        print(this_df).;
# #        print(this_df["QALYs"])
#         for j in this_df["cost"].values:
#             inc_cost = this_df["cost"].values[0] - this_df["cost"].values[1]
#             inc_cost_list.append(inc_cost)
#         for k in this_df["QALYs"].values:
#             inc_eff = this_df["QALYs"].values[0] - this_df["QALYs"].values[1]
#             inc_eff_list.append(inc_eff)
# #        if inc_cost/inc_eff > 100000:
#     return inc_cost_list, inc_eff_list
# =============================================================================
    
def run_comp_psa(run_times, strat_list, WTP):
#   Outputs a scatter plot for each strategy comparison for each gene
    perc_dict, result_dict = run_best_psa(run_times, strat_list, WTP)
#    print(inc_cost_list)
    icer_plot(result_dict["inc_cost"], result_dict["inc_eff"], 'g')
        
#        launch_plot()
    return
        

def save_dict(dictionary, filename):
    '''
    Saves the winning strategy psa results
    '''
    w = csv.writer(open(ps.psa/filename, "w"))
    for key, val in dictionary.items():
        w.writerow([key, val])
    return
            
def load_strat_dict(filename):
    strat_df = pd.read_csv(filename).dropna()
    strat_df = strat_df.set_index("Strategy")
    strat_dict = strat_df.to_dict()
    return strat_dict
    



#strat_dict = run_psa(1000, ps.strat_list, 100000) 
# make the result into a graph

def psa_plot(strat_dict, save=False):
    
    dom_dict = {x:y for x, y in strat_dict["Value"].items() if y != 0}
    print(dom_dict)
#    plt.title("PSA Optimal Strategy Percentage")
    plt.xlabel("Strategy")
    plt.ylabel("% of runs optimal")
    plt.ylim(0, 1.1)
    plt.bar(*zip(*dom_dict.items()))
    if save == True:
        file_name = 'PSA_all_combined.png'
        plt.savefig(ps.psa/file_name, dpi=400)
    plt.show()
    return

#def psa_plot_gene_h(strat_dict):
#    
#    dom_dict = {x:y for x, y in strat_dict.items() if y != 0}
#    fig, ((plt1, plt2), (plt3, plt4)) = plt.subplots(2, 2, figsize = (17, 14))
#    plt_list = [plt1, plt2, plt3, plt4]
##    plt.suptitle('Optimal Strategy Percentage of PSA',
##                 fontsize = 20, y = 0.94)
#    i = 0
#    for gene in ps.genes:
#        gene_dict = {x.strip(gene)[1:]:y for x, y in dom_dict.items() if gene in x}
#        print(gene_dict)
#        plt_list[i].set_title(gene, fontsize=14)
#        plt_list[i].set_xlabel("Strategy", fontsize=14)
#        plt_list[i].set_ylabel("% of runs optimal", fontsize=14)
#        if gene == "MLH1":
#            c = 'g'
#        elif gene == "MSH2":
#            c = 'r'
#        elif gene == "MSH6":
#            c = 'k'
#        elif gene == "PMS2":
#            c = 'm'
#        for item in (plt_list[i].get_xticklabels() + plt_list[i].get_yticklabels()):
#            item.set_fontsize(12)
#        plt_list[i].bar(*zip(*gene_dict.items()), color=c)
#        plt_list[i].set_ylim(0, 1.1)
#        
#        i += 1
#    file_name = 'PSA_combined_gender_composite_1.png'
#    plt.savefig(ps.psa/file_name, dpi=500)
#    plt.show()
#    return

def psa_plot_from_filename(filename):
    strat_dict = load_strat_dict(filename)["Value"]
#    psa_plot_gene_h(strat_dict)    
    return strat_dict

def psa():
    strat_dict = run_psa(10, ps.strats, 100000)
    icer_plot(strat_dict)
    return


def psa_best():
    WTP = 100000
    run_times = 1000
    strat_dict, result_dict = run_best_psa(run_times, ps.strats, WTP)
    icer_plot(result_dict["inc_cost"], result_dict["inc_eff"], "g", WTP, save=True)
    save_dict(result_dict, "psa_folf_vs_gen_results.csv")
    save_dict(strat_dict, "psa_pct_results.csv")
    psa_plot(strat_dict, save=True)
#    psa_plot_gene_h(strat_dict)
    return 
        
def get_acceptability_perc(df_list, strat_list, WTP, run_times):
# inputs are run_times and strat_list is from the strategies of df of run_CEA
    result_dict = {"dom_strat":[], "inc_cost":[], "inc_eff":[]}
    zeros = [0] * len(strat_list)
    strat_count = dict(zip(strat_list, zeros))
#    print(strat_count)
    print("in run_psa...")
    for dom_df in df_list:
        df = dom_to_fnl(dom_df, WTP)
        strat_count = count_best_strat(df, strat_count, WTP)
        inc_cost, inc_eff, dom_strat = raw_ce_psa(dom_df, "folf", "ag", WTP)
        result_dict["dom_strat"].append(dom_strat)
        result_dict["inc_cost"].append(inc_cost)
        result_dict["inc_eff"].append(inc_eff)
#            print(strat_count)
    strat_perc_dict = percent_ce(strat_count, run_times)
    return strat_perc_dict

def run_acceptability_curve(run_times, strat_list):
    # run 1000 icers
    df_list = [run_psa_icers(dom_strats=True) for i in range(run_times)]
    # save runs
    with open("psa_runs.txt", "wb") as output:
        pickle.dump(df_list, output)
    result_dict = {}
    for WTP in range(0, 210000, 10000):
        strat_perc_dist = get_acceptability_perc(df_list, strat_list, WTP, run_times)
        result_dict.update({WTP:strat_perc_dist})
    return result_dict

def run_acpt_curve_load(filename, run_times, strat_list):
    pickle_file = open(filename, "rb")
    df_list = pickle.load(pickle_file)
    result_dict = {}
    for WTP in range(0, 200000, 50000):
        strat_perc_dist = get_acceptability_perc(df_list, strat_list, WTP, run_times)
        result_dict.update({WTP:strat_perc_dist})
    return result_dict


def plot_acceptability_curve(result_dict):
    plt.xlabel("Willingness-to-pay")
    plt.ylabel("% Cost-effective")
    x_values = list(result_dict.keys())
    df = pd.DataFrame.from_dict({i: result_dict[i] for i in result_dict.keys()},
                                 orient="index")
    print(df)
    for i in range(df.shape[1]):
        y_values = list(df.iloc[:,i])
        y_smooth = gaussian_filter1d(y_values, sigma=2)
        plt.plot(x_values, y_smooth)
    plt.legend(["FOLFIRINOX", "G-nP", "Natural History"])
    plt.xticks([0, 50000, 100000, 150000, 200000])
    plt.show()
    return
    

def run_accept_curve():
    run_times = 1000
    results_dict = run_acceptability_curve(run_times, ps.strats)
    plot_acceptability_curve(results_dict)
    return

def run_accept_curve_load(filename):
    run_times = 1000
    result_dict = run_acpt_curve_load(filename, run_times, ps.strats)
    plot_acceptability_curve(result_dict)
    return

    
def psa_graph():
    return
    
    

def psa_bar_graph():
    return