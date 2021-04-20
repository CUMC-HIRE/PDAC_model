# main 

import panc_presets as ps
import simulator_functions as sim
import markov_plot as mplt
import markov_icer as mic
import owsa as ow
import probability_functions as pf
import numpy as np
from copy import copy
import csv


def main(run_mode):
    
    assert run_mode in ps.RUN_MODES, "invalid run_mode given"
    
    if run_mode == "MARKOV":
        for strat in ps.strat_dfs:
            arm = ps.arm(strat)
            t_matrix, D_matrix = sim.run_markov(arm, ps.RUN_TIME, 
                                            ps.LT_cmb, 
                                            arm.CHEMO_connects, arm.chemo_len)
            mplt.OS_plot(D_matrix, ps.death_states, strat)
            mplt.PFS_plot(D_matrix, ps.disease_states, strat)
        
    elif run_mode == "CALIBRATE":
        
        LT = ps.LT_cmb
        for strat in ps.strats_dfs:
#            print(strat)
            arm = ps.arm(strat)
            if arm.chemo_type == "NAT_HIST":
                ans = sim.calibrate_markov_single(arm, ps.RUN_TIME, LT, 3.6, 1)
            else:
                ans = sim.calibrate_markov(arm, ps.RUN_TIME, LT, 6, 1, 
                                           arm.CHEMO_connects, arm.chemo_len)
            mplt.OS_plot(ans[0][1], ps.death_states, arm)
            mplt.PFS_plot(ans[0][1], ps.disease_states, arm)
            np.save(ps.matrices/(arm.chemo_name + "_D_matrix_test.npy"), ans[0][1])
            np.save(ps.matrices/(arm.chemo_name + "_t_matrix_test.npy"), ans[0][0])
            print(str(strat) + ":" + str(ans[0][1][ps.surg_mnth][1]))
            
            
#            np.savetxt(ps.dump/(arm.chemo_name + "_D_matrix_test.csv"), ans[0][1], delimiter=",")
#            print(arm.chemo_name + ": " + str(ans[0][2]))
        
            
    elif run_mode == "ICER":
        ce_icers = mic.run_icers_w_dom(ps.strats_dfs, ps.ALL_STATES, ps.D_matrix_dict_edit)
#        ce_icers.to_csv(ps.icer/"panc_icers_coe.csv")
        return ce_icers
    # string ICER and CALIBRATE into one function
    elif run_mode == "POST_ICER":
        post_D_matrix_dict = ps.D_matrix_dict_edit.copy()
        for key in post_D_matrix_dict.keys():
            # updates each D_matrix in the copy of the dictionary to the
            # D_matrix without the chemo months
            d_matrix = post_D_matrix_dict[key]
            post_D_matrix_dict.update({key:d_matrix[ps.post_surg_mnth:]})
        ce_icers = mic.run_icers_w_dom(ps.strats_dfs, ps.ALL_STATES, post_D_matrix_dict)
        ce_icers.to_csv(ps.icer/"panc_post_surg_icers_v3.csv")
        return ce_icers
    
    elif run_mode == "PRE_ICER":
        post_D_matrix_dict = ps.D_matrix_dict_edit.copy()
        for key in post_D_matrix_dict.keys():
            # updates each D_matrix in the copy of the dictionary to the
            # D_matrix without the chemo months
            d_matrix = post_D_matrix_dict[key]
            post_D_matrix_dict.update({key:d_matrix[:ps.post_surg_mnth]})
        ce_icers = mic.run_icers_w_dom(ps.strats_dfs, ps.ALL_STATES, post_D_matrix_dict)
        ce_icers.to_csv(ps.icer/"panc_pre_surg_icers_v3.csv")
        return ce_icers
            
#        
    elif run_mode == "PLOT":
        mplt.comp_plot_OS(ps.D_matrix_dict_edit, ps.death_states)
        mplt.comp_plot_PFS(ps.D_matrix_dict_edit, ps.disease_states)
        for key in ps.D_matrix_dict.keys():
            mplt.OS_plot(ps.D_matrix_dict[key], ps.death_states, ps.arm(ps.arm_to_df[key]))
            mplt.PFS_plot(ps.D_matrix_dict[key], ps.disease_states, ps.arm(ps.arm_to_df[key]))
            
        OS_srvl_dict = mplt.all_median_survivals(ps.D_matrix_dict, ps.death_states)
        PFS_srvl_dict = mplt.all_median_survivals(ps.D_matrix_dict, ps.disease_states)
        return OS_srvl_dict, PFS_srvl_dict
    
    elif run_mode == "OWSA":
        ow.save_owsa_results()
        return
    
    elif run_mode == "MORB":
        for morb in ps.morb_table_dict.keys():
            D_matrix_dict = sim.update_matrix_dicts(ps.strats_dfs, ps.t_matrix_dict_edit, morb)
            ce_icers = mic.run_icers(ps.strats_dfs, ps.ALL_STATES, D_matrix_dict)
            clin_icers = mic.add_clinical_factors(ce_icers, D_mat_dict=D_matrix_dict)
            print(clin_icers)
            filename = "panc_icers_" + morb + "_comorb_v3.csv"
            clin_icers.to_csv(ps.icer/filename)
        return clin_icers
    
        
        
    return
    
    
        
        