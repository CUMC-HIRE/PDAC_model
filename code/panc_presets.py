# panc_presets

# authors: Myles Ingram and Brianna Lauren

import pathlib as pl
import pandas as pd
import probability_functions as pf
import numpy as np
#import data_manipulation as dm

# path to folders
src = pl.Path.cwd().parent
data_repo = src/"data"
dump = src/"dump"
icer = dump/"icer"
matrices = dump/"matrices"
graphs = dump/"graphs"
owsa = dump/"owsa"
psa = dump/"psa"


# model goes until death
CYCLE_LENGTH = 1/12 # units = months
START_AGE = 60
END_AGE = 72

RUN_TIME = int((END_AGE-START_AGE)/CYCLE_LENGTH)

discount = .97 # discount rate of .3

ALL_STATES = {0: "neoadj",
              1: "resect",
              2: "upfront_resect",
              3: "R0",
              4: "R1",
              5: "adj_chemo",
              6: "palliative",
              7: "2nd_line",
              8: "non_recur",
              9: "recurrence",
              10: "cancer_death",
              11: "all_cause",
              12: "comp_death"
              }

# list of connections to be modified during the owsa
OWSA_CONNECT = {0: [0, 1, 6, 7, 10], # 8 cycles
                1: [3, 4, 12],
                3: [8, 9],
                4: [8, 9],
                5: [],
                6: [6, 10],
                7: [7, 10],
                8: [8, 9],
                9: [9, 10]
                }


death_states = [10, 11, 12]

disease_states = [9, 10, 11, 12]

life_states = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

surg_mnth = 7
status_mnth = 8
post_surg_mnth = 9


WTP = 150000

strats = ["folf", "ag", "natural history"]#, "agc"]

neoadj_list = ["folf", "ag"]

D_matrix_dict = {"folf": np.load(matrices/"FOLF_D_matrix_use.npy"), 
                 "ag": np.load(matrices/"AG_D_matrix_use.npy"),
                 "natural history": np.load(matrices/"Natural History_D_matrix_use.npy")} 
#                 "gem": np.load(dump/"matrices"/"GEM_D_matrix_use.npy")} 
#                 "pexg": np.load(dump/"PEXG_D_matrix_use.npy")}
#                 "gemcap": np.loadtx("gemcap_D_matrix.npy")}

D_matrix_dict_agc = {"folf": np.load(matrices/"FOLF_D_matrix_use.npy"), 
                 "ag": np.load(matrices/"AG_D_matrix_use.npy"),
                 "natural history": np.load(matrices/"Natural History_D_matrix_use.npy"),
                 "agc": np.load(matrices/"AGC_D_matrix_test.npy")} 
#                 "gem": np.load(dump/"matrices"/"GEM_D_matrix_use.npy")} 
#                 "pexg": np.load(dump/"PEXG_D_matrix_use.npy")}
#                 "gemcap": np.loadtx("gemcap_D_matrix.npy")}

D_matrix_dict_edit = {"folf": np.load(matrices/"FOLF_D_matrix_edit.npy"), 
                 "ag": np.load(matrices/"AG_D_matrix_edit.npy"),
                 "agc": np.load(matrices/"AGC_D_matrix_edit.npy"),
                 "natural history": np.load(matrices/"Natural History_D_matrix_edit.npy")}
#                 
# =============================================================================
# D_matrix_dict = {"folf": np.load(dump/"matrices"/"FOLF_D_matrix_use.npy"), 
#                  "ag": np.load(dump/"matrices"/"AG_D_matrix_use.npy"),
#                  "natural history": np.load(dump/"matrices"/"Natural History_D_matrix_use.npy")} 
# #                 "gem": np.load(dump/"matrices"/"GEM_D_matrix_use.npy")} 
# #                 "pexg": np.load(dump/"PEXG_D_matrix_use.npy")}
# #                 "gemcap": np.loadtx("gemcap_D_matrix.npy")}
# =============================================================================

t_matrix_dict = {"folf": np.load(matrices/"FOLF_t_matrix_use.npy"),
                 "ag": np.load(matrices/"AG_t_matrix_use.npy"),
                 "natural history": np.load(matrices/"Natural History_t_matrix_use.npy")}
#                 "gem": np.load(dump/"matrices"/"GEM_t_matrix_use.npy")}
#                 "pexg": np.load(dump/"PEXG_D_matrix_use.npy")}
##                 "gemcap": np.loadtx("gemcap_D_matrix.npy")}

t_matrix_dict_agc = {"folf": np.load(matrices/"FOLF_t_matrix_use.npy"),
                 "ag": np.load(matrices/"AG_t_matrix_use.npy"),
                 "natural history": np.load(matrices/"Natural History_t_matrix_use.npy"),
                 "agc": np.load(matrices/"AGC_t_matrix_test.npy")}
#                 "gem": np.load(dump/"matrices"/"GEM_t_matrix_use.npy")}
#                 "pexg": np.load(dump/"PEXG_D_matrix_use.npy")}
##                 "gemcap": np.loadtx("gemcap_D_matrix.npy")}
t_matrix_dict_edit = {"folf": np.load(matrices/"FOLF_t_matrix_edit.npy"),
                 "ag": np.load(matrices/"AG_t_matrix_edit.npy"),
                 "natural history": np.load(matrices/"Natural History_t_matrix_edit.npy")}
#                 "agc": np.load(matrices/"AGC_t_matrix_edit.npy")}

RUN_MODES = ["CALIBRATE", "MARKOV", "PLOT", "ICER", "POST_ICER", "PRE_ICER",
             "OWSA", "MORB"]



bc_table = pd.read_csv(icer/'panc_icers_bc_v2.csv', index_col="arm")   

# ARM PARAMETERS
params = pd.ExcelFile(data_repo/"model_params"/'params.xlsx')

gen_params = params.parse("General Parameters", index_col="Parameters")

# costs for all arms
costs = params.parse("Costs", index_col="Parameters")
# utilities for all arms
util = params.parse("Utilities", index_col="Parameters")
# Folfirinox Arm params
folf = params.parse("FOLFIRINOX", index_col="Parameters")
# GEM arm params
gem = params.parse("Adj chemo (GEM)", index_col="Parameters")
# GEMCAP arm params
gemcap = params.parse("Adj chemo (GEMCAP)", index_col="Parameters")
# GEMABX (Gem Nab-Paclitaxel arm params
ag = params.parse("GemAbx", index_col="Parameters")
# PAXG arm params
agc = params.parse("GemAbxCis", index_col="Parameters")

pexg = params.parse("PEXG", index_col="Parameters")
# natural history arm params
nat_hist = params.parse("Nat_Hist", index_col="Parameters")

folf_OS = data_repo/"survival_curves"/"FOLF_OS_Dhir.csv" # Dhir
folf_PFS = data_repo/"survival_curves"/"FOLF_PFS.csv" # Michelakos
gem_OS = data_repo/"survival_curves"/"GEM_OS.csv" # ESPAC 4
gem_PFS = data_repo/"survival_curves"/"GEM_PFS.csv" # ESPAC 4
#gemcap_OS = data_repo/"survival_curves"/"GEMCAP_OS.csv" # ESPAC 4
#gemcap_PFS = data_repo/"survival_curves"/"GEMCAP_PFS.csv" # ESPAC 4
ag_OS = data_repo/"survival_curves"/"GEMABX_OS.csv" # Dhir
ag_PFS = data_repo/"survival_curves"/"GEMABX_PFS.csv" # Ielpo
nh_OS = data_repo/"survival_curves"/"NAT_HIST_OS_Shapiro.csv" #Shapiro
#pexg_OS = data_repo/"survival_curves"/"PEXG_OS.csv"
#pexg_PFS = data_repo/"survival_curves"/"PEXG_PFS.csv"


arm_to_df = {"folf": folf,
             "ag": ag,
             "natural history": nat_hist}
#             "agc": agc}
#             "gem": gem,
#             "pexg": pexg}

strats_dfs = [folf, ag, nat_hist]#, agc]#, gem] #, gemcap]
strats_dfs_cis = [agc]


class arm:
    '''
    This class contains the parameters used for each arm of the model. The arm
    is given a name on input and that name determines what the parameters of 
    the model will be.
    '''
    START_AGE = 60
    
    # all the states in all of the arms of the model
    states = ALL_STATES
    
    gen_params = gen_params
    
    def __init__(self, param_df, gen_params=gen_params):
        # Input: strategy-specific paramete dataframe
        # to make more readable and scalable have parameter dataframe as input and 
        # have a value in the table be neoadjuvant or adjuvant
#        if name == "folf":
#            param_df = folf
#            self.chemo_len = param_df.loc["Complete cycles", "Value"]
#            self.second_line = param_df.loc["Second-line survival, months", "Value"]
#        elif name == "gem":
#            param_df = gem
#
#        elif name == "gemcap":
#            param_df = gemcap
#            self.second_line = param_df.loc["Second-line survival, months", "Value"]
#            self.chemo_len = param_df.loc["Chemotherapy cycle length (months)", "Value"]
            
        # should remove from class to conserve space
        self.chemo_name = param_df.loc["Chemo Name", "Value"]
        self.chemo_type = param_df.loc["Chemo Type", "Value"]
        self.second_line = param_df.loc["Second-line survival, months", "Value"]
        self.chemo_len = param_df.loc["Chemotherapy cycle length (months)", "Value"]
        self.cycle_len = param_df.loc["Complete cycles", "Value"]
        self.pdac_death = pf.prob_conv(gen_params.loc["PDAC mortality", "Value"], 1, 1/60)
        self.surg_death = gen_params.loc["Surgical mortality", "Value"]
        self.panc_fist = gen_params.loc["Pancreatic fistula", "Value"]
        self.second_line_prob = pf.rate_to_prob(0.65, 1/2)
        self.second_line_srvl = pf.rate_to_prob(pf.exp_half_life_rate(self.second_line), 1)
        self.pall_care_prob = pf.rate_to_prob(0.35, 1/2)
        self.pall_srvl = pf.rate_to_prob(pf.exp_half_life_rate(4), 1)
        self.nat_hist_srvl = pf.rate_to_prob(pf.exp_half_life_rate(6), 1)
        self.dropout_rate = pf.prob_conv(param_df.loc["Dropout rate", "Value"],1, 1/self.chemo_len)

        self.tox_rate = pf.prob_conv(param_df.loc["Toxicity rate", "Value"], 1, 1/self.chemo_len)
        self.surg_compl_rate = param_df.loc["Surgical complication rate", "Value"]
        self.panc_fist_rate = param_df.loc["Pancreatic fistula rate", "Value"]
        self.R0_rate = param_df.loc["R0 rate", "Value"]
        self.pdac_recur = pf.rate_to_prob(param_df.loc["PDAC recurrence", "Value"], 1/12)
        self.hosp_tox = param_df.loc["Hospitalization for toxicity", "Value"]
        self.lymph_pos = param_df.loc["Lymph node positivity", "Value"]
        self.R0_recur_srvl = pf.rate_to_prob(pf.exp_half_life_rate(param_df.loc["R0 recurrence survival", "Value"]), 1)
        self.R1_recur_srvl = pf.rate_to_prob(pf.exp_half_life_rate(param_df.loc["R1 recurrence survival", "Value"]), 1)
        self.N0_recur_srvl = pf.rate_to_prob(pf.exp_half_life_rate(param_df.loc["N0 recurrence survival", "Value"]), 1)
        self.N1_recur_srvl = pf.rate_to_prob(pf.exp_half_life_rate(param_df.loc["N1 recurrence survival", "Value"]), 1)
        # costs
        self.resect_cost = gen_params.loc["Resection surgery cost", "Value"] 
        self.cap_cost = gen_params.loc["Capecitabine/radiation per month", "Value"]
        self.pall_cost = gen_params.loc["Palliative care cost", "Value"]/4 
        self.chemo_hosp_cost = gen_params.loc["Chemoradiation hospitalization costs", "Value"] * self.hosp_tox
        self.pdac_cost = gen_params.loc["PDAC costs per month (inpatient)", "Value"]
        self.admin_cost = param_df.loc["Administration cost per month", "Value"]
        self.chemo_cost = param_df.loc["Chemotherapy cost per cycle", "Value"] + self.admin_cost
#        self.tox_cost = param_df.loc["Toxicity cost per cycle", "Value"] * self.tox_rate * self.cycle_len 
        self.tox_cost = param_df.loc["Toxicity cost per month", "Value"] 
        self.second_line_cost = param_df.loc["Second-line cost per month", "Value"] 
        self.pdac_screen = 1525 # cost of endoscopic ultrasound from et al 
        # utility
        self.pdac_util = gen_params.loc["Progressive disease", "Value"]
        self.pall_util = gen_params.loc["Palliative care", "Value"]
        self.surgery_util = gen_params.loc["Recovery from surgery", "Value"]
        self.chemo_disutil = param_df.loc["Chemotherapy disutility", "Value"] * self.cycle_len / self.cycle_len
        self.tox_disutil = param_df.loc["Chemotherapy toxicity disutility", "Value"] * self.tox_rate / self.cycle_len 
        # should be in probability functions
        def recur_srvl(R0, R0_D, R1_D, N0, N0_D, N1_D):
            '''
            Uses law of total probability to determine the survival rate of
            recurrence
            
            Input: R0, R0_survival, R1_survival, N0, N0_survival, N1_survival
            
            Output: Recurrence survival probability per month
            '''
            recur_srvl = R0_D*R0 + R1_D*(1-R0) + N0_D*N0 + N1_D*(1-N0)
            return recur_srvl
        
        self.recur_srvl = recur_srvl(self.R0_rate, self.R0_recur_srvl, 
                                     self.R1_recur_srvl, self.lymph_pos, 
                                     self.N0_recur_srvl, self.N1_recur_srvl)
        # first state is the from state and the second state is the to state
        if self.chemo_type == "Neoadjuvant":
            self.params_dict = {("neoadj", "palliative"): self.dropout_rate * self.pall_care_prob, #.01909, #self.dropout_rate * self.pall_care_prob, #calibrated estimate
               ("neoadj", "2nd_line"): self.dropout_rate * self.second_line_prob,
               ("neoadj", "cancer_death"): self.pdac_death,
               ("resect", "R0"): self.R0_rate,
#               ("resect", "recurrence"): self.pdac_recur,
               ("resect", "comp_death"): self.surg_death,
               ("palliative", "cancer_death"): .057631, #self.pall_srvl # calibrated estimate
               ("2nd_line", "cancer_death"): self.second_line_srvl,
               ("recurrence", "cancer_death"): self.recur_srvl,
               ("non_recur", "recurrence"): self.pdac_recur,
               ("R0","recurrence"): self.pdac_recur,
               ("R1", "recurrence"): self.pdac_recur,
               }
            # state is the state at which the (cost, utility) is applied
            if self.chemo_name == "FOLF":
# =============================================================================
#                 self.cu_dict = {"neoadj": ((self.tox_cost + self.chemo_hosp_cost) , 
#                                      (.8 + (self.chemo_disutil + self.tox_disutil))),
# #                "resect": (self.resect_cost, self.surgery_util),
# #                "R0": (0, .8),
# #                "R1": (0, .8),
#                 "palliative": (self.pall_cost, self.pall_util),
#                 "2nd_line": (self.second_line_cost, self.pdac_util)#,
# #                "non_recur": (self.pdac_screen, .8), # check this 
# #                "recurrence": (self.pdac_cost, self.pdac_util)
#                 }
# =============================================================================
                self.cu_dict = {"neoadj": ((self.chemo_cost + self.tox_cost + self.chemo_hosp_cost + self.cap_cost) , 
                                     (.8 + (self.chemo_disutil + self.tox_disutil))),
                "resect": (self.resect_cost, self.surgery_util),
                "R0": (0, .8),
                "R1": (0, .8),
                "palliative": (self.pall_cost, self.pall_util),
                "2nd_line": (self.second_line_cost, self.pdac_util),
                "non_recur": (self.pdac_screen, .8), # check this 
                "recurrence": (self.pdac_cost, self.pdac_util)
                    }
            else:
# =============================================================================
#                 self.cu_dict = {"neoadj": ((self.tox_cost + self.chemo_hosp_cost) , 
#                                      (.8 + (self.chemo_disutil + self.tox_disutil))),
# #                "resect": (self.resect_cost, self.surgery_util),
# #                "R0": (0, .8),
# #                "R1": (0, .8),
#                 "palliative": (self.pall_cost, self.pall_util),
#                 "2nd_line": (self.second_line_cost, self.pdac_util)#,
# #                "non_recur": (self.pdac_screen, .8), # check this 
# #                "recurrence": (self.pdac_cost, self.pdac_util)
#                 }
# =============================================================================
                self.cu_dict = {"neoadj": ((self.chemo_cost + self.tox_cost + self.chemo_hosp_cost) , 
                                     (.8 + (self.chemo_disutil + self.tox_disutil))),
                "resect": (self.resect_cost, self.surgery_util),
                "R0": (0, .8),
                "R1": (0, .8),
                "palliative": (self.pall_cost, self.pall_util),
                "2nd_line": (self.second_line_cost, self.pdac_util),
                "non_recur": (self.pdac_screen, .8), # check this 
                "recurrence": (self.pdac_cost, self.pdac_util)
                }
                
            self.CONNECTIVITY = {0: [1, 6, 7, 10, 11], # 8 cycles
                     1: [3, 4,11, 12],
                     2: [],
                     3: [8, 9, 11],
                     4: [8, 9, 11],
                     5: [],
                     6: [6, 10, 11],
                     7: [7, 10, 11],
                     8: [8, 9, 11],
                     9: [9, 10, 11],
                     10: [10],
                     11: [11],
                     12: [12]
                     }
            self.CHEMO_connects = {0: [0, 6, 7, 10, 11], 1: [], 2: [], 3:[], 
                                   4: [], 5: [], 6: [6, 10, 11], 
                                   7: [7, 10, 11], 8: [], 9: [], 
                                   10: [10], 11: [11], 
                                   12: []}
            self.calib_params = [("neoadj", "cancer_death"),("recurrence", "cancer_death"),
                ("neoadj", "2nd_line"), ("neoadj", "palliative")
                , ("non_recur", "recurrence")]
#            ("palliative", "cancer_death")
            # put these in a different class
            if self.chemo_name == "FOLF":
                self.OS_target = folf_OS
                self.OS_nodes = [11, 23]
                self.PFS_target = folf_PFS
                self.PFS_nodes = [11, 23]
            elif self.chemo_name == "AG":
                self.OS_target = ag_OS
                self.OS_nodes = [11, 23]
                self.PFS_target = ag_PFS
                self.PFS_nodes = [11, 23]
            elif self.chemo_name == "AGC":
                self.OS_target = folf_OS
                self.OS_nodes = [11, 23]
                self.PFS_target = folf_PFS
                self.PFS_nodes = [11, 23]

        elif self.chemo_type == "NAT_HIST":
            self.params_dict = {("neoadj", "cancer_death"): self.pdac_death, # estimate to be calibrated
                                ("neoadj", "recurrence"): .05, # arbitrary starting estimate. Calibrated through model
                                ("recurrence", "cancer_death"): self.nat_hist_srvl # Shapiro 2016 et al
               }
# =============================================================================
#             self.cu_dict = {"adj_chemo": ((self.tox_cost + self.chemo_hosp_cost), 
#                                      (.8 + (self.chemo_disutil + self.tox_disutil))),
# #                "upfront_resect": (self.resect_cost, self.surgery_util),
# #                "R0": (0, .8),
# #                "R1": (0, .8),
#                 "palliative": (self.pall_cost, self.pall_util),
#                 "2nd_line": (self.second_line_cost, self.pdac_util)#,
# #                "non_recur": (self.pdac_screen, .8),
# #                "recurrence": (self.pdac_cost, self.pdac_util)
#                 }
# =============================================================================
            self.cu_dict = {"neoadj": (0, .8),
                "recurrence": (self.pall_cost, self.pall_util)
                }
            
            self.CONNECTIVITY = {0: [0, 9, 10, 11], # 8 cycles
                     1: [],
                     2: [],
                     3: [],
                     4: [],
                     5: [],
                     6: [],
                     7: [],
                     8: [],
                     9: [9, 10, 11],
                     10: [10],
                     11: [11],
                     12: []
                     }
            self.calib_params =  [("recurrence", "cancer_death"), ("neoadj", "recurrence"), 
                                  ("neoadj", "cancer_death")]
            self.OS_target = nh_OS
            self.OS_nodes = [12, 30]


            
# =============================================================================
# class arm_costs:
#     
#     def __init__(self, name):
#         
# class arm_util:
#     
#     def __init__(self, name)
# =============================================================================
                
# list of arms
folf_arm = arm(folf)
ag_arm = arm(ag)
gem_arm = arm(gem)

arm_list = [folf_arm, ag_arm, gem_arm]

# extracting probabilites from life tables for male
xl = pd.ExcelFile(data_repo/"morb_tables"/'Weighted_lifetable.xlsx')
lifeTable = xl.parse('Sheet1')
LT_prob_male = lifeTable.iloc[:, 2]
LT_rate_male = -(np.log(1-LT_prob_male))


# use ac mortality data beginning at starting age:
startIndex = int(START_AGE)
cut_ac_mortality = LT_prob_male[startIndex:]

all_cause_OS = data_repo/"morb_tables"/'all_cause_survival.npy'

LT_prob_female = lifeTable.iloc[:, 1]
LT_rate_female = -(np.log(1-LT_prob_female))

# turn into months

# use ac mortality data beginning at starting age:
cut_ac_mortality_fm = LT_prob_female[startIndex:]

LT_cmb = (cut_ac_mortality + cut_ac_mortality_fm)/2

xl_morb = pd.ExcelFile(data_repo/"morb_tables"/"20200219_comorb_life_table.xlsx")
comorb_table = xl_morb.parse("LifeTablesWithCancer")

no_fem = comorb_table.iloc[:35, 5]
no_male = comorb_table.iloc[2300:2335, 5]
low_fem = comorb_table.iloc[35:70, 5]
low_male = comorb_table.iloc[2335:2370, 5]
mod_fem = comorb_table.iloc[70:105, 5]
mod_male = comorb_table.iloc[2370:2405, 5]
sev_fem = comorb_table.iloc[105:140, 5]
sev_male = comorb_table.iloc[2405:2440, 5]

# first 6 years of life table
LT_start = LT_cmb[:6]
LT_no_cmb = pd.concat([LT_start, (no_fem + no_male.values)/2], ignore_index=True)
LT_low_cmb = pd.concat([LT_start, (low_fem + low_male.values)/2])
LT_mod_cmb = pd.concat([LT_start, (mod_fem + mod_male.values)/2])
LT_sev_cmb = pd.concat([LT_start, (sev_fem + sev_male.values)/2])


LT_no_cmb.index = range(60, 101)
LT_low_cmb.index = range(60, 101)
LT_mod_cmb.index = range(60, 101)
LT_sev_cmb.index = range(60, 101)

morb_table_dict = {"no": LT_no_cmb,
                         "low": LT_low_cmb,
                         "mod": LT_mod_cmb,
                         "sev": LT_sev_cmb
        }
# owsa cost names
# costs in the parameter tables
cost_pname_list = folf.index[18:22]
cost_gname_list = gen_params.index[6:11]
util_pname_list = folf.index[25:]
util_gname_list = gen_params.index[13:]


            
    
        


