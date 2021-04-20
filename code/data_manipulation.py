# functions that are useful for data manipulation
import numpy as np
import pandas as pd

def csv_to_lists(path):
    df = pd.read_csv(path)
    x_data = df.iloc[:,0]
    y_data = df.iloc[:,1]

	# convert from data frame to list
    x = x_data.values.tolist()
    y = y_data.values.tolist()
    
    return x, y


def frange(start, stop, step):
     i = start
     number_list = list()
     while i < stop:
         number_list.append(i)
         i = i + step
#         print(i)
     return number_list
         

def interp_KM(x, path, step):
# =============================================================================
#     Take range of values (x) to interpolate across using the x and y values
#     from the KM curves (xd and yd, respectively)
# =============================================================================
    xp, fp = csv_to_lists(path)
    y_interp = np.interp(frange(0, x, step), xp, fp)
    y_interp[0] = 100
    return y_interp

def combine_KM(x, path_1, path_2, step, percent_1):
# =============================================================================
#     Combines two KM curves into one curve. Percent_1 is the percentage of 
#     path_1 values that will make up the new curve. Represent percent_1 as 
#     decimal 
# =============================================================================
    percent_2 = 1 - percent_1
    y_interp_1 = interp_KM(x, path_1, step)
    y_interp_2 = interp_KM(x, path_2, step)
    
    y_new = y_interp_1 * percent_1 + y_interp_2 * percent_2
    x_new = x
    
    return x_new, y_new



def excel_to_lists(path, sheet_name):
	df = pd.read_excel(path, sheet_name)
	x_data = df.iloc[:,0]
	y_data = df.iloc[:,1]

	# convert from data frame to list
	x = x_data.values.tolist()
	y = y_data.values.tolist()

	return x, y

def flip(dictionary):
    new_dict = dict((v, k) for k, v in dictionary.items())
    return new_dict

 
def selection(df, keywords, column):
    df_new = df.loc[df[column].isin(keywords)]
    return df_new


def exclusion(df, keywords, column):
    df_new = df.loc[~df[column].isin(keywords)]
    return df_new


def keyword_search(df, array, column):
    keywords = []
    for i in array:
        for j in df[column]:
            print(j)
            if i.casefold() in str(j).casefold():
                keywords.append(j)
    df_new = selection(df, keywords, column)
    return df_new

        
        
