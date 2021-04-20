# Comorbidity Life Tables
This README was created by Brianna Lauren (brianna.lauren130@gmail.com) on February 20, 2020.
It explains the derivation of the comorbidity-specific life tables, based on the methods of Amy Knudsen and Claudia Seguin at MGH. 
The script is in comorb_lt.R

---
## Data
The README of the file 20200219_comorb_life_table.xlsx explains the sources of the data. To summarize, comorbidity specific cancer-free life tables were obtained from Lansdorp-Vogelaar et al (2014) via Amy Knudsen. Cancer mortality rates were obtained from the CDC Wonder database. 

---
## Code
* The working directory and file name are set at the top
* Data is read in as data tables (DT) from the spreadsheet for cancer-free mortality and cancer mortality.
* The `DT[i, j, by]` syntax is used to calculate total cancer mortality rates by age and gender.
* The ICD-10 codes specific to the cancer of interest (i.e., pancreatic cancer) are read in from the spreadsheet.
* The ICD-10 codes and the `DT[i, j, by]` syntax is used to calculate mortality rates specific to pancreatic cancer by age and gender.
* A for loop runs through each column of the cancer-free mortality table
    * The column names are used to determine the cohort age (`substr` first two characters), gender, and comorbidity. `grepl(pattern, x)` returns TRUE is the pattern is contained within x.
    * While age is less than or equal to 100: cancer mortality is added to the cancer-free mortality (`AllMort`) and pancreatic cancer mortality is then subtracted to get the mortality rate that includes all cancers except pancreatic cancer (`ExcludePancMort`).
    * Lists append values at each loop, which are then used to create and save a data frame.
* Output: data frame with gender, cohort age, comorbidity, age, all-cause mortality including cancer, and all-cause mortality including all cancers except pancreatic cancer. 