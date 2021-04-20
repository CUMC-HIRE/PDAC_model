# This script generates comorbidity-specific life tables that include all cancer mortality rates 
# except pancreatic cancer
# Requirements: Excel spreadsheet with the following sheets:
  # CancerFreeLifeTables: Each column represents a different gender, age, and comorbidity specific 
    # cancer-free annual mortality rate (from Lansdorp-Vogelaar et al, 2014)
  # CancerMort_Ages66-84: mortality rates for all malignant neoplasms for single ages 66-84
  # CancerMort_Ages85plus: mortality rates for all malignant neoplasms for 85+ years.
  # Panc_CauseOfDeathCodes: ICD-10 codes for pancreatic cancer

library("readxl")
library("data.table")
setwd("C:/Users/bnl2108/OneDrive for Business/Pancreatic cancer/model/data/morb_tables")
filename_lt = "20200219_comorb_life_table.xlsx"
filename_out = "panc_comorb_lt.csv"

# Get data for cancer-free survival by comorbidity level
CancerFreeLT <- read_excel(filename_lt, sheet = "CancerFreeLifeTables")
CancerFreeLT <- within(CancerFreeLT, rm("...1"))
# Get data for all cancer mortality
DT_CancerMort66_84 <- data.table(read_excel(filename_lt, sheet = "CancerMort_Ages66-84"))
DT_CancerMort85plus <- data.table(read_excel(filename_lt, sheet = "CancerMort_Ages85plus"))
# Calculate cancer motality rate by age and gender
DT_CancerMortRate66_84 <- DT_CancerMort66_84[, .(MortRate = sum(Crude.Rate)), by = .(Age = Single.Year.Ages.Code, Gender)]
DT_CancerMortRate85plus <- DT_CancerMort85plus[, .(MortRate = sum(Crude.Rate)), by = Gender]
Panc_Codes <- read_excel(filename_lt, sheet = "Panc_CauseOfDeathCodes", col_types = "text")
Panc_Codes <- Panc_Codes$'ICD 10 Cause of Death Codes'
# Calculate pancreatic cancer mortality rate by age and gender
DT_PancMortRate66_84 <- DT_CancerMort66_84[Cause.of.death.Code %in% Panc_Codes, .(MortRate = sum(Crude.Rate)), by = .(Age = Single.Year.Ages.Code, Gender)]
DT_PancMortRate85plus <- DT_CancerMort85plus[Cause.of.death.Code %in% Panc_Codes, .(MortRate = sum(Crude.Rate)), by = Gender]

l_Gender = c()
l_CohortAge = c()
l_Comorb = c()
l_Age = c()
l_AllMort = c()
l_ExcludePancMort = c()
for (i in names(CancerFreeLT)) {
  # Get data from column name
  CohortAge = as.numeric(substr(i, start = 1, stop = 2))
  age = CohortAge
  # Gender
  if (grepl("Women", i)) {gender = "Female"} 
  else if (grepl("Men", i)) {gender = "Male"}
  else {print("ERROR: invalid gender label")}
  # Comorbidity
  if (grepl("No", i)) {Comorb = "No"}
  else if (grepl("Low", i)) {Comorb = "Low"}
  else if (grepl("Moderate", i)) {Comorb = "Moderate"}
  else if (grepl("Severe", i)) {Comorb = "Severe"}
  else {print("ERROR: invalid comorbidity label")}
  
  while (age <= 100) {
    # Add to lists  
    l_CohortAge = c(l_CohortAge, CohortAge)
    l_Age = c(l_Age, age)
    l_Gender = c(l_Gender, gender)
    l_Comorb = c(l_Comorb, Comorb)
    
    # Calculate mortality rates
    i_age = age - CohortAge + 1
    CancerFreeMort = CancerFreeLT[[i]][i_age]
    if (age < 85) {
      CancerMort = DT_CancerMortRate66_84[Age == age & Gender == gender, MortRate]
      PancMort = DT_PancMortRate66_84[Age == age & Gender == gender, MortRate]
    } else {
      CancerMort = DT_CancerMortRate85plus[Gender == gender, MortRate]
      PancMort = DT_PancMortRate85plus[Gender == gender, MortRate]
    }
    AllMort = CancerFreeMort + CancerMort
    ExcludePancMort = AllMort - PancMort
    l_AllMort = c(l_AllMort, AllMort)
    l_ExcludePancMort = c(l_ExcludePancMort, ExcludePancMort)
    age = age + 1
  }
}
# Create dataframe with lists and export
DF_CancerLT = data.frame(Gender = l_Gender,
                         CohortAge = l_CohortAge,
                         Comorbidity = l_Comorb,
                         Age = l_Age,
                         AllCauseWithAllCancersIncluded = l_AllMort,
                         AllCauseExcludingPancCancers = l_ExcludePancMort)
write.csv(DF_CancerLT, filename_out)
