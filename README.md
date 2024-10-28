# Kuskokwim_Chinook Run Reconstruction
This repository contains Kuskokwim River Chinook salmon run reconstruction and SR analyses. 

## Folder Structure
Folder is created by year that contains several subfolders: Data, Model, and Outputs.  **Note  numbers and names of folders and sub folders may differ among years to reflect changes in methodologies and data**
### Year folder:  Contains each year's data, model, and outputs
#### Data folder: contains data files used.
		* Kusko_Chinook_data_lookup.csv  : This file lists abbrabiations usde data file and definition
		* Kusko_Chinook_RR_Input_yyyy.csv : This file contains harvest, escapement data used for the model 
		* Kusko_Chinook_RR_Age_yyyy.csv  : This file contains Harvest and Escapement scale age data 
		
#### Model folder: contains models used.
  * R_functions: Contains r functions to process
    * ADMBtoR.R : Set of functions reading and writing list data to text data 
    * Create_age_data.R: Read Age data 
	
	* Outputs folder: contains model results 

## How to update the model 
* Step 1: Copy and paste ENTIRE previous year's folder and rename to current year 
* Step 2: Update Data *.csv file 
* Outputs folder: contains model results 

## How to update the repository
1. Copy the repository via RStudio
2. **Pull**  to Update your project folder to mach the repository 
3. **Commit**  After updating your files and folders 
4. **Push**  to update the repository

## How to run the model 

## TMB (Template Model Builder)
Required packages: Rtool, TMB, RTMB, Reshape2.  
Rtool, TMB, RTMB:  Needed to run TMB model
Reshape2:  Used to transpose (long to wide, wide to long) data.

