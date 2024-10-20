ReadMe - Living with Uncertainty: Using Multi-Model Large Ensembles to Assess Emperor Penguin Extinction Risk for the IUCN Red List2024

All the code used to produce data for this project can be found in the Codes\_EP folder.

Folder **Codes\_EP** contains:

- **data\_marika** : sea ice data
- **Codes\_CMR** : Codes to run CMR model
- **Codes\_IPM** : Codes to run IPM model
- **Codes\_SAT** : Codes to run SAT model
- **Metapop\_allmodel.m** : Program to run the metapopulation model  
- **Figures** : Codes to make figures of the results

# **data\_marika**

Contains the sea ice and temperature data (nc files) given by Marika Holland and the program to standardize them and transform them into mat files.

- Climate data: Sea ice concentration data and temperature data for each climate model. 500 or 850 refers to the surface covered by the data around each colony (500km2 and 850km2).

  For example: CESM1-CAM5\_colony\_sic\_poly500.nc file for CESM1 model for 500km2 around each colony, tas\_regavh\_NH\_CESM1-CAM5.nc for temperature data.

- **Obtain\_sic\_data.m** to standardize and create the mat files (mat files in /Codes\_CMR/mat\_data/)

# **Codes\_CMR**

Contains the code for the CMR model.

- **Main\_proj\_Stef\_all\_model\_wvec.m** 

Code for the CMR model, computes one txt file for each climate ensemble. The code works for all climate models that we had.

There is also the code to project the mean of ensembles, and the code that starts from 2009.

- Folder **mat\_data** with env data formatted (obtained with /Codes\_EP/data\_marika/obtain\_sic\_data.m)
  - Folder **Not\_std** contains the data before standardization (modification of variance and mean to fit observations)
  - Other files are standardized. One mat file per climate model, contains ENS, format: ENS(ens).SIC(year, season, colony) or SICa
  - **Info\_models.mat** contains year start and end, number of ensembles for each climate model


- Folder **DataOBS\_SAT\_BILGECAN** with satellite data from Bilgecan 
- folder **SIC\_obs\_allcol** contains one csv file per colony (66), SIC observations (not anomalies) for each season between 1979-2021.


- **Explore\_ext\_event.m** 

File to adapt extreme event scenario (the one where extreme events increase with SIC declines) to observations. Find the proportion of extreme events that lead to a breeding failure. Save the proportion in proportion\_ext\_event.mat

GENERATION LENGTH

- **Proj\_pop\_glb\_nomov2.m**

To project the population without movement without density dependance on 300 years with SIC observations for each year to obtain stable population (**Astable\_SICobs.mat**, stable matrix) to apply generation time formula.

- **generation\_length.m** 

to calculate generation time with Bienvenu and legendre method, uses the file Astable\_SICobs.mat (the other Astable are for projected SIC)

- **empe\_sitesNewNB.csv** contains information on colonies

# **Codes\_IPM**

ENVIRONMENTAL DATA FOR 66 COLONIES

- We have one csv file for each colony in IPM\_R/All\_colonies/env\_data\_all\_col.csv

ex: SST\_forecast\_1920-2100\_col1.csv

The files SST\_forecast\_col1.csv are for the period 1979-2020 (exactly the same data as in files from 1920 to 2100).


IPM MODEL

- Final model in **main\_proj\_IPM\_final.m** 

Model runs with environmental data from the CESM-LE2, SSP370 (SST, wind and vwind). We keep NOW anomalies at zero because we don’t have forecasts.

Model runs for yearstart = 1921 or yearstart = 2009 to 2100.

Save one txt file per colony and ensemble, containing nsimulation (100) projections on nt (nt, nsim). Files are saved in Codes\_EP/GR\_results\_tot/GR\_results\_IPM/CESM2 or Codes\_EP/GR\_results\_tot/2009/ GR\_results\_IPM/CESM2

Section “Mean of Ensembles” to study internal variability, calculates mean of 50 ensembles and runs the model, saves files in Codes\_EP/GR\_results\_tot/2009/mean\_ensemble/GR\_results\_IPM/CESM2

# **Codes\_SAT**

- **globalSAT\_new\_july.m** principal file to produce results, also code to project mean of ensembles
- **params\_chains.csv** and **eps\_chains.csv** : posteriors given by Bilgecan Sen
- **SIC\_NOW\_2009-2018.csv** file with observations for each colony (as a decadal mean)
- **Obtain\_SIC\_proj\_data.m** file to format the climate data
- **SIC\_proj\_trans\_total.mat** are the standardized data (with bilgecan method to have same mean and variance as observations)


# **Metapop\_allmodel.m** 
is here to apply the metapopulation model to both CMR and IPM model, for each climate model. Also uses dispersion.m



# **Figures**

- **env\_data** : to produce figures to visualize environmental data (SIC, wind, vwind, SST)

- **Fig\_CMR\_climate\_model**: Figure with climate models for the CMR model, with extreme event scenarios.

- **Fig\_CMR\_climate\_scenario**: Figure with climate scenarios for CMR model and the internal variability

- **Fig\_ecomod\_mature**: Figure with three ecological models, as well as figure with eco-ensemble (three models gathered, percent change)

- **Fig\_Regions**: Plot projections of different regions (CCMLAR regions and genetic regions)

- **Results\_Pext\_critE**: calculate probability of extinction under criteria E for the IUCN: proba to go under a threshold of 10 or 100 mature individuals in 2073.

- **Table\_Pext\_colonies**: calculate probability of extinction of colonies under criteria A of the IUCN: calculate proba to go under three thresholds corresponding to IUCN statuses (Vulnerable, Endangered, Critically endangered)

- **Table\_Pext\_poptot\_model**: calculate probability of extinction under criteria A of the IUCN: calculate proba to go under three thresholds corresponding to IUCN statuses (Vulnerable, Endangered, Critically endangered). Probabilities for total population for ecological models, climate models, climate scenarios, extreme events, internal variability.

- **Table\_Pext\_regions**: calculate probability of extinction of regions under criteria A of the IUCN: calculate proba to go under three thresholds corresponding to IUCN statuses (Vulnerable, Endangered, Critically endangered)

- **Maps\_python** contains the code you can use on python (**Map\_colonies\_status.py**) and maps for 2073 and 2100, with SIC calculated as a mean of 10. **Data\_for\_AliceP** are files that Kristen gave me to plot the ice.






