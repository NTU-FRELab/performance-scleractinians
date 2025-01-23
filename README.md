# Data from: Energy acquisition and allocation strategies in scleractinian corals: insights from intraspecific trait variability

This README.txt file was generated on 2025-01-22 by Stéphane De Palmas [depalmasstephane@gmail.com](mailto: depalmasstephane@gmail.com). Knit DXXXXX.Rmd to generate html and extract R code. (Latest update on 2025/01/23)

1. **Author Information**
   *   First author
        *   Name: De Palmas, Stéphane 
        *   Institution 1: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
	*   Institution 2: Faculty of Science, University of the Ryukyus, Okinawa 903-0213, Japan
   *   Second author
	*   Name: Chen, Qi
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
   *   Third author
	*   Name: Guerbet, Arnaud
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
   *   Fourth author
	*   Name: Hsieh, Yunli Eric
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
   *   Fifth author
	*   Name: Hsu, Tsai-Hsuan Tony 
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
   *   Sixth author
	*   Name: Lin, Yuting Vicky 
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
   *   Seventh author
	*   Name: Sturaro, Nicolas
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
   *   Eighth author
	*   Name: Wang, Pei-Ling
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
   *   Nineth & corresponding author
	*   Name: Denis, Vianney
	*   Institution: Institute of Oceanography, National Taiwan University, Taipei 10617, Taiwan
	*   Email: [vianneydenis@ntu.edu.tw](mailto:vianneydenis@ntu.edu.tw)
2. **Date of data collection**: 2017
3. **Geographic location of data collection**: Taiwan, West Pacific
4. **Funding sources that supported the collection of the data**: National Science and Technology Council (NSTC), National Taiwan University
5. **Recommended citation for this dataset**
	De Palmas S, Chen Q, Guerbet A, Hsieh YE, Hsu THT, Lin YV, Sturaro N, Wang PL, Denis V (2025), Data from: Energy acquisition and allocation strategies in scleractinian corals: insights from intraspecific trait variability, Dryad, Dataset, [https://doi.org/10.5061/dryad.xwdbrv1hw](https://doi.org/10.5061/dryad.xwdbrv1hw)



## DESCRIPTION OF THE DATA AND FILE STRUCTURE

### DATA & FILE OVERVIEW

1. **Description of dataset**
    These data were generated to investigate scleractinian coral strategies in energy acquisition and allocation in the light of intraspecific variation. 
2. **File List**

   File 1: DePalmas_dataset_strategy_data.csv: Excel file containing coral trait information for 7 species and environmental information prior to collection.
     
### METHODOLOGICAL INFORMATION

A detailed description of data acquisition and processing can be found in the published manuscript in Coral Reefs (https://doi.org/10.1007/s00338-025-02617-w).

#### DATA-SPECIFIC INFORMATION
##### **DePalmas_dataset_strategy_data.csv**

1. Number of variables/columns: 75
2. Number of cases/rows: 197
3. Missing data codes: NA
4. Variable List
    column A - codeP 
    column B - date
    column C - spe
    column D - sp_check    
    column E - region
    column F - site
    column G - depth
    column H - season 
    column I - P
    column J - I
    column K - Volume
    column L - Weight
    column M - Surface
    column N - Vi
    column O - V.bio
    column P - W.fil
    column Q - W.fil.tis
    column R - W.fil.tis.bur
    columns S to X - countS.1, countS.2, countS.3, countS.4, countS.5, countS.6
    columns Y to AD - countC.1, countC.2, countC.3, countC.4, countC.5, countC.6
    column AE - V.chl.slurry
    column AF - V.extraction.T
    column AG - V.extraction.S
    columns AH to AK - A750T, A664T, A647T, A630T
    columns AL to AO - A750S, A664S, A647S, A630S
    column AP - SDS1V
    column AQ - V.prot.slurry
    column AR - STC
    column AS to AW - STCX4, STCX3, STCX2, STCX1, STCE
    column AX - abs.prot.blank
    column AY to BA - abs.prot.zoox1, abs.prot.zoox2, abs.prot.zoox3
    column BB to BD - abs.prot.host1, abs.prot.host2, abs.prot.host3
    column BE - STC-2(for host)
    column BF to BJ - STCX4b, STCX3b, STCX2b, STCX1b, STCEb    
    column BK - abs.prot.blank.b: absorption reading on the blank used for correction
    column BL to BO - Hdelta15N, Hdelta13C, Zdelta15N, Zdelta13C
    column BP to BS - H_N%, H_C%, Z_N%, Z_C% 
    column BT to BW - SST_mean, SST_sd, light_mean, light_sd 
    
5. Descriptions
The descriptions of each variables are as below:
    codeP: sample identification code 
    date: sampling date
    spe: species name
    sp_check: species identity checked under dissecting microscope (0: no, 1: yes)
    region: region where specimens were collected (NT: North Taiwan, GI: Green Island)
    site: name of the locality where specimens were collected
    depth: depth at which specimens were collected
    season: season when the specimens were collected 
    P: physiological investigations (0: no, 1: yes)
    I: isotopic investigations (0: no, 1: yes)
    Volume: volume of the skeletal fragment in cubic millimeter
    Weight: weight of the skeletal fragment in grams
    Surface: surface of the skeletal fragment in square millimeter
    Vi: Initial volume of slurry
    V.bio: volume from Vi used for biomass measurements
    W.fil: weight of the filter used for biomass calculation
    W.fil.tis: weight of the filter with addition of the slurry
    W.fil.tis.bur: weight of the burned filter after addition of the slurry
    countS.1, countS.2, countS.3, countS.4, countS.5, countS.6: 6 replicates of Symbiodiniaceae counting 
    countC.1, countC.2, countC.3, countC.4, countC.5, countC.6: 6 replicates of cnidocytes counting
    V.chl.slurry: volume from Vi used for chlorophyll measurements
    V.extraction.T: final volume of chlorophyll extraction for tissue
    V.extraction.S: final volume of chlorophyll extraction for skeleton
    A750T, A664T, A647T, A630T: wavelengths reading on spectrophotometer for tissue 
    A750S, A664S, A647S, A630S: wavelengths reading on spectrophotometer for skeleton
    SDS1V: volume from the Symbiodiniaceae pellet resuspended and used for protein measurements
    V.prot.slurry: volume from Vi used to isolate animal tissue for protein measurements
    STC: date of preparation of the plate for protein measurements
    STCX4, STCX3, STCX2, STCX1, STCE: Polynomial coefficients for assessing protein concentration 
    abs.prot.blank: absorption reading on the blank used for correction
    abs.prot.zoox1, abs.prot.zoox2, abs.prot.zoox3: absorbance reading on Symbiodiniaceae protein triplicates
    abs.prot.host1, abs.prot.host2, abs.prot.host3: absorbance reading on host protein triplicates 
    STC-2(for host): date of preparation of the plate for protein measurements (when Symbiodiniaceae and host protein were not measured simultaneously)
    STCX4b, STCX3b, STCX2b, STCX1b, STCEb: Polynomial coefficients for assessing protein concentration
    abs.prot.blank.b: absorption reading on the blank used for correction
    Hdelta15N, Hdelta13C, Zdelta15N, Zdelta13C: carbon and nitrogen isotopic ratios for host and zooxanthellae (H stands for host and Z stands for zooxanthellae)
    H_N%, H_C%, Z_N%, Z_C% : Nitrogen and carbon mass for host and Symbiodiniaceae used for C:N ratio calculation
    SST_mean, SST_sd, light_mean, light_sd: Sea surface temperature mean and standard deviation, light mean and standard deviation. All environmental parameters were calculated based on one month of data prior to the day of collection. 
    
6. **Other relevant information**:
    Host protein was measured with a different curve and blank (host and Symbiodiniaceae on different plates) for 13 specimens. Columns BD to BK have information for these 13 specimens.
