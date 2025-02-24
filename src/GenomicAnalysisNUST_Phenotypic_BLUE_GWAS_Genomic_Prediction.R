################################################################################
###
### This script will run statistical analyses of the phenotypic data from the
### Northern Uniform Soybean Tests (1993-2020) in three steps:
### 1. Phenotypic analysis Best Linear Unbiased Estimates (BLUE) via ASReml-R
### 2. Genome-wide Association Studies via rrBLUP::GWAS
### 3. Genomic Prediction via rrBLUP::mixed.solve
###
### Author: Cleiton Wartha
### > sessionInfo()
### R version 4.3.3 (2024-02-29 ucrt)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows 10 x64 (build 19045)
### 
### Matrix products: default
### locale:
### [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
### [5] LC_TIME=English_United States.utf8    
### time zone: America/Chicago
### tzcode source: internal
### attached base packages:
###   [1] stats     graphics  grDevices utils     datasets  methods   base     
### other attached packages:
### boot_1.3-29       vcfR_1.15.0       rrBLUP_4.6.3      asreml_4.2.0.257  Matrix_1.6-5      lattice_0.22-5    data.table_1.15.4 lubridate_1.9.3   forcats_1.0.0    
### stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2       readr_2.1.4       tidyr_1.3.0       tibble_3.2.1      ggplot2_3.5.1     tidyverse_2.0.0 
################################################################################

#Install packages and load libraries
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('rrBLUP')) install.packages('rrBLUP'); library(rrBLUP)
if (!require('asreml')) install.packages('asreml'); library(asreml)
if (!require('lattice')) install.packages('lattice'); library(lattice)
if (!require('data.table')) install.packages('data.table'); library(data.table)
if (!require('vcfR')) install.packages('vcfR'); library(vcfR)
if (!require('boot')) install.packages('boot'); library(boot)
options(scipen = 999) #prevent the triggering of scientific notation for integers

# Load in raw phenotypic data and respective metadata
NUST_raw_tall <- read_csv("data/Phenotype_Measures_Final_Master_NUST_1993_2020_years28.csv.gz") #raw phenotypic data in analyzable format
metadata <- read_csv("data/NUST_Metadata_1989_2020.csv", col_types = cols()) #metadata contains CV, LSD, Reps, Rows/Plot, and RowSp

# Create function to hard code renaming of experimental strains with cultivar name
clean_germplasm_id <- function(id) {
  id %>%
    str_replace_all("ORC_", "ORC") %>%
    str_replace_all("OAC_", "OAC") %>%
    str_replace_all("OT99-5", "ACDundas") %>%
    str_replace_all("OT99-4", "ACRodeo") %>%
    str_replace_all("HC94-1065", "Apex") %>%
    str_replace_all("E93001", "Apollo") %>%
    str_replace_all("ND02-2367", "Ashtabula") %>%
    str_replace_all("C1875", "Athow") %>%
    str_replace_all("ORC9004", "Ayr") %>%
    str_replace_all("ND95-931", "Barnes") %>%
    str_replace_all("ORC9205", "Blackjack21") %>%
    str_replace_all("C1804", "Bronson") %>%
    str_replace_all("Ky85-09073", "Calhoun") %>%
    str_replace_all("Ky88-5037", "CF461") %>%
    str_replace_all("HC85-6724", "Charleston") %>%
    str_replace_all("Charleston", "Charleston") %>%
    str_replace_all("LN90-4129", "Cisne") %>%
    str_replace_all("AConradBC", "Conrad94") %>%
    str_replace_all("HC89-2232", "Croton3.9") %>%
    str_replace_all("ND91-2735", "Daksoy") %>%
    str_replace_all("ND91-2330", "Danatto") %>%
    str_replace_all("HS93-3779", "Darby") %>%
    str_replace_all("SD02-22", "Davison") %>%
    str_replace_all("HS91-4523", "Defiance") %>%
    str_replace_all("HS1-3717", "Dennison") %>%
    str_replace_all("SD02-833", "Deuel") %>%
    str_replace_all("OAC05-30", "DH530") %>%
    str_replace_all("OAC06-14", "DH614") %>%
    str_replace_all("LN92-10507", "Dwight") %>%
    str_replace_all("M86-1322", "Faribault") %>%
    str_replace_all("U91-3607", "Fillmore") %>%
    str_replace_all("HF91-078", "Flint") %>%
    str_replace_all("M87-1621", "Freeborn") %>%
    str_replace_all("M87-731", "Glacier") %>%
    str_replace_all("M87-642", "Granite") %>%
    str_replace_all("ND09-5604", "Henson") %>%
    str_replace_all("A92-525014", "IA1006") %>%
    str_replace_all("A95-483010", "IA1008") %>%
    str_replace_all("A95-583021", "IA1009") %>%
    str_replace_all("A04-442030", "IA1019") %>%
    str_replace_all("A02-136032", "IA1021") %>%
    str_replace_all("A04-543037", "IA1022") %>%
    str_replace_all("A09-754003", "IA1026") %>%
    str_replace_all("IA2007BC", "IA2007R") %>%
    str_replace_all("IA2008BC", "IA2008R") %>%
    str_replace_all("A91-607052", "IA2021") %>%
    str_replace_all("A91-607024", "IA2022") %>%
    str_replace_all("A93-553034", "IA2036") %>%
    str_replace_all("A94-674017", "IA2038") %>%
    str_replace_all("A94-673009", "IA2039") %>%
    str_replace_all("A96-494018", "IA2050") %>%
    str_replace_all("A96-591033", "IA2051") %>%
    str_replace_all("A96-591046", "IA2052") %>%
    str_replace_all("A00-711013", "IA2068") %>%
    str_replace_all("A03-946036", "IA2072") %>%
    str_replace_all("A05-114010", "IA2077") %>%
    str_replace_all("A05-114019", "IA2078") %>%
    str_replace_all("A05-213013", "IA2079") %>%
    str_replace_all("A04-543012", "IA2093") %>%
    str_replace_all("A04-545020", "IA2094") %>%
    str_replace_all("A08-248043", "IA2102") %>%
    str_replace_all("A07-427027", "IA2105") %>%
    str_replace_all("A10-556015", "IA2109") %>%
    str_replace_all("AM90-211003", "IA3003") %>%
    str_replace_all("A91-701007", "IA3004") %>%
    str_replace_all("A92-726034", "IA3005") %>%
    str_replace_all("A94-774021", "IA3010") %>%
    str_replace_all("A97-973002", "IA3014") %>%
    str_replace_all("A99-315026", "IA3023") %>%
    str_replace_all("A03-946005", "IA3024") %>%
    str_replace_all("A03-946006", "IA3025") %>%
    str_replace_all("A04-642023", "IA3028") %>%
    str_replace_all("A06-911034", "IA3048") %>%
    str_replace_all("A07-626002", "IA3052") %>%
    str_replace_all("A05-312025", "IA4004") %>%
    str_replace_all("LN94-10527", "Ina") %>%
    str_replace_all("LN88-10534", "Iroquois") %>%
    str_replace_all("AKenwoodBC", "Kenwood94") %>%
    str_replace_all("K1191", "KS4694") %>%
    str_replace_all("ND98-2252", "LaMoure") %>%
    str_replace_all("LD00-2817P", "LD00-2817") %>%
    str_replace_all("LD00-2817W", "LD00-2817") %>%
    str_replace_all("LN89-295", "Macon") %>%
    str_replace_all("K1207", "Magellan") %>%
    str_replace_all("AMarcusBC", "Marcus95") %>%
    str_replace_all("M98-101137", "MN0101") %>%
    str_replace_all("ND(M)89-111", "MN0301") %>%
    str_replace_all("M92-185003", "MN0304") %>%
    str_replace_all("M06-289001", "MN0310CN") %>%
    str_replace_all("M91-821", "MN0901") %>%
    str_replace_all("M95-265009", "MN1005") %>%
    str_replace_all("M95-258296", "MN1006CN") %>%
    str_replace_all("M89-936", "MN1301") %>%
    str_replace_all("M91-1137", "MN1801") %>%
    str_replace_all("Md96-5722", "Monocacy") %>%
    str_replace_all("LS87-1615", "Mustang") %>%
    str_replace_all("ND12-21598", "ND17009GT") %>%
    str_replace_all("ND99-2621C", "ND99_2621") %>%
    str_replace_all("ND10-3464", "NDBenson") %>%
    str_replace_all("ND10-3067", "NDStutsman") %>%
    str_replace_all("U95-2418", "NE1900") %>%
    str_replace_all("U93-2412", "NE3297") %>%
    str_replace_all("U94-2306", "NE3399") %>%
    str_replace_all("U95-3231", "NE3400") %>%
    str_replace_all("OAC97-06", "OACStratford") %>%
    str_replace_all("OAC06-32", "OACDrayton") %>%
    str_replace_all("U91-3610", "Odell") %>%
    str_replace_all("HS90-3513", "Ohio_FG2") %>%
    str_replace_all("HF03-534", "OHS202") %>%
    str_replace_all("HS3-2840", "OHS303") %>%
    str_replace_all("HS4-3189", "OHS305") %>%
    str_replace_all("E91031", "Olympus") %>%
    str_replace_all("LN92-10855", "Pana") %>%
    str_replace_all("LN86-3357", "Piatt") %>%
    str_replace_all("ORC9901", "PRO28-53") %>%
    str_replace_all("C1832", "Probst") %>%
    str_replace_all("ORC9008", "RCAT_Columbus") %>%
    str_replace_all("ORC2002", "RCATLogan") %>%
    str_replace_all("ORC9902", "RCATWildcat") %>%
    str_replace_all("LN92-10725", "Rend") %>%
    str_replace_all("ND95-931CRR", "RG405RR") %>%
    str_replace_all("LS05-3229", "Saluki4411") %>%
    str_replace_all("HS88-4908", "Sandusky") %>%
    str_replace_all("ND96-1593", "Sargent") %>%
    str_replace_all("LN90-4187", "Savoy") %>%
    str_replace_all("ND01-3906", "Sheyenne") %>%
    str_replace_all("E98076", "Skylla") %>%
    str_replace_all("L83-3804", "Spry") %>%
    str_replace_all("HC94-421", "Stout") %>%
    str_replace_all("HS3-2324", "Streeter") %>%
    str_replace_all("HC89-2170", "Stressland") %>%
    str_replace_all("SD(M)92-1357", "Stride") %>%
    str_replace_all("HC94-422", "Sturdie") %>%
    str_replace_all("HS5-3417", "Summit") %>%
    str_replace_all("SD(M)92-1233", "Surge") %>%
    str_replace_all("HF93-082", "Tiffin") %>%
    str_replace_all("E93147", "Titan") %>%
    str_replace_all("ND90-2624", "Traill") %>%
    str_replace_all("HC90-196", "Troll") %>%
    str_replace_all("SD93-522", "Turner") %>%
    str_replace_all("HS88-4905", "Vertex") %>%
    str_replace_all("ND96-8929", "Walsh") %>%
    str_replace_all("ORC9404", "Westag97") %>%
    str_replace_all("HF02-0005", "Wyandot") %>%
    str_replace_all("LN87-2112", "Yale") %>%
    str_replace_all("OAC14-05C", "OAC14-05") %>%
    str_replace_all("-", "_") %>%
    str_replace_all("\\(", "_") %>%
    str_replace_all("\\)", "_") %>%
    str_replace_all("ORC ", "ORC") %>%
    str_replace_all("OAC ", "OAC") %>%
    str_replace_all("\\.", "_") %>%
    str_replace_all(" ", "")
}

NUST_raw_tall <- NUST_raw_tall %>%
  mutate(GermplasmId = clean_germplasm_id(GermplasmId)) %>% #use hard-coded function above
  mutate(Year = as.numeric(str_extract(Experiment, "\\d+$")),
         Concatenate = paste(Year, Experiment, Location, GermplasmId, sep = "")) %>%
  select(Concatenate, Year, everything(), -Entry, -Field, -Rep, -Plot, -Range, -Row) #remove unnecessary columns

# Remove Yield data from locations with CV greater than 15 (549 cases)
meta_CV_15_pls <- metadata %>%
  mutate(Value = as.numeric(Value)) %>%
  filter(Meta_data_type == "C.V." & Value >= 15)

#Loop through all the locations with CV > 15 to be NA with base r
i = 1
for(i in 1:nrow(meta_CV_15_pls)){
  
  NUST_raw_tall[ which(   NUST_raw_tall$Experiment == meta_CV_15_pls$Test_Year[i] 
                         & NUST_raw_tall$Location   == meta_CV_15_pls$Location2[i] 
                         & NUST_raw_tall$Phenotype  == "YieldBuA"), ]$Value  <-  rep(NA, length(NUST_raw_tall[ which(   NUST_raw_tall$Experiment == meta_CV_15_pls$Test_Year[i] 
                                                                                                                      & NUST_raw_tall$Location   == meta_CV_15_pls$Location2[i] 
                                                                                                                      & NUST_raw_tall$Phenotype  == "YieldBuA"), ]$Value ))
  
  NUST_raw_tall[ which(   NUST_raw_tall$Experiment == meta_CV_15_pls$Test_Year[i] 
                         & NUST_raw_tall$Location   == meta_CV_15_pls$Location2[i] 
                         & NUST_raw_tall$Phenotype  == "YieldRank"), ]$Value  <-  rep(NA, length(NUST_raw_tall[ which(   NUST_raw_tall$Experiment == meta_CV_15_pls$Test_Year[i] 
                                                                                                                       & NUST_raw_tall$Location   == meta_CV_15_pls$Location2[i] 
                                                                                                                       & NUST_raw_tall$Phenotype  == "YieldRank"), ]$Value ))
}

# Modify Experiment column to be a MG column
NUST_raw_tall <- NUST_raw_tall %>%
  drop_na(Value) %>% #drop the NA values with CV > 15
  mutate(Experiment = str_replace_all(Experiment, c("_.*" = "", "PT" = "MG", "UT" = "MG", "A" = "", "B" = "", "R" = "")),
         MG = case_when(
           Experiment == "MG00" ~ -1,
           Experiment == "MG0" ~ 0,
           Experiment == "MGIV" ~ 4,
           Experiment == "MGIII" ~ 3,
           Experiment == "MGII" ~ 2,
           Experiment == "MGI" ~ 1,
           TRUE ~ NA_real_),
         GermplasmId = paste(GermplasmId, Experiment, sep = "_")) %>%
  select(Concatenate, Year, Experiment, MG, everything())

# Convert data from tall to wide format
NUST_data <- NUST_raw_tall %>%
  arrange(Phenotype) %>% #organize traits in alphabetical order
  pivot_wider(names_from = Phenotype, values_from = Value)

# Remove phenotypes that will not be used in downstream analysis
NUST_data <- NUST_data %>%
  select(-c(BSR, BSR.Incid, BSR.Sev, BSRPlant, BSRStem, BTSa, Emergence, FE, FELS,
  HardSeed, Mottle, P.SB, PhytoRot, PSB, RootRot, RootRotRace25, SDSI., SDSDI,
  SDSDS, SDSR6,  SDSR6Date, SDSRank, SDSRDate, SDSS, SDSTest, SeedGerm,Stand, StemCanker,
  SMV, Frogeye, LINOLEIC, LINOLENIC, OLEIC, PALMITIC, STEARIC, DescriptiveCode, SCL,
  GreenStem, SDS, YieldRank)) %>% #remove these last three traits because they were not used in the GWAS/GP analysis
  mutate(across(c(Chlorosis, Height, Lodging, Maturity, Oil, Protein, SeedQuality, SeedSize, Shattering, YieldBuA), as.numeric)) #convert chr to num

# Calculate number of locations at which each line tested
test_table_year <- NUST_data %>%
  count(GermplasmId, Year)

count_locations <- test_table_year %>%
  group_by(GermplasmId) %>%
  summarise(Location_Count = sum(n))
# Plot histogram with base R - median = 13 locations
hist(count_locations$Location_Count,freq = TRUE, xlim = c(0,100),ylim = c(0,1000),
    breaks = 1000, col = "green", ylab = "Count of experimental strains", xlab = "Number of Environments", 
    main = "Number of environments at which experimental strains are tested")

################################################################################
### Phenotypic Analysis of Best Linear Unbiased Estimates across locations
### Analysis enabled via proprietary software ASReml-R
### Maturity was use as a covariate for most of the complex traits
###

blue.long <-c() #Create empty object to store results from loop

# Subset data on trait while running through each trait, remove NA rows
traits <- c("Chlorosis", "Height", "Lodging", "Maturity", "Oil", "Protein", "SeedQuality", "SeedSize", "Shattering", "YieldBuA")

###Loop over each trait
startTime <- Sys.time()
for (trait in traits) {
  NUST_trait <- trait

  #First prepare the input data
  if (NUST_trait %in% c("Chlorosis", "Maturity", "Shattering")) {
    pheno.temp <- NUST_data %>%
      select(Concatenate, Year, Experiment, MG, Location, GermplasmId, !!sym(NUST_trait)) %>%
      drop_na(!!sym(NUST_trait))
 
  } else {
    pheno.temp <- NUST_data %>%
      select(Concatenate, Year, Experiment, MG, Location, GermplasmId, !!sym(NUST_trait), Maturity) %>%
      drop_na(!!sym(NUST_trait), Maturity)
  }
  
  pheno.temp <- pheno.temp %>%
    mutate(across(c(Year, Location, GermplasmId, MG), as.factor))
  
  # Run ASRemlR Model
  #Simpler model for the three traits Chlorosis, Maturity and Shattering
  if (NUST_trait %in% c("Chlorosis", "Maturity", "Shattering")) {
    
    blue.fit <- asreml(fixed = get(NUST_trait) ~ Year + Location + Year * Location + GermplasmId,
                             data = pheno.temp,
                             workspace = "8gb")
  } else {
    pheno.temp <- pheno.temp %>%
      mutate(Maturity_cent = scale(Maturity, center = TRUE, scale = FALSE))
    #All remaining traits entail an adjustment using the Maturity date in Julian Date as a covariate to adjust for maturity differences
    blue.fit <- asreml(fixed = get(NUST_trait) ~ Year + Location + Year * Location + GermplasmId + Maturity_cent, #centralized Maturity as a covariate
                             data = pheno.temp,
                             workspace = "8gb")
  }
  # Get the BLUEs
  blue.asr <- as.data.frame(summary(blue.fit, coef = TRUE)$coef.fixed) # BLUEs
  int  <- blue.asr[rownames(blue.asr)%in% c("(Intercept)"), "solution"]  #retrieve intercept value ASReml-R 4.2 notation
  #Create tidy long format file to store the results from all traits
  blue.asr <- blue.asr %>% rownames_to_column()  %>% 
    filter(str_detect(rowname, "GermplasmId_")) %>%  #filter rownames with strains
    mutate(GermplasmId = str_replace(rowname, "GermplasmId_", "")) %>% 
    mutate(Value = solution + int) %>%
    rename(std.error= 'std error') %>%
    select(GermplasmId, Value, std.error) %>%
    mutate(Phenotype = NUST_trait) %>%
    relocate(c("Phenotype"), .before = Value)
  
  blue.long <- rbind(blue.long, blue.asr) #store results from each loop
  
  cat("\nDone with", NUST_trait,"\n") #print out completed trait for progress tracking
}
endTime <- Sys.time() #calculate end time of BLUE calculation
print(endTime - startTime) # print recorded Time difference of 1.676401 hours

#Convert BLUE file from long to wide format
BLUE.dat <- blue.long %>% select(-std.error) %>% #remove standard error not needed for downstream analysis
  pivot_wider(names_from = Phenotype, values_from = Value) %>%
  mutate(GermplasmId2 = GermplasmId) %>% #duplicate column to split the MG information
  separate(GermplasmId2, into = c("GermplasmId_name", "MG_factor"), sep = "_MG") %>%
  mutate(MG_value2 = paste0("MG", MG_factor),
         MG_value2 = case_when(
           MG_value2 == "MG00" ~ "-1",
           MG_value2 == "MG0" ~ "0",
           MG_value2 == "MGIV" ~ "4",
           MG_value2 == "MGIII" ~ "3",
           MG_value2 == "MGII" ~ "2",
           MG_value2 == "MGI" ~ "1",
           TRUE ~ MG_value2),
         MG = as.numeric(MG_value2),
         MG_factor= paste0("MG", MG_factor)) %>%
  select( -MG_value2) %>%
  relocate(c("MG","MG_factor", "GermplasmId_name"), .after = GermplasmId)

#save data frame with the output from the BLUE phenotypic analysis for downstream GWAS and Genomic Prediction
#write.csv(x= BLUE.dat, file = "results/BLUE_NUST_1993_2020_GermplasmId_Output_Table.csv", row.names = F)

###If you do not have a ASReml-R license, uncomment and read the BLUE output file below to continue with the analysis
#BLUE.dat <- read.csv("results/BLUE_NUST_1993_2020_GermplasmId_Output_Table.csv")

################################################################################
##### Formatting the compiled data on phenotypic descriptor traits          ####
#### Flower and pubescence color to be included in the GWAS                 ####
################################################################################
flwr.pub <- read.csv("data/NUST_Flower_Pubescence_Color_Input_Data.csv") #Load file

#Formatting from chr to numeric
flwr.pub <- flwr.pub %>%
  mutate(GermplasmId = clean_germplasm_id(GermplasmId), #Clean GermplasmId in flwr_pub using the same function defined above to 
  #convert experimental strain names into cultivar names and other minor formatting
         Flower = case_when(
           Flower == "P" ~ 1,
           Flower == "H" ~ NA,
           Flower == "NA" ~ NA,
           Flower == "W" ~ 0,
           TRUE ~ as.numeric(Flower)
         ),
         Pubescence = case_when(
           Pubescence == "T" ~ 1,
           Pubescence == "H" ~ NA,
           Pubescence == "NA" ~ NA,
           Pubescence == "G" ~ 0,
           TRUE ~ as.numeric(Pubescence)
         )) %>%
  distinct() %>%
  arrange(GermplasmId)
#------------------------------------------------------------------#
#    Formatting Flower Pub data and add with other Pheno data      #
#------------------------------------------------------------------#
pheno.dat <- BLUE.dat %>%
  left_join(., flwr.pub, by = c("GermplasmId_name"= "GermplasmId")) %>%
  arrange(GermplasmId)

#Rename column names function
rename_colnames <- function(name) {
  name %>%
    str_replace_all("SeedQuality", "Seed.Quality") %>%
    str_replace_all("SeedSize", "Seed.Size") %>%
    str_replace_all("GreenStem", "Green.Stem") %>%
    str_replace_all("YieldBuA", "Seed.Yield") %>%
    str_replace_all("YieldRank", "Yield.Rank")
}
#execute function
names(pheno.dat) <- rename_colnames(names(pheno.dat))

#------------------------------------------------------------------#
#    Formatting of the genotypic data                              #
#------------------------------------------------------------------#

######################################################################
###                         VCFToSimple                            ###
###   a function to generate simple numeric -101 format from vcf   ###
######################################################################
#function to convert vcf into simple numeric -101 format from Sushan Ru
# need vcfR library
# input: vcf file path and name
VCFToSimple <- function(infile){
  vcf <- read.vcfR(infile, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = F)
  fix <- as_tibble(getFIX(vcf))
  
  gt2 <- matrix(NA, ncol = ncol(gt), nrow = nrow(gt))
  colnames(gt2) <- colnames(gt)
  gt2[(gt == "0/0")|(gt == "0|0")] <- -1  #AA call
  gt2[(gt == "1/1")|(gt == "1|1")] <- 1   #BB call
  gt2[(gt == "0/1")|(gt == "1/0")|(gt == "0|1")|(gt == "1|0")] <- 0  #AB call
  gt2[(gt == "./.")|(gt == ".|.")] <- NA
  gt2 <- as_tibble(gt2) %>%
    mutate(SNPID = fix$ID)
  
  gt.simple <- fix %>%
    select(ID, CHROM, POS, REF, ALT) %>%
    rename(SNPID=ID, Chrom=CHROM, Position=POS) %>%
    left_join(gt2, by = 'SNPID')
  
  return(gt.simple)
}
#Define path of the input VCF file
## The marker data has been already filtered and imputed using TASSEL 5.0 (TASSEL - Trait Analysis by aSSociation, Evolution and Linkage)
#KNN imputation was used - there are still 15081 NA values
infile <- "data/NUST_geno_Wm82.a2_v1_Filt_KNNimp.vcf.gz" #original compressed file
#infile <- "data/NUST_Geno_Wm82.a2_v1_Raw.vcf" #Vishnu's EVA corrected
geno.imp <- VCFToSimple(infile) #execute function

#Note for rrBLUP: columns 1-3: SNPID, Chromosome, Position; Columns 4-end: genotype data by experimental strain -101
names(geno.imp) <- clean_germplasm_id(names(geno.imp)) #convert germplasmId names formatting from experimental strains to cultivar names

#Housekeeping - remove REF and ALT alleles not needed in downstream function
geno.imp <- geno.imp %>% select(-c(REF, ALT))

### Duplicate the marker data to experimental strains that were tested in more than one MG experiment
# Find and keep both duplicates
duplicates <- pheno.dat %>%
  filter(duplicated(GermplasmId_name) | duplicated(GermplasmId_name, fromLast = TRUE)) %>%
  arrange(GermplasmId_name, MG) %>%
  rename(GermplasmId_2 = GermplasmId)
# Remove all duplicate rows from main data set
pheno.dat.reformat <- pheno.dat %>%
  filter(!GermplasmId %in% duplicates$GermplasmId_2) %>%
  rename(GermplasmId_2 = GermplasmId_name) %>%
  relocate(GermplasmId_2, .before = MG) %>%
  relocate(GermplasmId, .after = MG_factor) %>%
  rename(GermplasmId_name = GermplasmId)

# Add reformatted duplicate line set back into full data set
pheno.dat.final <- bind_rows(pheno.dat.reformat, duplicates) %>%
  arrange(GermplasmId_2, MG) %>%
  rename(GermplasmId = GermplasmId_2) %>% #rename germplasmId column
  select(-GermplasmId_name) %>% #delete secondary GermplasmId_name column
  data.frame()
# Format genotypic data to include name_MG# for lines tested in >1 MG
temp.geno1 <- t(geno.imp)
rownames <- as.matrix(row.names(temp.geno1))
temp.geno2 <- cbind(rownames, temp.geno1)
colnames(temp.geno2) <- temp.geno2[1,]
colnames(temp.geno2)[1] <- "GermplasmId"

# Make temporary table with names of duplicate lines
duplicates.sub <- duplicates %>%
  select(GermplasmId_2, GermplasmId_name) %>%
  rename(GermplasmId = GermplasmId_name)

# Convert to data frames for following join step
temp.geno2 <- as.data.frame(temp.geno2)

# Make join table
duplicates.sub.geno <- left_join(duplicates.sub, temp.geno2, by = c("GermplasmId"= "GermplasmId")) %>%
  filter(complete.cases(.[ , 3])) #column 3 has the first marker

# Remove lines in genotype table if present in duplicate table
temp.geno2 <- temp.geno2 %>% filter(!GermplasmId %in% duplicates.sub.geno$GermplasmId)
# Adding germplasmid_2 to genotype table
temp.geno2 <- temp.geno2 %>%
  add_column(GermplasmId_2 = temp.geno2$GermplasmId, .before = "GermplasmId")
# Bind duplicates.sub.geno to temp geno table
temp.geno3 <- bind_rows(temp.geno2, duplicates.sub.geno)
# Transpose and set as data.frame to reset the table, drop row 2 since can't keep row
# Rename column line name headers using GermplasmId_2 names
temp.geno3 <- t(temp.geno3) %>% as.data.frame()
colnames(temp.geno3) <- temp.geno3[1, ] #rename the columns to ensure all GermplasmId_2 stay
colnames(temp.geno3)[1] <- "SNPID" #rename first column
temp.geno3 <- temp.geno3[-c(1:2), ] #delete the first two rows because the GermplasmId names were transferred to colnames

#Sort the GermplasmIds and convert chr to num
geno.dat <- temp.geno3 %>% select(SNPID, Chrom, Position, sort(names(.))) %>%
  mutate(across(-c(SNPID), as.numeric)) #convert chr to num all variables except SNPID
#Delete unnecessary intermediate files
temp.geno1 <- NULL; temp.geno2 <- NULL; temp.geno3 <- NULL; 

################################################################################
#------------------------------------------------------------------------------#
#     Loop through all traits for rrBLUP::GWAS() analysis                      #
#------------------------------------------------------------------------------#
################################################################################
# Initialize an empty data object to store all GWAS results
GWAS.out <- tibble(SNPID = geno.dat$SNPID, Chrom = geno.dat$Chrom, Position = geno.dat$Position)
# Subset data on trait while running through each trait, remove NA rows
traits.GWAS <- c("Chlorosis","Flower", "Height", "Lodging", "Maturity", "Oil", "Protein",  "Pubescence", "Seed.Quality", "Seed.Size", "Seed.Yield", "Shattering")
###Loop over each trait
startTime <- Sys.time()
for (trait in traits.GWAS) {
  GWAS.trait <- trait
# Loop through each trait column starting from the 5th column
# for (k in 5:ncol(Final_phenotype_data)) {
#   # Reset genotype data each time
#   genotype_data2 <- Final_genotype_table
#   
#   # Get name of trait from column heading
#   trait_name <- colnames(Final_phenotype_data)[k]
  
  # Make working data set for trait with just experimental strain, MG_factor, and trait data
  pheno.temp <- pheno.dat.final %>%
    select(GermplasmId, MG_factor, !!sym(GWAS.trait)) %>%
    drop_na(!!sym(GWAS.trait)) %>% #Remove NAs from temp table from both Trait column AND Maturity date column
  distinct(GermplasmId, .keep_all = TRUE)

  common_ges <- intersect(pheno.temp$GermplasmId, colnames(geno.dat)) #identify the experimental strains with both marker and pheno data
  geno <- geno.dat %>% select(SNPID, Chrom, Position, all_of(common_ges)) #select genotypes with marker and pheno data
  pheno <- pheno.temp %>% filter(GermplasmId %in% names(geno)) #remove phenotypes with no marker data
  #Run the rrBLUP::GWAS() function
  GWAS.fit <- GWAS(pheno = pheno, geno = geno, fixed=c("MG_factor"), n.PC=3, min.MAF=0.05, n.core=1, P3D=TRUE, plot=TRUE)
 # Add trait name to the GWAS result
  #GWAS.fit <-  GWAS.fit %>% mutate(Phenotype = GWAS.trait)
  # Combine the results
  GWAS.out <- left_join(GWAS.out, GWAS.fit)
  # Record status
  cat("\nDone with", GWAS.trait,"\n") #print out completed trait for progress tracking
}
endTime <- Sys.time() #calculate end time of BLUE calculation
print(endTime - startTime) # print recorded Time difference of 2.43 hours
write.csv(GWAS.out, "results/GWAS_Output_Table_All_traits.csv", row.names = F)
###If you do not want to wait for the rrBLUP::GWAS function to run and want to
#visualize the results, uncomment and read the GWAS output file below
#GWAS.out <- read.csv("results/GWAS_Output_Table_All_traits.csv")

################################################################################
#### Genomic Prediction using rrBLUP::mixed.solve()
### 10-fold cross validation repeated 100 times
### Predictive abilities compiled across all genotypes and within each maturity group
################################################################################

#Address the minimal NA values that were not imputed by TASSEL KNNi
### NAs are not allowed in the rrBLUP::mixed.solve()
### Naive imputation with replacing NA with mean of the marker
geno.mat <- geno.dat %>%
  select(-SNPID, -Chrom, -Position) %>% #select genotypes with marker and pheno data
  t() %>% data.frame() %>% # transpose to have markers in columns to create a numeric matrix (n x m) for mixed.solve
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>% #replace NA with the mean
  as.matrix() #convert the data.frame to matrix

# Assign value for the number of cross validations to perform
FOLD <- 10
Reps <- 100 #time consuming analysis with 100 replicates 

# Initialize an empty data frame to store all correlation results
cor.tab <- tibble(Rep.number = character(), Type = character(), Trait = character(), PredAbility = numeric())

common_ges <- intersect(pheno.dat.final$GermplasmId, rownames(geno.mat)) #identify common genotypes with marker and pheno data
PredVal <- pheno.dat.final %>% select(GermplasmId:Seed.Yield) %>%
          filter(GermplasmId %in% common_ges) %>%
          select(-MG) #data frame to store all the predicted values

#PredVal <- data.frame(GermplasmId = rownames(geno.mat)) #store all the predicted values
PredVal.comb <- tibble()
# Subset data on trait while running through each trait, remove NA rows
traits.GP <- c("Chlorosis", "Height", "Lodging", "Maturity", "Oil", "Protein", "Seed.Quality", "Seed.Size",  "Seed.Yield", "Shattering")

startTime <- Sys.time()
# Loop through all traits for analysis
for (P in 1:Reps) {
  for (trait in traits.GP) {
    GP.trait <- trait #obtain chr of trait name

    # Make working data set for trait with just experimental strain, MG_factor, and trait data
    pheno.temp <- pheno.dat.final %>%
      select(GermplasmId, MG_factor, !!sym(GP.trait)) %>%
      drop_na(!!sym(GP.trait)) %>% #Remove NAs from temp table from both Trait column AND Maturity date column
      distinct(GermplasmId, .keep_all = TRUE)
    
    common_ges <- intersect(pheno.temp$GermplasmId, rownames(geno.mat)) #identify the experimental strains with both marker and pheno data
    geno <- geno.mat[common_ges, ] #subset the matrix using row indexing
    #geno <- geno.dat %>% select(all_of(common_ges))#select genotypes with marker and pheno data
    pheno <- pheno.temp %>% filter(GermplasmId %in% rownames(geno)) #remove phenotypes with no marker data
    
    # Create design matrix (n x p) X for fixed effects to account for maturity group
    X <- model.matrix(~ 0 + pheno$MG_factor)
    
    # Shuffle data
    set.seed(P) #random number generation for reproducible results
    shuf.ind <- sample(nrow(pheno)) #random sample to shuffle data
    pheno.shuf <- pheno[shuf.ind, ] %>%
      mutate(!!paste0("Pred.", GP.trait) := NA) #create a column to store the GEBVs from each fold
    geno.shuf <- geno[shuf.ind, ]
    X.shuf <- X[shuf.ind, ]

    # Fold Cross Validation loop
    N.drop <- floor(nrow(pheno.shuf) / FOLD) #define the number of genotypes to drop according to the fold - using floor instead of round to ensure they will fit in the data.frame
    
    for (M in 1:FOLD) {
      row.drop <- ((M - 1) * N.drop + 1):(M * N.drop) #define the specific rows that will have the pheno and geno data removed from the GP model
      y <- pheno.shuf[-row.drop, ] #define phenotypic data
      Z <- geno.shuf[-row.drop, ]
      X.fold <- as.matrix(X.shuf[-row.drop, ])
      
      # Estimate marker effects using rrBLUP package and make predictions in test set
      rrModel <- mixed.solve(y = y[[GP.trait]], X = X.fold, Z = Z) #
      
      # Calculate genomic estimated breeding values
      mrkEff <- as.matrix(rrModel$u)
      rrGebv <- geno.shuf %*% mrkEff
      rrGebv <- as.data.frame(rrGebv)
      
      # Add GEBVs to table for dropped rows only
      # Assign values to the new variable for dropped rows only
      new.var.name <- paste0("Pred.", GP.trait)
      pheno.shuf <- pheno.shuf %>%
        mutate(!!new.var.name := replace(!!sym(new.var.name), row_number() %in% row.drop, rrGebv[row.drop, 1]))
        
      cat("\nDone with fold: ", M, "\n")
    }
    # Calculate correlation across all genotypes
    cor.all <- cor(pheno.shuf[[GP.trait]], pheno.shuf %>% pull(!!sym(new.var.name)), use = "complete.obs")
    cor.tab <- cor.tab %>%
      add_row(Rep.number = paste(P), Type = "All", Trait = GP.trait, PredAbility = cor.all)
    
    # Calculate correlations within MG sets
    for (MG in unique(pheno.shuf$MG)) {
      cor.MG <- cor(pheno.shuf %>% filter(MG_factor == MG) %>% pull(!!sym(GP.trait)),
                    pheno.shuf %>% filter(MG_factor == MG) %>% pull(!!sym(new.var.name)),
                    use = "complete.obs")
      cor.tab <- cor.tab %>%
        add_row(Rep.number = paste(P), Type = paste0(MG), Trait = GP.trait, PredAbility = cor.MG)
    }
    
    # Store all the predicted values in the long format
    PredValRep <- pheno.shuf %>%
      select(GermplasmId, MG_factor, !!sym(new.var.name)) %>%
      rename(!!paste0("Pred.", GP.trait) := !!sym(new.var.name)) %>%
      mutate(Rep.number = P)
    
    PredVal.comb <- bind_rows(PredVal.comb, PredValRep)
    
    cat("\nFinishing trait: ", GP.trait, "\n")
    }
    cat("\nFinishing Rep: ", P, "\n")
    write.csv(cor.tab, "results/Across_MG_CrossValidation_combined_cor_table.csv", row.names = F) #save the pred abilities from the intermediate runs
  }
  
# Write combined predicted values to a single output file
write.csv(PredVal.comb, "results/Across_MG_CrossValidation_combined_pred_values.csv", row.names = F)
  
# Write combined correlation results to a single output file
write.csv(cor.tab, "results/Across_MG_CrossValidation_combined_cor_table.csv", row.names = F)
endTime <- Sys.time()
cat("Total time taken: ", endTime - startTime, "\n")

#Create summary table of the correlation coefficients - predictive abilities
cor.Smry <- cor.tab %>% group_by(Trait, Type) %>%
  select(-Rep.number) %>%
  summarise_each(funs(mean,min,max,sd,se=sd(.)/sqrt(n())))


################################################################################
#Bootstrap the correlation coefficients to investigate their significance
#' Generalized bootstrapper by J Neyhart
#' 
#' @description Computes bootstrap resamples of a statistic
#' 
#' @param x A numeric vector.
#' @param y NULL (default) or a numeric vector for statistics that require a second
#' numeric vector
#' @param fun A character vector of the desired statistic function. For instance,
#' to calculate the variance of values in the numeric vector, pass "var".
#' @param boot.reps The number of bootstrapping replications.
#' @param alpha The significance level for calculating a confidence interval.
#' 
#' @import boot
#' 
#' @export
#' 
bootstrap <- function(x, y = NULL, fun, boot.reps = 1000, alpha = 0.05) {
  
  # Error handling
  boot.reps <- as.integer(boot.reps)
  
  # Error if the function does not exist
  stopifnot(exists(fun))
  # Otherwise get the function
  fun_torun <- get(fun)
  
  # Prob must be between 0 and 1
  alpha_check <- alpha > 0 | alpha < 1
  
  if (!alpha_check)
    stop("'alpha' must be between 0 and 1.")
  
  # Combine the data
  mat <- cbind(x, y)
  # How many columns
  n_col <- ncol(mat)
  
  
  ## Define a function for bootstrapping a single vector
  boot_fun_vec <- function(data, i) {
    # Use the index i to sample the data
    data_i <- data[i,,drop = FALSE]
    # Execute the function
    fun_torun(data_i[,1])
    
  }
  
  ## Define a function for bootstrapping a two vectors
  boot_fun_mat <- function(data, i) {
    # Use the index i to sample the data
    data_i <- data[i,]
    # Execute the function
    fun_torun(data_i[,1], data_i[,2])
    
  }
  
  # Calculate the base statistic
  base_stat <- if (n_col > 1) fun_torun(mat[,1], mat[,2]) else fun_torun(mat[,1])
  
  # If the base is not NA, proceed
  if (!is.na(base_stat)) {
    
    # Perform the bootstrapping
    if (n_col > 1) {
      boot_results <- boot::boot(data = mat, statistic = boot_fun_mat, R = boot.reps)
      
    } else {
      boot_results <- boot::boot(data = mat, statistic = boot_fun_vec, R = boot.reps)
      
    }
    
    # Standard error
    se <- sd(boot_results$t)
    # Bias
    bias <- mean(boot_results$t) - base_stat
    
    
    # Confidence interval
    ci_upper <- quantile(boot_results$t, 1 - (alpha / 2))
    ci_lower <- quantile(boot_results$t, (alpha / 2))
    
  } else {
    
    se <- bias <- ci_lower <- ci_upper <- NA
    
  }
  
  # Assemble list and return
  data.frame(statistic = fun, base = base_stat, se = se, bias = bias,
             ci_lower = ci_lower, ci_upper = ci_upper, row.names = NULL)
}

#Define parameters for the bootstrapping analysis
boot_reps <- 100000
alpha <- 0.05
set.seed(2025)

#PredAbil significance via bootstrapping analysis
predAbSig <- cor.tab %>% ##input data in the long format
  group_by(Trait, Type) %>% 
  do(cbind(bootstrap(x = .$PredAbility, y = NULL, fun = "mean", boot.reps = boot_reps, alpha = alpha), n_reps = length(.$PredAbility))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
  ungroup() %>%
  mutate(base_rounded = round(base, 2),
         base_signif = paste(base_rounded, annotation, sep= "")) %>%
  mutate(parameter = "mu") %>%
  relocate(parameter, .before = statistic)
