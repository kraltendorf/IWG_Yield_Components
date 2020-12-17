# IWG_Yield_Components

**Yield components analysis in intermediate wheatgrass spaced plants using correlation and structural equation modeling analyses. Here I include the five R scripts required to do this analysis and to create the associated figures.**

### Script 1 - Preparing and Analyzing Phenotypic Data.R
First we take raw data downloaded from the IWG databse, format it, and select the traits necessary for this analysis. Then, we run two types of linear models: 1) across all years and environments; and 2) within unique environments. We also caculate and output ANOVA tables, and estimated marginal means (emmeans) for genets.

*Required files:*
1. NAM_Data.txt # this file is the master phenotypic data available from the IWG database and available upon request (and will be made public after the entire dissertation is published)
2. backbone.csv # this allows you to create a balanced dataset to run the linear model

### Script 2 - Calculating Heritability and Correlations Between Traits.R
Here we use the emmeans to create a formatted correlation table, calculate broad and narrow sense heritabilities, calculate and output a trait summary table (min, mean, max, sd), and create Figure 3 (boxplots of phenotypic values). 

### Script 3 - Structural Equation Modeling.R
Using standardized emmeans, we test the initial model against the data, calculate and evaluate modification indices, and determine final models. We extract the path coefficients from the output and create Figure 4 and the Supplementary Tables S9 - 12. 

### Script 4 - Weather Data.R
Creating the weather and precipitation figure using data from NOAA and The Land Institute Weather Station. 
