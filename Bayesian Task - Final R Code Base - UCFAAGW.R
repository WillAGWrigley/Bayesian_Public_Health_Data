# Geog0115 - Spatial Bayesian Modelling - risk of mortality due to Suicide: 

# Install & import relevant libraries and packages: 

library("SpatialEpi")
library("INLA")
library("sf")
library("tmap")
library("spdep")
library("raster")
library("shiny")
library("shinydashboard")


# Read relevant datasets 

SuicideDataNST <- read.csv("Will Data - UK Suicide Fatalities NST.csv", header = TRUE, sep = ",")
SuicideDataST <- read.csv("Will Data - UK Suicide Fatalities ST.csv", header = TRUE, sep=, ",")

laboundaries <- read_sf("England Local Authority Shapefile.shp")
regionboundaries <- read_sf("England Regions Shapefile.shp")

# Create an empty map for England's local authorities


tm_shape(regionboundaries) + tm_polygons(alpha = 0, border.alpha = 0.3) +
  tm_shape(regionboundaries) + tm_polygons(alpha = 0, border.col = "black") +
  tm_compass(position = c ("right", "top")) + tm_scale_bar(position = c ("right", "bottom"))





# Data Preperation -- in order to estimate the risk of fatlity due to suicide
# at a region level in England we will need to create a data frame 
# that contains both number of fatalities, the denominators (i.e., reference
# population size) and the number of expected casualties. 

# Aggregate the data 

aggregate_suicide_DF <- aggregate(list(SuicideDataNST$Casulties, SuicideDataNST$Population), FUN = sum, by = list(SuicideDataNST$Code, SuicideDataNST$Name))
names(aggregate_suicide_DF) <-c("regcode", "regname", "fatalities", "population")

# Reassign labelling to columns
names(SuicideDataNST) <-c("regcode", "regname", "fatalities", "population")

# Calculate the expected number of fatalities 
# Sort original DF:
sort_Suicide <- SuicideDataNST[order(SuicideDataNST$regcode, SuicideDataNST$regname),]

expected_suicide <- expected(population = sort_Suicide$population, cases = sort_Suicide$fatalities, n.strata = 1)

# Match the expected counts to their corresponding Local Authority

aggregate_suicide_DF$expected_suicide <- expected_suicide[match(aggregate_suicide_DF$regcode, unique(sort_Suicide$regcode))]

# View the outputs of our aggregated data 

head(aggregate_suicide_DF, n=10)

# Dropping weird column name accidentally created during code iteration: 

keeps <- c("regcode", "regname", "fatalities", "population", "expected_suicide")
aggregate_suicide_DF[keeps]

head(aggregate_suicide_DF, n=10)

# Merge the data frames together - allowing adding the GIS data for 
# regional boundaries:

suicide_analysis_data_nontemp <- merge(regionboundaries, aggregate_suicide_DF, by.x = "code", by.y="regcode", all.x=TRUE)

names(suicide_analysis_data_nontemp)

# Create an adjacency matrix -- i.e. neighbour matrix 

# g --> we need this to be able to create the model: 

#adjacencyMatrix <- poly2nb(suicide_analysis_data_nontemp)

#nb2INLA("adjacencyObject.adj", adjacencyMatrix)

g <- inla.read.graph(filename = "adjacencyObject.adj")



#  Model Formation# 

# We need to specify the likelihood function which is typically defining the 
#statistical model needed for quantifying the risks of fatality due to suicide in the UK. 

suicide_analysis_data_nontemp$uI <- 1:nrow(suicide_analysis_data_nontemp)

# This type of model is known as the Besag-York-Mollie (BYM) model. 
# It is commonly referred to as the Spatial Conditional Autoregressive Model (CAR).

# Indexing the variable for the spatial random effects: 

# specify the id number for the random effect


#We specified the likelihood function as Poisson distribution 

# We specified the likelihood function as Poisson distribution because we are dealing with count data as an outcome.

# By default, we use priors of the loggamma distribution 

# specify priors
prior <- list(
  prec = list(prior = "loggamma", param = c(1, 0.0005)), # set prior for spatial random effects
  phi = list(prior = "loggamma", param = c(1, 0.0005))   # set prior for spatial random effects
)

# Write the formula with the outcome fatalities on the left side.

# defining out formula - fatalities a column in the dataframe
suicide_formula = fatalities ~ 1 + f(uI, model = "bym2", graph = g, hyper = prior)

# Bayesian inference using INLA - fitting the formula 
 
suicide_riskestimates <- inla(formula = suicide_formula, family = "poisson", data = suicide_analysis_data_nontemp, E = expected_suicide, control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE))
summary(suicide_riskestimates)



# intepreting overall risks fixed effects 

# Mapping Relative Risks --> exploring area specific risks --> extract risk rations and store:

suicide_riskratio <- suicide_riskestimates$summary.fitted.values
head(suicide_riskratio, n=10)

#The important estimates here are the posterior mean and 95% CrI results from the risk ratio object: 

suicide_analysis_data_nontemp$RR <- suicide_riskratio[, "mean"]       # Relative risk
suicide_analysis_data_nontemp$LL <- suicide_riskratio[, "0.025quant"] # Lower credibility limit
suicide_analysis_data_nontemp$UL <- suicide_riskratio[, "0.975quant"] # Upper credibility limit

# Here we can view the total RR for our dataset

head(suicide_analysis_data_nontemp, n = 10)

# Visualisation of the areas with low and high risks as well as show a map that
# indicates which areas are statistically significant: 

# Understand the lower and higher value of RR

summary(suicide_analysis_data_nontemp$RR)
#min (0.8041), (1.1549)


## Exceedence Probabilities:

# Set the threshold value
threshold <- 1.0  # Adjust the threshold as per your requirements

# Calculate the exceedance probabilities
suicide_analysis_data_nontemp$exceedance_prob <- ifelse(suicide_analysis_data_nontemp$RR > threshold, 1, 0)

# Calculate the deviation from the threshold
suicide_analysis_data_nontemp$deviation <- suicide_analysis_data_nontemp$RR - threshold



# Create risk categories for legend for Map

SuicideRiskCategorylist <- c("0.01 to 0.10", "0.11 to 0.25", "0.26 to 0.50", "0.51 to 0.75", "0.76 to 0.99", "1.00 (null value)",
                      ">1.00 to 1.10", "1.11 to 1.25", "1.26 to 1.50", "1.51 to 1.75", "1.76 to 2.00", "2.01 to 3.00",
                      "Above 3.00")

# Create the colours for the above categories - from extreme blues to extreme reds
SuicideRRPalette <- c("#33a6fe","#65bafe","#98cffe","#cbe6fe","#dfeffe","#fef9f9","#fed5d5","#feb1b1","#fe8e8e","#fe6a6a","#fe4646","#fe2424","#fe0000")

# Now generate categories
suicide_analysis_data_nontemp$RelativeRiskCat <- NA
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR>= 0 & suicide_analysis_data_nontemp$RR <= 0.10] <- -5
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 0.10 & suicide_analysis_data_nontemp$RR <= 0.25] <- -4
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 0.25 & suicide_analysis_data_nontemp$RR <= 0.50] <- -3
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 0.50 & suicide_analysis_data_nontemp$RR <= 0.75] <- -2
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 0.75 & suicide_analysis_data_nontemp$RR < 1] <- -1
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR == 1] <- 0
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 1.00 & suicide_analysis_data_nontemp$RR <= 1.10] <- 1
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 1.10 & suicide_analysis_data_nontemp$RR <= 1.25] <- 2
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 1.25 & suicide_analysis_data_nontemp$RR <= 1.50] <- 3
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 1.50 & suicide_analysis_data_nontemp$RR <= 1.75] <- 4
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 1.75 & suicide_analysis_data_nontemp$RR <= 2.00] <- 5
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 2.00 & suicide_analysis_data_nontemp$RR <= 3.00] <- 6
suicide_analysis_data_nontemp$RelativeRiskCat[suicide_analysis_data_nontemp$RR> 3.00 & suicide_analysis_data_nontemp$RR <= 250] <- 7

# Column creation for the significance categories:

suicide_analysis_data_nontemp$Significance <- NA
suicide_analysis_data_nontemp$Significance[suicide_analysis_data_nontemp$LL<1 & suicide_analysis_data_nontemp$UL>1] <- 0    # NOT SIGNIFICANT
suicide_analysis_data_nontemp$Significance[suicide_analysis_data_nontemp$LL==1 | suicide_analysis_data_nontemp$UL==1] <- 0  # NOT SIGNIFICANT
suicide_analysis_data_nontemp$Significance[suicide_analysis_data_nontemp$LL>1 & suicide_analysis_data_nontemp$UL>1] <- 1    # SIGNIFICANT INCREASE
suicide_analysis_data_nontemp$Significance[suicide_analysis_data_nontemp$LL<1 & suicide_analysis_data_nontemp$UL<1] <- -1   # SIGNIFICANT DECREASE


# Generate maps for relative risk and signifcant regions: 

## # # #  Generate categories based on the "RR" variable
categories <- cut(suicide_analysis_data_nontemp$RR,
                  breaks = c(0, 0.10, 0.25, 0.50, 0.75, 0.99, 1.00, 1.10, 1.25, 1.50, 1.75, 2.00, 3.00, Inf),
                  labels = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7),
                  include.lowest = TRUE)

## # # # #  Add the categories as a column to the dataset
suicide_analysis_data_nontemp$RelativeRiskCat <- categories


# Plot the map of relative risk
tm_shape(suicide_analysis_data_nontemp) + 
  tm_fill("RelativeRiskCat", style = "cat", title = "Relative Risk", palette = SuicideRRPalette, labels = SuicideRiskCategorylist) +
  tm_shape(regionboundaries) + tm_polygons(alpha = 0.05, border.col = "white", border.lwd = 0.05) + tm_text("name", size = 0.58, fontface = "bold", col = "black") +
  tm_layout(frame = FALSE, legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7,legend.position = c("right", "top")) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))


# Map of significance of regions:
tm_shape(suicide_analysis_data_nontemp) + 
  tm_fill("Significance", style = "cat", title = "Significance Categories", palette = c("white", "#fed5d5", "#fe0000"), labels = c("Significantly low", "Not Significant", "Significantly high")) +
  tm_shape(regionboundaries) + tm_polygons(alpha = 0.10, border.col = "white", border.lwd = 0.05) + tm_text("name", size = 0.58, fontface = "bold", col = "black") +
  tm_layout(frame = FALSE, legend.position = c("right", "bottom"), legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))


# Plot the map of deviation from threshold - chose to show exceedance as proximity from threshold 
tm_shape(suicide_analysis_data_nontemp) +
  tm_fill("deviation", style = "cont", title = "Exceedance Probability", palette = c("white", "#fed5d5", "#fe0000"), labels = c("Below Threshold", "Threshold", "Above Threshold")) +
  tm_shape(regionboundaries) + tm_polygons(alpha = 0.05, border.col = "white", border.lwd = 0.05) + tm_text("name", size = 0.58, fontface = "bold", col = "black") +
  tm_layout(frame = FALSE, legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7, legend.position = c("right", "top")) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))


#  RSHINY -- building a dashboard

# UI
# UI
ui <- fluidPage(
  dashboardPage(
    dashboardHeader(title = "Suicide Analysis Dashboard"),
    dashboardSidebar(),
    dashboardBody(
      # First plot: Significance of regions
      fluidRow(
        box(
          width = 6,
          height = 500,
          title = "Significance of Regions",
          solidHeader = TRUE,
          status = "primary",
          tmapOutput("significance_plot")
        )
      ),
      # Second plot: Relative Risk
      fluidRow(
        box(
          width = 6,
          height = 500,
          title = "Relative Risk",
          solidHeader = TRUE,
          status = "primary",
          tmapOutput("relative_risk_plot")
        )
      ),
      # Third plot: Exceedance Probability
      fluidRow(
        box(
          width = 6,
          height = 500,
          title = "Exceedance Probability",
          solidHeader = TRUE,
          status = "primary",
          tmapOutput("exceedance_chart")
        )
      )
    )
  )
)

# Server
server <- function(input, output) {
  output$significance_plot <- renderTmap({
    # Map of significance of regions
    tm_shape(suicide_analysis_data_nontemp) + 
      tm_fill("Significance", style = "cat", title = "Significance Categories",
              palette = c("white", "#fed5d5", "#fe0000"),
              labels = c("Significantly low", "Not Significant", "Significantly high")) +
      tm_shape(regionboundaries) +
      tm_polygons(alpha = 0.10, border.col = "white", border.lwd = 0.05) +
      tm_text("name", size = 0.58, fontface = "bold", col = "black") +
      tm_layout(frame = FALSE, legend.position = c("right", "bottom"),
                legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7) +
      tm_compass(position = c("right", "top")) +
      tm_scale_bar(position = c("right", "bottom"))
  })
  
  output$relative_risk_plot <- renderTmap({
    # Map of relative risk
    tm_shape(suicide_analysis_data_nontemp) + 
      tm_fill("RelativeRiskCat", style = "cat", title = "Relative Risk",
              palette = SuicideRRPalette, labels = SuicideRiskCategorylist) +
      tm_shape(regionboundaries) +
      tm_polygons(alpha = 0.05, border.col = "white", border.lwd = 0.05) +
      tm_text("name", size = 0.58, fontface = "bold", col = "black") +
      tm_layout(frame = FALSE, legend.outside = TRUE,
                legend.title.size = 0.8, legend.text.size = 0.7,
                legend.position = c("right", "top")) +
      tm_compass(position = c("right", "top")) +
      tm_scale_bar(position = c("right", "bottom"))
  })
  
  output$exceedance_chart <- renderTmap({
    # Exceedance probability chart logic
    tm_shape(suicide_analysis_data_nontemp) +
      tm_fill("deviation", style = "cont", title = "Exceedance Probability",
              palette = c("white", "#fed5d5", "#fe0000"),
              labels = c("Below Threshold", "Threshold", "Above Threshold")) +
      tm_shape(regionboundaries) +
      tm_polygons(alpha = 0.05, border.col = "white", border.lwd = 0.05) +
      tm_text("name", size = 0.58, fontface = "bold", col = "black") +
      tm_layout(frame = FALSE, legend.outside = TRUE,
                legend.title.size = 0.8, legend.text.size = 0.7,
                legend.position = c("right", "top")) +
      tm_compass(position = c("right", "top")) +
      tm_scale_bar(position = c("right", "bottom"))
  })
}

# Hey Presto - run the Shiny App :) 
shinyApp(ui, server)



