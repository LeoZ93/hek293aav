##################################################################
##################### qp with Error propagation ##################
##################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(robustbase)
library(readxl)
library(Deriv)

#################
## Data Import ##
#################

hek_data <- read_excel("data/hek/biomass/hek_data.xlsx")

#####################
## Biomass fitting ##
#####################

growth_results <- data.frame(state = character(),
                             timepoint = numeric(),
                             bm = numeric(),
                             err = numeric(),
                             err_perc = numeric(),
                             mu = numeric(),
                             mu_error = numeric(),
                             mu_err_perc = numeric(),
                             stringsAsFactors = FALSE)

# extract a row by index to perform analysis for one metabolite
## for example 4-Hydroxy-Prolin in row 6
df1 <- as.data.frame(as.numeric(hek_data[1, 2:81]))
df1 <- t(df1)
df1 <- df1*1000

# reshape columns
column_sequences <- list(c(1:4,9:12,17:20,25:28,33:36),
                         c(5:8,13:16,21:24,29:32,37:40),
                         c(41:44,49:52,57:60,65:68,73:76),
                         c(45:48,53:56,61:64,69:72,77:80))

df_list <- lapply(column_sequences, function(seq) {
  df_subset <- df1[, seq]
  names(df_subset) <- paste0("V", 1:length(seq))
  df_subset
})

df1 <- bind_rows(df_list)

times <- data.frame(t(c(0,0,0,0,4,4,4,4,24,24,24,24,48,48,48,48,72,72,72,72)))
names(times) <- names(df1)  # Ensure column names match

# Add the new row to the dataframe
df1 <- rbind(df1, times)

logistic <- function(x, a, b, h) {
  a / (1 + exp(-b * (x - h)))
}

deriv_bm <- Deriv(logistic, "x")

parderiv_a_bm <- Deriv(logistic, "a")
parderiv_b_bm <- Deriv(logistic, "b")
parderiv_h_bm <- Deriv(logistic, "h")
deriv_bm <- Deriv(logistic, "x")
parderiv_deriv_bm_a <- Deriv(deriv_bm, "a")
parderiv_deriv_bm_b <- Deriv(deriv_bm, "b")
parderiv_deriv_bm_h <- Deriv(deriv_bm, "h")

x_values <- c(4, 24, 48, 72)

# Assuming growth_results is initialized somewhere before
for(i in 1:4) {
  
  #  if (exists("a_bm")) {
  #    rm(model, a_bm, b_bm)
  #  }
  
  if(i == 1) {
    state_value <- "HP_TR"
  } 
  else if(i == 2) {
    state_value <- "HP_MO"
  } 
  else if(i == 3) {
    state_value <- "LP_TR"
  } 
  else if(i == 4) {
    state_value <- "LP_MO"
  } 
  
  
  df <- df1[c(i,5),]
  df <- data.frame(t(df))
  colnames(df) <- c("biomass", "time_point")
  
  # Exclude outlier rows
  #  df <- df[-c(1:4),]
  df$biomass <- as.numeric(df$biomass)
  df$time_point <- as.numeric(df$time_point)
  
  n = 1
  
  # Create a variable to keep track of whether the model fitting succeeded
  success = FALSE
  
  while(!success && n > -1){
    # Try to fit the model
    tryCatch({
      # Attempt to fit the model with the current value of r
      model <- nlrob(biomass ~ logistic(time_point, a, b, h), 
                     start = list(a = max(df$biomass), b = n, h = median(df$time_point)), 
                     data = df)
      
      # If the model fitting did not produce an error, set success to TRUE
      success = TRUE
    },
    error = function(e){
      # If an error occurred, print the error message and decrease r by 0.01
      #    print(paste("Failed with r =", r, ": ", e$message, sep = ""))
    },
    {n = n - 0.01})
  }
  
  # If the loop ended because r reached 0, print a message to inform the user
  if(n < -1){
    print("Failed to fit the model with all tried values of r.")
  }
  
  # If the model fitting succeeded, the model object will contain the fitted model
  if(success){
    print(paste("Successfully fit the model", i, " for state", state_value, " with b =", n))
    #    print(summary(model))
  }
  # Fit model
  smry <- summary(model)
  
  if(is.na(smry$coefficients[1,2])) {
    print("Model does not exist")
  } else {
    
    a_bm <- smry$coefficients[1,1]
    b_bm <- smry$coefficients[2,1]
    a_bm_err <- smry$coefficients[1,2]
    b_bm_err <- smry$coefficients[2,2]
    h <- smry$coefficients[3,1]
    h_err <- smry$coefficients[3,2]
    #  c <- smry$coefficients[4,1]
    
    
    #      df$predicted_biomass <- predict(model)
    
    # Plot
    #      ggplot(df, aes(x = time_point, y = biomass)) +
    #        geom_point() +
    #        geom_line(aes(y = predicted_biomass), color = "blue") +
    #        theme_minimal() +
    #        labs(x = "Time point", y = "Biomass",
    #             title = "Cell growth over time {state_value}",
    #             subtitle = "Observed data and fitted logistic growth curve")
    
    for (x in x_values) {
      
      bm_val <- logistic(x, a_bm, b_bm, h)
      bm_err <- sqrt( (parderiv_a_bm(x,a_bm,b_bm, h)*a_bm_err)^2 + ( parderiv_b_bm(x, a_bm, b_bm, h) * b_bm_err)^2 + (parderiv_h_bm(x, a_bm, b_bm, h) * h_err)^2 )
      
      deriv_bm_val <- deriv_bm(x, a_bm, b_bm, h)
      deriv_bm_err <- sqrt( ( (parderiv_deriv_bm_a(x,a_bm, b_bm, h)) * a_bm_err)^2 + ( (parderiv_deriv_bm_b(x,a_bm,b_bm, h)) * b_bm_err)^2 + ( parderiv_deriv_bm_h(x,a_bm,b_bm,h)*h_err )^2)
      
      mu_val <- deriv_bm_val / bm_val
      mu_err <- sqrt( ( (1 / bm_val) * deriv_bm_err )^2 + ( (-deriv_bm_val/bm_val^2) * bm_err )^2 )
      
      
      #      bm_val <- logistic(x, a_bm, b_bm, h)
      
      #      deriv_bm_val <- deriv_bm(x, a_bm, b_bm, h)
      
      #      mu_val <- deriv_bm_val / bm_val
      
      growth_results <- growth_results %>%
        bind_rows(data.frame(state = state_value,
                             timepoint = x,
                             bm = bm_val,
                             #  mu = mu_val
                             err = bm_err,
                             err_perc = bm_err / bm_val,
                             mu = mu_val,
                             mu_error = mu_err,
                             mu_err_perc = mu_err / mu_val
        ))
    }
  }
}


################
## qp fitting ##
################

Results_Metabolomics_Media <- read_excel("data/hek/metabolome/Results_Metabolomics_Media_1.xlsx")
# remove glutamine to take our values
Results_Metabolomics_Media <- Results_Metabolomics_Media[-59, ]

## Exponential
exponential <- function(x, a, r) {
  a * exp(r * x)
}

exp_deriv <- Deriv(exponential, "x")
exp_deriv_parderiv_a <- Deriv(exp_deriv, "a")
exp_deriv_parderiv_r <- Deriv(exp_deriv, "r")


## Quadratic
quadratic <- function(x, a, b, c) {
  a * x^2 + b * x + c
}

quad_deriv <- Deriv(quadratic, "x")
quad_deriv_parderiv_a <- Deriv(quad_deriv, "a")
quad_deriv_parderiv_b <- Deriv(quad_deriv, "b")
quad_deriv_parderiv_c <- Deriv(quad_deriv, "c")

## Cubic
cubic <- function(x,a,b,c,d) {
  a * x^3 + b * x^2 + c * x + d
}

cubic_deriv <- Deriv(cubic, "x")
cubic_deriv_parderiv_a <- Deriv(cubic_deriv, "a")
cubic_deriv_parderiv_b <- Deriv(cubic_deriv, "b")
cubic_deriv_parderiv_c <- Deriv(cubic_deriv, "c")
cubic_deriv_parderiv_d <- Deriv(cubic_deriv, "d")

###################
## Gln Depletion ##
###################

gln_depletion <- read_excel("data/hek/metabolome/gln_depletion.xlsx")

#df <- gln_depletion[c(1,3), c(2,3,6,7,10,11,14,15)]
#df <- gln_depletion[c(1,3), c(4,5,8,9,12,13,16,17)]
df <- gln_depletion[c(1,3), 2:17]
df <- data.frame(t(df))
colnames(df) <- c("concentration", "time_point")

# Convert columns to numeric
df$concentration <- as.numeric(df$concentration)*1000
df$time_point <- as.numeric(df$time_point)

model <- nlrob(concentration ~ exponential(time_point, a, r), 
               start = list(a = 6, r = 0.3), 
               data = df)

smry <- summary(model)
k_depl <- smry$coefficients[2,1]


#####

Results_Metabolomics_Media <- read_excel("data/hek/metabolome/Results_Metabolomics_Media_1.xlsx")
lac_ammonia_glc_fitting <- read_excel("data/hek/biomass/lac_ammonia_glc_fitting.xlsx")

# remove glutamine to take our values
Results_Metabolomics_Media <- Results_Metabolomics_Media[-59, ]

specific_rates <- data.frame(metabolite = character(),
                             state = character(),
                             time_point = numeric(),
                             qs = numeric(),
                             std_err = numeric(),
                             err_perc = numeric(),
                             ub = numeric(),
                             lb = numeric(),
                             stringsAsFactors = FALSE)


exp_analytes <- c("Arginine",
                  "Asparagine",
                  "Aspartic acid",
                  "Glucose",
                  "Glutamate",
                  "Glutamine",
                  "Histidine",
                  "Isoleucine",
                  "Leucine",
                  "Lysine",
                  "Methionine",
                  "Phenylalanine",
                  "Serine",
                  "Threonine",
                  "Tryptophan",
                  "Tyrosine",
                  "Valine",
                  "alpha-Ketoglutaric acid",
                  "Choline",
                  "Citric acid",
                  "Cystine",
                  "Deoxyuridine",
                  "Fumaric acid",
                  "Homoserine",
                  "Isocitric acid",
                  "Malic acid",
                  "4-Hydroxy-proline"
)

#exp_analytes = "Glutamine"
quad_analytes_1 <- c("Alanine")
quad_analytes_2 <- c("Pro", "Gly", "Lac")
cub_analytes <- c("Nh3")

name_suffixes <- c("HP_TR", "HP_MO", "LP_TR", "LP_MO")

#exp_analytes = c("Tyrosine")

for (analyte in exp_analytes) {
  
  extracted_row <- subset(Results_Metabolomics_Media, Results_Metabolomics_Media$Analyte == analyte)
  df1 <- extracted_row[, 4:83]
  
  column_sequences <- list(c(1:4,9:12,17:20,25:28,33:36),
                           c(5:8,13:16,21:24,29:32,37:40),
                           c(41:44,49:52,57:60,65:68,73:76),
                           c(45:48,53:56,61:64,69:72,77:80))
  
  df_list <- lapply(column_sequences, function(seq) {
    df_subset <- df1[, seq]
    names(df_subset) <- paste0("V", 1:length(seq))
    df_subset
  })
  
  df1 <- bind_rows(df_list)
  
  times <- data.frame(t(c(0,0,0,0,4,4,4,4,24,24,24,24,48,48,48,48,72,72,72,72)))
  names(times) <- names(df1)  # Ensure column names match
  
  # Add the new row to the dataframe
  df1 <- rbind(df1, times)
  
  for(i in 1:4){
    
    #    if (exists("a_met")) {
    #      rm(model, a_met, b_met)
    #    }
    
    # Extract data for current row
    df <- df1[c(i,5),]
    df <- data.frame(t(df))
    colnames(df) <- c("concentration", "time_point")
    df <- df[-c(1,2,3,4),]
    
    # Convert columns to numeric
    df$concentration <- as.numeric(df$concentration)/1000
    df$time_point <- as.numeric(df$time_point)
    
    df <- df[!is.na(df$concentration), ]
    
    #    if (analyte == "Glucose") {
    #      df[is.na(df)] <- 0.00001
    #    }
    
    file_name <- paste0("/home/users/lzehetner/",analyte,"_",name_suffixes[i], ".png")
    
    png(file_name, width = 1200, height = 700)
    
    ggplot(df, aes(x = time_point, y = concentration)) +
      geom_point(size = 10) +  # Increase the size of the dots
      geom_line(aes(y = predicted_concentration), color = "blue", size = 5) +
      theme_minimal() +
      labs(x = "Hours post transfection", y = "Concentration [mmol/l]",
           title = "",
           subtitle = "") +
      theme(
        text = element_text(size = 0),  # Increase the font size for all text elements
        axis.title = element_text(size = 50),  # Increase the font size for axis titles
        axis.text = element_text(size = 50)  # Increase the font size for axis labels
      )
    
    dev.off()
    
    r = 1
    
    # Create a variable to keep track of whether the model fitting succeeded
    success = FALSE
    
    while(!success && r >=-10 ){
      # Try to fit the model
      tryCatch({
        # Attempt to fit the model with the current value of r
        model <- nlrob(concentration ~ exponential(time_point, a, r), 
                       start = list(a = mean(df[c(1:4),1]), r = r), 
                       #                       method = "mtl",
                       data = df)
        
        # If the model fitting did not produce an error, set success to TRUE
        success = TRUE
      },
      error = function(e){
        # If an error occurred, print the error message and decrease r by 0.01
        #    print(paste("Failed with r =", r, ": ", e$message, sep = ""))
      },
      {r = r - 0.001})
    }
    
    # If the loop ended because r reached 0, print a message to inform the user
    
    
    # If the model fitting succeeded, the model object will contain the fitted model
    #    if(model$status == "converged"){
    #      print(paste("Successfully fit the model", i, "for analyte", analyte, "with r =", r))
    #      print(summary(model))
    #    }
    
    if (c(analyte == "Aspartic acid" && c(i == 3 | i == 4)) | c(analyte == "4-Hydroxy-proline") ) {
      
      print("Fitting LP Aspartic acid")
      
      r = 1
      
      # Create a variable to keep track of whether the model fitting succeeded
      success = FALSE
      
      while(!success && r >=-1 ){
        # Try to fit the model
        tryCatch({
          # Attempt to fit the model with the current value of r
          model <- nls(concentration ~ exponential(time_point, a, r), 
                       start = list(a = mean(df[1,1]), r = r), 
                       #                       method = "mtl",
                       data = df)
          
          # If the model fitting did not produce an error, set success to TRUE
          success = TRUE
        },
        error = function(e){
          # If an error occurred, print the error message and decrease r by 0.01
          #    print(paste("Failed with r =", r, ": ", e$message, sep = ""))
        },
        {r = r - 0.001})
      }
    }
    
    df$predicted_concentration <- predict(model)
    
    file_name <- paste0("/home/users/lzehetner/",analyte,"_",name_suffixes[i], ".png")
    
    png(file_name, width = 1200, height = 700)
    
    ggplot(df, aes(x = time_point, y = concentration)) +
      geom_point() +
      geom_line(aes(y = predicted_concentration), color = "blue") +
      theme_minimal() +
      labs(x = "Time point", y = "Concentration",
           title = "Concentration over time",
           subtitle = "Observed data and fitted exponential growth curve")
    
    dev.off()
    
    smry <- summary(model)
    
    if (!is.na(smry$coefficients[1,2])){
      print(paste("Successfully fit the model", i, "for analyte", analyte, "with r =", r))
      #      print(summary(model))
    }
    
    a_met <- smry$coefficients[1,1]
    b_met <- smry$coefficients[2,1]
    a_met_err <- smry$coefficients[1,2]
    b_met_err <- smry$coefficients[2,2]
    
    tps <- c(4,24,48,72)
    
    for (x in tps) {
      if (i == 1) {
        state_value = "HP_TR"
        if (x == 4) {
          bm_val <- growth_results[1,3]
          bm_err <- growth_results[1,4]
        }
        if (x == 24) {
          bm_val <- growth_results[2,3]
          bm_err <- growth_results[2,4]
        }
        if (x == 48) {
          bm_val <- growth_results[3,3]
          bm_err <- growth_results[3,4]
        }
        if (x == 72) {
          bm_val <- growth_results[4,3]
          bm_err <- growth_results[4,4]
        }
      }
      if (i == 2) {
        state_value = "HP_MO"
        if (x == 4) {
          bm_val <- growth_results[5,3]
          bm_err <- growth_results[5,4]
        }
        if (x == 24) {
          bm_val <- growth_results[6,3]
          bm_err <- growth_results[6,4]
        }
        if (x == 48) {
          bm_val <- growth_results[7,3]
          bm_err <- growth_results[7,4]
        }
        if (x == 72) {
          bm_val <- growth_results[8,3]
          bm_err <- growth_results[8,4]
        }
      }
      if (i == 3) {
        state_value = "LP_TR"
        if (x == 4) {
          bm_val <- growth_results[9,3]
          bm_err <- growth_results[9,4]
        }
        if (x == 24) {
          bm_val <- growth_results[10,3]
          bm_err <- growth_results[10,4]
        }
        if (x == 48) {
          bm_val <- growth_results[11,3]
          bm_err <- growth_results[11,4]
        }
        if (x == 72) {
          bm_val <- growth_results[12,3]
          bm_err <- growth_results[12,4]
        }
      }
      if (i == 4) {
        state_value = "LP_MO"
        if (x == 4) {
          bm_val <- growth_results[13,3]
          bm_err <- growth_results[13,4]
        }
        if (x == 24) {
          bm_val <- growth_results[14,3]
          bm_err <- growth_results[14,4]
        }
        if (x == 48) {
          bm_val <- growth_results[15,3]
          bm_err <- growth_results[15,4]
        }
        if (x == 72) {
          bm_val <- growth_results[16,3]
          bm_err <- growth_results[16,4]
        }
      }
      
      if (analyte == "Glutamine") {
        deriv_val_1 <- exp_deriv(x, a_met, b_met)
        if (x == 4) {
          deriv_val <- deriv_val_1 - k_depl * mean(df[1:4,1])
        }
        if (x == 24) {
          deriv_val <- deriv_val_1 - k_depl * mean(df[5:8,1])
        }
        if (x == 48) {
          deriv_val <- deriv_val_1 - k_depl * mean(df[9:12,1])
        }
        if (x == 72) {
          deriv_val <- deriv_val_1 - k_depl * mean(df[13:16,1])
        }
        
        
        deriv_err <- sqrt((exp_deriv_parderiv_a(x, a_met, b_met) * a_met_err)^2 + (exp_deriv_parderiv_r(x, a_met, b_met) * b_met_err)^2)
        
        qs_val <- deriv_val/bm_val
        qs_err <- sqrt((1/bm_val*deriv_err)^2 + ((-deriv_val/bm_val^2)*bm_err)^2)
        
        specific_rates <- specific_rates %>%
          bind_rows(data.frame(metabolite = analyte,
                               state = state_value,
                               time_point = x,
                               qs = qs_val,
                               std_err = qs_err,
                               err_perc = qs_err / qs_val,
                               ub = qs_val + qs_err,
                               lb = qs_val - qs_err
                               #                             std_err = qs_err,
                               #                             err_perc = qs_err / qs_val
          ))
      }
      if (analyte != "Glutamine") {
        deriv_val <- exp_deriv(x, a_met, b_met)
        deriv_err <- sqrt((exp_deriv_parderiv_a(x, a_met, b_met) * a_met_err)^2 + (exp_deriv_parderiv_r(x, a_met, b_met) * b_met_err)^2)
        
        qs_val <- deriv_val/bm_val
        qs_err <- sqrt((1/bm_val*deriv_err)^2 + ((-deriv_val/bm_val^2)*bm_err)^2)
        
        specific_rates <- specific_rates %>%
          bind_rows(data.frame(metabolite = analyte,
                               state = state_value,
                               time_point = x,
                               qs = qs_val,
                               std_err = qs_err,
                               err_perc = qs_err / qs_val,
                               ub = qs_val + qs_err,
                               lb = qs_val - qs_err
                               #                             std_err = qs_err,
                               #                             err_perc = qs_err / qs_val
          ))
      }
    }
  }
}

for (analyte in quad_analytes_1) {
  
  extracted_row <- subset(Results_Metabolomics_Media, Results_Metabolomics_Media$Analyte == analyte)
  df1 <- extracted_row[, 4:83]
  
  column_sequences <- list(c(1:4,9:12,17:20,25:28,33:36),
                           c(5:8,13:16,21:24,29:32,37:40),
                           c(41:44,49:52,57:60,65:68,73:76),
                           c(45:48,53:56,61:64,69:72,77:80))
  
  df_list <- lapply(column_sequences, function(seq) {
    df_subset <- df1[, seq]
    names(df_subset) <- paste0("V", 1:length(seq))
    df_subset
  })
  
  df1 <- bind_rows(df_list)
  
  times <- data.frame(t(c(0,0,0,0,4,4,4,4,24,24,24,24,48,48,48,48,72,72,72,72)))
  names(times) <- names(df1)  # Ensure column names match
  
  # Add the new row to the dataframe
  df1 <- rbind(df1, times)
  
  for(i in 2:2){
    
    #    if (exists("a_met")) {
    #      rm(model, a_met, b_met, c_met)
    #    }
    # Extract data for current row
    df <- df1[c(i,5),]
    df <- data.frame(t(df))
    colnames(df) <- c("concentration", "time_point")
    df <- df[-c(1,2,3,4),]
    
    # Convert columns to numeric
    df$concentration <- as.numeric(df$concentration)/1000
    df$time_point <- as.numeric(df$time_point)

    #    df <- df[!is.na(df$concentration), ]
    
    model <- try(nlrob(concentration ~ quadratic(time_point, a, b, c), 
                       start = list(a = -0.1, b = -0.1, c = 0.1), 
                       data = df)
    )
    
    df$predicted_concentration <- predict(model)
    
    smry <- summary(model)
    
    #    if(model$status != "converged") {
    
    #      model <- try(nls(concentration ~ quadratic(time_point, a, b, c), 
    #                       start = list(a = 0.1, b = 0.1, c = 0.1), 
    #                       data = df))
    #    }
    
    a_met <- smry$coefficients[1,1]
    b_met <- smry$coefficients[2,1]
    c_met <- smry$coefficients[3,1]
    a_met_err <- smry$coefficients[1,2]
    b_met_err <- smry$coefficients[2,2]
    c_met_err <- smry$coefficients[3,2]
    
    tps <- c(4,24,48,72)
    
    for (x in tps) {
      if (i == 1) {
        state_value = "HP_TR"
        if (x == 4) {
          bm_val <- growth_results[1,3]
          bm_err <- growth_results[1,4]
        }
        if (x == 24) {
          bm_val <- growth_results[2,3]
          bm_err <- growth_results[2,4]
        }
        if (x == 48) {
          bm_val <- growth_results[3,3]
          bm_err <- growth_results[3,4]
        }
        if (x == 72) {
          bm_val <- growth_results[4,3]
          bm_err <- growth_results[4,4]
        }
      }
      if (i == 2) {
        state_value = "HP_MO"
        if (x == 4) {
          bm_val <- growth_results[5,3]
          bm_err <- growth_results[5,4]
        }
        if (x == 24) {
          bm_val <- growth_results[6,3]
          bm_err <- growth_results[6,4]
        }
        if (x == 48) {
          bm_val <- growth_results[7,3]
          bm_err <- growth_results[7,4]
        }
        if (x == 72) {
          bm_val <- growth_results[8,3]
          bm_err <- growth_results[8,4]
        }
      }
      if (i == 3) {
        state_value = "LP_TR"
        if (x == 4) {
          bm_val <- growth_results[9,3]
          bm_err <- growth_results[9,4]
        }
        if (x == 24) {
          bm_val <- growth_results[10,3]
          bm_err <- growth_results[10,4]
        }
        if (x == 48) {
          bm_val <- growth_results[11,3]
          bm_err <- growth_results[11,4]
        }
        if (x == 72) {
          bm_val <- growth_results[12,3]
          bm_err <- growth_results[12,4]
        }
      }
      if (i == 4) {
        state_value = "LP_MO"
        if (x == 4) {
          bm_val <- growth_results[13,3]
          bm_err <- growth_results[13,4]
        }
        if (x == 24) {
          bm_val <- growth_results[14,3]
          bm_err <- growth_results[14,4]
        }
        if (x == 48) {
          bm_val <- growth_results[15,3]
          bm_err <- growth_results[15,4]
        }
        if (x == 72) {
          bm_val <- growth_results[16,3]
          bm_err <- growth_results[16,4]
        }
      }
      
      deriv_val <- quad_deriv(x, a_met, b_met, c_met)
      deriv_err <- sqrt((quad_deriv_parderiv_a(x, a_met, b_met, c_met) * a_met_err)^2 + (quad_deriv_parderiv_b(x, a_met, b_met, c_met) * b_met_err)^2 + (quad_deriv_parderiv_c(x, a_met, b_met, c_met)*c_met_err)^2 )
      
      qs_val <- deriv_val/bm_val
      qs_err <- sqrt((1/bm_val*deriv_err)^2 + ((-deriv_val/bm_val^2)*bm_err)^2)
      
      specific_rates <- specific_rates %>%
        bind_rows(data.frame(metabolite = analyte,
                             state = state_value,
                             time_point = x,
                             qs = qs_val,
                             std_err = qs_err,
                             err_perc = qs_err / qs_val,
                             ub = qs_val + qs_err,
                             lb = qs_val - qs_err
                             #                             std_err = qs_err,
                             #                             err_perc = qs_err / qs_val
        ))
    }
  }
}

for (analyte in quad_analytes_2) {
  
  extracted_row <- subset(lac_ammonia_glc_fitting,  lac_ammonia_glc_fitting[, 1]== analyte)
  df1 <- extracted_row[, 2:161]
  
  # df1 <- as.data.frame(Results_Metabolomics_Media[14, 4:83])
  
  # reshape columns
  column_sequences <- list(c(1:4, 9:12, 17:20, 25:28, 33:36, 41:44, 49:52, 57:60, 65:68, 73:76),
                           c(5:8, 13:16, 21:24, 29:32, 37:40, 45:48, 53:56, 61:64, 69:72, 77:80),
                           c(81:84, 89:92, 97:100, 105:108, 113:116, 121:124, 129:132, 137:140, 145:148, 153:156),
                           c(85:88, 93:96, 101:104, 109:112, 117:120, 125:128, 133:136, 141:144, 149:152, 157:160)
  )
  
  df_list <- lapply(column_sequences, function(seq) {
    df_subset <- df1[, seq]
    names(df_subset) <- paste0("V", 1:length(seq))
    df_subset
  })
  
  df1 <- bind_rows(df_list)
  
  times <- data.frame(t(c(0,0,0,0,4,4,4,4,21,21,21,21,24,24,24,24,27,27,27,27,45,45,45,45,48,48,48,48,51,51,51,51,69,69,69,69,72,72,72,72)))
  names(times) <- names(df1)  # Ensure column names match
  
  # Add the new row to the dataframe
  df1 <- rbind(df1, times)
  
  for(i in 1:4){
    
    if (exists("a_met")) {
      rm(model, a_met, b_met, c_met)
    }
    
    # Extract data for current row
    df <- df1[c(i,5),]
    df <- data.frame(t(df))
    colnames(df) <- c("concentration", "time_point")
    df <- df[-c(1,2,3,4),]
    
    # Convert columns to numeric
    df$concentration <- as.numeric(df$concentration)/1000
    df$time_point <- as.numeric(df$time_point)
    df <- df[!is.na(df$concentration), ]
    
    model <- try(nlrob(concentration ~ quadratic(time_point, a, b, c), 
                       start = list(a = -1, b = -1, c = 0.5), 
                       data = df)
    )
    
    df$predicted_concentration <- predict(model)
    
    
    #    if(model$status != "converged") {
    
    #      model <- try(nls(concentration ~ quadratic(time_point, a, b, c), 
    #                       start = list(a = -0.1, b = -0.1, c = 0.1), 
    #                       data = df))
    #    }
    
    
    smry <- summary(model)
    a_met <- smry$coefficients[1,1]
    b_met <- smry$coefficients[2,1]
    c_met <- smry$coefficients[3,1]
    a_met_err <- smry$coefficients[1,2]
    b_met_err <- smry$coefficients[2,2]
    c_met_err <- smry$coefficients[3,2]
    
    tps <- c(4,24,48,72)
    
    for (x in tps) {
      if (i == 1) {
        state_value = "HP_TR"
        if (x == 4) {
          bm_val <- growth_results[1,3]
          bm_err <- growth_results[1,4]
        }
        if (x == 24) {
          bm_val <- growth_results[2,3]
          bm_err <- growth_results[2,4]
        }
        if (x == 48) {
          bm_val <- growth_results[3,3]
          bm_err <- growth_results[3,4]
        }
        if (x == 72) {
          bm_val <- growth_results[4,3]
          bm_err <- growth_results[4,4]
        }
      }
      if (i == 2) {
        state_value = "HP_MO"
        if (x == 4) {
          bm_val <- growth_results[5,3]
          bm_err <- growth_results[5,4]
        }
        if (x == 24) {
          bm_val <- growth_results[6,3]
          bm_err <- growth_results[6,4]
        }
        if (x == 48) {
          bm_val <- growth_results[7,3]
          bm_err <- growth_results[7,4]
        }
        if (x == 72) {
          bm_val <- growth_results[8,3]
          bm_err <- growth_results[8,4]
        }
      }
      if (i == 3) {
        state_value = "LP_TR"
        if (x == 4) {
          bm_val <- growth_results[9,3]
          bm_err <- growth_results[9,4]
        }
        if (x == 24) {
          bm_val <- growth_results[10,3]
          bm_err <- growth_results[10,4]
        }
        if (x == 48) {
          bm_val <- growth_results[11,3]
          bm_err <- growth_results[11,4]
        }
        if (x == 72) {
          bm_val <- growth_results[12,3]
          bm_err <- growth_results[12,4]
        }
      }
      if (i == 4) {
        state_value = "LP_MO"
        if (x == 4) {
          bm_val <- growth_results[13,3]
          bm_err <- growth_results[13,4]
        }
        if (x == 24) {
          bm_val <- growth_results[14,3]
          bm_err <- growth_results[14,4]
        }
        if (x == 48) {
          bm_val <- growth_results[15,3]
          bm_err <- growth_results[15,4]
        }
        if (x == 72) {
          bm_val <- growth_results[16,3]
          bm_err <- growth_results[16,4]
        }
      }
      
      deriv_val <- quad_deriv(x, a_met, b_met, c_met)
      deriv_err <- sqrt((quad_deriv_parderiv_a(x, a_met, b_met, c_met) * a_met_err)^2 + (quad_deriv_parderiv_b(x, a_met, b_met, c_met) * b_met_err)^2 + (quad_deriv_parderiv_c(x, a_met, b_met, c_met)*c_met_err)^2 )
      
      qs_val <- deriv_val/bm_val
      qs_err <- sqrt((1/bm_val*deriv_err)^2 + ((-deriv_val/bm_val^2)*bm_err)^2)
      
      specific_rates <- specific_rates %>%
        bind_rows(data.frame(metabolite = analyte,
                             state = state_value,
                             time_point = x,
                             qs = qs_val,
                             std_err = qs_err,
                             err_perc = qs_err / qs_val,
                             ub = qs_val + qs_err,
                             lb = qs_val - qs_err
                             #                             std_err = qs_err,
                             #                             err_perc = qs_err / qs_val
        ))
    }
  }
}

for (analyte in cub_analytes) {
  extracted_row <- subset(lac_ammonia_glc_fitting,  lac_ammonia_glc_fitting[, 1]== analyte)
  df1 <- extracted_row[, 2:161]
  
  # df1 <- as.data.frame(Results_Metabolomics_Media[14, 4:83])
  
  # reshape columns
  column_sequences <- list(c(1:4, 9:12, 17:20, 25:28, 33:36, 41:44, 49:52, 57:60, 65:68, 73:76),
                           c(5:8, 13:16, 21:24, 29:32, 37:40, 45:48, 53:56, 61:64, 69:72, 77:80),
                           c(81:84, 89:92, 97:100, 105:108, 113:116, 121:124, 129:132, 137:140, 145:148, 153:156),
                           c(85:88, 93:96, 101:104, 109:112, 117:120, 125:128, 133:136, 141:144, 149:152, 157:160)
  )
  
  df_list <- lapply(column_sequences, function(seq) {
    df_subset <- df1[, seq]
    names(df_subset) <- paste0("V", 1:length(seq))
    df_subset
  })
  
  df1 <- bind_rows(df_list)
  
  times <- data.frame(t(c(0,0,0,0,4,4,4,4,21,21,21,21,24,24,24,24,27,27,27,27,45,45,45,45,48,48,48,48,51,51,51,51,69,69,69,69,72,72,72,72)))
  names(times) <- names(df1)  # Ensure column names match
  
  # Add the new row to the dataframe
  df1 <- rbind(df1, times)
  
  for(i in 1:4){
    
    #    if (exists("a_met")) {
    #      rm(model, a_met, b_met, c_met)
    #    }
    
    # Extract data for current row
    df <- df1[c(i,5),]
    df <- data.frame(t(df))
    colnames(df) <- c("concentration", "time_point")
    df <- df[-c(1,2,3,4),]
    #    if (analyte == "Nh3") {
    #      if (i == 1) {
    #        df <- df[-4,]
    #      }
    #    }
    
    
    # Convert columns to numeric
    df$concentration <- as.numeric(df$concentration)/1000
    df$time_point <- as.numeric(df$time_point)
    
    df <- df[!is.na(df$concentration), ]
    
    model <- nlrob(concentration ~ cubic(time_point, a, b, c, d), 
                   start = list(a = -1, b = -1, c = 1, d = 1), 
                   data = df)
    
    df$predicted_concentration <- predict(model)
    
    #    if(model$status != "converged") {
    
    #      model <- try(nls(concentration ~ cubic(time_point, a, b, c, d), 
    #                       start = list(a = -1, b = -1, c = 1, d = 1), 
    #                       data = df))
    #    }
    
    smry <- summary(model)
    a_met <- smry$coefficients[1,1]
    b_met <- smry$coefficients[2,1]
    c_met <- smry$coefficients[3,1]
    d_met <- smry$coefficients[4,1]
    a_met_err <- smry$coefficients[1,2]
    b_met_err <- smry$coefficients[2,2]
    c_met_err <- smry$coefficients[3,2]
    d_met_err <- smry$coefficients[4,2]
    
    tps <- c(4,24,48,72)
    
    for (x in tps) {
      if (i == 1) {
        state_value = "HP_TR"
        if (x == 4) {
          bm_val <- growth_results[1,3]
          bm_err <- growth_results[1,4]
        }
        if (x == 24) {
          bm_val <- growth_results[2,3]
          bm_err <- growth_results[2,4]
        }
        if (x == 48) {
          bm_val <- growth_results[3,3]
          bm_err <- growth_results[3,4]
        }
        if (x == 72) {
          bm_val <- growth_results[4,3]
          bm_err <- growth_results[4,4]
        }
      }
      if (i == 2) {
        state_value = "HP_MO"
        if (x == 4) {
          bm_val <- growth_results[5,3]
          bm_err <- growth_results[5,4]
        }
        if (x == 24) {
          bm_val <- growth_results[6,3]
          bm_err <- growth_results[6,4]
        }
        if (x == 48) {
          bm_val <- growth_results[7,3]
          bm_err <- growth_results[7,4]
        }
        if (x == 72) {
          bm_val <- growth_results[8,3]
          bm_err <- growth_results[8,4]
        }
      }
      if (i == 3) {
        state_value = "LP_TR"
        if (x == 4) {
          bm_val <- growth_results[9,3]
          bm_err <- growth_results[9,4]
        }
        if (x == 24) {
          bm_val <- growth_results[10,3]
          bm_err <- growth_results[10,4]
        }
        if (x == 48) {
          bm_val <- growth_results[11,3]
          bm_err <- growth_results[11,4]
        }
        if (x == 72) {
          bm_val <- growth_results[12,3]
          bm_err <- growth_results[12,4]
        }
      }
      if (i == 4) {
        state_value = "LP_MO"
        if (x == 4) {
          bm_val <- growth_results[13,3]
          bm_err <- growth_results[13,4]
        }
        if (x == 24) {
          bm_val <- growth_results[14,3]
          bm_err <- growth_results[14,4]
        }
        if (x == 48) {
          bm_val <- growth_results[15,3]
          bm_err <- growth_results[15,4]
        }
        if (x == 72) {
          bm_val <- growth_results[16,3]
          bm_err <- growth_results[16,4]
        }
      }
      
      deriv_val_1 <- cubic_deriv(x, a_met, b_met, c_met, d_met)
      
      if (x == 4) {
        deriv_val <- deriv_val_1 + k_depl * mean(df[1:4,1])
      }
      if (x == 24) {
        deriv_val <- deriv_val_1 + k_depl * mean(df[9:12,1])
      }
      if (x == 48) {
        deriv_val <- deriv_val_1 + k_depl * mean(df[21:24,1])
      }
      if (x == 72) {
        deriv_val <- deriv_val_1 + k_depl * mean(df[33:36,1])
      }
      
      deriv_err <- sqrt((cubic_deriv_parderiv_a(x, a_met, b_met, c_met, d_met) * a_met_err)^2 + (cubic_deriv_parderiv_b(x, a_met, b_met, c_met, d_met) * b_met_err)^2 + (cubic_deriv_parderiv_c(x, a_met, b_met, c_met, d_met)*c_met_err)^2 )
      
      qs_val <- deriv_val/bm_val
      qs_err <- sqrt((1/bm_val*deriv_err)^2 + ((-deriv_val/bm_val^2)*bm_err)^2)
      
      specific_rates <- specific_rates %>%
        bind_rows(data.frame(metabolite = analyte,
                             state = state_value,
                             time_point = x,
                             qs = qs_val,
                             std_err = qs_err,
                             err_perc = qs_err / qs_val,
                             ub = qs_val + qs_err,
                             lb = qs_val - qs_err
                             #                             std_err = qs_err,
                             #                             err_perc = qs_err / qs_val
        ))
    }
  }
}



write_csv(specific_rates, "/home/users/lzehetner/data/hek/specific_exchange_rates.csv")
write_csv(growth_results, "/home/users/lzehetner/data/hek/growth_results.csv")

