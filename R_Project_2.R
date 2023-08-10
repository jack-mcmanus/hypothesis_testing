library(ggplot2)
library(dbplyr)
library(rmarkdown)
library(shiny)
library(ISwR)
library(MASS)
library(Sleuth3)
library(tidyr)

##################################################
##################################################

#                 Part 1

##################################################
##################################################

####################################

#                Part I a)

####################################

generate_bernoulli <- function(sample_size, pr, num_dataset){
  
  #initialize array where each row corresponds to 1 dataset
  #e.g. row 5 column 4 corresponds to the 4th draw in the 5th simulation
  array_of_samples <- array(dim = c(num_dataset, sample_size))
  
  for (i in 1:num_dataset){
    array_of_samples[i,] <- rbinom(sample_size, 1, prob=pr)
    
  }
  return(array_of_samples)
}
b_means <- function(data, n){
  #Get number of rows
  rows <- dim(data)[1]
  
  #Initialize array where each row corresponds to 1 dataset
  #e.g. row 5 column 4 corresponds to the 4th draw in the 5th simulation
  array_of_means <- array(dim = c(rows, 1))
  for (i in 1:rows){
    array_of_means[i,] <- (1/(n-1)) * sum(data[i,])
  }

  return(array_of_means)
}

two_sided_z_test <- function(null, alpha, sample_size, prob, num_dataset, test_type){
  #Generate 1000 datasets
  data <-generate_bernoulli(sample_size,prob,num_dataset)
  
  #Compute sample means
  x_bar <- b_means(data, sample_size)
  
  
  #Compute s.d. under assumption population follows Bernoulli distribution with
  #population mean = null = 0.25
  sd <- sqrt(null*(1-null))

  #Compute standard error
  se <- sd/sqrt(sample_size)
  
  #Divide alpha by 2 for two-tailed test
  alpha2 <- alpha/2
  
  #Find z-alpha at given level of alpha = 0.01
  z_alpha <- qnorm(alpha2, lower.tail=FALSE)
  
  #Under the null, we would reject the null if the sample mean < lmb
  #Under the null, we would reject the null if the sample mean > umb
  lmb <- null - (z_alpha * se)
  umb <- null + (z_alpha * se)
  
  
  #Standardized sample mean
  st_means <- (x_bar-null)/se
  
  
  #Initialize array to store reject/accept values
  reject_accept <- array(dim = c(num_dataset, 1))
  
  counter <-0
  #Reject/Fail to reject sample means
  for(i in 1:num_dataset){
    if(x_bar[i,] < lmb || x_bar[i,] > umb){
      reject_accept[i,] <- 1
      counter <- counter + 1
    }else{
      reject_accept[i,] <- 0
    }
  }
  
  #Calculate the proportion of times we reject the null hypothesis
  reject_proportion <- counter/num_dataset

  #Confirm Type I error rate = alpha
  TIER <- pnorm(-z_alpha) + 1 - pnorm(z_alpha)
  
  #calculate Type I error rate
  
  #Under the alternative, we would reject the null if the standarized sample mean < z4b
  #Under the alternative, we would reject the null if the standardized sample mean > z4b2
  z4B <- -z_alpha - ((prob-null)/se)
  z4b2 <- z_alpha - ((prob-null)/se)
  
  power <- pnorm(z4B) + 1 - pnorm(z4b2)
  
  #Compute power
  TIIER <- 1-power
  
  #Combine data, sample means, z-test statistics, p-value, and accept/reject into a single arrray
  combined_array <- cbind(x_bar, st_means, lmb, umb, reject_proportion, TIER, TIIER, power)
  
  #Label the combined array
  #Generate row labels
  row_labels <- c()
  for (i in 1:num_dataset) {
    row_labels[i] <- paste("Dataset ", i)
  }
  
  #Create column labels
  column_labels <- c()
  column_labels[1] <- "Sample Means"
  column_labels[2] <- "Standardized Means"
  column_labels[3] <- "Lower Bound of Acceptance Region"
  column_labels[4] <- "Upper Bound of Acceptance Region"
  column_labels[5] <- "Proportion of Rejections"
  column_labels[6] <- "Type I ER"
  column_labels[7] <- "Old Reject Proportion"
  column_labels[8] <- "Power"
  
  #Assign row and column labels to the array
  rownames(combined_array) <- row_labels
  colnames(combined_array) <- column_labels
  
  #Convert to dataframe
  combined_dataframe <- as.data.frame(combined_array)
  
  return(combined_dataframe)
  
}

one_sided_z_test <- function(null, alpha, sample_size, prob, num_dataset, test_type){
  #Generate 1000 datasets
  data <-generate_bernoulli(sample_size,prob,num_dataset)
  
  #Compute sample means
  x_bar <- b_means(data, sample_size)
  
  
  #Compute s.d. under assumption population follows Bernoulli distribution with
  #population mean = null = 0.25
  sd <- sqrt(null*(1-null))
  
  #Compute standard error
  se <- sd/sqrt(sample_size)
  
  #Find z-alpha at given level of alpha = 0.01
  z_alpha <- qnorm(alpha, lower.tail=TRUE)

  #Under the null, we would reject the null if the sample mean < mb
  mb <- null + (z_alpha * se)
  print(mb)
  
  #Standardized sample mean
  st_means <- (x_bar-null)/se
  
  #Initialize array to store reject/accept values
  reject_accept <- array(dim = c(num_dataset, 1))
  
  counter <-0
  #Reject/Fail to reject sample means
  for(i in 1:num_dataset){
    if(x_bar[i,] <= mb){
      reject_accept[i,] <- 1
      counter <- counter + 1
    }else{
      reject_accept[i,] <- 0
    }
  }
  
  #Calculate the proportion of times we reject the null hypothesis
  reject_proportion <- counter/num_dataset
  
 
  #Calculate Type I error rate
  TIER <- pnorm(z_alpha)
  
  #calculate type 2 error rate
  #Under the alternative, we would reject the null if the standarized sample mean < z4b
  #Under the alternative, we would reject the null if the standardized sample mean > z4b
  z4b <- z_alpha - ((prob-null)/se)

  #Calculate power
  power <- pnorm(z4b)
  
  #Compute Type II Error
  TIIER <- 1-power
  
  
  #Combine data, sample means, z-test statistics, p-value, and accept/reject into a single arrray
  combined_array <- cbind(x_bar, st_means, mb, reject_proportion, TIER, TIIER, power)
  
  #Generate row labels
  row_labels <- c()
  for (i in 1:num_dataset) {
    row_labels[i] <- paste("Dataset ", i)
  }
  
  #Label the combined array
  column_labels <- c()
  column_labels[1] <- "Sample Means"
  column_labels[2] <- "Standardized Means"
  column_labels[3] <- "Bound of Acceptance Region"
  column_labels[4] <- "Proportion of Rejections"
  column_labels[5] <- "Type I ER"
  column_labels[6] <- "Type II ER"
  column_labels[7] <- "Power"
  
  #Assign row and column labels to the array
  rownames(combined_array) <- row_labels
  colnames(combined_array) <- column_labels
  
  #Convert to dataframe
  combined_dataframe <- as.data.frame(combined_array)

  return(combined_dataframe)
  
}


#This function pulls the data out of my samples 
extract_info <-function(d, n, p){
  arr <- array(dim = c(1,6))
  arr[1,1] <- p
  arr[1,2] <- n
  arr[1,3] <- d[1,ncol(d)-3]
  arr[1,4] <- d[1,ncol(d)-2]
  arr[1,5] <- d[1,ncol(d)-1]
  arr[1,6] <- d[1,ncol(d)]
  df <- data.frame(arr)
  return(df)
}


graph_results <- function(d, tt){
    df_list <- split(d, rep(1:5, each = 4))

    #Create a function to generate a plot with a dynamic title
    plot_func <- function(sub_df) {
      # Extract the p_value for the current plot
      p_val <- sub_df[1, 1]

      ggplot(data = sub_df, aes(x = `Num Samples`, y = `Proportion of Rejections`)) +
        geom_point(aes(color = "Proportion of Rejections", label = "Proportion of Rejections"), size = 3, alpha = 0.5) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
        geom_point(aes(x = `Num Samples`, y = `Power`, color = "Power", label = "Power"), size = 3, alpha = 0.5) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
        scale_color_manual(values = c("red", "blue")) +
        labs(title = paste(tt, "p =", p_val),
             x="Number of Samples", y="Value", color = "Lines") +
        theme(plot.title = element_text(hjust=0.5))
    }

    #Use lapply to apply the function to each subset of data
    plot_list <- lapply(df_list, plot_func)

    #Print the plots
    print(plot_list)

    return()
}

n <- c(100, 200, 300, 400) #sample sizes
p = c(0.1, 0.2, 0.25, 0.3, 0.4)#probability of success
datasets = 1000 #number of datasets to generate
null <- 0.25 #p-naught
alpha <- 0.01 #alpha level for two-sided z-test
tt <- "Two-Sided, "

df_info <- data.frame(matrix(nrow = 20, ncol = 6))
colnames(df_info) <- c("P", "Num Samples", "Proportion of Rejections", "Type I ER", "Type II ER", "Power")
counter <- 1
for(j in 1:5){
  for(i in 1:4){
    k <- two_sided_z_test(null, alpha, n[i], p[j], datasets, tt)
    d <- extract_info(k,n[i],p[j])
    df_info[counter,] <- d
    counter <- counter +1
  }
}
options(scipen = 999)
df_info
graph_results(df_info, tt)

####################################

#                Part I b)

####################################

#Define Variables for Bernoulli pulls
n <- c(100, 200, 300, 400) #sample sizes
p = c(0.1, 0.2, 0.25, 0.3, 0.4)#probability of success
datasets = 1000 #number of datasets to generate
null <- 0.25 #p-naught
alpha <- 0.01 #alpha level for z-test
tt <- "One-Sided, "


df_info <- data.frame(matrix(nrow = 20, ncol = 6))
colnames(df_info) <- c("P", "Num Samples", "Proportion of Rejections", "Type I ER", "Type II ER", "Power")
counter <- 1
for(j in 1:5){
  for(i in 1:4){
    k <- one_sided_z_test(null, alpha, n[i], p[j], datasets, tt)
    d <- extract_info(k,n[i],p[j])
    df_info[counter,] <- d
    counter <- counter +1
  }
}
options(scipen = 999)
df_info
graph_results(df_info, tt)


##################################################
##################################################

#                 Part 2

##################################################
##################################################

#Step 1: Read in data
# Construct the path to the CSV file in the same directory as the script
csv_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "numbats.csv")

# Load the CSV data
bad_data <- read.csv(csv_file)

#Step 2: Remove rows where month = N/A
data <- bad_data[!is.na(bad_data$month) & bad_data$dryandra == TRUE,]

#Step 3: Create summer/fall and winter/spring datasets for conveienece
sf_data <- data[data$month %in% c("Dec", "Jan", "Feb", "Mar", "Apr", "May"),]
ws_data <- data[data$month %in% c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov"),]

#Step 3: Define Hypothesis Test
#We need to use the Two-Sample T-Test
alpha <-0.01

#NOTE: Assumption of my test - we observed numbats in either seasons
#in the years 1906, 1954, 2007, 2014, 2015, 2017, 2019, 2020, 2021, 2022, 2023.
#Hence I used the # of years we observed a numbat in either season as the number
#of total seasons for the purpose of calculating the average number of numbats
#per season
num_seasons <- 11

#Compute counts of sightings per season
season_counts <- array(dim = c(num_seasons, 2))
years <- c(1906, 1954, 2007, 2014, 2015, 2017, 2019, 2020, 2021, 2022, 2023)

#Count the number of rows that contain each year
sfcounts <- sapply(years, function(year) sum(sf_data$year == year))
wscounts <- sapply(years, function(year) sum(ws_data$year == year))

#Combine the counts into a matrix
combined_counts <- cbind(sfcounts, wscounts)

#Create a new dataset with the counts and row names
obs_per_season <- data.frame(combined_counts)

#Label Dataframe
rownames(obs_per_season) <- c(1906, 1954, 2007, 2014, 2015, 2017, 2019, 2020, 2021, 2022, 2023)
colnames(obs_per_season) <- c("SF_Obs", "WS_Obs")

#Compute average number of sightings per season
u1 <- mean(obs_per_season$SF_Obs) #average number of numbat sightings in the summer/fall
u2 <- mean(obs_per_season$WS_Obs) #average number of numbat sightings in the winter/spring


#Calculate the sample standard deviation of s1 and s2
s1 <- sd(obs_per_season$SF_Obs)
s2 <- sd(obs_per_season$WS_Obs)

#Calculate pooled estimator
sp = sqrt(((num_seasons-1)*s1^2 + (num_seasons-1)*s2^2)/(num_seasons+num_seasons-2))

#Calculate observed t-stat
tobs = (u1 - u2 - 0)/(sp*sqrt((1/num_seasons) + (1/num_seasons)))

#Get t_alpha
t_alpha <- qt(0.01, 2*num_seasons-2)

#Reject or fail to reject
if(tobs < t_alpha){
  print("We reject the null hypothesis in favor of the alternative at alpha of 0.01.")
  
}else{
  print("We fail to reject the null hypothesis at alpha of 0.01.")
}

#Create and print useful summary table
sum_vect <- c(u1, u2, s1, s2, sp, tobs, t_alpha)
summary_table <- data.frame(sum_vect)
rownames(summary_table) <- c("SF Mean", "WS Mean", "SF Sample s.d", "WS Sample s.d.", "Pooled Sample s.d.", "tobs", "t_alpha")
colnames(summary_table) <- c("Value")
summary_table











