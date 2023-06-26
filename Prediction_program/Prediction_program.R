# Load required libraries
library(rJava)
library(rjson)
library(RCurl)
library(randomForest)

# Function to convert SMILES to PaDEL descriptors
smiles_to_descriptors <- function(smiles) {
  # Set up the URL and parameters for the PaDEL web service
  url <- "http://padel.nus.edu.sg/cgi-bin/web.cgi"
  params <- list(smiles = smiles, selection = "true", nbits = "1024", fingerprints = "fingerprinter")
  
  # Send request to the PaDEL web service and get the response
  response <- postForm(url, .params = params)
  
  # Extract the descriptor data from the response
  json_data <- fromJSON(response)
  descriptors <- json_data$data
  
  # Convert descriptors to a numeric vector
  descriptors <- as.numeric(descriptors)
  
  return(descriptors)
}

# Example SMILES strings
smiles_strings <- c("SMILES", "SMILES")

# Calculate descriptors for each SMILES string
descriptors <- sapply(smiles_strings, smiles_to_descriptors)

# Load the pre-trained RF model
load("random_forest_model.RDa")

# Predict using the RF model
predictions <- predict(random_forest_model, as.data.frame(t(descriptors)))

# Print the predictions
for (i in 1:length(smiles_strings)) {
  cat("SMILES:", smiles_strings[i], "\n")
  cat("Prediction:", predictions[i], "\n\n")
}
