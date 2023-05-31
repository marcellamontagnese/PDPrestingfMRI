# This code creates a chord diagram of the NBS subnetwork.

# Load the necessary libraries.
library(readxl)
library(circlize)

# Read the data.
df <- read_excel("dataset_NBS.xlsx")

# Transform the input data into an adjacency matrix.
adjacencyData <- with(df, table(origin, destination))
# Create the chord diagram.
chordDiagram(adjacencyData, transparency = 0.4)
