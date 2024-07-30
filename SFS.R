# Function to read alignment from a text file
read_alignment <- function(file) {
    # Read the file line by line
    lines <- readLines(file)
    
    # Split each line into individual characters to form the alignment matrix
    alignment <- do.call(rbind, strsplit(lines, ""))
    
    return(alignment)
}

# Function to compute site frequency spectrum (SFS)
compute_SFS <- function(alignment) {
    # Number of sequences
    num_sequences <- nrow(alignment)
    
    # Initialize an empty vector to store counts for each frequency
    sfs <- numeric(num_sequences - 1)
    
    # Loop through each column in the alignment
    for (i in 1:ncol(alignment)) {
        # Count the number of unique characters in the column
        counts <- table(alignment[, i])
        
        # Remove gaps "-" from the counts
        counts <- counts[!names(counts) %in% "-"]
        
        # Increment the count for each frequency in the SFS
        for (count in counts) {
            if (count > 1 && count <= num_sequences) {
                sfs[count - 1] <- sfs[count - 1] + 1
            }
        }
    }
    
    return(sfs)
}

# Path to the alignment file, change with your file in matrix array format, you can use the py scripts in the repository
file <- "alignment.txt"

# Read the alignment from the file
alignment <- read_alignment(file)
print(alignment)
# Compute the SFS
SFS <- compute_SFS(alignment)

# Output the SFS
print(SFS)
