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

# Path to the alignment file
file <- "alignment.txt"

# Read the alignment from the file
alignment <- read_alignment(file)
print(alignment)
# Compute the SFS
SFS <- compute_SFS(alignment)

# Output the SFS
print(SFS)

# The function for computeing TajimasD values
TajimasD <- function(sfs){
    #For the computed SFS
    n <- length(sfs) + 1
    ss <- sum(sfs)
    
    a1 <- sum(1 / seq_len(n-1))
    a2 <- sum(1 / seq_len(n-1)^2)
    
    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2 / a1^2
    
    e1 <- c1 / a1
    e2 <- c2 / (a1^2 + a2)
    
    Vd <- e1 * ss + e2 * ss * (ss - 1) 
    
    theta_pi <- sum(2 * seq_len(n-1) * (n - seq_len(n-1)) * sfs)/(n*(n-1))
    theta_w <- ss / a1
    res <- (theta_pi - theta_w) / sqrt(Vd)
    return(res)
}

#Compute TajimaD
TajimasD(SFS)

