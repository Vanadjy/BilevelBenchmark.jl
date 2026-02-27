#!/bin/bash

echo "Running Julia script..."
# Run the Julia script. Replace myscript.jl with your file name.
echo "Julia script finished."

# Define the path to the Julia script
JULIA_SCRIPT="./numerics/ProfilesHub.jl"
ARGUMENT1="log" # or "dec"
ARGUMENT2="All_Final" # or "All_All" or "All_Backward"

# Execute the Julia script
julia $JULIA_SCRIPT $ARGUMENT1 $ARGUMENT2