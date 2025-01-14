#!/bin/bash

# Iterate through all subdirectories recursively
find . -type d | while read -r dir; do
    # Create a .gitkeep file in each directory if it doesn't already exist
    if [ ! -f "$dir/.gitkeep" ]; then
        touch "$dir/.gitkeep"
        echo "Created .gitkeep in $dir"
    fi
done

echo "Finished creating .gitkeep files in all subdirectories."