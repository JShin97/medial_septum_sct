#!/bin/bash

# ==============================================================================
# SCRIPT CONFIGURATION
# ==============================================================================

# 1. Define the target directory where the processed files will be placed.
TARGET_DIR="/u/scratch/j/jshin/medial_septum_sct/GEOSubmission" 

# 2. List of the three base file names to move.
FILES_TO_MOVE=(
    "features.tsv.gz"
    "barcodes.tsv.gz"
    "matrix.mtx.gz"
)

# Define the range of directories to process (from 1 to 24)
START_NUM=2
END_NUM=24

# The static part of the nested file path within each TriSeptum directory.
NESTED_PATH="outs/filtered_feature_bc_matrix"

# ==============================================================================
# EXECUTION
# ==============================================================================

echo "Starting scRNA-seq file consolidation..."

# Create the target directory if it does not exist
mkdir -p "$TARGET_DIR"

if [ $? -ne 0 ]; then
    echo "Error: Could not create target directory '$TARGET_DIR'. Exiting."
    exit 1
fi

echo "Target directory created/verified: $TARGET_DIR"

cd /u/scratch/j/jshin/medial_septum_sct/Data/Raw/Counts

# Loop through the directory indices (from START_NUM to END_NUM)
for i in $(seq $START_NUM $END_NUM); do
    DIR_NAME="TriSeptum-$i"
    SOURCE_BASE_DIR="$DIR_NAME/$NESTED_PATH"

    echo "--- Processing $DIR_NAME (Directory $i of $END_NUM) ---"

    # Check if the main subdirectory exists before proceeding
    if [ ! -d "$DIR_NAME" ]; then
        echo "Warning: Directory '$DIR_NAME' not found. Skipping."
        continue
    fi

    # Loop through each of the three specific files
    ALL_FILES_MOVED=true
    for FILE_BASE in "${FILES_TO_MOVE[@]}"; do
        SOURCE_FILE="$SOURCE_BASE_DIR/$FILE_BASE"
        NEW_FILE_NAME="${DIR_NAME}_${FILE_BASE}"
        DEST_FILE="$TARGET_DIR/$NEW_FILE_NAME"

        # Check if the source file exists in the expected location
        if [ -f "$SOURCE_FILE" ]; then
            # Move and rename the file in one command
            mv "$SOURCE_FILE" "$DEST_FILE"
            echo "  -> Moved: $FILE_BASE to $NEW_FILE_NAME"
        else
            echo "  -> Warning: File not found at '$SOURCE_FILE'. This file was skipped."
            ALL_FILES_MOVED=false
        fi
    done

    if $ALL_FILES_MOVED; then
        echo "Status: All 3 files successfully moved for $DIR_NAME."
    else
        echo "Status: Finished processing $DIR_NAME with missing files."
    fi

done

echo "==================================================="
echo "Consolidation complete. All available files are in '$TARGET_DIR'."
