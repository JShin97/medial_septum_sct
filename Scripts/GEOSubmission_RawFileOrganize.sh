CRIPT CONFIGURATION
# ==============================================================================

# 1. Define the target directory where the processed files will be placed.
TARGET_DIR="/u/scratch/j/jshin/medial_septum_sct/GEOSubmission"

# Define the range of directories to process (from 1 to 24)
START_NUM=2
END_NUM=24

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

cd /u/scratch/j/jshin/medial_septum_sct/Data/FASTQ

for i in $(seq $START_NUM $END_NUM); do
	DIR_NAME="XY-TriSeptum-$i"

	echo "--- Processing $DIR_NAME (Directory $i of $END_NUM) ---"

	# Use 'find' to locate all files matching the pattern, then execute 'mv'
	# -type f: Ensures only regular files are matched (not directories)
	# -name "*.fastq.gz": Specifies the file extension pattern to match
	# -exec mv {} "$TARGET_DIR" \;: Executes the 'mv' command on each found file ({})
	#                             and moves it to the target directory.
	find "$DIR_NAME" -type f -iname "*.fastq.gz" -exec mv {} "$TARGET_DIR" \;

    # Check the exit status of the find/mv operation for this specific directory
    if [ $? -eq 0 ]; then
        echo "Successfully moved all found '.fastq.gz' files from $DIR_NAME."
    else
        # Note: 'find' will often return 0 even if no files were found.
        # A non-zero return usually indicates a permission or serious I/O error.
        echo "An error occurred during the move operation for $DIR_NAME. Check permissions."
    fi

done

echo "==================================================="
echo "All processing steps finished. Review the output for any warnings."
echo "Final consolidated files are in: $TARGET_DIR"
