#!/bin/bash

# --------------------------
# Configuration
# --------------------------
INPUT_FILE="runs_on_tape.txt"    # Your file with full /mss/... paths
EMAIL="aleksandr.bulgakov1999@gmail.com"                         # Optional: your_email@jlab.org
TO_STAGE="to_stage.txt"
LOG="jcache_staging.log"

# --------------------------
# 1. Verify input
# --------------------------
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: File $INPUT_FILE not found."
    exit 1
fi

# --------------------------
# 2. Check which files are NOT in /cache/
# --------------------------
echo "Checking cache status..."
> "$TO_STAGE"

while read -r mss_file; do
    cache_file="${mss_file/mss/cache}"
    if [ ! -f "$cache_file" ]; then
        echo "$mss_file" >> "$TO_STAGE"
    fi
done < "$INPUT_FILE"

# --------------------------
# 3. Run jcache get only for missing files
# --------------------------
num_to_stage=$(wc -l < "$TO_STAGE")
if [ "$num_to_stage" -eq 0 ]; then
    echo "✅ All files are already in /cache/. Nothing to stage."
    rm "$TO_STAGE"
    exit 0
fi

echo "⏳ Staging $num_to_stage files with jcache..."
while read -r filepath; do
    if [ -n "$EMAIL" ]; then
        jcache get "$filepath" -e "$EMAIL" >> "$LOG" 2>&1
    else
        jcache get "$filepath" >> "$LOG" 2>&1
    fi
done < "$TO_STAGE"

echo "✅ Done staging. See log: $LOG"
rm "$TO_STAGE"
