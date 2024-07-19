#!/usr/bin/bash

# A text file for storing unique ids.
IDFILE='data/ids.txt'

# Create an empty text file.
touch ${IDFILE}

function generate_ids_from_reads() {
    READDIR=${1:-data/reads}
    # List all files in the reads directory.
    ls -l ${READDIR} | \
    # Extract the last column containing filenames.
    awk '{print $9}' | \
    # Remove file format suffix.
    sed 's/_R[12].fq//' | \
    # Remove duplicate sample names.
    uniq | \
    # Strip out empty spaces.
    grep "\S"
}

# Extract unique read ids.
if [ -f 'data/design.csv' ]; then
    # Retrieve from design file
    tail -n +2 data/design.csv | csvcut -c 1 > ${IDFILE}
else
    # Retrieve for names of read files.
    generate_ids_from_reads > ${IDFILE}
fi