#!/bin/bash
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

# Manual parsing for long options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_data) input_data="$2"; shift ;;
        --labels) labels="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check that required arguments are provided
if [[ -z "$input_data" || -z "$labels" ]]; then
    echo "Usage: $0 --input_data INPUT --labels OUTPUT"
    exit 1
fi

# Determine if the input file is gzipped
if file "$input_data" | grep -q 'gzip compressed'; then
    reader="zcat"
else
    reader="cat"
fi

# Run the pipeline and write output
$reader "$input_data" | awk 'BEGIN {OFS="\t"} {if ($1 ~ /neg/) {print 0} else {print 1}}' > "$labels"
