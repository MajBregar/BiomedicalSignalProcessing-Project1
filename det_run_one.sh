#!/bin/bash

record_id="${1:-100}"
database_dir="${2:-converted_mitbih}"
record_path="$database_dir/$record_id"

# Sampling rate from header BEFORE cd
samp_rate=$(sampfreq "$record_path")
output_plots="true"

#run from src folder
cd "$(dirname "$0")/src" || exit 1
octave --no-gui --quiet --eval \
    "Run_Detector(\"../$record_path\", $samp_rate, $output_plots)"

