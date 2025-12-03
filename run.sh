#!/bin/bash

record_id="100"
converted_record_path="../converted_mitbih/$record_id"
samp_rate=$(sampfreq "converted_mitbih/$record_id")
output_plots="true"

cd "$(dirname "$0")/src" || exit 1

output=$(octave --no-gui --quiet --eval "Run_Detector(\"$converted_record_path\", $samp_rate, $output_plots)")
echo "$output"
