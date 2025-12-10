# BiomedicalSignalProcessing-Project1

A QRS complex detector, written in Octave, implementing the algorithm described in the paper:
Ligtenberg, A., & Kunt, M. (1983). *A robust-digital QRS-detection algorithm for arrhythmia monitoring*.  
Computers and Biomedical Research, 16(3), 273–286.  


## Directory Structure
All record files must be placed in their dedicated `converted_*` database directories.  
This includes the following file types:

- `.hea`
- `.dat`
- `.atr`
- `m.hea`
- `m.mat`

## Running the Detector
To run the detector on all records:

1. Open `det_run_all.sh`  
2. Update the paths inside the script to match your directory layout.

## Running Evaluation
To evaluate the detector outputs:

1. Open `det_eval_all.sh`  
2. Update the paths inside the script as needed.

Both scripts expect the database directories (`converted_*`) to contain the relevant record files.

## QRS directories
The directories `mitbih_improved_qrs` and `mitbih_original_qrs` contain all generated binary annotation files of our improved and original algorithm implementations.