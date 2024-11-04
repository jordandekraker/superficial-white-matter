#!/bin/bash

# Define the path to the Python script
SCRIPT_PATH="/home/bic/igoodall/Desktop/blur_calc/run_blur.py"

# Check if the script exists
if [ ! -f "$SCRIPT_PATH" ]; then
  echo "Error: $SCRIPT_PATH not found!"
  exit 1
fi
# Run the Python script
/data/mica1/03_projects/ian/anaconda3/envs/zbrains/bin/python3 "$SCRIPT_PATH"