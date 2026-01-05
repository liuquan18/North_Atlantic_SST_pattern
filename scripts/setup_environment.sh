#!/bin/bash
# Setup script for Levante HPC system
# This script creates the conda environment for the North Atlantic SST pattern analysis

set -e  # Exit on error

echo "======================================"
echo "Setting up North Atlantic SST Pattern Analysis Environment"
echo "======================================"

# Check if we're on Levante
if [[ $(hostname) == *"levante"* ]]; then
    echo "Running on Levante..."
    
    # Load conda module (adjust if needed based on Levante's module system)
    module load python3/unstable
    module load openmpi
else
    echo "Warning: Not running on Levante. Proceeding anyway..."
fi

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo ""
echo "Project root: $PROJECT_ROOT"
echo ""

# Check if environment already exists
if conda env list | grep -q "nalt"; then
    echo "Environment 'nalt' already exists."
    read -p "Do you want to remove and recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n nalt -y
    else
        echo "Keeping existing environment. Use 'conda activate nalt' to activate."
        exit 0
    fi
fi

# Create conda environment from environment.yml
echo ""
echo "Creating conda environment from environment.yml..."
conda env create -f "$PROJECT_ROOT/environment.yml"

echo ""
echo "======================================"
echo "Setup complete!"
echo "======================================"
echo ""
echo "To activate the environment, run:"
echo "  conda activate north_atlantic_sst"
echo ""
echo "To deactivate, run:"
echo "  conda deactivate"
echo ""
