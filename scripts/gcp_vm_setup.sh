#!/usr/bin/env bash
# Set up a CPU-only Linux VM (Debian or Ubuntu, x86-64) to run the Track 2 analysis. No GPU, no secrets, no open ports.
#
#   bash scripts/gcp_vm_setup.sh                 # installs system packages, a virtual environment and the pinned dependencies
#   source .venv/bin/activate && export PYTHONPATH=src
#   python -m track2.run_all --list
#
# Suggested VM: 8 vCPU and 32 GB of memory (the knowledge graph and the network null model are the memory-heavy steps), 100 GB of disk.
# Nothing here creates cloud resources: create the VM yourself and run this script on it.
set -euo pipefail

sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-pip git build-essential curl unzip
python3 --version

python3 -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements-track2.lock
# optional (docking stage only): pip install -r requirements-docking.txt

export PYTHONPATH=src
python -m pytest -q
python -m track2.run_all --list
echo "Ready. Downloads (ChEMBL, STRING, Reactome, DepMap, ...) are made by the stages and recorded with their checksums."
