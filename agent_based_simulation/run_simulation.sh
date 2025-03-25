#!/bin/bash
#SBATCH -J agent_based
#SBATCH -o sim.out
#SBATCH -e sim.err
#SBATCH -c 32
#SBATCH -t 3-00:00:00
#SBATCH --mem=32G

cargo run --release