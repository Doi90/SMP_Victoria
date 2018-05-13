#!/bin/bash

for taxon in {1..5}

  sbatch ~/../../data/cephfs/punim0440/scripts/slurm/benefit_master.slurm $taxon