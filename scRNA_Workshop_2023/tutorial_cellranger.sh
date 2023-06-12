#!/bin/bash -l

# Batch script to run an OpenMP threaded job under SGE.

# Request wallclock time (format hours:minutes:seconds).
#$ -l h_rt=06:00:00

# Request RAM for each core/thread 
# (must be an integer followed by M, G, or T)
#$ -l mem=4G

# Request TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set the name of the job.
#$ -N tutorial_cellranger

# Request 8 cores.
#$ -pe smp 8

# Set the working directory to somewhere in your scratch space.
# Replace "<your_UCL_id>" with your UCL user ID
#$ -wd /home/regmgbe/Scratch/Tutorial

# Run the application.
module load cellranger/6.0.1

cellranger count --id=heart \
                   --transcriptome=refdata-gex-mm10-2020-A \
                   --fastqs=/home/regmgbe/Scratch/Tutorial/fastqs/heart_1k_v3_fastqs \
                   --sample=heart_1k_v3 \
                   --expect-cells=1000 \
                   --localcores=8 \
                   --localmem=32

cellranger count --id=neuron \
                   --transcriptome=refdata-gex-mm10-2020-A \
                   --fastqs=/home/regmgbe/Scratch/Tutorial/fastqs/neuron_1k_v3_fastqs \
                   --sample=neuron_1k_v3 \
                   --expect-cells=1000 \
                   --localcores=8 \
                   --localmem=32
