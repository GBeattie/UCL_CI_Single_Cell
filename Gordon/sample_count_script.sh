# Batch script to run a serial job on Legion with the upgraded
# software stack under SGE.
# 1. Force bash as the executing shell.
#$ -S /bin/bash
# 2. Request 2 hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=48:00:00
## 48 hours is the max, but if samples aren't finished running simply delete the _lock file in the run folder and run again, cellranger will stop from where it left off.
# 3. Request 1 gigabyte(s) of RAM
#$ -l mem=8G
# 4. Request 10 gigabyte(s) of TMPDIR space (default is 10 GB)
#$ -l tmpfs=20G
# 5. Set the name of the job.
#$ -N count_log
#$ -pe smp 12
# 6. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
#$ -wd /home/##USERID/Scratch/ ## set working directory
# 7. Your work *must* be done in $TMPDIR
#cd $TMPDIR

# 8. Run the application.

module load cellranger/6.0.1

## Variables
FASTQS= ## directory where all your fastq files are (put all in one folder, not in sub-directories)
REF=/home/regmgbe/Programmes/refdata-cellranger-mm10-3.0.0 ## path to folder containing desired cellranger reference, mm10 shown here
SAMPLES=$(ls $FASTQS | grep -oE '.*_S' | uniq | sed 's/.\{2\}$//') ## This is a grep for sample IDs, will work if fastq names formatted properly (usually are)

for s in $SAMPLES; do
	echo; echo $s; echo
	CMD="cellranger count --id=$s --sample=$s --fastqs=$FASTQS --transcriptome=$REF --localmem=96 --localcores=12"
	echo; echo $CMD; echo
	$CMD
done
exit
