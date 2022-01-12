#!/usr/bin/bash
#input file name in command line
#SBATCH --job-name=ncbi-down      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fanzhen@ufl.edu    # Where to send mail.  Set this to your email address
#SBATCH --account=seonghee
#SBATCH --qos=seonghee-b
#SBATCH --ntasks=1                  # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=4            # Number of cores per MPI task
#SBATCH --nodes=1                    # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=1         # Maximum number of tasks on each node
#SBATCH --mem=10gb          # Memory (i.e. RAM) per processor
#SBATCH --time=3-00:00:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=mpi_test_%j.log     # Path to the standard 

file=$1
ml sra

cat $file | while read line; do 
	fastq-dump --split-files $line -O raw_read
done
