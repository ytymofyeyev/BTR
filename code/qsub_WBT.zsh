#!/bin/zsh
#
# request Bourne shell as shell for job
#$ -S /bin/zsh

#$ -N jobs

#$ -cwd

#$ -o log
#$ -e log

#$ -t 1-50

# $ -m be
# $ -M ytymofye@its.jnj.com

echo $(hostname) $SGE_TASK_ID $date

cat exeWBT.R | R --vanilla -q

exit 0;
