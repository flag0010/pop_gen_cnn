#!/bin/bash -l
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=6:00:00
#PBS -e error.txt
#PBS -o out.txt

cd /home/mcgaughs/flag0010/FILET/



python extract.data.set.py
