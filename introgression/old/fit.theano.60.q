#!/bin/bash -l
#PBS -l nodes=1:ppn=12,mem=16gb,walltime=7:00:00
#PBS -e error.txt
#PBS -o out.txt

cd /home/mcgaughs/flag0010/FILET/

python fit.theano.60.py
