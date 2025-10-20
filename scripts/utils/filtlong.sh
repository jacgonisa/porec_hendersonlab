#!/bin/bash

#Need to import the singularity module 
#Run this in the nextflow environment 
/home/jg2070/Filtlong/bin/filtlong --min_length 1000 --min_mean_q 90 $1 > $2
