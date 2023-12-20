#!/bin/bash

for i in {1..15}
do
	#Rscript
	Rscript meth_stat_part.R meth_p_$i.RData dmc_p_$i.RData &  
done	
