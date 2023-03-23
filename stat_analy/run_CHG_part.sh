#!/bin/bash     

for i in {1..16}
do 
	#Rscript
	Rscript meth_stat_part.R part_data/CHG_p_$i.RData part_data/CHG_dmc_$i.RData &  
	sleep 2s 
done 	