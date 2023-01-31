#!/bin/bash
types=("all")
sizes=(20 40 60 80 100)
cats=("common" "disease-gene" "disease-specific" "normal-specific" "undefined")
for size in ${sizes[@]};
do
for type in ${types[@]};
do
#sbatch --export=size=$size,type=$type /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_separation.sh
for cat in ${cats[@]};
do
sbatch --export=size=$size,type=$type,cat=$cat /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_go_enrichment.sh
sbatch --export=size=$size,type=$type,cat=$cat /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_go_enrichment_coexpr.sh
#sbatch --export=size=$size,type=$type,cat=$cat /home/j.aguirreplans/Projects/Scipher/SampleSize/scripts/calculate_proximity_multiple_drugs.sh
done
done
done
