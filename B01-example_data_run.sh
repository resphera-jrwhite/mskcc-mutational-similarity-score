#!/bin/bash
echo "Begin Example Runs...";
for i in `ls supporting-materials/example_data/*Tumor*`
do 
	for j in `ls supporting-materials/example_data/*`
	do
		sleep 0.5;
		printf "\n\n**********************************************************\n\nComparing $i to $j...\n";
		./A01-mutational-profile-similarity.pl -i $i -j $j -p ASX -t LUAD
	done
done
printf "\n\n**********************************************************\n\nExample Runs Complete...\n";
exit 0 

