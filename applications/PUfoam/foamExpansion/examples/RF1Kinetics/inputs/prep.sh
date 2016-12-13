#!/bin/bash

# this is for post-processing after using probeLocation

# remove the lines starting with #
for file in $(ls *); 
	do     
		sed '/^#/d' $file > $file.txt; 
	done
# rename the files
# for file in *.txt;
# 	do
# 		mv "$file" "${file/.txt/v1.txt}"
# 	done
# clean up directory
ls > list
egrep -v '*txt|prep.sh' list > list2
mv list2 list
rm -rf $(<list)
rm -f prep.sh.txt
echo "Done!"

