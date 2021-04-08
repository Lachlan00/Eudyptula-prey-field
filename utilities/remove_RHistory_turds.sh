#! /bin/bash
# This script removes the '.Rhistory' and '.Rapp.history' file turds that appear
# Clean RHistory
echo 'Looking for ".Rhistory" file turds'
lines=$(find . -name "*.Rhistory" -type f | wc -l)
if [ $lines -eq 0 ]
then
	echo 'No turds to clean up!'
else
	find . -name "*.Rhistory" -type f
	# prompt
	read -p "Delete these turds? (y/n) " -n 1 -r
	echo
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
	    find . -name "*.Rhistory" -type f -delete
	    echo 'Files deleted'
	fi
fi
# Clean Rapp.history
echo 'Looking for ".Rapp.history" file turds'
lines=$(find . -name "*.Rapp.history" -type f | wc -l)
if [ $lines -eq 0 ]
then
	echo 'No turds to clean up!'
else
	find . -name "*.Rapp.history" -type f
	# prompt
	read -p "Delete these turds? (y/n) " -n 1 -r
	echo
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
	    find . -name "*.Rapp.history" -type f -delete
	    echo 'Files deleted'
	fi
fi