#! /bin/bash
# This script removes the MacOS '._' file turds that appear
echo 'Looking for "._*" file turds'
lines=$(find . -name "._*" -type f | wc -l)
if [ $lines -eq 0 ]
then
	echo 'No turds to clean up!'
else
	find . -name "._*" -type f
	# prompt
	read -p "Delete these turds? (y/n) " -n 1 -r
	echo
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
	    find . -name "._*" -type f -delete
	    echo 'Files deleted'
	fi
fi