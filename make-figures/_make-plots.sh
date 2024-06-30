for script in $(ls *.py); do
	echo "python $script"
    python $script 2>&1 /dev/null
done
