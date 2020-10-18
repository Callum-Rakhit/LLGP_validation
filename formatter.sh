# Get information about total reads and % ROI covered, plus headers

for i in $(ls /home/callum/LLGP_validation/LLGP_metrics_files/*.metrics.exoncoverage); 
do grep -A 2 "# COVERAGE BINS" $i | tail -1 >> exoncoverage.summary;
done

for i in $(ls /home/callum/LLGP_validation/LLGP_metrics_files_plus/*.exoncoverage); 
do grep -A 2 "# COVERAGE BINS" $i | tail -1 >> exoncoverage.summaryplus;
done

sed -i -r 's/\S+//19' exoncoverage.summaryplus 

for i in $(ls /home/callum/LLGP_validation/LLGP_metrics_files/*.metrics.READS); 
do grep -A 2 "# BASIC FILTERING METRICS" $i | tail -1 >> READS.summary;
done

# Convert names to lower case
find . -name '*.*' -exec sh -c 'a=$(echo "$0" | sed -r "s/([^.]*)\$/\L\1/"); [ "$a" != "$0" ] && mv "$0" "$a"' {} \;

for i in $(ls /home/callum/LLGP_validation/LLGP_metrics_files_plus/*.reads); 
do grep -A 2 "# BASIC FILTERING METRICS" $i | tail -1 >> READS.summaryplus;
done

# Get the duplication rates from the .json
for i in $(ls /home/callum/LLGP_validation/LLGP_metrics_files/*.json); 
do grep -A 1 "PERCENT_DUPLICATION" $i | awk '{print $2}' | cut -d \, -f 1 | awk -F '"TAS"' '{print $1}' | awk 'NF' >> duplication.summary;
done

# Get the duplication rates from the .json and subsamples
for i in $(ls /home/callum/LLGP_validation/LLGP_metrics_files_plus/*.json); 
do grep -A 1 "PERCENT_DUPLICATION" $i | awk '{print $2}' | cut -d \, -f 1 | awk -F '"TAS"' '{print $1}' | awk 'NF' >> duplication.summaryplus;
done
