#!/bin/bash

## provide path to snappy analysis folder and truthfile.tab 
## columns of truth file [1] "expected" name [2] chr [3] pos [4] pos-end(optional) [5] ref [6] alt [7] VAF [8] VD [9] /sample_analysis_folder 
## [10]snappy_id(optional) [11] unique_variant_ID (e.g. chr_pos_sampleID)
## If working on with vcf files from elsewhere (i.e. analysis not on machine)then use this in analysis folder 
## e.g. ls > names \\ for line in names \\ do \\ cp  ${line}/default/variants/*VD.filtered.vcf.gz ${line}.vcf.gz\\ done 
## put these names into truthfile.tab and point to vcf directly in bcftools below (just use ${path}/${fileName})

path=$1
truth=$2
file=$(cat $truth)
break=$(echo "xxxxxxxxxxxxxxxxxxxxx")

## csnappy aligned bcftools version 
## if changing state version in validation and keep consistant through validation work. pr-cds01 version = 1.7(in $PATH already)

# bcftools="/home/alexsmith/biotools/ngs_software/samtools/bcftools-1.5/bcftools" 
# tabix="/home/alexsmith/biotools/ngs_software/samtools/htslib-1.5/tabix"

mkdir truth_extracted_folder
OLDIFS=$IFS; IFS=$'\n'; for line in $file 
do 
    fileName=$(echo $line | awk 'BEGIN{FS="\t"; OFS="\t"} {print $9}' -)    
    echo $fileName
    ID=$(echo $line | awk 'BEGIN{FS="\t"; OFS="\t"} {print $11}' -)
    base=$(basename  $fileName .vcf) 
    echo $base
    variantPos=$(echo $line | awk 'BEGIN{FS="\t"; OFS="\t"} {print $2,$3,$4}' -)
    echo $variantPos >  tempfile2
    # bcftools query   -H -f 'actual\t%CHROM\t%POS\tend\t%REF\t%ALT\t[%AF\t%VD\t%SAMPLE]\n' -R tempfile2  ${path}/${fileName}/default/variants/*VD.vcf.gz --output ${fileName}${ID}"filtered_by_truthPos" 
    bcftools query   -H -f 'actual\t%CHROM\t%POS\tend\t%REF\t%ALT\t[%AD\t%DP\t%SAMPLE\t%GT]\t%FILTER\t%QD\n' -R tempfile2  ${path}/${fileName}.filtered.vcf.gz --output ${fileName}${ID}"filtered_by_truthPos"
    bcftools filter  -T tempfile2  ${path}/${fileName}.filtered.vcf.gz  --output ${fileName}${ID}"filtered_by_truthPos.vcf" -O z
    tabix  ${fileName}${ID}"filtered_by_truthPos.vcf" 
    out=$(cat ${fileName}${ID}"filtered_by_truthPos") 
    out2="$out"$'\n'"$line"$'\n'"$break" 
    echo "$out2" >${fileName}${ID}"filtered_by_truthPos_format"
    rm tempfile2
    mv ${fileName}${ID}"filtered_by_truthPos" truth_extracted_folder/
    mv ${fileName}${ID}"filtered_by_truthPos_format" truth_extracted_folder/
    mv ${fileName}${ID}"filtered_by_truthPos.vcf"* truth_extracted_folder/
# ~/biotools/samtools_bcf_htslib_1_5/htslib-1.5/tabix  $base"filtered_truthset.vcf"
done	
cat truth_extracted_folder/*filtered_by_truthPos_format > truth_extracted_folder/cat_all_extracted_by_truth
bcftools merge -m none --force-samples truth_extracted_folder/*filtered_by_truthPos.vcf > truth_merged_truth_extracted.vcf
