#!/bin/bash
## this script is for testing sampled down experiment for number of variant supporting reads of all variants in outputed snappy vcf. edited and updated nov27 2019. This is used for comparing depth to vaf as well as VD to VAF, and for GATK piupelines, depth vs halotype caller outpuuteed PASS variants. Adjust data extracted from VCF as appropriate. For R plotting, Just point R to correct colums of output. 
## provide merged and filtered vcf from subsampled experiments....merging all vcfs and then isec against origianl concordant variant calls to remove background(see lod_routine.readme) could replace with merge at higher snappy cut variant filtration cut-offs and simply filter/extract by coordinates of truth set. However, tend to miss some of low VAF variants dur to subsampling error.
## uses snappy bcftools/samtools/bedtools versions
## uses variant filter config in snappy analysis with very low thresholds to capture all initially ( so nothing removed to to subsampling error) 
##BEWARE multiallieic calls which pulling out corresponding depth from bedgraph as bedtools will only pull for once per position. Isec filter will remove any aditional multi-allieic calls as variant calling becomes less reliable in low subsampled samples.ted depth as depth for same repeated position in bed file   
##original non-sampled vcf output of snappy (with regular variant filtration) giving expected values for  AF (in column 6 here)<original_calls> and bed file for these calls too <original_bed>. Need to: sort -n -k1 -k2 on these files. Also filter by original_filtered.bed Such sample specific events can be removed at the end manually if under thesholds of 3% VAF and <6VD and/or not found in original sample  calls. 

## Below is for  GATK halotype caller pipeline.

original_concordant_vcf=$1 # selected_filtered_calls_from_panel_being_tested_non_subsampled
merged_filtered_vcf=$2 # filtered_merged_subsampled_calls_isec_with_original_non_subsampled_calls_above
bedgraph_dir=$3 # bedgraphs of subsampled experiments with same basename as sample as named in vcf sample tag
outname=$4
# bgzip ${merged_filtered_vcf}  
# tabix ${merged_filtered_vcf}.gz
bcftools query ${merged_filtered_vcf}  -l  > sample_names.list
samples=$(cat sample_names.list)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\t[%PL\t]\t[%DP\t]\t[%GT\t]\n'  ${original_concordant_vcf} | sort -n -k1 -k2 -> original_calls
# awk 'BEGIN{FS="\t"; OFS="\t"} {sum=$10+$11+$12+$13} {print sum / 4}' original_calls> expectedVAF ## expected vaf is average of repeated experiments
awk 'BEGIN{FS="\t"; OFS="\t"} {print $9}' original_calls> original_representative_depth # to order with in plotting
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $2+1}' original_calls | sort -n -k1 -k2 -> original.bed
echo $samples
# samplestest="700_HD798"
mkdir depth
for line in $samples
do 
    echo $line
    bcftools view -s $line -O z ${merged_filtered_vcf}  >  ${line}.vcf.gz
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%PL]\t[%DP]\t[%GT]\n'  ${line}.vcf.gz > ${line}bcftools_temp.txt
    awk -v var="$line" 'BEGIN{FS="\t"; OFS="\t"} {print $0, var}' ${line}bcftools_temp.txt  | tail -n +1 - | sort -n -k1 -k2 -> ${line}bcftools_filter2.txt
    ls ${bedgraph_dir}${line}.coverage.bedgraph > $line.name
    bedgraph=$(cat  $line.name)    
    ###generate depth at sites of interest to use to calculate expected VD by multiplying VAF from base sample with depth extracted here
    bedtools intersect -wao -b ${bedgraph} -a original.bed > depth/$line.bedgraph_depth
    awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$2,$7}' depth/$line.bedgraph_depth | sort -n -k1 -k2  > depth/${line}expected.depth
    paste original_representative_depth depth/${line}expected.depth ${line}bcftools_filter2.txt > ${line}bcftools_filtered_temp.txt
    awk 'BEGIN{FS="\t"; OFS="\t"} {print $0, $4/$1}' ${line}bcftools_filtered_temp.txt > ${line}bcftools_filtered_final.txt
    rm ${line}bcftools_filter2.txt
    rm ${line}bcftools_*temp.txt
    rm ${line}.name
    rm ${line}.vcf.gz
done

#####cat all outputs, repalce '.' with \t , and where 0 in VD column added expected VD in new column
cat *filtered_final.txt > ${outname}_sampled_down_all_out_tidy_data
awk 'BEGIN{FS="\t"; OFS="\t"} {if ($11=="."){$15=$4} print}' ${outname}_sampled_down_all_out_tidy_data  > ${outname}_sampled_down_all_out_tidy_data_format1 
sed -e s'/\.\s/\t/g' ${outname}_sampled_down_all_out_tidy_data_format1 > ${outname}_sampled_down_all_out_tidy_data_format2

# awk 'BEGIN{FS="\t"; OFS="\t"} {if ($9=="." || $9<6 || $10<0.02){$14=$13} print}' HD798_sampled_down_all_out_tidy_data  | awk 'BEGIN{FS="\t"; OFS="\t"} {if ($9>=6&&$10>=0.02){$15=$9} print}' -  | awk 'BEGIN{FS="\t"; OFS="\t"} {if ($9=="." || $10<0.02 || $9<6){$16=$4} print}' -  | awk 'BEGIN{FS="\t"; OFS="\t"} {if ($10>=0.02&&$9>=6){$17=$11} print}' -  > HD798_sampled_down_all_out_tidy_data_format1
# sed -e 's/\.\s/\t/g' HD798_sampled_down_all_out_tidy_data_format1 > HD798_sampled_down_all_out_tidy_data_format2
