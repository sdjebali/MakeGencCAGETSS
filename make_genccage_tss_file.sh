#!/bin/bash

# make_genccage_tss_file.sh
# makes a merged tss file from gencode tss and cage clusters, providing information about the original elements
# each merged tss is made of

# This script takes as input:
#############################
# - a stranded gencode tss file in gff format (already filtered for quality of the tss and with no header)
# - a stranded cage cluster file in bed format (already filtered for quality and with no header)
# This script provides as output in the working directory and indexed by the inputs' names:
# - a genccage tss file in gff format with information about the initial elements it is made of and in particular
#  (contrarily to previous approached) the coord of the initial cage clusters from which each genccage element is made
# - some statistics in standard output and some log information in standard error 


# Check that all the necesary arguments are present, otherwise exit
###################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: make_genccage_tss_file.sh genctss.gff cageclus.bed [keep_interm_files] >&2
    echo where: >&2
    echo - genctss.gff is a gff version 2 file of Gencode TSS obtained with the script make_TSS_file_from_annotation_with_confidence_better.sh >&2
    echo - cageclus.bed is a bed file of cage clusters \(\for example HMM TSS predictions\) >&2
    echo - keep_interm_files is an optional non empty string for keeping intermediate files >&2
    echo Will output a gff version 2 file of gencode + cage tss clusters with information about initial elements merged >&2
    echo as well as some statistics \in standard output and some log information \in standard error >&2
    echo "" >&2
    exit 1
fi


# Variables from input
######################
genctss=$1
cageclus=$2
keep=$3
b=`basename ${genctss%.gff}`
b2=${b%.gtf} 
c=`basename ${cageclus%.bed}`

# Programs
##########
GFF2GFF=Awk/gff2gff.awk
MERGE=mergeBed
INTERSECT=intersectBed
OVERLAP=bin/overlap

# 0. Make a gff file out of the cage cluster bed file
#####################################################
echo I am making a gff version 2 file for the cage clusters >&2
nb1=`wc -l $cageclus | awk '{print $1}'`
echo Initial number of cage clusters is $nb1
awk '{print $1, "rik", "cageclus", $2+1, $3, ".", $6, ".", "id", $1"_"$2"_"$3"_"$6}' $cageclus | awk -f $GFF2GFF > $c.gff
echo done >&2

# 1. Extend the gencode tss by 50bp on each side
################################################
echo I am extending the gencode tss by 50bp on each side >&2
nb2=`wc -l $genctss | awk '{print $1}'`
echo Initial number of gencode tss is $nb2
awk -v ext=50 '{$4=$4-ext; $5=$5+ext; $3="TSSExt"ext; print $0}' $genctss | awk '{if($4<1){$4=1}print $0}' | awk -f $GFF2GFF > $b2\_ext50eachside.gff
echo done >&2

# 2. Merge the gencode tss extended in a stranded way
#####################################################
echo I am merging the gencode tss extended \in a stranded way >&2
awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5, ".", ".", $7}' $b2\_ext50eachside.gff | $MERGE -i stdin -s -n | awk 'BEGIN{OFS="\t"}{$6=$5; $5="."; print}' > $b2\_ext50eachside_merged.bed
nb3=`wc -l $b2\_ext50eachside_merged.bed | awk '{print $1}'`
echo Number of gencode tss clusters after stranded merge is $nb3
echo done >&2

# 3. Merge the gencode tss clusters with the cage clusters
##########################################################
echo I am merging the gencode tss clusters with the cage clusters >&2
cat $cageclus $b2\_ext50eachside_merged.bed | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6}' | $MERGE -i stdin -s -n | awk 'BEGIN{OFS="\t"}{$7=$6; $6=$5; $5="."; print}' > $c\_$b2\_ext50eachside_merged.bed
nb4=`wc -l $c\_$b2\_ext50eachside_merged.bed | awk '{print $1}'`
echo Number of \(genctss + cage\) clusters after merging is $nb4
echo done >&2

# 4. Add many pieces of information to the file
###############################################
# a. annotation type
#####################
# CAGEOnly
# GencOnly
# GencCAGE
# b. gencode tss and genes each cluster comes from
##################################################
# c. cage cluster identifier (as coordinates)
#############################################
echo I am adding many pieces of information to the merged gencode tss + cage cluster file >&2
# a. Annotation class
#####################
echo "    1. the annotation class" >&2
$INTERSECT -a $c\_$b2\_ext50eachside_merged.bed -b $b2\_ext50eachside_merged.bed -s -wa -c > $c\_$b2\_ext50eachside_merged_overgencodenb.bed
$INTERSECT -a $c\_$b2\_ext50eachside_merged_overgencodenb.bed -b $cageclus -s -wa -c > $c\_$b2\_ext50eachside_merged_over_$b2\clusters_nb_over_cagetssclusters_nb.bed
awk '{print $1, "rikcrg", "tss", $2+1, $3, ".", $6, ".", "class:", (($7!=0) ? (($8!=0) ? "GencCAGE" : "GencOnly") : "CAGEOnly")}' $c\_$b2\_ext50eachside_merged_over_$b2\clusters_nb_over_cagetssclusters_nb.bed | awk -f $GFF2GFF > $c\_$b2\_ext50eachside_merged_over_$b2\clusters_nb_over_cagetssclusters_nb.gff
echo "Number of (genc tss + cage) clusters of each class" 
awk '{print $NF}' $c\_$b2\_ext50eachside_merged_over_$b2\clusters_nb_over_cagetssclusters_nb.gff | sort | uniq -c | sort -k1,1nr 
echo "    1. done" >&2
# b. add the info of genctss and gene
#####################################
echo "    2. the Gencode tss idenfier and gene" >&2
# Overlap to get tss coord
##########################
$OVERLAP $c\_$b2\_ext50eachside_merged_over_$b2\clusters_nb_over_cagetssclusters_nb.gff $b2\_ext50eachside.gff -st 1 -m -1 -i 2 -f genctssclus -v | awk '{$(NF-3)=""; $(NF-2)=""; print}' | awk -f $GFF2GFF > $c\_with_$b2\_ext50eachside_merged_genctsscoord.gff
# Overlap to get genes of tss
#############################
$OVERLAP $c\_with_$b2\_ext50eachside_merged_genctsscoord.gff $b2\_ext50eachside.gff -st 1 -m 10 -i 2 -f genctssclus -nr -v | awk '{$(NF-3)=""; $(NF-2)=""; print}' | awk -f $GFF2GFF > $c\_with_$b2\_ext50eachside_merged_genctsscoord_gnlist.gff
echo "    2. done" >&2
# c. add the info of initial cage cluster coordinates
#####################################################
echo "    3. the initial cage clusters' coordinates" >&2
echo "    3. done" >&2
$OVERLAP $c\_with_$b2\_ext50eachside_merged_genctsscoord_gnlist.gff $c.gff -st 1 -m 10 -i 2 -f cageclus -v | awk '{$(NF-3)=""; $(NF-2)=""; print}' | awk -f $GFF2GFF > $c\_with_$b2\_ext50eachside_merged_genctsscoord_gnlist_cageclusid.gff
echo done >&2

# 5. Clean
###########
echo I am cleaning unless I am asked not to do so >&2
if [ ! -n "$3" ]
then
rm $c.gff $b2\_ext50eachside.gff $b2\_ext50eachside_merged.bed $c\_$b2\_ext50eachside_merged.bed $c\_$b2\_ext50eachside_merged_overgencodenb.bed $c\_$b2\_ext50eachside_merged_over_$b2\clusters_nb_over_cagetssclusters_nb.bed $c\_$b2\_ext50eachside_merged_over_$b2\clusters_nb_over_cagetssclusters_nb.gff $c\_with_$b2\_ext50eachside_merged_genctsscoord.gff $c\_with_$b2\_ext50eachside_merged_genctsscoord_gnlist.gff
fi
gzip $c\_with_$b2\_ext50eachside_merged_genctsscoord_gnlist_cageclusid.gff
echo done >&2

