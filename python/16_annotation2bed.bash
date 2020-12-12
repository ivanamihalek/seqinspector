#!/bin/bash


# Part of seqinspector toy NGS QC pipeline
# Ivana Mihalek,  2020

# Quick adn dirty reduction of the file size for the annotation of regions
# by the genes they cover. Requires hg19.ncbiRefSeq.gtf file
# available from https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/


cut -f 1,4,5,9 ../annotation/hg19.ncbiRefSeq.gtf > ../annotation/hg19.ncbiRefSeq.bed
# ncbi first, to keep the annotation column
~/third/bedtools/bedtools intersect  -a  ../annotation/hg19.ncbiRefSeq.bed -b  ../pileup/coverage/merged_target_regions.bed \
    > ../annotation/target_regions_annotated.bed
rm -f ../annotation/hg19.ncbiRefSeq.bed
