#!/usr/bin/perl

#
# This script is designed to be used in combination with others, see https://github.com/mw55309/RNAfreak
#
# The script takes a SAM output file from htseq-count, and pipes the non-uniquely mapped
# reads through samtools and BEDTools to create "multi-map groups" - that is, groups of
# genes that multi-mapped reads uniquely map to.  It then counts the number of MMGs 
# created and outputs to STDOUT
#
# Run without any arguments to see usage
#
# NB if you have not used Ensembl genes in the pipeline, you need to change "XF:Z:ENS" below
#

use strict;
use POSIX;

# if not in your $PATH, edit the variables below
my $samtools = "samtools";
my $bedtools = "intersectBed";

# we should have two arguments, a SAM file 
# and a BED file
unless (@ARGV==3) {
	warn "Usage: perl count_from_bed.pl <sam file> <bed file> <fasta file>\n";
	exit;
}

# get SAM file from the command line
my $sam = shift;
unless (-f $sam) {
	warn "$sam is not a file\n";
	exit;
}

# get BED file from the command line
my $bed = shift;
unless (-f $bed) {
	warn "$bed is not a file\n";
	exit;
}

# get FASTA file from the command line
my $fasta = shift;
unless (-f $fasta && -f "$fasta.fai") {
	warn "$fasta or its index does not exist\n";
	exit;
}

# figure out if it's zipped
my $cat = "cat";
if ($sam =~ m/\.gz/) {
	$cat = "zcat";
}

# anything matching XF:Z:ENS has successfully mapped a read to a gene
# so we can ignore them.  Everything else is fair game
# BEDTools expects a BAM file so we need to pipe through samtools
open(IN, "$cat $sam \| grep -v \"XF:Z:ENS\" \| $samtools view -S -b -T $fasta - \| $bedtools -abam - -bed -wb -b $bed |");

# variables to hold the results
my $r = undef;
my $g - undef;

while(<IN>) {

	# split on whitespace
	my @d = split(/\s+/);

	# get the read ID
	my $read = $d[3];

	# get the hit ID 
	my $hit  = $d[15];

	# hit ID from BED file is assumed to be exon_id.gene_id
	my ($exon,$gene) = split(/\./, $hit);
	
	# count that gene for that read
	# only once using hash refs
	$g->{$read}->{$gene}++;
}
close IN;

# go through $g and count
# instances of gene groups
my %c;
while(my($read,$hr) = each %{$g}) {
	
	# get gene list for this read
	my @genes = sort keys %{$g->{$read}};

	# create a key
	my $genekey = join(",", @genes);

	# count that key in a hash
	$c{$genekey}++;
}

# print out group counts
while(my($group,$count) = each %c) {
	print "$group\t$count\n";
}	

