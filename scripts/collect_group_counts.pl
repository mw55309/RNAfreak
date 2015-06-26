#!/usr/bin/perl

#
# This script is designed to be used in combination with others, see https://github.com/mw55309/RNAfreak
#
# The script takes a read count limit, a proportion and a list of cts files output from count_from_bed.pl
# A group has to be present at greater than the read count limit in more than the proportion of files
# to be in the final result.  Otherwise, the script simply takes the cts files and creates a matrix
# of counts for each gene group
#
# Run without any arguments to see usage
#
# NB if you have not used Ensembl genes in the pipeline, you need to change "XF:Z:ENS" below
#

use strict;
use POSIX;

unless (@ARGV>0) {
	warn "Usage: perl collect_group_counts.pl <read count limit> <proportion of datasets> <cts file 1> <cts file 2> ... <cts file n>\n";
	exit;
}

# this is the number of reads a group must have in
# each experiment to remain in the dataset
my $num = shift;

# check we have an argument
unless (defined $num) {
	warn "Please provide a number of reads as the first argument\n";
	exit;
}

# check it is a number
if ($num =~ m/\D/) {
	warn "$num does not look like a number, quitting\n";
	exit;
}

# this is the proportion (fraction) of experiments
# that must meet the criteria for a group to remain
# in the dataset
my $prop = shift;
unless (ceil($prop) == 1) {
	warn "$prop does not look like a fraction\n";
	exit;
}

# everything else is a SAM file
my @cts = @ARGV;
unless (scalar(@cts) > 0) {
	warn "please provide a list of counts files\n";
	exit;
}

# variable to store the results
my $g;

# iterate over CTS files
foreach my $c (@cts) {

	# read file
	open(IN, $c);
	while(<IN>) {
		chomp();

		# first column is group, second column is counts
		my ($g1,$ct) = split(/\t/);

		# split and sort groups to be sure
		my @g = split(/,/, $g1);
		my $gkey = join(",", sort @g);

		# store the count for this group in this file
		$g->{$gkey}->{$c} = $ct;
	}
	close IN;
}

# calculate the limit - that is
# the number of files a group has to
# be present in to remain in the final
# results
my $limit = int(scalar(@cts) * $prop);

# initialize group numbering system
my $gid = 1;

# print out headers
print "GROUP\tMEMBERS\t", join("\t", @cts), "\n";

# iterate over hash ref
while(my($gp,$hr) = each %{$g}) {

	# variable to store the number of files
	# a group is present in
	my $count = 0;

	# start to build the output line, starting
	# with group ID and group members
	my $line = "G" . $gid . "\t" . $gp;

	# foreach file we measured
	foreach my $c (@cts) {

		# if this group was in this file...
		my $value = 0;
		if (exists $hr->{$c}) {
			$value = $hr->{$c};
		}

		# add it to the line anyway
		$line .= "\t" . $value;

		# add 1 to the count if the group is present
		# at greater than $num counts
		if ($value >= $num) {
			$count++;
		}

	}

	# only if the group is present in greater than
	# $limit files do we print out the results
	if ($count >= $limit) {
		print $line, "\n";
		$gid++;
	}	

}
