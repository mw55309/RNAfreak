#!/usr/bin/perl

#
# This script is designed to be used in combination with others, see https://github.com/mw55309/RNAfreak
# 
# The script takes a matrix of group counts from collect_group_counts.pl and merges
# all groups that are completely contained within a larger group
#
# Run without any arguments to see usage
#

use strict;

# we should have one argument, a text matrix
# of group counts
unless (@ARGV==1) {
	warn "Usage: perl collapse_groups.pl <group counts file>\n";
	exit;
}

# open file and extract titles
my $gcounts = shift;
my @titles;
open(IN, "$gcounts");
while(<IN>) {
	chomp();
	@titles = split(/\t/, $_);
	last;
}

# $g will store the counts, indexed by
# the group membership column (column 2)
my $g;
while(<IN>) {
	chomp();
	my @d = split(/\t/);

	# $d[1] is the group membership i.e. a list of genes
	# $g->{...} is a hashref, keyed on the group membership
	# @{$g->{...}} "casts" this to an array type
	@{$g->{$d[1]}} = @d[2..$#d];

}

# sort the groups on length
# note here string length is a proxy for group size
# and relies on Ensembl identifiers all being the
# same length
my @groups = sort {length($b) <=> length($a)} keys %{$g};

# hash to hold indices of groups that have been merged
my %skip;

# variable to hold the merged counts
my $mg;

# go through the groups, from largest to smallest
for(my $i=0;$i<@groups;$i++) {

	# next if this group has already been
	# merged
	next if (exists $skip{"G".$i});

	# $g1 is the current group
	my $g1 = $groups[$i];

	# set the up the merged counts starting at the level of the current group
	$mg->{$groups[$i]} = $g->{$groups[$i]};

	# iterate through all groups below this one
	# in the list
	for (my $j=($i+1);$j<@groups;$j++) {
		
		# next if this group has already been
		# merged
		next if (exists $skip{"G".$j});

		# $g2 is the group being investigated
		# to see if it is completely contained in $g1
		my $g2 = $groups[$j];

		# can we merge these two groups?
		my $merge = &can_merge($g1,$g2);

		if ($merge == 1) {
			# yes merge it
			
			# add the values for $g2 to the values
			# for group $i in $mg
			my @add = @{$g->{$groups[$j]}};
			for ($a=0;$a<@add;$a++) {
				$mg->{$groups[$i]}->[$a] += $add[$a];
			}
			
			# group $j has been merged, skip in future
			$skip{"G".$j}++;
		} 
	}
}

# print out the titles and the merged groups
print join("\t", @titles), "\n";
my $mgc = 1;
foreach my $key (keys %{$mg}) {
	print "MG", $mgc++, "\t", $key, "\t", join("\t", @{$mg->{$key}}), "\n";
}





sub can_merge {

	# the groups
	my $g1 = shift;
	my $g2 = shift;

	# hashes to hold members
	my %h1;
	my %h2;

	# store the group members in the hashes
	if (length($g1) > length($g2)) {
		foreach my $element (split(/,/, $g1)) {
			$h1{$element}++;
		}

		foreach my $element (split(/,/, $g2)) {
			$h2{$element}++;
		}
	} else {
		foreach my $element (split(/,/, $g2)) {
			$h1{$element}++;
		}

		foreach my $element (split(/,/, $g1)) {
			$h2{$element}++;
		}
	}

	# go through the biggest group
	# and delete corresponding members from the 
	# smallest group
	foreach my $key (keys %h1) {
		if (exists $h2{$key}) {
			delete $h2{$key};
		}
	}

	# if there is something left in %h2, we cannot merge
	if (scalar(keys %h2) > 0) {
		return 0;
	} else {
		# otherwise we can merge
		return 1;
	}

}

