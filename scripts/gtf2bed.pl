#!/usr/bin/perl

my $gtf = shift;
unless (-f $gtf) {
	warn "Please provide a GTF file\n";
	exit;
}

# find exon features
open(IN, "cat $gtf \| awk '\$3==\"exon\"' |");

# iterate over each exon
while(<IN>) {
	chomp();
	my @data = split(/\t/);

	# it's for BEDTools, which is zero-based, half open
	# see https://code.google.com/p/bedtools/wiki/FAQ#What_does_zero-based,_half-open_mean?
	# "In other words, BEDTools interprets the .start. column as being 1 basepair higher than what is represented in the file"
	
	# print the chromosome
	print $data[0], "\t";

	# print the start -1
	print $data[3] - 1, "\t";

	# print the end
	print $data[4], "\t";

	# create a name for the feature
	my $eid = "no_exon";
        if (m/exon_id "(\S+)";/) {
                $eid = $1;
        }

	my $gid = "no_gene";
        if (m/gene_id "(\S+)";/) {
                $gid = $1;
        }

	# print out the name
	print "$eid.$gid", "\t";

	# print no score
	print "-", "\t";

	# print the strand
	print $data[6], "\n";

}
close IN;


