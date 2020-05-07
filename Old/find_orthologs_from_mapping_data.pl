#!/usr/bin/perl -w

use warnings;
use strict;

	if ( scalar (@ARGV) != 2 ) { die "USAGE: $0 <*map.og output> <species list>\n\n"; }
# The *map.og file look like this

# # clid  key     type    length  start   end     rawscore        normscore       evalue
# 2       7220:0018f3     0       412     2736    3147    16468   25.2229 0
# 2       7245:00139d     0       412     2732    3143    16468   25.2229 0
# 2       7227:0003ca     0       1088    1       1088    16424.1 24.8624 0
# 2       7244:001675     0       1888    5113    7000    15937.6 31.2135 0
# 2       7260:000e2e     0       1888    5362    7249    15937.6 31.2135 0
# 2       69004:002e77    0       1217    1066    2282    15738.5 23.9515 0
# ...

# the species list is one species name per line and can either be the taxon ID
# or the 5-letter name.

# load the species names into a hash
my %species_list = ();
open ( SP_LIST, $ARGV[1] ) or die;

while ( my $line = <SP_LIST> ) {
	if ( $line =~ /^#/ ) { next; }
	chomp ($line);
	$species_list{$line} = 1;
}
close (SP_LIST);

# start reading the map.og file
open ( OG_MAP, $ARGV[0] ) or die;

my $line = <OG_MAP>;

while ( $line =~ /^#/ ) {
	$line = <OG_MAP>;
}

my @sc_all = ();    # single-copy in all species
my @pr_all = ();    # present-in-all species
my @pr_maj = ();    # present in the majority (all except one or two species)
my @pr_oth = ();    # present in 'other' species, less than the majority
my @oth = ();       # orthologs in species other than those in %species_list

# the following will keep the gene names found in a cluster
my %genes_found = ();

my %fh =();

# initialize each array within the above hash
foreach my $species ( sort keys %species_list) {
	$genes_found{$species} = ();
}

while ( $line ) {
	if ( $line=~/^#/ ) { last; }
	chomp ($line);

	my @f = split (/\s/, $line);

	my $cl_name = $f[0];
	
	my %counters =();
	foreach my $species ( sort keys %species_list ) {
		$counters{$species} = 0;
	}

	# the following will hold the cluster members of the current cluster only, as
	# they're shown in the original file.
	# By keeping the original record, you can print the gene name of the species of
	# interest, if it belongs in "single copy ortholog, present in all 4 flies", 
	# "present in 3 flies", "present in 2 flies", etc.
	# In purely technical terms, keeping the original lines in the following array
	# makes them available in later parts of the "while" loop.
	my @cluster_members = ();
	
	# similar to the above, but keeps only the gene name and also groups genes
	# according to the species they belong
	my %species_genes_cl = ();

	while ( $line && ( $f[0] eq $cl_name ) ) {
		# get the species name
		my ( $curr_species ) = ( $f[1] =~ /^(\S+):\S+$/ );
		
		if ( $species_list{$curr_species} ) {
			$counters{$curr_species}++;
			push ( @{$genes_found{$curr_species}}, $f[1] );
			push ( @cluster_members, $line );
			push ( @{$species_genes_cl{$curr_species}}, $f[1] );
		}
		
		$line = <OG_MAP>;

		if ( !$line ) { last; } # end of file

		chomp ($line);
		@f = split ( /\s/, $line );
	}
	# count the number of species having a single gene copy in this OG
	my $sc_sum = 0;
	foreach my $species ( sort keys %species_list ) {
		if ( $counters{$species} == 1 ) {
			$sc_sum++;
		}
	}
	
	# count the number of species present in this OG
	my $pr_sum = 0;
	foreach my $species ( sort keys %species_list ) {
		if ( $counters{$species} >= 1 ) {
			$pr_sum++;
		}
	}
	
	if ( $sc_sum == scalar ( keys %species_list ) ) { # if all species have a single-copy in this OG
		# put all species counts in a single string
		my $all_counts = "$cl_name";
		foreach my $species ( sort keys %species_list ) {
			$all_counts .= "\t$counters{$species}";
		}
		push ( @sc_all, $all_counts );
		
		# print the gene names of each SC cluster so that you can
		# later retrieve the amino acid sequences (for tree building)
		open ( SC_OUT, ">$cl_name.names" ) or die;
		foreach my $cluster_member (sort @cluster_members) {
			my @ff = split ( /\s/, $cluster_member );
			print SC_OUT "$ff[1]\n";
		}
		close (SC_OUT);
		#print STDERR "SC\t$dmela_cnt\n";
		
		# print names of genes for each species that belong to this category
		foreach my $species ( sort keys %species_list ) {
			open ( OUT, ">>$species.with_orthologs.details.list" ) or die;
			foreach my $tmp ( @{$species_genes_cl{$species}} ) {
					# this will loop through all the genes for $species
					print OUT "$tmp\tSC_ALL\n";
			}
			close (OUT);
		}

	}
	elsif ( $pr_sum == scalar ( keys %species_list ) ) {
		# put all species counts in a single string
		my $all_counts = "$cl_name";
		foreach my $species ( sort keys %species_list ) {
			$all_counts .= "\t$counters{$species}";
		}
		push ( @pr_all, $all_counts );
	
		# print names of genes for each species that belong to this category
		foreach my $species ( sort keys %species_list ) {
			open ( OUT, ">>$species.with_orthologs.details.list" ) or die;
			foreach my $tmp ( @{$species_genes_cl{$species}} ) {
					# this will loop through all the genes for $species
					print OUT "$tmp\tPR_ALL\n";
			}
			close (OUT);
		}
	}
	elsif ( $pr_sum == scalar ( keys %species_list ) - 1 || $pr_sum == scalar ( keys %species_list ) - 2 ) {
		# put all species counts in a single string
		my $all_counts = "$cl_name";
		foreach my $species ( sort keys %species_list ) {
			$all_counts .= "\t$counters{$species}";
		}
		push ( @pr_maj, $all_counts );

		# print names of genes for each species that belong to this category
		foreach my $species ( sort keys %species_list ) {
			open ( OUT, ">>$species.with_orthologs.details.list" ) or die;
			foreach my $tmp ( @{$species_genes_cl{$species}} ) {
					# this will loop through all the genes for $species
					print OUT "$tmp\tPR_MAJ\n";
			}
			close (OUT);
		}
	}
	elsif ( $pr_sum >= 2 && $pr_sum <= scalar ( keys %species_list ) - 3 ) {
		# put all species counts in a single string
		my $all_counts = "$cl_name";
		foreach my $species ( sort keys %species_list ) {
			$all_counts .= "\t$counters{$species}";
		}
		push ( @pr_oth, $all_counts );
	
		# print names of genes for each species that belong to this category
		foreach my $species ( sort keys %species_list ) {
			open ( OUT, ">>$species.with_orthologs.details.list" ) or die;
			foreach my $tmp ( @{$species_genes_cl{$species}} ) {
					# this will loop through all the genes for $species
					print OUT "$tmp\tPATCHY\n";
			}
			close (OUT);
		}
	}
	else {
		# put all species counts in a single string
		my $all_counts = "$cl_name";
		foreach my $species ( sort keys %species_list ) {
			$all_counts .= "\t$counters{$species}";
		}
		push ( @oth, $all_counts );

		# print names of genes for each species that belong to this category
		foreach my $species ( sort keys %species_list ) {
			open ( OUT, ">>$species.with_orthologs.details.list" ) or die;
			foreach my $tmp ( @{$species_genes_cl{$species}} ) {
					# this will loop through all the genes for $species
					print OUT "$tmp\tOTHERS\n";
			}
			close (OUT);
		}
	}
}

#print scalar (@sc), "\n";
#print scalar (@pr), "\n";
#print scalar (@pr8_9), "\n";
#print scalar (@pr2_7), "\n";
#print scalar (@oth), "\n";

# First, print the gene names that are found in clusters

foreach my $species ( sort keys %species_list ) {
	open ( FOUND, ">$species.with_orthologs.list" ) or die;
	foreach my $gene ( @{$genes_found{$species}} ) {
		print FOUND "$gene\n";
	}
}


# Second, loop through each array and count the number of genes for each species
my %counters =();
print "\t"; # this is for the header
foreach my $species ( sort keys %species_list ) {
	print "\t$species"; # header
	$counters{$species} = 0;
}
print "\n"; # header

foreach my $tmp (@sc_all) {
	my @f = split ( /\t/, $tmp );

	my @sorted_species_list = sort keys %species_list;
	for ( my $i = 0; $i < scalar (@sorted_species_list); $i++ ) {
		$counters{$sorted_species_list[$i]} += $f[$i+1];
	}
}
print "Single copy";
foreach my $species ( sort keys %species_list ) {
	print "\t$counters{$species}";
}
print "\n";





foreach my $species ( sort keys %species_list ) {
	$counters{$species} = 0;
}

foreach my $tmp (@pr_all) {
	my @f = split ( /\t/, $tmp );

	my @sorted_species_list = sort keys %species_list;
	for ( my $i = 0; $i < scalar (@sorted_species_list); $i++ ) {
		$counters{$sorted_species_list[$i]} += $f[$i+1];
	}
}
print "Present in all";
foreach my $species ( sort keys %species_list ) {
	print "\t$counters{$species}";
}
print "\n";






foreach my $species ( sort keys %species_list ) {
	$counters{$species} = 0;
}

foreach my $tmp (@pr_maj) {
	my @f = split ( /\t/, $tmp );

	my @sorted_species_list = sort keys %species_list;
	for ( my $i = 0; $i < scalar (@sorted_species_list); $i++ ) {
		$counters{$sorted_species_list[$i]} += $f[$i+1];
	}
}
print "Present in maj";
foreach my $species ( sort keys %species_list ) {
	print "\t$counters{$species}";
}
print "\n";






foreach my $species ( sort keys %species_list ) {
	$counters{$species} = 0;
}

foreach my $tmp (@pr_oth) {
	my @f = split ( /\t/, $tmp );

	my @sorted_species_list = sort keys %species_list;
	for ( my $i = 0; $i < scalar (@sorted_species_list); $i++ ) {
		$counters{$sorted_species_list[$i]} += $f[$i+1];
	}
}
print "Present in >=2";
foreach my $species ( sort keys %species_list ) {
	print "\t$counters{$species}";
}
print "\n";




foreach my $species ( sort keys %species_list ) {
	$counters{$species} = 0;
}

foreach my $tmp (@oth) {
	my @f = split ( /\t/, $tmp );

	my @sorted_species_list = sort keys %species_list;
	for ( my $i = 0; $i < scalar (@sorted_species_list); $i++ ) {
		$counters{$sorted_species_list[$i]} += $f[$i+1];
	}
}
print "Orth with other";
foreach my $species ( sort keys %species_list ) {
	print "\t$counters{$species}";
}
print "\n";

