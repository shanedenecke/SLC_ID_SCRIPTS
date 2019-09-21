#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use Bio::SeqIO;
use File::Find::Rule;

my %opts = ();
getopts('p:b:t:o:', \%opts);
# p: directory suffix
# b: average bootstrap support (BS)
# t: tree certainty (TC)

#print STDERR (keys%opts), "\n";

if( scalar ( keys%opts ) != 3 ){
	my @temp = split ( /\//, $0 );
	die	
"
\nUSAGE: $temp[-1] -p <directory suffix> (-b|-t) <threshold value> -o <output directory>
This script gets trees whose BS or TC is above the threshold in sub-directories
that are located in the current directory and whose names begin with the suffix

";
}

my $output_dir = $opts{o};
$output_dir =~ s/\/$//g; # remove the final slash, if it's added
mkdir ( $output_dir, 0755 ) or die "Cannot create output directory\n";

my $directory_suffix = $opts{p};

# get all subdirectories
my $pwd = ".";
my $find_rule = File::Find::Rule->new();
# set the rules
$find_rule->maxdepth(1);
$find_rule->directory;
# apply the rules
my @sub_dirs = $find_rule->in($pwd);

#print STDERR join ( "\n", @sub_dirs );

# this hash will hold the concatenated sequences. Sequence will grow as new
# sequence will be coming from genes with strong phylogenetic signal
my %concat = ();

# this hash will hold the individual gene trees (with bootstrap support values)
my %trees = ();

foreach my $sub_dir (@sub_dirs) {
	unless ( $sub_dir =~ /$directory_suffix$/ ) {
		# ignore sub-directories not ending with the given suffix
		next;
	}
	#print STDERR "Entering: $sub_dir\n";
	chdir ( $sub_dir );

	# See whether the gene used for this tree has a strong phylogenetic signal
	my $has_strong_signal = &has_strong_signal();

	#print STDERR $has_strong_signal, "\n";

	if ( $has_strong_signal == 1 ) {
		#print STDERR "  It has strong signal\n";
		my $aln_file = `ls *aln.trimmed`;
		my $in_fasta = Bio::SeqIO->new(-file => $aln_file, -format => 'Fasta') or die;
		
		# push all sequences into an array so that you can sort the sequences
		# alphabetically. This way you'll have the same order in every file.
		# Of course, the exact same species should exist in EVERY tree you examine!
		#my @sequences = ();
		
		while ( my $seq = $in_fasta->next_seq() ) {
			#print STDERR $seq->display_id, "\n";
			#print STDERR $seq->desc,"\n";
			#print STDERR $seq->seq,"\n";
			# display_id should only contain the 5-letter species name

			# initialize the hash
			if ( !exists ($concat{$seq->display_id}) ) {
				$concat{$seq->display_id} = $seq->seq;
				#print STDERR "    Initialization\n";
			}
			# ...or simply grow the sequence of a given species
			else {
				$concat{$seq->display_id} .= $seq->seq;
				#print STDERR "    Extension\n";
			}
		}
		
		# Also take the tree with the BS values
		my $tree_file = `ls RAxML_bipartitions.*.renamed.aln.trimmed.phy.tre`;
		my ( $pub_og_id ) = ( $tree_file =~ /RAxML_bipartitions.(\S+).renamed.aln.trimmed.phy.tre/ );
		open ( TREE, $tree_file ) or die "Can't find the tree file in $sub_dir/";
		# there's just one line in this file; the newick-formatted tree
		my $nw_tree = <TREE>;
		chomp ($nw_tree);
		$trees{$pub_og_id} = $nw_tree;
	}
	
	# move back to "home" directory
	#print STDERR "Exiting: $sub_dir\n\n";
	chdir "../";
}

# print the concatenated alignment (fasta) file
open ( CONCAT, ">$output_dir/concat_genes.fs.aln.trimmed" ) or die;
foreach my $concat_seq (keys%concat) {
	print CONCAT ">$concat_seq\n";
	print CONCAT $concat{$concat_seq}, "\n";
}
close (CONCAT);

# also print two more files: (a) the OG IDs, (b) the trees from each OG.
open ( OG, ">$output_dir/og_ids.txt" ) or die;
open ( TREES, ">$output_dir/og_trees.tre" ) or die;

foreach my $og_id ( keys%trees ) {
	print OG "$og_id\n";
	print TREES $trees{$og_id}, "\n";
}
close (OG);
close (TREES);


sub has_strong_signal {
	my $return_value = 0;
	# if filtering is done based on BS
	if ( $opts{b} ) {
		if ( -e "RAxML_bipartitions.average_bootstrap_support.txt" ) {
			open ( BS, "RAxML_bipartitions.average_bootstrap_support.txt" ) or die;
			while ( my $line = <BS> ) {
				chomp ($line);
				if ( $line =~ /: (\S+)$/ ) {
					my $bs = $1;
					#print STDERR "BS: $bs\n";
					if ( $bs >= $opts{b} ) {
						$return_value = 1;
					}
				}
			}
		}
	}
	
	# if filtering is done based on TC
	if ( $opts{t} ) {
		if ( -e "RAxML_info.T1" ) {
			open ( TC, "RAxML_info.T1" ) or die;
			while ( my $line = <TC> ) {
				chomp ($line);
				if ( $line =~ /Relative tree certainty for this tree: (\S+)$/ ) {
					my $tc = $1;
					if ( $tc >= $opts{t} ) {
						$return_value = 1;
					}
				}
			}
		}
	}
	return $return_value;
}
