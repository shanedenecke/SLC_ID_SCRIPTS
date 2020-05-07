#!/usr/bin/perl -ws

#use Bio::SeqIO;

if ( scalar (@ARGV) == 0 ) {
	my @tmp = split ( /\//, $0 );
	die "USAGE: $tmp[0] <input fasta file> [<outgroup>] [<CPUs>]\n\n";
}

my $in_file = $ARGV[0];
chomp ($in_file);

my $PWD = `pwd`;
chomp ($PWD);

my $OUTGROUP; # if it's left undef, then no outgroup will be used in RAxML
if ( $ARGV[1] ) {
	$OUTGROUP = $ARGV[1];
}

my $CPUs = 16;
if ( $ARGV[2] ) {
	$CPUs = $ARGV[2];
}

#my $SPECIES_TO_USE = "DMELA,DMOJA,MDEST,AGAMB,CQUIN,HMELP,BMORI,TCAST,NVITR,APISU,DPULE";

# 1) run prinseq to change the headers
#print STDERR "Running prinseq-lite...\n";
#system ( "prinseq-lite.pl -fasta $in_file -rm_header -seq_id prot_ -out_good stdout > $in_file.renamed" );
#print STDERR " Done!\n\n";

# 1) The following is OrthoDB-specific; if you don't have orthologous groups
#    from OrthoDB as input, then you can use prinseq (above), in order to
#    get unique sequence IDs that will also be phylip-compatible
#print STDERR "Compacting headers (5-letter species name only) and keeping only certain species...\n";

#open ( RENAMED, ">$in_file.renamed" ) or die;

#my $in_fasta  = Bio::SeqIO->new(-file => $in_file, -format => 'Fasta') or die;

#my $count = 0;

#while ( my $seq = $in_fasta->next_seq() ) {
#	#print STDERR $seq, "\n";
#	print STDERR $seq->name(),"\n";
#	#my ( $sp_name ) = ( $seq->desc() =~ /.+ ([A-Z]{5})$/ );
#	my ( $sp_name ) = ( $seq->desc() =~ /([A-Z]{5}):/ );
#	if ( $SPECIES_TO_USE =~ /$sp_name/ ) {
#		print RENAMED ">$sp_name\n";
#		print RENAMED $seq->seq(),"\n";
#		
#		$count++;
#	}
#}
#close (RENAMED);

#my $asdf = <STDIN>;

# checking if all the species were found
#my @tmp = split ( /,/, $SPECIES_TO_USE );
#if ( scalar (@tmp) != $count ) {
#	print STDERR "### Not all species in the list were found in the fasta file #####\n";
#}

# 1) reduce headers
print STDERR "Compacting headers (5-letter species name only) and keeping only certain species...\n";
system ( "perl -pe 's/:.+\$//g;' $in_file > $in_file.renamed" );
print STDERR "Done!\n\n";

# 2) run muscle
print STDERR "Running muscle...\n";
system ( "muscle -in $in_file.renamed -out $in_file.renamed.aln" );
print STDERR "Done!\n\n";

# 3) run trimAl
print STDERR "Running trimAl...\n";
#system ( "trimal -in $in_file.renamed.aln -out $in_file.renamed.aln.trimmed -fasta -w 3 -gt 0.95 -st 0.01 -htmlout $in_file.renamed.aln.trimmed.html" );
system ( "trimal -in $in_file.renamed.aln -out $in_file.renamed.aln.trimmed -fasta -automated1" );
print STDERR "Done!\n\n";

# 4) convert to phylip format
print STDERR "Running seqret...\n";
system ( "../../convert.sh $in_file.renamed.aln.trimmed > $in_file.renamed.aln.trimmed.phy" );
print STDERR "Done!\n\n";

# 5) run RAxML
print STDERR "Running RAxML...\n";
if ( $OUTGROUP ) {
	system ( "raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $CPUs -m PROTGAMMAAUTO -s $in_file.renamed.aln.trimmed.phy -w $PWD -n $in_file.renamed.aln.trimmed.phy.tre -o $OUTGROUP" );
}
else {
	system ( "raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 100 -T $CPUs -m PROTGAMMAAUTO -s $in_file.renamed.aln.trimmed.phy -w $PWD -n $in_file.renamed.aln.trimmed.phy.tre" );
}
print STDERR "Done!\n\n";

# 6) Calculate average bootstrap support
print STDERR "Calculating average bootstrap support...\n";
my $avg_boot_sup = `nw_labels RAxML_bipartitions.$in_file.renamed.aln.trimmed.phy.tre -L | awk '{sum+=\$1; lines++} END {print sum/lines/100}'`;
open ( AVG_BOOT_SUP, ">RAxML_bipartitions.average_bootstrap_support.txt" ) or die;
print AVG_BOOT_SUP "average bootstrap_support (x100): ", $avg_boot_sup, "\n";
close (AVG_BOOT_SUP);
print STDERR "Done!\n\n";

# 7) Calculate tree certainty
print STDERR "Calculating tree certainty...\n";
system ( "raxmlHPC-AVX -L MRE -z RAxML_bootstrap.$in_file.renamed.aln.trimmed.phy.tre -m GTRCAT -n T1" );
print STDERR "Done!\n";
# The tree certainty value is stored in the "RAxML_info.T1" file; it's in
# the line "Relative tree certainty for this tree: <number from 0 to 1>"
