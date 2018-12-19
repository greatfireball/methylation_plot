#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Bio::SeqIO;

use Log::Log4perl;

# Initialize Logger
my $log_conf = q(
   log4perl.rootLogger              = DEBUG, SCREEN
   log4perl.appender.SCREEN         = Log::Log4perl::Appender::Screen
   log4perl.appender.SCREEN.stderr  = 0
   log4perl.appender.SCREEN.layout  = Log::Log4perl::Layout::PatternLayout
   log4perl.appender.SCREEN.layout.ConversionPattern = [%d] %p %m %n
);
Log::Log4perl::init(\$log_conf);
my $log = Log::Log4perl->get_logger();

use version 0.77; our $VERSION = version->declare("v0.1");

use Term::ProgressBar;

my $man = 0;
my $help = 0;
my $version = 0;

my ($gff, $methylfile, $output, $fasta);

my $region_len        = 2000;    # default value from Xiang et al. (2010)
my $num_bins_per_part = 20;      # default value from Xiang et al. (2010)

my $bins              = {};      # final output
my @samples           = ();      # later used for sample naming and tracking of wrong input lines

GetOptions(
    'help|?'         => \$help,
    man              => \$man,
    'version'        => \$version,
    'gff=s'          => \$gff,
    'methfile=s'     => \$methylfile,
    'output=s'       => \$output,
    'regionlength=i' => \$region_len,
    'fasta=s'        => \$fasta,
    ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

if ($version)
{
    print $VERSION,"\n";
    exit;
}

unless (defined $gff && -e $gff)
{
    $log->logdie("Please specify a GFF3 file with gene annotations via --gff parameter");
}

unless (defined $methylfile && -e $methylfile)
{
    $log->logdie("Please specify a methylation file via --methfile parameter");
}

unless (defined $output)
{
    $log->logdie("Please specify an output file via --output parameter");
}

if (-e $output)
{
    $log->logdie("Output file exists, please specify a non existent file for output");
}

unless (defined $fasta)
{
    $log->warn("No fasta file was defined via --fasta paramert. Therefore, border testing is only possible for upstream region!")
}

if (defined $fasta && ! -e $fasta)
{
    $log->logdie("The specified fasta file can not be accessed");
}

# Import sequence information
my $contig_length = {};

if ($fasta)
{
    my $seqio = Bio::SeqIO->new(-file => $fasta);
    while (my $seq = $seqio->next_seq())
    {
	$log->logdie("Contig ".$seq->id()." already known!") if (exists $contig_length->{$seq->id()});

	$contig_length->{$seq->id()} = $seq->length();
    }

    $log->info("Imported ".int(keys %{$contig_length})." contigs");
}

# open GFF and import genes
my $annotation_data = {};
my $num_genes = 0;
my $num_genes_without_enough_distance = 0;

open(FH, "<", $gff) || $log->logdie("Unable to open input GFF: $!");
while (<FH>)
{
    next if (/^#/);

    chomp;
    my ($chr, undef, $type, $start, $stop, undef, $strand, undef) = split("\t", $_);

    if (defined $fasta && ! exists $contig_length->{$chr})
    {
	$log->logdie("The chromosome $chr was not found in input fasta");
    }

    next unless ($type eq "gene");

    # check if the up-/downstream region crosses the contig borders, if possible
    my $valid_gene = 1;
    my $upstream   = $start-$region_len;
    my $downstream = $stop+$region_len;
    if ($upstream<=0 || (defined $fasta && $downstream>$contig_length->{$chr}))
    {
	$num_genes_without_enough_distance++;
	$valid_gene = 0;
	$upstream = 1 if ($upstream <= 0);
	$downstream = $contig_length->{$chr} if ($downstream>$contig_length->{$chr});
    }

    push(@{$annotation_data->{$chr}{$strand}}, {start => $start, stop => $stop, valid => $valid_gene, upstream => $upstream, downstream => $downstream });
    $num_genes++;

}

close(FH) || $log->logdie("Unable to close input GFF: $!");

$log->info("Imported $num_genes genes for ".int(keys %{$annotation_data})." contigs, but $num_genes_without_enough_distance genes have to short distance to contig border");

# Import methylation information
my $methylation = {};
my $skipped_missing_annotation = 0;
my $skipped_missing_valid_gene = 0;
my $imported_methylation = 0;

# prepare ProgressBar
my $filesize = -s $methylfile;
my $progress = Term::ProgressBar->new({name => 'Methylation import', count => $filesize, remove => 1, ETA => 'linear'});
$progress->minor(0);
my $next_update = 0;
open(FH, "<", $methylfile) || $log->logdie("Unable to open input methylation file: $!");
while (<FH>)
{
    $next_update = $progress->update(tell(FH)) if tell(FH) >= $next_update;

    next if (/^#/);

    chomp;

    my ($chr, $pos, $strand, @methylated) = split(/\s+/, $_);

    # store the information about the number of samples if not already set
    unless (@samples)
    {
	@samples = map { 'sample_'.$_ } (1..int(@methylated));
    }

    # check the number of samples and die if not correct
    unless (@samples == @methylated)
    {
	$log->logdie(sprintf("Wrong number if samples data points in input line %d (expected %d, but found %d)", $., int(@samples), int(@methylated)));
    }

    unless (exists $annotation_data->{$chr})
    {
	$skipped_missing_annotation++;
	next;
    }

    # get all genes for the position:
    my @genes = ();
    my @strands = ();
    if ($strand eq ".")
    {
	@strands = keys %{$annotation_data->{$chr}};
    } else {
	@strands = ($strand);
    }

    foreach my $anno_strand (@strands)
    {
	my @subset = grep {$_->{valid} && $_->{upstream} <= $pos && $_->{downstream} >= $pos} (@{$annotation_data->{$chr}{$anno_strand}});
	push(@genes, @subset);
    }

    unless (@genes)
    {
	$skipped_missing_valid_gene++;
	next;
    }

    foreach my $gene (@genes)
    {
	update_bins($bins, $gene, \@methylated, $region_len, $pos, $num_bins_per_part);
    }
    $imported_methylation++;
}
close(FH) || $log->logdie("Unable to close input methylation file: $!");
$progress->update($filesize) if $filesize >= $next_update;
$log->info("Imported methylation status for $imported_methylation. Skipped $skipped_missing_annotation positions, due to missing annotation on contig. Skipped $skipped_missing_valid_gene positions, due to lack of valid genes. ");

sub update_bins
{
    my ($bins, $gene, $methylated, $region_len, $pos, $num_bins_per_part) = @_;

    my ($bin, $part);

    my $bin_size_upstream_downstream = $region_len / $num_bins_per_part;
    my $bin_size_gene = abs($gene->{start}-$gene->{stop}+1) / $num_bins_per_part;

    # upstream?
    if ($pos < $gene->{start})
    {
	$part = "upstream";
	$bin  = int(($pos-$gene->{upstream})/$bin_size_upstream_downstream);
    }
    # downstream ?
    elsif ($pos > $gene->{stop})
    {
	$part = "downstream";
	$bin  = int(($pos-$gene->{stop})/$bin_size_upstream_downstream);
    }
    # should be in gene
    else
    {
	$part = "gene";
	$bin = int(($pos-$gene->{start})/$bin_size_gene);
    }

    for(my $i=0; $i<@{$methylated}; $i++)
    {
	$bins->{$part}[$bin]{methylated}[$i]+=$methylated->[$i];
	$bins->{$part}[$bin]{total}[$i]++;
    }
}

__END__

=head1 methylation_plot.pl

C<methylation_plot.pl> - Simple program for data preparation for plots according to Xiang et al. (2010)

=head1 SYNOPSIS

    methylation_plot.pl [options] --gff GFF-File --methyfile Methylation-File --output output.csv
     Options:
       --gff             GFF3 with gene annotations
       --methfile        tab-seperated file with methylation information
       --output          name of the output file
       --regionlength    Length of up-/downstream region (default 2k)
       --fasta           Sequence file to enable border checking for up-/downstream
       --help            brief help message
       --man             full documentation
       --version         prints the current program version

=head1 PARAMETERS/OPTIONS

=over 8

=item B<--gff>

Specify the gff3 file which will be used for locations of genes.

=item B<--methfile>

Specify the file containing the information about methylated cysteins in the following format:

     chromosome[tab]position[tab]strand

=item B<--output>

Where should the tab seperated output be stored.

=item B<--regionlength>

Specifies the length of up-/downstream region to check.

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut
