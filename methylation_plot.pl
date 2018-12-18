#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

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

my $man = 0;
my $help = 0;
my $version = 0;

my ($gff, $methylfile, $output);

GetOptions(
    'help|?'     => \$help,
    man          => \$man,
    'version'    => \$version,
    'gff=s'      => \$gff,
    'methfile=s' => \$methylfile,
    'output=s'   => \$output
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

__END__

=head1 methylation_plot.pl

C<methylation_plot.pl> - Simple program to plot data according to Xiang et al. (2010)

=head1 SYNOPSIS

    methylation_plot.pl [options] --gff GFF-File --methyfile Methylation-File --output output.pdf
     Options:
       --gff             GFF3 with gene annotations
       --methfile        tab-seperated file with methylation information
       --output          name of the output file
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

Where should all the pdf output be stored.

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut
