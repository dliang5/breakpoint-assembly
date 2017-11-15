#!/usr/bin/perl
use strict;

my $usage = "Use: $0 input.fastqn";

if (scalar @ARGV != 1) { die $usage; }
if (! (open (IN, "<$ARGV[0]"))) { die "Can't open $ARGV[0]: $!n"; }
if (-e "$ARGV[0].fasta") { die "$ARGV[0].fasta extists - please remove or renamen"; }
if (-e "$ARGV[0].qual") { die "$ARGV[0].qual extists - please remove or renamen"; }
if (! (open (FASTA, ">$ARGV[0].fasta"))) { die "Can't write to $ARGV[0].fasta: $!n"; }
if (! (open (QUAL, ">$ARGV[0].qual"))) { die "Can't write to $ARGV[0].qual: $!n"; }

my $id = "";
my $seq = "";
my $qual = "";
my $line = ;
if ($line !~ /@/) { die "$ARGV[0] does not look like a FASTQ file (it should start with @)n"; }
while ($line =~ s/^@//o) {
  chomp ($line);
  my $id = $line;
  $line = ;
  while ($line !~ s/^+//o) {
    chomp ($line);
    $seq .= $line;
    $line = ; }
  chomp ($line);
  if (($line =~ /S/o) && ($line ne $id)) { die "ID of $id not followed by same identifier for quality scoresn"; }
  $line = ;
  while (length ($qual) < length ($seq)) {
    chomp ($line);
    $qual .= $line;
    $line = ; }
  print FASTA ">$idn$seqn";
  print QUAL ">$idn";
  my @a = split //o, $qual;
  foreach my $i (@a) {
    my $q = 10 * log (1 + 10 ** ((ord ($i) - 64) / 10)) / log (10);
    printf QUAL "%1.0f ", $q; }
  print QUAL "n";
  $id = "";
  $seq = "";
  $qual = ""; }