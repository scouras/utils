#!/usr/bin/perl
BEGIN { push @INC, "$ENV{CODE_BASE}/scouras" }
use warnings;
use strict;
use Phenix;
use Utility;
use Data::Dumper;
use Date::Manip;
use IPC::Open2;

sub usage { 
  print "$0 <title> <batch 1> [<batch N>]\n";
  exit;
}

my $TITLE = $ARGV[0] or die &usage;
my @BATCHES = ();

foreach my $arg ( @ARGV[1..$#ARGV] ) { 
  if ( -e $arg ) { push @BATCHES, $arg } 
  else { 
    die "Couldn't find batch file, '$arg'.\n";
  }
}




