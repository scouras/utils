#!/usr/bin/perl
use warnings;
use strict;


##############################################################################
# Legacy Utility Functions
# Alexander D. Scouras
# alexscouras@gmail.com
# 2004-2013
#-----------------------------------------------------------------------------
# These are functions from Utility.pm that were either superceeded, 
# deprecated, or simply are no longer used. But that I hadn't decided to 
# outright delete for some reason yet. Basically because they might turn out
# to be useful for some other legacy code, and I want to be able to quickly
# restore them if it comes up (without checking github archives).
##############################################################################

#################################################### ADJUST RESIDUES
# Only used by ~/phd/code/scouras/sasa_sum.
# Shifts a number by an offset, subject to min/max constraints.  
# Poor man's way of shifting a residue range.  

sub Adjust_Residues {
  my $residue = $_[0];
  my $offset  = $_[1];
  my $start   = $_[2];
  my $finish  = $_[3];
  $residue += $offset;
  if ( $residue < $start  ) { $residue = $start }
  if ( $residue > $finish ) { $residue = $finish }
  return $residue;
}














