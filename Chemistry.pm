#!/usr/bin/perl

use warnings;
use strict;





my $ph1 = $ARGV[0];
my $ph2 = $ARGV[1];
my $ph3 = $ARGV[2];
my $v   = $ARGV[3];

printf ( "V1\@pH$ph1: %.4f\tV2\@pH$ph2: %.4f\n\n", mix_two_buffers_to_target_ph($ph1, $ph2, $ph3, $v));





sub calc_buffer_needed_for_target_ph {

  my $current_ph      = $_[0];
  my $target_ph       = $_[1];
  my $target_volume   = convert_volumes_to_numbers($_[2]);
  my $buffer          = $_[3];

}




sub mix_two_buffers_to_target_ph { 
  my $ph1     = $_[0];
  my $ph2     = $_[1];
  my $ph3     = $_[2];
  my $volume  = $_[3];

  my $mol1 = 10**(-$ph1);
  my $mol2 = 10**(-$ph2);
  my $mol3 = 10**(-$ph3);
  my $v    = $volume;
  #my $v    = convert_volume_to_numbers($volume);

  my $X = ($v*($mol3-$mol2))/($mol1-$mol2);
  #my $X = ($v*($mol3+$mol2))/($mol1+$mol2);
  my $Y = ($v-$X);

  return ($X, $Y);

}


sub convert_ph_volume_to_mols { 
  my $ph = $_[0];
  my $volume = $_[1];
  return 10**(-$ph) * $volume;
}



sub convert_volume_to_numbers { 

  my $volume = $_[0];

  if ( $volume =~  /l/ ) { return $volume }
  if ( $volume =~ /ml/ ) { return $volume / 1,000 }
  if ( $volume =~ /ul/ ) { return $volume / 1,000,000 }
  if ( $volume =~ /nl/ ) { return $volume / 1,000,000,000 }
  
  return $volume;
}


sub convert_volumes_to_human_readable { 

  my $volume = $_[0];

  if ( $volume > 1    ) { return sprintf "%.4f l",  $volume }
  if ( $volume > 1e-3 ) { return sprintf "%.4f ml", $volume * 1,000 }
  if ( $volume > 1e-6 ) { return sprintf "%.4f ul", $volume * 1,000,000 } 
  if ( $volume > 1e-9 ) { return sprintf "%.4f nl", $volume * 1,000,000,000 } 

  return $volume;
}
