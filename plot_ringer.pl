#!/usr/bin/perl
BEGIN { push @INC, "$ENV{CODE_BASE}/scouras" }
use warnings;
use strict;
use Utility;
use Data::Dumper;


### TODO:
# Find colorblind friendly palatte
# Make lines thicker, at least in legend. Probably longer too.
# lmargin too wide in analysis mode, too thin in paper mode
#


sub usage {
  print "$0 <conf file> \n"
      . "<-compare>     compare electron density\n"
      . "<-dihedral>    plot dihedrals \n"
      . "<-res=1,2,3>   limit to residues \n"
      . "<-chi=3>       maximum chi angle\n"
      . "<-ala>         include alanines for reference\n"
      . "<-den=X>       don't plot residues below peak density X\n"
      . "<-nomatrix>    skip plotting comparison matrix\n"
      . ($_[0] || "\n");
  exit;
}

#my $PLOT_WIDTH_COMP   = 8000;
#my $PLOT_HEIGHT_COMP  =  300;
my $PLOT_WIDTH_COMP   = 4000;
my $PLOT_HEIGHT_COMP  =  500;


my $STYLE = 'presentation';
#my $STYLE = 'paper';
#my $STYLE = 'analysis';


my $COLORSCHEME = 'grayscale';
#my $COLORSCHEME = 'colors';


my @USEFUL_STATISTICS = qw(
  Min Max 
  Mean SD 
  Deg_Over_1.0sig Deg_Over_0.3sig
  Area Area_0.0s Area_0.3s
);

my $PLOT_USEFUL  = 1;
my $PLOT_FILTERS = 0;
my %FILTERS = ();
my @FILTERS = qw(rmsd);
#my @FILTERS = qw(scale_sd n05 p0sigma);
my @OLD_AND_DEAD = qw(scale mean n10 norm shift mmsd mpsd p3sigma);

my %AA = %Utility::AA;

###############################################################################
#                                                            Parse Command Line
###############################################################################

my $CONF              = (grep { not /^-/    } @ARGV)[0] or &usage;
my $DO_DH_PLOTS       = (grep { /-dih/      } @ARGV) ? 1 : 0;
my $INCLUDE_ALA       = (grep { /-ala/      } @ARGV) ? 1 : 0;
my $DO_COMPARISONS    = (grep { /-comp/     } @ARGV) ? 1 : $DO_DH_PLOTS;
my $DO_MATRIX_PLOT    = (grep { /-nomatrix/ } @ARGV) ? 0 : 1;

my $SHOW_LEGEND       = (grep { /-nolegend/ } @ARGV) ? 0 : 1;
my $SHOW_TITLE        = (grep { /-notitle/  } @ARGV) ? 0 : 1;
my $SHOW_LABELS       = (grep { /-nolabel/  } @ARGV) ? 0 : 1;
my $SHOW_AXES         = (grep { /-noaxes/   } @ARGV) ? 0 : 1;

my $DENSITY_FILTER    = (grep { /-den=/     } @ARGV)[0];
my $FILTER_RES        = (grep { /-res=/     } @ARGV)[0];
my $CHI_MAX           = (grep { /-chi=/     } @ARGV)[0];
my $VERBOSE           = (grep { /-v/        } @ARGV) ? 1 : 0;

#============================================================ Process Arguments

#----------------- Selected Residues
my @FILTER_RES;
my @FILTER_FLAG;
if ( defined $FILTER_RES ) { 
  $FILTER_RES =~ /^-res=(.*)$/ 
    or usage("Bad residue range '$FILTER_RES'.\n");
  @FILTER_RES = Utility::Parse_Range($1);
  map { $FILTER_FLAG[$_] = 1 } @FILTER_RES;
}

#----------------- Chi Angle Limit
if ( defined $CHI_MAX ) { 
  $CHI_MAX =~ /^-chi=(\d+)$/ 
    or usage("Bad chi range '$CHI_MAX'.\n");
  $CHI_MAX = int($1);
} else { 
  $CHI_MAX = 4;
}

#----------------- Density Filter
if ( defined $DENSITY_FILTER ) { 
  $DENSITY_FILTER =~ /^-den=([\d\.\-]+)$/ 
    or usage("Bad density '$DENSITY_FILTER'.\n");
  $DENSITY_FILTER = $1;
}

#----------------- Include Alanines for Density Reference
if ( $INCLUDE_ALA ) { 
  $AA{'ALA'}{sc} = ['chi1'];
  $AA{'ala'}{sc} = ['chi1'];
}

#===================================================== Parse Configuration File

print "Reading data files.\n";
my ($TITLE, $OUTDIR, @PDBS) = parse_conf_file ( $CONF );
if ( not $TITLE ) { $TITLE = 'ringer_out' }
if ( not defined $OUTDIR ) { 
  $OUTDIR = Utility::Clean_Filename($TITLE);
}
if ( not ( -e $OUTDIR and -d $OUTDIR ) ) { 
  mkdir $OUTDIR or die "Couldn't open ringer output directory, $OUTDIR. $!";
}
chdir $OUTDIR or die "Couldn't chir to output driectory, $OUTDIR. $!";

print "Output Directory: $OUTDIR\n";



&initialize_plots;




##############################################################################
#                         Extract consensus sequence, residue numbers, binning
##############################################################################
print "Validating consensus sequences\n";
my @IDS = ();
my @SEQ = ();
my @CHI = ();
my $temp_max = Utility::Max ( 
                  map { $_->{'length'} + $_->{'offset'} } 
                  @PDBS
                );


#================================================= Determine Consensus Sequence
foreach my $k ( 0..$temp_max ) { 
  my %seq = ();
  my $chi_max =  0;
  # Count residue occurances at each location
  foreach my $i ( 0..$#PDBS ) { 
    my %pdb_i = %{$PDBS[$i]};
    my $k_i   = $k + $pdb_i{'offset'};
    my $res   = $pdb_i{'seq'}[$k_i];
    my $chi   = $pdb_i{'chi'}[$k_i];
    next if not defined $res;

    $seq{$res}++;
    if ( $chi > $chi_max ) { $chi_max = $chi };
  }
  # Set consensus to most popular
  $SEQ[$k] = (sort { -($seq{$a}<=>$seq{$b}) } keys %seq)[0];
  $CHI[$k] = $chi_max;
  if ( defined $CHI_MAX and $CHI_MAX < $CHI[$k] ) { $CHI[$k] = $CHI_MAX }

  #print "Concensus for $k is $SEQ[$k].\n";

  # Double check for disagreements just to notes them.
  foreach my $i ( 0..$#PDBS ) { 
    my %pdb_i = %{$PDBS[$i]};
    my $k_i   = $k + $pdb_i{'offset'};
    my $res   = $pdb_i{'seq'}[$k_i];
    next if not defined $res;

    if ( $res ne $SEQ[$k] ) { 
      printf "Sequence Disagreement: %10s %10s Res%4i %4s vs. %4s\n",
             $pdb_i{'name'},
             $pdb_i{'desc'},
             $k,
             $SEQ[$k],
             $res,
             ;
    }
  }
}

#================= Get Sequence Extent and Filter Residues
my @DEG = @{$PDBS[0]{'deg'}};
my @SEQ_LIST = ();
my ($SEQ_MIN, $SEQ_MAX) = (grep { defined $SEQ[$_] } 0..$#SEQ)[0,-1];
print "Min, Max: $SEQ_MIN, $SEQ_MAX\n";
if ( @FILTER_FLAG ) { @SEQ_LIST = grep { $FILTER_FLAG[$_] } 0..$#SEQ } 
else                { @SEQ_LIST = grep { defined $SEQ[$_] } 0..$#SEQ }


#================= Remove Residues With No Chi Angles
my @SEQ_LIST_FINAL = ();
foreach my $i ( 0..$#SEQ_LIST ) { 
  my $k = $SEQ_LIST[$i];
  next if not defined $SEQ[$k];
  next if $CHI[$k] < 1;
  push @SEQ_LIST_FINAL, $k;
}
@SEQ_LIST = @SEQ_LIST_FINAL;

#================= Fake Density Data for Skipped Residues
my $FAKE_DENSITY = [ (-10) x scalar @DEG ];





##############################################################################
#                                                       PRECOMPUTE COMPARISONS
##############################################################################


##### Precompute common arrays


foreach my $i ( 0..$#PDBS ) { 
  my $pdb = $PDBS[$i];
  my @seq = @{$pdb->{'seq'}};
  my @chi = @{$pdb->{'chi'}};
  my @den = @{$pdb->{'den'}};
  print "Precomputing Electron Density Filters $i $pdb->{'name'}\n";
  foreach my $k ( @SEQ_LIST ) { 
    
    # Get the sequence
    my $ki = $k + $pdb->{'offset'};

    # Skip if misaligned
    next if not defined $seq[$ki];
    # Skip if not in filter
    #next if not $FILTER_FLAG[$ki];
   
    foreach my $c ( 1..$CHI[$k]) { 
      
      next if not defined $den[$ki][$c];
      my @m   = @{$den[$ki][$c]};
      my @n;

      #==== Basic Statistics
      my ($mean, $sd) = Utility::Mean_and_Stddev(\@m);
      my $min = Utility::Min(\@m);
      my $max = Utility::Max(\@m);

      if ( defined $DENSITY_FILTER and $max < $DENSITY_FILTER ) { 
        #print "Residue $k $seq[$k]  Chi $c failed density test ($max < $DENSITY_FILTER)\n";
        #$den[$ki][$c] = $FAKE_DENSITY;
        #@m = @$FAKE_DENSITY;
        $den[$ki][$c] = undef;
        next;
      }

      if ( 0 ) { 
        printf "  K %-4i  RES %-4s CHI %-4i MEAN %10.6f   SD %10.6f %s\n",
                  $k, $seq[$k], $c, $mean, $sd,
                  (join ',', map { sprintf "%10.6f", $_ } @m[0, 60, 120, 180, 240, 300]),
                  ;
      }


      #==== Filter Peak Electron Density Below 3
 
      # Normalize Area to 1
      @n = @m;
      Utility::Normalize_Histogram(\@n);
      $pdb->{'norm'}[$k][$c] = [@n];

      # Normalize Min to 0 and then Sum to 1
      @n = @m;
      Utility::Histogram_Shift_Min_to_X(\@n);
      $pdb->{'shift'}[$k][$c] = [@n];
      Utility::Normalize_Histogram(\@n);
      $pdb->{'scale'}[$k][$c] = [@n];


   
      my $threshold;
      #===== Noise over 3Sigma
      $threshold = 0.3;
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      if ( scalar grep { $_ != 0 } @n ) { 
        Utility::Histogram_Shift_Min_to_X(\@n);
        Utility::Normalize_Histogram(\@n);
        $pdb->{'p3sigma'}[$k][$c] = [@n];
      }

      $threshold = 0.0;
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      if ( scalar grep { $_ != 0 } @n ) { 
        Utility::Histogram_Shift_Min_to_X(\@n);
        Utility::Normalize_Histogram(\@n);
        $pdb->{'p0sigma'}[$k][$c] = [@n];
      }
       
      #===== MMSD 
      $threshold = $mean-$sd; 
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      Utility::Histogram_Shift_Min_to_X(\@n);
      Utility::Normalize_Histogram(\@n);
      $pdb->{'mmsd'}[$k][$c] = [@n];

      #===== MPSD
      $threshold = $mean+$sd; 
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      Utility::Histogram_Shift_Min_to_X(\@n);
      Utility::Normalize_Histogram(\@n);
      $pdb->{'mpsd'}[$k][$c] = [@n];

      #===== MEAN
      $threshold = $mean; 
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      Utility::Histogram_Shift_Min_to_X(\@n);
      Utility::Normalize_Histogram(\@n);
      $pdb->{'mean'}[$k][$c] = [@n];

      #===== Scale SD
      @m = @{$pdb->{'scale'}[$k][$c]};
      ($mean, $sd) = Utility::Mean_and_Stddev(\@m);
      $threshold = $sd;
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      Utility::Histogram_Shift_Min_to_X(\@n);
      Utility::Normalize_Histogram(\@n);
      $pdb->{'scale_sd'}[$k][$c] = [@n];


      #===== Noise 5%
      @n = sort { $a <=> $b } @m;
      my $noise = 0.05;
      my $total = 0.00;
      for my $i ( 0..$#n ) { 
        $total += $n[$i];
        if ( $total > $noise ) { 
          $threshold = $n[$i];
          last;
        }
      }
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      Utility::Histogram_Shift_Min_to_X(\@n);
      Utility::Normalize_Histogram(\@n);
      $pdb->{'n05'}[$k][$c] = [@n];

      #===== Noise 10%
      @n = sort { $a <=> $b } @m;
      $noise = 0.10;
      $total = 0.00;
      for my $i ( 0..$#n ) { 
        $total += $n[$i];
        if ( $total > $noise ) { 
          $threshold = $n[$i];
          last;
        }
      }
      @n = map { $_ < $threshold ? 0 : $_-$threshold } @m;
      Utility::Histogram_Shift_Min_to_X(\@n);
      Utility::Normalize_Histogram(\@n);
      $pdb->{'n10'}[$k][$c] = [@n];



    }
  }
}


##############################################################################
#                                                    DIHEDRAL ANGLE STATISTICS
##############################################################################


foreach my $i ( 0..$#PDBS ) { 
  my $pdb = $PDBS[$i];
  my @seq = @{$pdb->{'seq'}};
  my @den = @{$pdb->{'den'}};
  print "Analyzing Dihedral Angles $i $pdb->{'name'}\n";
  foreach my $k ( @SEQ_LIST ) { 
    
    # Get the sequence
    my $ki   = $k + $pdb->{'offset'};
    my $res = $seq[$ki];
    # Skip if misaligned
    next if not defined $res;
   
    foreach my $c ( 1..$CHI[$k]) { 

      next if not defined $den[$ki][$c];
      my @m   = @{$den[$ki][$c]};
      my @n;

      #==== Basic Statistics
      my ($mean, $sd) = Utility::Mean_and_Stddev(\@m);
      my $min     = Utility::Min(\@m);
      my $max     = Utility::Max(\@m);
      my $range   = $max - $min;
      
      my $sum_reg = Utility::Sum(\@m);
      my $sum_abs = Utility::Sum( map { abs($_) } @m );
      my $sum_pos = Utility::Sum( grep { $_ > 0 } @m );
      my $sum_one = Utility::Sum( grep { $_ > 1 } @m );

      my $fg_mean = (scalar grep { $_>$mean } @m) / (scalar @m) * 100.0;
      my $fg_sd   = (scalar grep { $_>$sd   } @m) / (scalar @m) * 100.0;
      my $fg_zero = (scalar grep { $_>0     } @m) / (scalar @m) * 100.0;
      my $fg_one  = (scalar grep { $_>1     } @m) / (scalar @m) * 100.0;
      
      my $dg_mean = (scalar grep { $_>$mean } @m) / (scalar @m) * 360.0;
      my $dg_sd   = (scalar grep { $_>$sd   } @m) / (scalar @m) * 360.0;
      my $dg_zero = (scalar grep { $_>0     } @m) / (scalar @m) * 360.0;
      my $dg_one  = (scalar grep { $_>1     } @m) / (scalar @m) * 360.0;
      my $dg_03s  = (scalar grep { $_>0.3   } @m) / (scalar @m) * 360.0;

      if ( 0 ) { 
        printf "  K %-4i  RES %-4s CHI %-4i MEAN %10.6f   SD %10.6f %s\n",
                  $k, $res, $c, $mean, $sd,
                  (join ',', map { sprintf "%10.6f", $_ } @m[0, 60, 120, 180, 240, 300]),
                  ;
      }


      @n = @{$pdb->{'scale_sd'}[$k][$c]};
      my ($ssd_mean, $ssd_sd) = Utility::Mean_and_Stddev(\@n);
      my $ssd_min = Utility::Min(\@n);
      my $ssd_max = Utility::Max(\@n);
      my $ssd_fg_mean = (scalar grep { $_>$ssd_mean } @n) / (scalar @n) * 100.0;
      my $ssd_fg_sd   = (scalar grep { $_>$ssd_sd   } @n) / (scalar @n) * 100.0;
      my $ssd_fg_zero = (scalar grep { $_>0         } @n) / (scalar @n) * 100.0;


      #----------- Signal cannot be better than previous Chi angles
      my $ISigI = $max;
      if ( $c > 1 ) { 
        my $ISigI_Last = $pdb->{'dh'}{'ISigI'}[$ki][$c-1];
        if ( $ISigI_Last < $ISigI ) { 
          $ISigI = $ISigI_Last;
        }
      }
      $pdb->{'dh'}{'ISigI'          }[$ki][$c] = $ISigI;

      $pdb->{'dh'}{'Mean'           }[$ki][$c] = $mean;
      $pdb->{'dh'}{'Min'            }[$ki][$c] = $min;
      $pdb->{'dh'}{'Max'            }[$ki][$c] = $max;
      $pdb->{'dh'}{'SD'             }[$ki][$c] = $sd;
      $pdb->{'dh'}{'Deg_Over_1.0sig'}[$ki][$c] = $dg_one;
      $pdb->{'dh'}{'Deg_Over_0.3sig'}[$ki][$c] = $dg_03s;



      $pdb->{'dh'}{'raw_mean'     }[$ki][$c] = $mean;
      $pdb->{'dh'}{'raw_sd'       }[$ki][$c] = $sd;
      $pdb->{'dh'}{'raw_min'      }[$ki][$c] = $min;
      $pdb->{'dh'}{'raw_max'      }[$ki][$c] = $max;
      $pdb->{'dh'}{'raw_range'    }[$ki][$c] = $range;
      $pdb->{'dh'}{'raw_sum_reg'  }[$ki][$c] = $sum_reg;
      $pdb->{'dh'}{'raw_sum_pos'  }[$ki][$c] = $sum_pos;
      $pdb->{'dh'}{'raw_sum_abs'  }[$ki][$c] = $sum_abs;
      $pdb->{'dh'}{'raw_sum_one'  }[$ki][$c] = $sum_one;
      $pdb->{'dh'}{'raw_fg_mean'  }[$ki][$c] = $fg_mean;
      $pdb->{'dh'}{'raw_fg_sd  '  }[$ki][$c] = $fg_sd;
      $pdb->{'dh'}{'raw_fg_zero'  }[$ki][$c] = $fg_zero;
      $pdb->{'dh'}{'raw_fg_one'   }[$ki][$c] = $fg_one;

      $pdb->{'dh'}{'ssd_mean'    }[$ki][$c] = $ssd_mean;
      $pdb->{'dh'}{'ssd_sd'      }[$ki][$c] = $ssd_sd;
      $pdb->{'dh'}{'ssd_max'     }[$ki][$c] = $ssd_max;
      $pdb->{'dh'}{'ssd_fg_mean' }[$ki][$c] = $ssd_fg_mean;
      $pdb->{'dh'}{'ssd_fg_sd'   }[$ki][$c] = $ssd_fg_sd;
      $pdb->{'dh'}{'ssd_fg_zero' }[$ki][$c] = $ssd_fg_zero;

      #===== Area Statistics
      #----- Sums of the bins of the histogram. Some thresholded above a value
      my $area    = 0;
      my $area0_0 = 0;
      my $area0_3 = 0;
      for my $m ( @m ) { 
        $area += $m;
        if ( $m > 0.0 ) { $area0_0 += $m }
        if ( $m > 0.3 ) { $area0_3 += $m }
      }
      $pdb->{'dh'}{'Area'            }[$ki][$c] = $area;
      $pdb->{'dh'}{'Area_0.0s'}[$ki][$c] = $area0_0;
      $pdb->{'dh'}{'Area_0.3s'}[$ki][$c] = $area0_3;
      #print "Area     = $area\n";
      #print "Area 0.0 = $area0_0\n";
      #print "Area 0.3 = $area0_3\n";
    }
  }
}


#my @KEYS_DH = keys %{$PDBS[0]{'dh'}};
#foreach my $key ( @KEYS_DH ) { 
#  print_dihedral_analysis(\@PDBS, \@SEQ, $SEQ_MIN, $SEQ_MAX, $key, -1);
#}



##############################################################################
#                                                          DENSITY COMPARISONS
##############################################################################


if ( $DO_COMPARISONS ) { 


print "Comparing Densities\n";
foreach my $i ( 0..$#PDBS ) { 
  my $pdb_i = $PDBS[$i];
  my @seq_i = @{$pdb_i->{'seq'}};
  my @den_i = @{$pdb_i->{'den'}};
  print "  I: $i $pdb_i->{'name'}\n";
  foreach my $j ( ($i+1)..$#PDBS ) { 
    my $pdb_j = $PDBS[$j];
    my @seq_j = @{$pdb_j->{'seq'}};
    my @den_j = @{$pdb_j->{'den'}};
    print "    J: $j $pdb_j->{'name'}\n";
 
    my $identical = 0;
    foreach my $k ( @SEQ_LIST ) { 
      
      # Get the sequence
      my $k_i   = $k + $pdb_i->{'offset'};
      my $k_j   = $k + $pdb_j->{'offset'};
      my $res_i = $seq_i[$k_i];
      my $res_j = $seq_j[$k_j];

      # Skip if misaligned
      next if not defined $res_i;
      next if not defined $res_j;
      #next if $res_i ne $res_j;
  
      $identical++;
     
      foreach my $c ( 1..$CHI[$k] ) { 
        
        # Skip if one residue is smaller
        next if $c > $pdb_i->{'chi'}[$k_i];
        next if $c > $pdb_j->{'chi'}[$k_j];
        next if not defined $den_i[$k_i][$c];
        next if not defined $den_j[$k_j][$c];

  
        my @m_i   = @{$den_i[$k_i][$c]};
        my @m_j   = @{$den_j[$k_j][$c]};
        my (@n_i, @n_j);
        my ($n_i, $n_j);

        
        # Correlation Coefficient
        $FILTERS{'cor' }[$i][$j][$k][$c] = Utility::Correlation(\@m_i, \@m_j); 
        $FILTERS{'rmsd'}[$i][$j][$k][$c] = Utility::RMSD(\@m_i, \@m_j);
      
        # Displacement/Similarity Filters
        #foreach my $f ( @FILTERS ) { 
        #  #print "$f\n";
        #  $n_i = $pdb_i->{$f}[$k_i][$c];
        #  $n_j = $pdb_j->{$f}[$k_j][$c];
        #  my $dis;
        #  if ( not defined $n_i or not defined $n_j ) { 
        #    $dis = 0;
        #  } else { 
        #    $dis = Utility::Displacement_Raw($n_i, $n_j);
        #  }
        #  $FILTERS{$f}[$i][$j][$k][$c] = 1 - ($dis/2.0);
        #}
          

      }
    }
    $identical /= scalar @SEQ_LIST;
    $IDS[$i][$j] = $identical;
    printf "      Identity: %.2f%%\n", $identical * 100;
  }
}


print_comparison ( \@PDBS, \@SEQ_LIST, \@SEQ, $SEQ_MIN, $SEQ_MAX, $FILTERS{'cor'}, "cor", $TITLE );
foreach my $f ( @FILTERS ) { 
  print_comparison ( \@PDBS, \@SEQ_LIST, \@SEQ, $SEQ_MIN, $SEQ_MAX, $FILTERS{$f}, $f, $TITLE );
} 

} # /END COMPARISONS

###############################################################################
#                                                                  Ringer Plots
###############################################################################


if ( $DO_DH_PLOTS ) { 


  print "Plotting Dihedral Distributions\n";

  foreach my $k ( @SEQ_LIST ) { 
    
    my $res = $SEQ[$k];
    next if not defined $res;
    
    foreach my $c ( 1..$CHI[$k] ) { 

      my $title  = sprintf "%s - %s %d - Chi%d", $TITLE, $res, $k, $c;
      if ( $STYLE eq 'paper' ) { 
        $title  = sprintf "%s%d {/Symbol c}%d", ucfirst(lc($res)), $k, $c;
        #$title  = sprintf "%s %d - Chi %d", ucfirst(lc($res)), $k, $c;
      }
      my $file   = sprintf "%s - %04d %s Chi%d", $TITLE, $k, $res, $c;
      my @arrays = ();
      my @labels = ();
      my @annos  = ();
      my @shifts = ();
      my @scales = ();
      my @norms  = ();
      my @filter = ();
      my @f_sd   = ();
      my @f_mmsd = ();
      my @f_mpsd = ();
      my @f_mean = ();

      print "Residue $k $res Chi$c\n";

      foreach my $i ( 0..$#PDBS ) { 
        my %pdb = %{$PDBS[$i]};
        my $off = $pdb{'offset'};
        #push @labels, "$pdb{'name'}/$pdb{'desc'}";
        #push @labels, $pdb{'name'};
        my ( $min, $max, $mean, $sd, $isigi ) = map { $pdb{dh}{$_}[$k][$c] || 0.0 } qw(Min Max Mean SD ISigI);
        if ( $STYLE ne 'paper' ) { 
          push @labels, sprintf "%s - %5.2f - %5.2f", $pdb{'name'}, $max, $isigi;
        } else { 
          push @labels, $pdb{'name'};
        }
        #push @labels, sprintf "%s - %5.2f/%5.2f - %5.2f/%5.2f - %5.2f", $pdb{'name'}, 
        #                                          $min, $max, $mean, $sd, ($mean/$sd);
        push @arrays, $pdb{'den'}[$k+$off][$c];
      }
      
      # Find worst ringer correlation and annotate it
      my $best_cor  = 0;
      my $best_i    = 0;
      my $best_j    = 0;
      my $worst_cor = 1;
      my $worst_i   = 0;
      my $worst_j   = 0;
      foreach my $i ( 0..$#PDBS ) { 
        foreach my $j ( $i+1..$#PDBS ) { 
          my $cor = $FILTERS{'cor'}[$i][$j][$k][$c] || $worst_cor;
          if ( $cor > $best_cor ) { 
            $best_cor = $cor;
            $best_i   = $i;
            $best_j   = $j;
          }
          if ( $cor < $worst_cor ) { 
            $worst_cor = $cor;
            $worst_i   = $i;
            $worst_j   = $j;
          }
        }
      }


      if ( scalar @PDBS > 1 ) { 
        push @annos, "Extremes: $PDBS[0]{'name'}-$PDBS[-1]{'name'}   " . 
                     join " ", 
                      #(sprintf "Cor: %5.2f", $FILTERS{'cor'}[0][-1][$k][$c] || 0) . 
                      (map { sprintf "%s: %5.2f", $_, ($FILTERS{$_}[0][-1][$k][$c] || 0) } 'cor', @FILTERS),
                     ;

        if ( scalar @PDBS > 2 ) { 
          push @annos, "Best: $PDBS[$best_i]{'name'}-$PDBS[$best_j]{'name'}   " . 
                       join " ", 
                        #(sprintf "Cor: %5.2f", $FILTERS{'cor'}[0][-1][$k][$c] || 0) . 
                        (map { sprintf "%s: %5.2f", $_, ($FILTERS{$_}[$best_i][$best_j][$k][$c] || 0) } 'cor', @FILTERS),
                       ;
                       #(sprintf "Cor: %5.2f", $FILTERS{'cor'}[$best_i][$best_j][$k][$c] || 0) . 
                       #(map { sprintf "%s: %5.2f", $_, ($FILTERS{$_}[$best_i][$best_j][$k][$c] || 0) } @FILTERS),
                       #"",
                       #;

          push @annos, "Worst: $PDBS[$worst_i]{'name'}-$PDBS[$worst_j]{'name'}   " . 
                       join " ", 
                        #(sprintf "Cor: %5.2f", $FILTERS{'cor'}[0][-1][$k][$c] || 0) . 
                        (map { sprintf "%s: %5.2f", $_, ($FILTERS{$_}[$worst_i][$worst_j][$k][$c] || 0) } 'cor', @FILTERS),
                       ;
                       #(sprintf "Cor:  %5.2f", $FILTERS{'cor'}[$worst_i][$worst_j][$k][$c] || 0) . 
                       #(map { sprintf "%s: %5.2f", $_, ($FILTERS{$_}[$worst_i][$worst_j][$k][$c] || 0) } @FILTERS),
                       #"",
                       #;
        }
      } else { 

        my @list = ();
        if ( $PLOT_USEFUL ) { 
          @list = @USEFUL_STATISTICS
        } else { 
          @list = sort keys %{$PDBS[0]{'dh'}};
        }
        push @annos,
          map { sprintf "%s: %5.2f", $_, $PDBS[0]{'dh'}{$_}[$k][$c] } @list
      }

      plot_dihedral_distributions ( 
        $title,
        $file,
        \@DEG,
        \@arrays,
        \@labels,
        \@annos,
        uc $res,
        $c,
      );

      foreach my $f ( @FILTERS ) { 

        if ( $PLOT_FILTERS ) { 

          plot_dihedral_distributions ( 
            $title . " - $f",
            $file  . " - $f",
            \@DEG,
            [ map { $PDBS[$_]{$f}[$k+$PDBS[$_]{'offset'}][$c] } 0..$#PDBS ],
            \@labels,
            \@annos,
            uc $res,
            $c,
          );
        }
      }
    }
  }
}   



##############################################################################
##############################################################################
#                                                                  SUBROUTINES
##############################################################################
##############################################################################

##############################################################################
#                                                 PRINT DIHEDRAL ANALYSIS FILE
##############################################################################


sub print_dihedral_analysis { 

  my $pdbs      = $_[0];
  my $seq       = $_[1];
  my $seq_min   = $_[2];
  my $seq_max   = $_[3];
  my $field     = $_[4];
  my $default   = $_[5];
  if ( not defined $default ) { $default = 'NA' }


  my $file_dh_all  = "dh.$field.all.dat";
  my $file_dh_chi1 = "dh.$field.chi1.dat";

  open DH_ALL,  ">$file_dh_all"  or die "Couldn't open file dh, '$file_dh_all'. $!";
  open DH_CHI1, ">$file_dh_chi1" or die "Couldn't open file dh, '$file_dh_chi1'. $!";
 
 
  #==== Print Header
  print DH_ALL "RES\tCHI";
  print DH_CHI1 "RES";
  foreach my $i ( 0..$#$pdbs ) { 
    my $pdb = $pdbs->[$i];
    my $name = $pdb->{'name'};
    print DH_ALL "\t$name";
    print DH_CHI1 "\t$name";
  } 
  print DH_ALL "\n";
  print DH_CHI1 "\n";


  #==== Print Data
  foreach my $k ( $seq_min..$seq_max ) { 
    my $res = $seq->[$k];
    #next if not defined $res;

    my $angles = scalar @{$Utility::AA{lc $res}{'sc'}};
    #foreach my $c ( 1..1 ) {
    foreach my $c ( 1..$angles ) {
      print DH_ALL "$k\t$c";
      print DH_CHI1 "$k" if $c == 1;
      foreach my $i ( 0..$#$pdbs ) { 
        my $pdb = $pdbs->[$i];
        my $name = $pdb->{'name'};
        my $ki = $k + $pdb->{'offset'};
        my $x = $pdb->{'dh'}{$field}[$ki][$c];
        if ( defined $x ) { $x = sprintf "%.4f", $x }
        else              { $x = $default }
        printf DH_ALL "\t$x";
        if ( $c == 1 ) { 
          printf DH_CHI1 "\t$x";
        }
      }
      print DH_ALL "\n";
      print DH_CHI1 "\n" if $c == 1;
    }
  }
  close DH_ALL  or die "Couldn't close file dh, '$file_dh_all'. $!";
  close DH_CHI1 or die "Couldn't close file dh, '$file_dh_chi1'. $!";
}






###############################################################################
#                                                     PRINT COMPARISON MATRICES
#------------------------------------------------------------------------------
# This subroutine prints several files for the comparison:
# 
#   - Matrix Comparison
#   matrix.COMPARE.SET.CHI.txt
#   matrix.COMPARE.SET.CHI.transposed.txt
# 
#   - Linear Comparison
#   comparison.COMPARISON.txt
#
#   - Pymol ColorSchemes
#   color.PDB1.PDB2.COMPARISON.ca.pymol
#   color.PDB1.PDB2.COMPARISON.chi.pymol
###############################################################################


sub print_comparison { 

  my $pdbs      = $_[0];
  my $seq_list  = $_[1];
  my $seq       = $_[2];
  my $seq_min   = $_[3];
  my $seq_max   = $_[4];
  my $data      = $_[5];
  my $compare   = $_[6];
  my $crystal   = $_[7];

  print_comparison_matrix(@_);
  print_comparison_linear(@_);
  print_comparison_pymol (@_);
  return;

}

##############################################################################
#                                                      PRINT COMPARISON LINEAR
#-----------------------------------------------------------------------------
# Going to split these out for easier comprehension
##############################################################################

sub print_comparison_linear { 

  my $pdbs      = $_[0];
  my $seq_list  = $_[1];
  my $seq       = $_[2];
  my $seq_min   = $_[3];
  my $seq_max   = $_[4];
  my $data      = $_[5];
  my $compare   = $_[6];
  my $crystal   = $_[7];


  my $baseline = 1.0;
  my $file_linear = "comparison.$compare.txt";
  open LIN, ">$file_linear" or die "Couldn't open file, '$file_linear'. $!";


  ########################################################### PRINT HEADER ROW
  print LIN "PDB1\tPDB2\tNAME1\tNAME2\tRESI\tRESN\tCHI\tCOMP\tISIGI1\tISIGI2\n";

  foreach my $i ( 0..$#$pdbs ) { #====================================== PDB I
    my $pdb_i   = $pdbs->[$i];
    my $name_i  = $pdb_i->{name};
    $name_i =~ s/\s/_/g;
    foreach my $j ( ($i+1)..$#$pdbs ) { #=============================== PDB J
      my $pdb_j   = $pdbs->[$j];
      my $name_j  = $pdb_j->{name};
      $name_j =~ s/\s/_/g;

      #========== Loop Over Residues
      foreach my $k ( @$seq_list ) {
      #foreach my $k ( $seq_min..$seq_max ) { 
        my $resn = $seq->[$k];
        my $angles = scalar @{$AA{lc $resn}{'sc'}};

        #======== Loop over Chi Angles
        foreach my $chi ( 1..$angles ) { 
          my $c = $data->[$i][$j][$k][$chi];
          next if not defined $c;

          my $isigi_i = $pdb_i->{dh}{ISigI}[$k][$chi] || 0.0;
          my $isigi_j = $pdb_j->{dh}{ISigI}[$k][$chi] || 0.0;

          printf LIN "%i\t%i\t%s\t%s\t%i\t%s\t%i\t%.4f\t%.2f\t%.2f\n",
                      $i,
                      $j,
                      $name_i,
                      $name_j,
                      $k,
                      $resn,
                      $chi,
                      $c,
                      $isigi_i,
                      $isigi_j,
                      ;                      
        }
      }
    }
  }
  close LIN or die "Couldn't close file, '$file_linear'. $!";
}



##############################################################################
#                                                      PRINT COMPARISON MATRIX
#-----------------------------------------------------------------------------
# Going to split these out for easier comprehension
##############################################################################

sub print_comparison_matrix { 

  my $pdbs      = $_[0];
  my $seq_list  = $_[1];
  my $seq       = $_[2];
  my $seq_min   = $_[3];
  my $seq_max   = $_[4];
  my $data      = $_[5];
  my $compare   = $_[6];
  my $crystal   = $_[7];

  my $baseline = 1.0;
  #=============== Setup Files
  my %mtx = ();
  my @mtx = ();

  my @sets = ();
  if ( scalar @$pdbs <= 1 ) { return }
  elsif ( scalar @$pdbs == 2 ) { @sets = ('step') }
  elsif ( scalar @$pdbs  > 2 ) { @sets = qw(all step first last) }

  for my $set ( @sets ) { 
    for my $chi ( qw(chi1 chix chiz) ) { 
      $mtx{$set}{$chi} = {
        file  => "matrix.$compare.$set.$chi.txt",
        trans => "matrix.$compare.$set.$chi.transposed.txt",
        set   => $set,
        chi   => $chi,
      };
      my $m = $mtx{$set}{$chi};
      # Open Transposed File Pointers
      open $m->{FPT}, ">$m->{trans}" 
        or die "Couldn't open matrix file, '$m->{trans}'. $!";
      push @mtx, $m;
    }
  }

  my %HEAD_WRITTEN = ();

  #=============== Group File Pointer Sets
  my @chi1 = grep { $_->{chi} eq 'chi1'  } @mtx;
  my @chix = grep { $_->{chi} eq 'chix'  } @mtx;
  my @chiz = grep { $_->{chi} eq 'chiz'  } @mtx;
  #my @all  = grep { $_->{set} eq 'all'   } @mtx;
  #my @r100 = grep { $_->{set} eq 'first' } @mtx;
  #my @r260 = grep { $_->{set} eq 'last'  } @mtx;
  #my @step = grep { $_->{set} eq 'step'  } @mtx;

  #=============== Print Header Start
  map { print {$_->{FPT}} 'CHI_ANGLE' } @chix;
  map { print {$_->{FPT}} 'RESIDUE'   } @chi1; 
  map { print {$_->{FPT}} 'RESIDUE'   } @chiz; 

  #=============== Print Header - Sequence RESI + Fractional Chi Angle
  foreach my $k ( @$seq_list ) { 
    my $resn = $seq->[$k];
    #next if $resn eq 'ALA' and not $INCLUDE_ALA;
    #map { print {$_->{FPT}} "\t$k.0" } @chix;
    #foreach my $m (@mtx) { 
    #  if ($m->{chi} eq 'chi1') { print {$_->{FPT}} "\t$k" }
    #  if ($m->{chi} eq 'chiz') { print {$_->{FPT}} "\t$k" }
    #}
    map { print {$_->{FPT}} "\t$k" } @chi1;
    map { print {$_->{FPT}} "\t$k" } @chiz;

    if ( not defined $resn ) { 
      $HEAD_WRITTEN{$k}{0} = 1;
      map { print {$_->{FPT}} "\t$k.0" } @chix;
      next;
    }

    my $angles = scalar @{$AA{lc $resn}{'sc'}};
    if ( $angles == 0 ) { $angles = 1 }
      
    foreach my $chi ( 1..$angles ) { 
      #print "Writing header for $k $resn $chi/$angles\n";
      $HEAD_WRITTEN{$k}{$chi} = 1;
      my $fraction = sprintf "\t%.4f", $k+(($chi-1)/$angles);
      map { print {$_->{FPT}} $fraction } @chix;
    }
  }
  map { print {$_->{FPT}} "\n" } @mtx;

  
  #============================================================ FIRST PDB LOOP
  foreach my $i ( 0..$#$pdbs ) { 
    my $name_i = $pdbs->[$i]{'name'};
    $name_i =~ s/\s/_/g;
    #========================================================== FIRST PDB LOOP
    foreach my $j ( ($i+1)..$#$pdbs ) { 
      my $name_j = $pdbs->[$j]{'name'};
      $name_j =~ s/\s/_/g;
       
      #========== Filter Which Matrix Rows/Columns Get Printed
      #           Print Row Headers (Comparison Titles)
      my $title = "$name_i\:\:$name_j";
      my @m = ();
      foreach my $m ( @mtx ) { 
        if    ( ($i == ($j-1)  ) and ($m->{set} eq 'step' ) ) { push @m, $m }
        elsif ( ($i ==     0   ) and ($m->{set} eq 'first') ) { push @m, $m }
        elsif ( ($j == $#$pdbs ) and ($m->{set} eq 'last' ) ) { push @m, $m }
        elsif (                      ($m->{set} eq 'all'  ) ) { push @m, $m }
      }
      map { print {$_->{FPT}} $title } @m;

      #========== Loop Over Residues
      foreach my $k ( @$seq_list ) {
        my $resn = $seq->[$k];

        #--- Default missing residue
        if ( (not defined $resn) or ($resn eq 'ALA') or ($resn eq 'GLY')) { 
        #if ( (not defined $resn) ) { 
          map { printf {$_->{FPT}} "\t%.4f", $baseline } @m;
          next;
        }

        #--- Default missing chi angles
        my $angles = scalar @{$AA{lc $resn}{'sc'}};
        if ( $angles == 0 ) { 
          map { printf {$_->{FPT}} "\t%.4f", $baseline } @m;
          next;
        }

        #======== Loop over Chi Angles
        foreach my $chi ( 1..$angles ) { 
          my $c = $data->[$i][$j][$k][$chi];
          if ( not defined $c ) { $c = $baseline }

          foreach my $m ( @m ) {  
            my $print = 0;
            if ( ($chi == 1      ) and ($m->{chi} eq 'chi1') ) { $print = 1 } 
            if ( ($chi == $angles) and ($m->{chi} eq 'chiz') ) { $print = 1 } 
            if (                       ($m->{chi} eq 'chix') ) { $print = 1 }
            if ( $print ) { 
              printf {$m->{FPT}} "\t%.4f", $c;
            }
          }
        }
      }
      map { print {$_->{FPT}} "\n" } @m;
    }
  }
  foreach my $m ( @mtx ) {
    close $m->{FPT} or die "Couldn't close file, '$m->{trans}'. $!";
    `cat $m->{trans} | transpose > $m->{file}`;
    unlink $m->{trans} or die "Couldn't delete file matrix transpose, '$m->{trans}'. $!";
  }

  foreach my $m ( @mtx ) { 
    if ( $DO_MATRIX_PLOT ) { 
      #next if $m->{set} eq 'all';
      plot_comparison_matrix($m->{file}, $crystal, $compare, $m->{set}, $m->{chi});
    }
  }
}

##############################################################################
#                                                       PRINT COMPARISON PYMOL
#-----------------------------------------------------------------------------
# Going to split these out for easier comprehension
##############################################################################

sub print_comparison_pymol { 

  my $pdbs      = $_[0];
  my $seq_list  = $_[1];
  my $seq       = $_[2];
  my $seq_min   = $_[3];
  my $seq_max   = $_[4];
  my $data      = $_[5];
  my $compare   = $_[6];
  my $crystal   = $_[7];

  my $baseline = 1.0;
  
  if ( not ( -e "colors" and -d "colors" ) ) { 
    mkdir "colors" or die "Couldn't mkdir 'colors'. $!";
  }
  #============================================================ FIRST PDB LOOP
  foreach my $i ( 0..$#$pdbs ) { 
    my $name_i = $pdbs->[$i]{'name'};
    $name_i =~ s/\s/_/g;
    #========================================================== FIRST PDB LOOP
    foreach my $j ( ($i+1)..$#$pdbs ) { 
      my $name_j = $pdbs->[$j]{'name'};
      $name_j =~ s/\s/_/g;
       
      #========== Prep PyMOL Colors 
      my $file_color_ca  = "color.$name_i.$name_j.$compare.ca.pymol";
      my $file_color_chi = "color.$name_i.$name_j.$compare.chi.pymol";
      open COL_CA,  ">colors/$file_color_ca"  or die "Couldn't open '$file_color_ca'. $!";
      open COL_CHI, ">colors/$file_color_chi" or die "Couldn't open '$file_color_chi'. $!";
      print COL_CA  "alter all, b=$baseline\n";
      print COL_CHI "alter all, b=$baseline\n";

      #========== Loop Over Residues
      foreach my $k ( @$seq_list ) {
      #foreach my $k ( $seq_min..$seq_max ) { 
        my $resn = $seq->[$k];
        my $angles = scalar @{$AA{lc $resn}{'sc'}};

        #======== Loop over Chi Angles
        foreach my $chi ( 1..$angles ) { 
          my $c = $data->[$i][$j][$k][$chi];
          next if not defined $c;


          if ( $chi == 1 ) { 
            print COL_CA  "alter resi $k, b=$c\n";
            print COL_CHI "alter resi $k, b=$c\n";
          } elsif ( $chi==2 ) { 
            print COL_CHI "alter resi $k and n. *g+*d+*e+*z, b=$c\n";
          } elsif ( $chi==3 ) { 
            print COL_CHI "alter resi $k and n. *d+*e+*z, b=$c\n";
          } elsif ( $chi==4 ) { 
            print COL_CHI "alter resi $k and n. *e+*z, b=$c\n";
          } 
        }
      }

      print COL_CA  "cmd.spectrum('b', 'rainbow_rev', minimum=0, maximum=1)\n";
      print COL_CHI "cmd.spectrum('b', 'rainbow_rev', minimum=0, maximum=1)\n";
      close COL_CA  or die "Couldn't close file, 'colors/$file_color_ca'. $!";
      close COL_CHI or die "Couldn't close file, 'colors/$file_color_chi'. $!";
    }
  }
     
}

##############################################################################
#                                                                  PLOT MATRIX
##############################################################################

sub plot_comparison_matrix { 
  my $file    = $_[0];
  my $crystal = $_[1];
  my $compare = $_[2];
  my $pairs   = $_[3];
  my $chis    = $_[4];

  #return;

  my $zmin = 0.5;
  my $zmax = 1.0; 
  my $reverse = " ";
  my $zticksize = 0.05;
  my $colorscheme = 'rainbow';
  if ( $compare eq 'rmsd' ) {
    $zmin = 0.3;
    $zmax = 1.0;
    $reverse = " ";
    #$reverse = "reversecolors ";
    $zticksize = 0.1;
    $colorscheme = 'white-red';
  }


  #Ringer Correlation Matrix - IPA Crystal C - Ref Pairs - Chi 1
  my $title = sprintf "%s - %s - %s Pairs - %s",
                      $crystal,
                      ucfirst $compare,
                      ucfirst $pairs,
                      ucfirst $chis,
                      ;

  my $cmd = "plot_matrix_categories "
          . "$file "
          . "title='$title' "
          . "height=$PLOT_HEIGHT_COMP "
          . "width=$PLOT_WIDTH_COMP "
          . (@FILTER_RES ? "catx " : "axisx ")
          . "caty "
          . "colorscheme=$colorscheme "
          . "zmin=$zmin "
          . "zmax=$zmax "
          . "zcolors=1 "
          . "zticksize=$zticksize "
          . "$reverse "
          . "xlabel=Residues "
          . "ylabel='' "
          . (sprintf "xtickinterval=%.0f ", ( scalar @FILTER_RES ? (scalar @FILTER_RES / 5.0) : 5.0 ))
          ;

  print "  $cmd\n";
     
  my $flags = { verbose => $VERBOSE, fatal => 1 }; 
  my %result = Utility::Try_Command($cmd, "Plot $title", $flags);
  if ( $result{'value'} ) { 
      Utility::Report_Error ( 
        "Command failed",
        \%result,
        $flags,
      );
  } 
        

}


##############################################################################
#                                                   PLOT DIHEDRAL DISTRIBUTION
##############################################################################


my (@COLORS_BLUE_RED, @COLORS_1, @COLORS_2, @COLORS_3, @COLORS_5, @COLORS_9, @COLORS_GRAY2, @COLORS_GRAYSCALE, @COLORS_MISC);
my (%PLOT_CFG_PRES, %PLOT_CFG_ANAL, %PLOT_CFG_PAPER);

my %DIHEDRAL_BOUNDARIES     = ();
my %ROTAMER_BOUNDARIES      = ();
my %ROTAMER_BOUNDARIES_TEXT = ();

sub initialize_plots { 

  #======================================================= DIHEDRAL BOUNDARIES
  %DIHEDRAL_BOUNDARIES = (
    TETRA   => { s => -180, d => 120, c => 180 },
    PLANAR  => { s => -225, d =>  90, c =>   0 },
    TRPX2   => { s => -240, d => 120, c =>   0 },
    PROLINE => { s => -180, d =>  30, c =>   0 },
  );

  %ROTAMER_BOUNDARIES_TEXT = (
    ALA => [ qw(TETRA)],
    ARG => [ qw(TETRA TETRA TETRA TETRA) ],
    ASN => [ qw(TETRA PLANAR) ],
    ASP => [ qw(TETRA PLANAR) ],
    CYS => [ qw(TETRA) ],
    GLN => [ qw(TETRA TETRA PLANAR) ],
    GLU => [ qw(TETRA TETRA PLANAR) ],
    GLY => [],
    HIS => [ qw(TETRA PLANAR) ],
    ILE => [ qw(TETRA TETRA) ],
    LEU => [ qw(TETRA TETRA) ],
    LYS => [ qw(TETRA TETRA TETRA TETRA) ],
    MET => [ qw(TETRA TETRA TETRA) ],
    PHE => [ qw(TETRA PLANAR) ],
    PRO => [ qw(PROLINE) ],
    SER => [ qw(TETRA) ],
    THR => [ qw(TETRA) ],
    TRP => [ qw(TETRA PLANAR) ],
    TYR => [ qw(TETRA PLANAR) ],
    VAL => [ qw(TETRA) ],
  );

  for my $res ( keys %ROTAMER_BOUNDARIES_TEXT ) { 
    foreach my $chi ( 0..$#{$ROTAMER_BOUNDARIES_TEXT{$res}} ) {
      my $name     = $ROTAMER_BOUNDARIES_TEXT{$res}[$chi];
      my $boundary = $DIHEDRAL_BOUNDARIES{$name};
      $ROTAMER_BOUNDARIES{$res}[$chi] = $boundary;
    }
  }

  #==================================================================== COLORS
  ##### Some pretty default color series for the plot lines.
  @COLORS_1 = ( 
    'rgb "black"',
  );

  @COLORS_2 = (
    #'rgb "black"',
    #'rgb "gray"',
    'rgb "dark-blue"',
    'rgb "dark-red"',
  );

  # For 5 plots
  @COLORS_5 = (
    'rgb "dark-violet"', 
    'rgb "blue"',
    'rgb "forest-green"',
    #'rgb "dark-green"',
    'rgb "gold"',
    #'rgb "orange"',
    'rgb "red"',
  );

  # For 9 plots
  @COLORS_9 = (
    'rgb "dark-violet"', 
    'rgb "dark-blue"',
    'rgb "blue"',
    'rgb "dark-cyan"',
    'rgb "forest-green"', 
    #'rgb "dark-green"', # not color blind friendly to dark-red
    'rgb "gold"',
    #'rgb "dark-goldenrod"',
    'rgb "orange"',
    'rgb "red"',
    'rgb "dark-red"',
  );

  # Miscellaneous color options
  @COLORS_MISC = (
    'rgb "black"',
    'rgb "red"',
    'rgb "orange"',
    'rgb "yellow"',
    'rgb "green"',
    'rgb "forest-green"',
    'rgb "chartreuse"', 
    'rgb "dark-green"',
    'rgb "dark-cyan"',
    'rgb "cyan"',
    'rgb "blue"',
    'rgb "dark-gray"',
    'rgb "brown"',
    'rgb "navy"',
  );


  @COLORS_BLUE_RED = (
    'rgb "blue"',
    'rgb "light-blue"',
    'rgb "magenta"',
    'rgb "light-red"',
    'rgb "red"',
  );


  @COLORS_GRAY2 = (
    'rgb "gray"',
    'rgb "black"',
  );

  @COLORS_GRAYSCALE = ( 
    'rgb "black"', 
    'rgb "grey20"', 
    'rgb "grey40"', 
    'rgb "grey60"', 
    'rgb "grey80"', 
  );


  #============================================================ GNUPLOT STYLES
  # Linetypes:
  #  -1: -------
  #   0: .......
  #   1: -------
  #   2: - - - - 
  #   3: .......
  #   4: ???
  #   5: -..-..-
  #
  #   But only LT -1 and 0 are supported. 



  %PLOT_CFG_PRES = (
    width_inches  =>    4,
    height_inches =>    3,
    width_pixels  => 1600,
    height_pixels => 1200,

    general_font  => "Helvetica,18",
    title_font    => "Helvetica,22",
    label_font    => "Helvetica,20",
    key_font      => "Helvetica,16",
    axis_font     => "Helvetica,20",
    anno_font     => "Helvetica,12",

    line_width    =>    6.0,

    lmargin       =>    5.0,
    rmargin       =>    1.0,
    tmargin       =>    2.0,
    bmargin       =>    8.0,

    border_sides  =>   15.0, # bits specify sides and are all on
    border_lw     =>    4.0,
    
    title_x_coord => "graph",
    title_y_coord => "screen",
    title_x       =>    0.50,   # Parameters for using set label
    title_y       =>    0.95,
    title_justify => "center",
    #title_x       =>    0.00,  # Parameters when I was using set title. 
    #title_y       =>   -0.50,  # Now using set label so I can left justify it

    #------------- X/Y Tick Marks & Gridlines
    xtics_lt      =>    0,
    xtics_lw      =>    3.0,
    xtics_lc      =>  'rgb "gray"',
    mxtics_lt     =>    0,
    mxtics_lw     =>    1.0,
    mxtics_lc     =>  'rgb "light-gray"',
    xtics_mx      =>    2,

    ytics_lt      =>    0,
    ytics_lw      =>    3.0,
    ytics_lc      =>  'rgb "gray"',
    ytics_mx      =>    2,

    #------------- Legend Properties
    key_x         =>    0.00, # x coordinate
    key_y         =>    0.00, # y coordinate
    key_spacing   =>    1.00, # vertical spacing
    key_samplen   =>    0.20, # length of color bar
    key_width     =>   -4.00, # proportions of something in there
    key_position  => "bottom left",
    key_justify   => "Left",
    key_orient    => "vertical",
    key_reverse   => "reverse",
    #key_gravity   => "bottom left Left reverse vertical",
    #key_gravity   => "center left Left reverse vertical",

    #------------- Annotations Properties
    anno_width    => -0.25,
    anno_height   =>  0.03,
    anno_x        =>  0.99,
    anno_y        =>  0.02,
    anno_columns  =>  7.00,
    anno_position => "right front",
    
    
    #------------- Other Properties
    cutoff        =>    0.3,  # density cutoff marker
    
    #------------- Labels
    xlabel_x      =>    0.00,
    xlabel_y      =>    0.20,
    ylabel_x      =>    1.40,
    ylabel_y      =>    0.00,
  );

  #=============== Analysis - Smaller than Normal?
  %PLOT_CFG_ANAL = %PLOT_CFG_PRES;
  $PLOT_CFG_ANAL{width_pixels } = 800;
  $PLOT_CFG_ANAL{height_pixels} = 600;


  #=============== Papers - Minimal Decorations - Very condensed
  %PLOT_CFG_PAPER = (
    width_inches  =>    4,
    height_inches =>    3,
    width_pixels  => 1600,
    height_pixels => 1200,

    general_font  => "Helvetica,32", # Every conceivable font
    title_font    => "Helvetica,36",
    label_font    => "Helvetica,28",
    key_font      => "Helvetica,24",
    axis_font     => "Helvetica,32",
    anno_font     => "Helvetica,18",

    line_width    =>   12.0,
    border_sides  =>   15.0, # Bits specify sides and are all on
    border_lw     =>    6.0,
    
    lmargin       =>    5.20, # Margins
    rmargin       =>    1.50,
    tmargin       =>    2.00,
    bmargin       =>    3.50,

    title_x_coord => "graph", # Title coordinate systems and positions
    title_y_coord => "screen",
    title_x       =>    0.000,  
    title_y       =>    0.950,
    title_justify => "left",

    key_x         =>    0.95, # x coordinate
    key_y         =>    0.94, # y coordinate
    key_spacing   =>    3.00, # vertical spacing
    key_samplen   =>    0.20, 
    key_width     =>    0.00,
    key_position  => "center right",  # Strangely not ignored
    key_justify   => "Right",         
    key_orient    => "horizontal",
    key_reverse   => "reverse",

    xtics_mx      =>    2, # Minor tick count between lables
    xtics_lt      =>    0,    
    xtics_lw      =>    5.0,
    xtics_lc      =>  'rgb "gray"',
    mxtics_lt     =>    0,
    mxtics_lw     =>    1.0,
    mxtics_lc     =>  'rgb "light-gray"',

    ytics_mx      =>    2, # Minor tick count between lables
    ytics_lt      =>    0,
    ytics_lw      =>    1.0,
    ytics_lc      =>  'rgb "gray"',

    cutoff        =>    0.3,

    xlabel_x      =>    0.00,
    xlabel_y      =>    0.20,
    ylabel_x      =>    1.00,
    ylabel_y      =>    0.00,
  );

 
}


################################################## PLOT DIHEDRAL DISTRIBUTIONS

sub plot_dihedral_distributions { 
  my $title     = $_[0];
  my $file      = $_[1];
  my $xaxis     = $_[2];
  my $arrays    = $_[3];
  my $labels    = $_[4];
  my $annos     = $_[5];
  my $resn      = $_[6];
  my $chi       = $_[7];

  #=========================================================== Initialize Plot

  #--------------- Determine Mode
  my %PC = ();
  if    ( $STYLE eq 'presentation' ) { %PC = %PLOT_CFG_PRES }
  elsif ( $STYLE eq 'analysis'     ) { %PC = %PLOT_CFG_ANAL }
  elsif ( $STYLE eq 'paper'        ) { %PC = %PLOT_CFG_PAPER }
  else { die "Unknown mode '$STYLE'.\n" }

  #my $image_dir = "images_dihedral";
  #if ( not -e $image_dir ) { mkdir $image_dir or die "Couldn't mkdir $image_dir. $!" }

  #--------------- Set Filenames and Clean Title and Labels
  $title =~ s/_/ /g;
  $title =~ s/\.pdb//g;
  $file  =~ s/_/ /g;
  $file  =~ s/\.pdb//g;
  
  my $DATA = "tmp_plot_dh.dat";
  my $GPLT = "tmp_plot_dh.gplt";
  #my $EPS   = "tmp_plot_dh.eps";
  my $PNG  = Utility::Clean_Filename($file) . ".png";
  my $EPS  = Utility::Clean_Filename($file) . ".eps";

  if ( not $SHOW_TITLE ) { 
    $title = "";
  }

  my $XLAB = "";
  my $YLAB = "";
  if ( $SHOW_LABELS ) {
    $XLAB = "Dihedral Angle ({\\312})";
    $YLAB = "Electron Density ({/Symbol s})";
  }

  #--------------- Determine Color Scheme
  my @COLORS = @COLORS_1;
  my $count = scalar @$arrays;
  if    ( $count == 1 ) { @COLORS = @COLORS_1 }
  elsif ( $count == 2 ) { @COLORS = @COLORS_2 } 
  elsif ( $count  > 5 ) { @COLORS = @COLORS_9 }
  elsif ( $count  > 2 ) { @COLORS = @COLORS_5 }

  #if    ( scalar @$arrays  > 5 ) { @COLORS = @COLORS_9 }
  #elsif ( scalar @$arrays  > 2 ) { @COLORS = @COLORS_5 }
  #elsif ( scalar @$arrays  > 1 ) { @COLORS = @COLORS_2 }

  #@COLORS = @COLORS_BLUE_RED;
  #@COLORS = @COLORS_GRAY2;
  #@COLORS = @COLORS_GRAYSCALE;


  #========================================================== Create Data File

  #--------------- Find 180 degrees in the axis
  my $midpoint = -1;
  for my $i ( 0..($#$xaxis-1) ) { 
    if ( $xaxis->[$i] < 180 and $xaxis->[$i+1] >= 180 ) { 
      $midpoint = $i+1;
      last;
    }
  }
  if ( $midpoint < 0 ) { die "WTF Midpoint wasn't set.\n" . Dumper($xaxis) }


  #--------------- Decide if 0 centered or 180 centered
  my @angle_sequence = ();
  my $dh_boundary = $ROTAMER_BOUNDARIES{$resn}[$chi-1] || $DIHEDRAL_BOUNDARIES{TETRA};
  my $rotate_degrees =   0;
  my $deg_start      =   0;
  my $deg_end        = 360;
  
  if ( $dh_boundary->{c} == 180 ) { 
    @angle_sequence = ( 0..$#$xaxis );
  } else { 
    @angle_sequence = ( $midpoint..$#$xaxis, 0..($midpoint-1) );
    $rotate_degrees = 1;
    $deg_start = -180;
    $deg_end   =  180;
  }
  push @angle_sequence, scalar @$xaxis;

  #--------------- Write the temporary data file
  open DAT, ">$DATA", or die "Couldn't open data file, '$DATA'. $!";
  #print "Length xaxis: " . (scalar @$xaxis) . "\n"
  #    . "Ang Seq[0]:   $angle_sequence[0]\n"
  #    . "Ang Seq[-1]:  $angle_sequence[-1]\n"
  #    ;
  
  foreach my $i ( @angle_sequence ) { 
    my ($index, $angle); 
    $index = $i;
    $angle = $xaxis->[$i];
    if ( $i > $#$xaxis ) { 
      $index = $angle_sequence[0];
      $angle = $deg_end;
    } else { 
      if ( $rotate_degrees and $angle >= 180 ) { $angle -= 360 }
    }

    printf DAT $angle;
    foreach my $array ( 0..$#$arrays ) { 
      printf DAT "\t%.6f", ($arrays->[$array][$index] || 0);
    }
    print DAT "\n";
  }
  close DAT or die "Couldn't lose data file, '$DATA'. $!";



  #=========================================================== Write Plot File
  # Enhanced Postscript Control Characters
  #   http://gnuplot.sourceforge.net/docs_4.2/node411.html
  #   http://www.math.u-bordeaux1.fr/~mleguebe/docs/gnuplot_liite34.pdf
  # ------------------------------------------------------
  # a^x         : superscripted x
  # a_x         : subscripted x
  # x@^x_y      : aligns superscript x over subscript y
  # &{space}    : inserts space of specified length
  # ~a{.8-}     : overprints '-' on 'a', vertically offset by .8 times the current fontsize
  # {/Symbol a} : Greek characters
  # {/243}      : Octal characters
  # The file "ps_guide.ps" in the /docs/psdoc subdirectory of the gnuplot 
  # source distribution contains more examples of the enhanced syntax
  #===========================================================================

#set key at screen $PC{key_x},$PC{key_y} center right Right reverse horizontal samplen $PC{key_samplen} spacing $PC{key_spacing} width $PC{key_width} font "$PC{key_font}"
#set key at screen 0.05,0.1 center left Left reverse horizontal samplen 0.1 spacing 1.0 width -4 font "Helvetica,14"
#set title  "$title"  offset -15.0,-0.25    font "Helvetica,24"
#set title  "$title"  offset 0.0,-0.25    font "Helvetica,22"


my $PLOT = <<__GNUPLOT__
#----------------- Output Format
set terminal postscript eps color solid "$PC{general_font}" enh round size $PC{width_inches} in, $PC{height_inches} in
#set terminal postscript eps color solid "Helvetica" 20 enh round size $PC{width_inches} in, $PC{height_inches} in
set output "$EPS"
set clip one
set clip two
set clip points

#----------------- Plot Margins
set lmargin $PC{lmargin}
set rmargin $PC{rmargin}
set tmargin $PC{tmargin}
set bmargin $PC{bmargin}

#----------------- Title and Labels
#set title  "$title"  \\
#           offset $PC{title_x},$PC{title_y} \\
#           font "$PC{title_font}"

set label "$title" \\
          at    $PC{title_x_coord} $PC{title_x}, \\
                $PC{title_y_coord} $PC{title_y}  \\
                $PC{title_justify} \\
          font "$PC{title_font}"
          

set xlabel "$XLAB"  offset $PC{xlabel_x},$PC{xlabel_y}  font "$PC{label_font}"
set ylabel "$YLAB"  offset $PC{ylabel_x},$PC{ylabel_y}  font "$PC{label_font}"

#----------------- Range and Ticks
set xrange [$deg_start:$deg_end]
set yrange []

set xtics $dh_boundary->{s},$dh_boundary->{d} font "$PC{axis_font}"
set mxtics $PC{xtics_mx}

set ytics 1.0 font "$PC{axis_font}"
set mytics $PC{ytics_mx}

set grid xtics mxtics ytics mytics \\
  lt $PC{xtics_lt} lw $PC{xtics_lw} lc $PC{xtics_lc}, \\
  lt $PC{mxtics_lt} lw $PC{mxtics_lw} lc $PC{mxtics_lc}

#----------------- Thick Border
set border $PC{border_sides} lw $PC{border_lw}

#----------------- Reference Lines
zero(x)   = 0.0
cutoff(x) = $PC{cutoff}

__GNUPLOT__
;

if ( $SHOW_LEGEND ) {
  $PLOT .= <<__GNUPLOT__
#----------------- Key
set key at screen   $PC{key_x},$PC{key_y}   \\
                    $PC{key_position}       \\
                    $PC{key_justify}        \\
                    $PC{key_orient}         \\
                    $PC{key_reverse}        \\
        samplen     $PC{key_samplen}        \\
        spacing     $PC{key_spacing}        \\
        width       $PC{key_width}          \\
        font       "$PC{key_font}"

__GNUPLOT__
;
} else { 
  $PLOT .= "set key off\n";
}


  ##### Add the annotations

  my @annos           = reverse(map { s/_/ /g; $_ } @$annos);
  my $anno_count      = scalar @$annos;
  #my $anno_bottom     =  0.02;
  #my $anno_left       =  0.97;
  #my $anno_height     =  0.03;
  #my $anno_width      = -0.25;
  #my $anno_column_max =  7.00;

  if ( $STYLE ne 'paper' ) { 
    foreach my $anno_i (0..$#annos) { 
      my $anno = $annos[$anno_i];
      my $column = POSIX::floor($anno_i/$PC{anno_columns});
      my $row    = $anno_i % $PC{anno_columns};
      $PLOT .= sprintf "set label '%s' at screen %.4f,%.4f %s font '%s'\n",
                        $anno,
                        ($column * $PC{anno_width }) + $PC{anno_x},
                        ($row    * $PC{anno_height}) + $PC{anno_y},
                        $PC{anno_position},
                        $PC{anno_font},
                        ; 
    }
  }
#  if ( $STYLE ne 'paper' ) { 
#    foreach my $anno_i (0..$#annos) { 
#      my $anno = $annos[$anno_i];
#      my $column = POSIX::floor($anno_i/$anno_column_max);
#      my $row    = $anno_i % $anno_column_max;
#      $PLOT .= sprintf "set label '%s' at screen %.4f,%.4f right front font '$PC{anno_font}'\n",
#                        $anno,
#                        ($column * $anno_width ) + $anno_left,
#                        ($row    * $anno_height) + $anno_bottom,
#                        ; 
#    }
#  }
    

  #=============== Add the complicated plot line
  my @plot_lines = ();

  if ( scalar @$arrays == 0 ) { die "Passed no arrays to plot\n" }
  if ( scalar @$labels == 0 ) { die "Passed no labels to plot\n" }

  foreach my $i ( 0..$#$arrays ) { 
    next if not defined $labels->[$i];
    my $j = $i+2;
    my $lab = $labels->[$i];
    my $color = "gray" . ($i+1) . "0";
    #$lab =~ s/IPA //;
    #$lab =~ s/K.*/K/;

    $lab =~ s/3PZW/WT/;
    $lab =~ s/DM.*1xASU/DM/;

    push @plot_lines, "'$DATA' u 1:$j title '$lab' with lines lt 1 lw $PC{line_width} lc $COLORS[$i]";
  }

  #=============== Add cutoff marker
  push @plot_lines, "cutoff(x) title '' with lines lw 1 lt rgb 'black'";
  push @plot_lines, "zero(x)   title '' with lines lw 2 lt rgb 'black'";

  $PLOT .= ("plot " . (join ",  \\\n     ", @plot_lines) . "\n");

  open GPLT, ">$GPLT" or die "Couldn't open gplt file, '$GPLT'. $!";
  print GPLT $PLOT;
  close GPLT or die "Couldn't close gplt file, '$GPLT'. $!";

  ##### Plot the File and Convert to PNG
  my $GNUPLOT = "/usr/local/bin/gnuplot $GPLT";
  #my $GNUPLOT = "/usr/bin/gnuplot $GPLT";
  print `$GNUPLOT`;

  my $CONVERT_CMD = "convert -alpha off "
                          . "-antialias "
                          . "-background white "
                          . "-density 600 "
                          . "$EPS "
                          . "-background white "
                          . "-resize $PC{width_pixels}x$PC{height_pixels} "
                          . "$PNG";
  print `$CONVERT_CMD`;

  if ( not -e $PNG ) { 
    die "Failed to plot file, '$PNG'.\n";
  }

}




##############################################################################
#                                                     PARSE CONFIGURATION FILE
##############################################################################


sub parse_conf_file { 

  my $conf = $_[0];
  if ( not -e $conf ) { die "Couldn't find configuration file, '$conf'.\n" }

  my @lines = split /\n/, `cat $conf`;
  my $title;
  my $struct_dir;
  my $out_dir;
  my @dirs;

  foreach my $line ( @lines ) { 
    $line = Utility::Clean_Line($line);
    next if not $line;
    if ( $line =~ /^TITLE\s+(.*)$/ ) { 
      $title = $1;
      next;
    }
    if ( $line =~ /^STRUCTURE_DIR\s+(.*)$/ ) { 
      $struct_dir = $1;
      if ( not -e $struct_dir ) { 
        die "Can't find struture directory, '$struct_dir'.\n";
      }
      next;
    }

    if ( $line =~ /^OUT_DIR\s+(.*)$/ ) { 
      $out_dir = $1;
      if ( not -e $out_dir ) { 
        mkdir $out_dir or die "Couldn't mkdir output dir, '$out_dir'. $!";
      }
      next;
    }

    if ( $line =~ /^PDB/ ) { 
      my ( $x, $name, $desc, $file, $chain, $offset ) = split /\s+/, $line;
      if ( not defined $file ) { next } 
      if ( not -e $file ) { $file = "$struct_dir/$file"; }
      if ( not -e $file ) { die "Couldn't find file, '$file'." }
      my ( $seq, $chi, $deg, $den ) = read_ringer_signal_plot($file, $chain);
      $name =~ s/_/ /g;
      $desc =~ s/_/ /g;
      push @dirs, {
        name    => $name,
        desc    => $desc,
        file    => $file,
        chain   => $chain,
        offset  => $offset,
        length  => (scalar @$seq),
        #length  => scalar (grep {defined} @$seq),
        seq     => $seq,
        chi     => $chi,
        deg     => $deg,
        den     => $den,
      };
      #print Dumper (\@dirs);
      #exit;
      next;
    }

    die "Unknown line, '$line'";
  }

  return $title, $out_dir, @dirs;
}


##############################################################################
#                                                 READ RINGER SIGNAL PLOT FILE
##############################################################################


sub read_ringer_signal_plot { 

  my $file  = $_[0];
  my $chain = $_[1];

  if ( not -e $file or -d $file ) { 
    if ( -e "$file/ringer.signal_plot.txt" ) { 
      $file .= "/ringer.signal_plot.txt";
    } else { 
      die "Can't find file $file\n";
    }
  }

  my @contents = map { [ split "\t", $_ ] } split "\n", `cat $file`;

  my $header = shift @contents;
  my @header = @$header;
  shift @header;

  my @matrix = Utility::Transpose(\@contents);

  my @SEQ = ();
  my @CHI = ();
  my @DEG = ();
  my @DEN = ();

  foreach my $i ( 0..$#header ) { 

    my $h = $header[$i];

    $h =~ /^(\w+)_(\w)_(\d+)_chi(\d)$/ 
      or die "Couldn't parse column header $h\n";

    my ( $resn, $ch, $resi, $chi ) = ( $1, $2, $3, $4 );
    next if ( $chain ne $ch );
    my @m = @{$matrix[$i]};
    if ( not @DEG ) { @DEG = map { /^(\d+) .*$/; $1 } @m }

    $SEQ[$resi] = $resn;
    $CHI[$resi] = scalar @{$Utility::AA{lc $resn}{'sc'}};
    $DEN[$resi][$chi] = [ map { /^\d+ (.*)$/; $1 } @m ];
    if ( $resn eq 'ALA' ) { $CHI[$resi] = 1 }
  }
  return ( \@SEQ, \@CHI, \@DEG, \@DEN );
}


