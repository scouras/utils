package RotLib::Rotamer;
BEGIN{
  push @INC,
            "/users/scouras/code",
            "/users/scouras/code/scouras",
            "/net/programs/perl",
#            "/net/programs/perl/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/",
            "/net/programs/perl/lib/perl5/site_perl/5.8.5/i386-linux-thread-multi",
}
#use 5.008006;
use strict;
use warnings;

use Data::Dumper;
use Utility qw(:all);
use PerlIO::gzip;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use RotLib ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);


# CVS VERSION INFORMATION (AUTOMAGICALLY GENERATED BY CVS)
'$Revision: 1.1.1.1 $'          =~ /^.Revision: (.*) \$/; our $VERSION    = $1;
'$Date: 2009/03/25 18:32:30 $'  =~ /^.Date: (.*) \$/;     our $CHECKED_IN = $1;
'$Author: scouras $'            =~ /^.Author: (.*) \$/;   our $AUTHOR     = $1;



my %CONFIG_ROT = (

  files => {
    # CONFIGURATION INPUT
    #dihedral_def_file => 'dihedral.dunbrack.def',
    dihedral_def_file         => 'dihedral.full.def',
    rotamer_def_static_file   => 'rotamers.static.def',
    rotamer_def_dynamic_file  => 'rotamers.dynamic.def',
    
    
    # DATA INPUT
    major_hist        => 'rotamer_major_hist.dat',
    major_time        => 'rotamer_major_t.dat',
    dihedral_examples => '*_t.dat',
    dihedral_time     => 'summary_t.dat',
    dihedral_time_old => 'residue_joined_t.dat.NOT',

    # PROGRAMS
    summary_program_file  => "code/scouras/join_chi_angles.pl",
  },

  write_cache         => 1,
  bin_width           => 20,
  overwrite_rotamers  => 0,
  overwrite_dihedrals => 0,
  metadata            => [qw(total_rotamers rotamer_names residue res_name sub_count)],

  readonly            => 0,

);


###################################### INITIALIZE ROTAMER CONVERSION

sub Initialize {

  my $class     = $_[0];
  my $CONFIG    = $_[1];
  my $AA        = $_[2];
  my $PROTEINS  = $_[3] || {};

  my $ROT = { %CONFIG_ROT };
  
  bless $ROT, $class;

  if    ( $CONFIG->{'rotamer_mode'} eq 'static' ) {
    $ROT->{'files'}{'rotamer_def_file'} = $ROT->{'files'}{'rotamer_def_static_file'}; 
    $ROT->{'files'}{'cache'} = 'summary_t.static.dat';
  } elsif ( $CONFIG->{'rotamer_mode'} eq 'dynamic' ) {
    $ROT->{'files'}{'rotamer_def_file'} = $ROT->{'files'}{'rotamer_def_dynamic_file'}; 
    $ROT->{'files'}{'cache'} = 'summary_t.dynamic.dat';
  } else { 
    die "Unknown rotamer mode '$CONFIG->{'rotamer_mode'}'."; 
  }
  
  $ROT->{'AA'} = $AA;
  $ROT->{'main'}                  = $CONFIG;
  $ROT->{'bin_count'}             = POSIX::floor(360 / $ROT->{'bin_width'});
  $ROT->{'overwrite_rotamers'}   += $CONFIG->{'overwrite_rotamers'};
  $ROT->{'overwrite_dihedrals'}  += $CONFIG->{'overwrite_dihedrals'};
  $ROT->{'overwrite'}             = $ROT->{'overwrite_rotamers'} + $ROT->{'overwrite_dihedrals'};
  $ROT->{'files'}{'rotamer_def'}  = $CONFIG->{'root_dir'} . "/" . $CONFIG_ROT{'files'}{'rotamer_def_file'};
  $ROT->{'files'}{'dihedral_def'} = $CONFIG->{'root_dir'} . "/" . $CONFIG_ROT{'files'}{'dihedral_def_file'};
  $ROT->{'files'}{'summary_program'} 
                                  = $CONFIG->{'home_dir'} . "/" . $CONFIG_ROT{'files'}{'summary_program_file'};
                                  
  $ROT->Import_Dihedral_Definitions ( );
  $ROT->Import_Rotamer_Definitions  ( );
  $ROT->Initialize_Major_Rotamers   ( );

  $CONFIG->{'ROTAMER_OBJECT'} = $ROT;

  return $ROT;

};



######################################## IMPORT DIHEDRAL DEFINITIONS
# Each line of a dihedral files list a residue name, a
# dihedral angle name, and the 4 atoms which define that dihedral 
# angle.  Angle dependence MUST be in increasing order, as there is
# no other reliable way to determine dependence without parsing the
# full bonding structure of a "residue".  Note that here, residue
# means anything molecule-fragment used as a residue in ILMM that 
# could have dihedral definitions made for it.  
#-------------------------------------------------------------------
#   ALA omega CA- C- N CA
#   ALA phi C- N CA C
#   ALA psi N CA C N+
#   #ALA chi1 N CA CB HB1
#   ARG omega CA- C- N CA
#   ARG phi C- N CA C
#   ARG psi N CA C N+
#   ARG chi1 N CA CB CG
#   ARG chi2 CA CB CG CD
#   ARG chi3 CB CG CD NE
#   ARG chi4 CG CD NE CZ
#-------------------------------------------------------------------

sub Import_Dihedral_Definitions {
  
  my $self  = $_[0];
  my $AA    = $self->{'AA'};
  my $file  = $self->{'files'}{'dihedral_def' };
  my $line;
  
  open FILE, $file or die "Couldn't open dihedral definition file, '$file'. $!";

  while ( defined ( $line = <FILE> ) ) {
    next if ( ( $line = Clean_Line($line) ) eq '' );

    my ( $res_name, $angle, @atoms ) = split /\s+/, $line;
    
    # check for 4 atoms in definition.
    if ( scalar @atoms != 4 ) {
      die "Dihedral definition must have exactly 4 atoms.  $file : $res_name has " 
        . scalar @atoms . ": " . ( join ', ', @atoms ) . "\n";
    }

    my $res = $AA->{$res_name};
    if ( not defined $res->{'angle_names'} ) {
      $res->{'angle_names'} = [];
      $res->{'rotamer_names'} = [];
      $res->{'rotamer_indecies'} = [];
    }
    push @{$res->{'angle_names'}}, $angle;
    $res->{'angle_indecies'}{$angle} = $#{$res->{'angle_names'}};
    $res->{'dihedrals'}{$angle} = 
      { 
        'index' => scalar @{$res->{'angle_names'}} - 1,
        'atoms' => \@atoms,
      };
    $res->{'initialized_dihedrals'} = 1;
  }

  close FILE or die "Couldn't close dihedral definition file, '$file'. $!";


}







######################################### IMPORT ROTAMER DEFINITIONS
# 
#-------------------------------------------------------------------
#      1	g+	30	90
#      2	t	-30	30
#      3	g-	-90	-30
#      1	g+	-150	-90
#      2	t	150	-150
#      3	g-	90	150
#      ASN 2
#      ASP 2
#      ASH 2
#      GLN 3
#      GLU 3
#      GLH 3
#-------------------------------------------------------------------

sub Import_Rotamer_Definitions {
  my $self  = $_[0];
  
  my $AA    = $self->{'AA'};
  my $file  = $self->{'files'}{'rotamer_def'};
  my $line;
  
  open FILE, $file or die "Couldn't open rotamer definition file,' $file'. $!";

  my @rotamer_bins = ();
  my @primes = ();
  my @confs = (); # for error checking that bin names are unique in a minor rotamer
  my $doing = '';
  my $minor_count = 0;
  my $symmetry = undef;
  my @starts = ();
  
  while ( defined ( $line = <FILE> ) ) {
    next if ( ( $line = Clean_Line($line) ) eq '' );

    #=================================================== TRANSITIONS
    
    #----- Transition to inputting minors
    if ( $line =~ /^\d/ ) {
      if ( $doing eq 'angle' ) {
        # reinitialize for another round of minor rotamer defs
        @rotamer_bins = ();
        @primes = ();
        @confs = ();
        @starts = ();
      }
      $doing = 'minor';
      
      
    #----- Transition to inputting angles
    } else { 
      if ( $doing eq 'minor' ) {
       
        @rotamer_bins = sort { $a->{'start'} <=> $b->{'start'} } @rotamer_bins;
        
        # count minor bins
        my @numbers = map { $_->{'number'} } @rotamer_bins;
        $minor_count = scalar @{Unique ( @numbers )};
        if ( $minor_count != Max ( @numbers ) ) {
          die "The minor rotamer definition does not have the bins indexed correctly on line:\n"
            . "Numbers: " . (join ', ', @numbers) . "\n"
            . "Minor Counts: (" . (Dumper($minor_count)) . "\n"
            . "$line\n";
        }

        # Some rotamers have 2 fold rotational symmetry because
        # you cannot tell some atoms apart in crystal structures
        # or because they are literally symmetric (valine, 
        # tyrosine).  If this is such a case, check that all the 
        # definitions are symmetric and note the symmetry.
        # We will in fact handle any fold rotational symmetry.
        $symmetry = undef;
        foreach my $i ( 1..$minor_count ) {
          my @minors = grep { $i == $_->{'number'} } @rotamer_bins;
          my $sub_count = scalar @minors;
          my $prime_bin = ( sort { ($a->{'start'}%360) <=> ($b->{'start'}%360) } @minors )[0];
          my $min_start = $prime_bin->{'start'};
          my $max_start = $prime_bin->{'stop'};
          $prime_bin->{'is_prime'} = 1;
          $primes[$prime_bin->{'number'}] = $prime_bin;

          if ( not defined $symmetry ) { $symmetry = $sub_count; }
          elsif ( $symmetry != $sub_count ) { 
            die "There is a problem in the rotamer definitions file on line $..  There is a rotamer "
              . "definition where some bins have $symmetry fold symmetry and other bins "
              . "have $sub_count fold symmetry.  They must all have the same symmetry to "
              . "make any sense. \n";
          }
          # add symmetry to reach minor.
          foreach ( @minors ) {
            $_->{'symmetry' } = $symmetry;
            $_->{'min_start'} = $min_start;
            $_->{'max_start'} = $max_start;
            $_->{'prime_bin'} = $prime_bin;
          }
        }

     
        # verify that the minor rotamer def is contiguous and non-overlapping
        my $first_start = -1;
        my $last_stop = -1;
        my $connected = 0;
        foreach my $rot ( @rotamer_bins ) {
          my ($number, $conformation, $start, $stop) = 
            ($rot->{'number'}, $rot->{'conformation'},
             $rot->{'start'},  $rot->{'stop'});
    
          # different start and stop
          die "Start ($start) and Stop ($stop) "
            . "angles must be different in rotamer definitions $!"
            if ($start == $stop);
       
          # Some checks that we completed 360 degrees
          if ( $first_start < 0 ) { $first_start = $start };
          if ( $connected ) { die "Minor rotamer wraps around more than a full circle. $!"; }
          if ( $stop == $first_start ) { $connected = 1; }
        
          # Check that there are no breaks between stops and starts
          if ( $last_stop >= 0  and $last_stop != $start) {
            die "Minor rotamer has a break between $last_stop and $start. $!";
          }

          $last_stop = $stop;
        }
        
        # not 360 degrees
        if ( not $connected ) {
          die "Rotamer has a break between $last_stop and $first_start. $!";
        }

        # cache starts
        @starts = map { $_->{'start'} } @rotamer_bins;
      }
      $doing = 'angle';
    }


    #==================================================== INPUT DATA
    
    # New Minor Rotamer Definition
    if ( $doing eq 'minor' ) { 
      my ( $number, $conf, $start, $stop ) = split /\s+/, $line;
      push @rotamer_bins, 
        { 'number'        => $number,
          'conformation'  => $conf,
          'start'         => $start % 360,
          'stop'          => $stop  % 360,
          'crosses_0'     => (( $start % 360 ) > ( $stop % 360 ) ),
          'crosses_180'   => (( $start % 360 < 180) and ($stop % 360 > 180)),
        };
      if ( not exists $confs[$number] ) {
        $confs[$number] = $conf;
      } elsif ( $confs[$number] ne $conf ) {
        die "Rotamer is not well defined: "
          . "multiple conformations exist "
          . "($conf, $confs[$number]) "
          . "for a single index $number. $!";
      }
    }

    # List of Angles Prior Minor Rotamer Applies To.
    elsif ( $doing eq 'angle' ) {
      my ( $res, @angles ) = split /\s+/, $line;
      if ( not @angles ) { die "Rotamer definition line assigns no angles for residue '$res': \n'$line'.\n $!"; }
      foreach my $angle ( @angles ) {
        my $fullname;
        if ( $angle =~ /^\d+$/) { $fullname = "chi$angle"; }
        else                    { $fullname = $angle;      }
        push @{$AA->{$res}{'rotamer_names'}}, $fullname;
        $AA->{$res}{'dihedrals'}{$fullname}{'count'     } = $minor_count;
        $AA->{$res}{'dihedrals'}{$fullname}{'symmetry'  } = $symmetry;
        $AA->{$res}{'dihedrals'}{$fullname}{'minors'    } = [@rotamer_bins];
        $AA->{$res}{'dihedrals'}{$fullname}{'primes'    } = [@primes];
        $AA->{$res}{'dihedrals'}{$fullname}{'starts'    } = [@starts];
        $AA->{$res}{'initialized_rotamers'} = 1;
      }
    }
  }
  close FILE or die "Couldn't close rotamer definition file, '$file'. $!";
}







#========================================= INITIALIZE MAJOR ROTAMERS

sub Initialize_Major_Rotamers {
  
  my $self  = $_[0];

  my $AA = $self->{'AA'};
  
  while ( my ( $res_name, $res ) = each %$AA ) {
    
    if ( not $res->{'initialized_dihedrals'} 
      or not $res->{'initialized_rotamers'} ) {
        $res->{'total_rotamers'} = 1;
        $res->{'angle_names'} = [];
        $res->{'rotamer_names'} = [];
        $res->{'rotamer_indecies'} = [];
        next;
    }
    next if $res->{'finalized_rotamers'};

    my $total_rotamers = 1;
    my $dihedrals       = $res->{'dihedrals'};
    my @all_angles      = @{$res->{'angle_names'}};
    # SORT INTO NORMAL ORDER BY @ALL_ANGLE_NAMES
    my @rotamer_names  = @{$res->{'rotamer_names'}};
    my @rotamer_indecies = map { $res->{'angle_indecies'}{$_} } @rotamer_names;
    $res->{'rotamer_indecies'} = \@rotamer_indecies;

    $res->{'quickstarts' } = [];
    $res->{'quicknumbers'} = [];
    foreach my $i ( reverse ( 0..$#rotamer_names ) ) {
      
      # initialize major rotamer calculations
      my $angle = $rotamer_names[$i];
      $dihedrals->{$angle}{'multiplier'} = $total_rotamers;
      $total_rotamers *= $dihedrals->{$angle}{'count'};
    
      # initialize quickstarts
      $res->{'quickstarts' }[$i] = [map { $_->{'start' } } @{$dihedrals->{$angle}{'minors'}}];
      $res->{'quicknumbers'}[$i] = [map { $_->{'number'} } @{$dihedrals->{$angle}{'minors'}}];
    }
    $res->{'total_rotamers'} = $total_rotamers;
    $res->{'quickmults'  } = [map { $dihedrals->{$_}{'multiplier'} } @rotamer_names ];

#    $res->{'chi_angles'}          = [grep { /^chi/ } @all_angles];
#    $res->{'chi_angle_indecies'}  = [grep { $all_angles[$_] =~ /^chi/   } 0..$#all_angles];
    $res->{'phi_index'}           = (grep { $all_angles[$_] =~ /^phi/   } 0..$#all_angles)[0];
    $res->{'psi_index'}           = (grep { $all_angles[$_] =~ /^psi/   } 0..$#all_angles)[0];
    $res->{'omega_index'}         = (grep { $all_angles[$_] =~ /^omega/ } 0..$#all_angles)[0];


    $res->{'major_to_minor'} = [];
    $res->{'minor_to_major'} = {};
    for my $major ( 0..$total_rotamers-1) {
      my @minors = grep { defined $_ and $_ and $_ > 0 } $self->Get_Rotamer_Minors ( $res_name, $major );
      $res->{'major_to_minor'}[$major] = \@minors;
      $res->{'minor_to_major'}{join '_', @minors} = $major;
    }
    
    my %conformations = ();
    foreach my $angle ( @rotamer_names ) {
      my $minors = $res->{'dihedrals'}{$angle}{'minors'};
      foreach my $minor ( 0..$#$minors ) {
        $conformations{$angle}[$minor] = $minors->[$minor]{'conformation'};
      }
    }
    $res->{'conformations'} = \%conformations;
    $res->{'finalized_rotamers'} = 1;
  }
  
  # Some aggregate data
  $self->{'highest_rotamer_count'} 
    = Max ( map { $_->{'total_rotamers'} } map { $AA->{$_} } @{$self->{'main'}{'residue_list'}});
  $self->{'highest_angle_count'} 
    = Max ( map { scalar @{$_->{'angle_names'}}} map { $AA->{$_} } @{$self->{'main'}{'residue_list'}} );
  $self->{'highest_rotamer_angle_count'} 
    = Max ( map { scalar grep { /^chi/ } @{$_->{'rotamer_names'}}} map { $AA->{$_} } @{$self->{'main'}{'residue_list'}} );

}


sub Lookup_Rotamer_Major_to_Minor {
  my $self     = $_[0];
  my $res      = $_[1];
  my $major    = $_[2];
  return @{$self->{'AA'}{$res}{'major_to_minor'}[$major]};
}

sub Lookup_Rotamer_Minor_to_Major {
  my $self     = $_[0];
  my $res      = $_[1];
  my $minors   = $_[2];
  return $self->{'AA'}{$res}{'minor_to_major'}{join '_', @$minors};
}


#==================================== GET ROTAMER NUMBER FROM ANGLES
# for each chi angle, compare the chi angle to that chi's rotamers, 
# then multiply that minor index by it's multiplier and add to the  
# major index sum

sub Get_Rotamer_Number {
  my $self      = $_[0];
  my $res       = $_[1];
  my $angles    = $_[2];
  
  my $AA        = $self->{'AA'};
  my $major     = 0;
  my @minors;

  my %res = %{$AA->{$res}};
  my ($angle_name, $angle, %dihedral, $dihedral, $starts, $minor, $i);
  foreach my $index ( 0..$#$angles ) {
    $angle_name = $res{'angle_names'}[$index];
    next if not $angle_name =~ /^chi\d+$/;
    
    $angle = $angles->[$index] % 360;
    %dihedral = %{$res{'dihedrals'}{$angle_name}};
    $starts = $dihedral{'starts'};
        
    $minor = $#$starts; # default last, though usually replaced
    
    # find the minor bin after the one containing this angle
    foreach $i ( 0..$#$starts ) {
      if ( $angle <= $starts->[$i] ) { 
        $minor = $i - 1;
        last;
      }
    }
    $major += ( $dihedral{'minors'}[$minor]{'number'} - 1 ) * 
                $dihedral{'multiplier'};
    
    $minors[$index] = $minor;
  }
  return $major, [grep {defined $_} @minors];
}


#============================= GET ROTAMER MINORS FROM ROTAMER MAJOR

sub Get_Rotamer_Minors {
  
  my $self     = $_[0];
  my $res_name = $_[1];
  my $major    = $_[2];

  my $AA = $self->{'AA'};
  my @minors = ();
  my %res = %{$AA->{$res_name}};
  foreach my $index (0..$#{$res{'rotamer_names'}} ) {
    my $angle_name = $res{'rotamer_names'}[$index];
    my $mult = $res{'dihedrals'}{$angle_name}{'multiplier'};
    my $div  = POSIX::floor( $major / $mult );
    $minors[$index] = $div + 1;
    $major -= $div * $mult;
  }
  if ( wantarray ) { return @minors; }
  else             { return \@minors; }
}

sub Get_Rotamer_Minor_Names {

  my $self = $_[0];
  my $res_name = $_[1];
  my $minors = $_[2];

  my $AA = $self->{'AA'};
  my %res = %{$AA->{$res_name}};
  my @names = ();

  foreach my $index (0..$#{$res{'rotamer_names'}} ) {
    my $angle_name = $res{'rotamer_names'}[$index];
    my $m = (grep { $_->{'number'} == ($minors->[$index]+1) } 
                  @{$res{'dihedrals'}{$angle_name}{'minors'}}) [0];
    push @names, $m->{'conformation'};
  }
  if ( wantarray ) { return @names; }
  else             { return \@names; }
}


sub Get_Rotamer_Minor_Name_Options {
  my $self = $_[0];
  my $res_name = $_[1];

  my $AA = $self->{'AA'};
  my %res = %{$AA->{$res_name}};
  my @angle_names = @{$res{'rotamer_names'}};
  my @names = ();
  foreach my $angle ( 0..$#angle_names ) { 
    my @m = @{$res{'dihedrals'}{$angle_names[$angle]}{'minors'}};
    my @indecies = sort { $a <=> $b } Utility::Unique ( map { $_->{'number'} } @m );
    foreach my $i ( @indecies ) { 
      my $m = (grep { $_->{'number'} == $i } @m )[0];
      $names[$angle][$i-1] = $m->{'conformation'};
    }
  }
  if ( wantarray ) { return @names }
  else             { return \@names }
}



####################################################################
#                       DISTRIBUTION INITIALIZATION AND MANIPULATION
####################################################################





#========================================= EXTRACT ROTAMER INTERVALS

sub Extract_Rotamer_Intervals {
  my $self      = $_[0];
  my $RES       = $_[1];
  my $intervals = $_[2];

  my %NEW_RES = 
    map { $_ => $RES->{$_} } 
      qw(res_name headers rotamer_names total_rotamers residue files angle_names populated for);

  my $data      = $RES->{'data'};
  my @new_data  = ();
  my $i_time    = $RES->{'headers'}{'time'};

  for my $i ( 0..$#$data ) {

    my $time = $data->[$i][$i_time];
    next if $time != int($time); # TODO: Skip high resolution data for now.
    last if   Past_Range ( $intervals, $time );
    next if not In_Range ( $intervals, $time );
    push @new_data, [@{$data->[$i]}];
  }
  
  $NEW_RES{'timepoints' } = scalar @new_data;
  $NEW_RES{'data'       } = \@new_data;
  
  $self->Generate_Histograms ( \%NEW_RES );

  return \%NEW_RES;
}

####################################################################
#                                                        IMPORT DATA
####################################################################

################################################## READ ROTAMER DATA
# Determines if a cache file exists.  If not, reads the raw dihedral
# data and processes it, then saves the cache file.

sub Read_Rotamer_Data {

  my $self        = $_[0];
  my $intervals   = $_[1];
  my $res_name    = $_[2];
  my $res_dir     = $_[3];
  
  my $RES;
  
  # Try Cache First
  if ( not $self->{'overwrite'} ) {
    $RES = $self->Read_Rotamer_File ( $intervals, $res_name, $res_dir );
    if ( $RES and $self->{'main'}{'just_make_cache_rotamer_files'} ) { return 1 }
  }
  

#  # try to piece together a summary file from other files
#  if ( not $self->{'overwrite'} and not $RES ) {
#    $RES = $self->Aggregate_and_Read_Rotamer_File ( $intervals, $res_name, $res_dir );
#  }
    
  # if overwrite or if the data doesn't exist, then calculate the 
  # rotamers from the raw dihedral data
  if ( not $RES ) {

    my $RES_FULL = $self->Read_Dihedral_Files ( $res_name, $res_dir );
    if ( not $RES_FULL ) {
      $RES_FULL = $self->Read_Dihedral_Files ( $res_name, $res_dir, 1 );
    }
    if ( $self->{'write_cache'} ) {
      $self->Write_Rotamer_Cache_Files ( $RES_FULL, $res_dir );
    }
    $RES = $self->Extract_Rotamer_Intervals ( $RES_FULL, $intervals );
  }
  
  return $RES if $RES->{'populated'};

  return;

  die "Could not generate rotamer data for unknown reason.\n";
  
}



################################## AGGREGATE AND IMPORT ROTAMER FILE
sub Aggregate_and_Read_Rotamer_File {
  my $self        = $_[0];
  my $intervals   = $_[1];
  my $res_name    = $_[2];
  my $res_dir     = $_[3];

  $self->Run_Aggregator ( $res_dir ) or return 0;
  return $self->Read_Rotamer_File ( $intervals, $res_name, $res_dir );
}

sub Run_Aggregator { 
  my $self      = $_[0];
  my $res_dir   = $_[1];
  
  my $program = $self->{'files'}{'summary_program'};
  my $cmd = "$program $res_dir";
  # TODO: Test if successful
  `$cmd`;
  return 1;
}
 

################################################ GENERATE HISTOGRAMS

sub Generate_Histograms {
  my $self    = $_[0];
  my $RES     = $_[1];

  $self->Generate_Dihedral_Histograms ( $RES );
  $self->Generate_Rotamer_Histograms  ( $RES );

}

################################################ IMPORT ROTAMER FILE
# Reads a rotamer file, that is, a cache file that combines both
# dihedral and rotamer data at each time point.

sub Read_Rotamer_File {
  
  my $self        = $_[0];
  my $intervals   = $_[1];
  my $res_name    = $_[2];
  my $res_dir     = $_[3];

  my (@majors, @minors, @raw, @times);
  
 
  # Find the data file
  my $file = "$res_dir/$self->{'files'}{'cache'}";
  if ( not -e $file ) { 
    if ( -e "$file.gz" ) { 
      $file = "$file.gz";
    } else { 
      return;
    }
  };

  if ( $self->{'main'}{'just_make_cache_rotamer_files'} ) { return 1 }

  #print "Reading cache file $file\n";

  my $IN;
  open $IN, "<gzip(autopop)", $file or die "Couldn't open file, '$file'. $!";

  #my $IN = Open_File ( $file );
  #open IN, $file or die "Couldn't open rotamer data file, '$file'. $!";
  
  # for looking up which column is which
  my $headers = lc <$IN>;
  my @headers = split /\s+/, $headers;
  my %headers = ();
  my $i = 0;
  foreach my $h (@headers) {
    $h =~ s/#//g;
    $headers{$h} = $i++;
  }
  my $AA                = $self->{'AA'};
  my $res               = $AA->{$res_name};
  my @angle_names       = grep { $headers{$_} } @{$res->{'angle_names'}};
  my @rotamer_names     = @{$res->{'rotamer_names'}};
  my @angle_indecies    = map { $headers{$_} } @angle_names;
  my @rotamer_indecies  = map { $headers{$_} } @rotamer_names;

#  print "\nHeaders: " . join (', ', @headers) . "\n";
  return if ( not exists $headers{"rotamer_major"} );
  foreach my $angle ( @{$res->{'rotamer_names'}} ) {
#    print "Looking for angle $angle\n";
    return if not exists $headers{$angle};
#    print "Found $angle\n";
    return if not exists $headers{"rotamer_$angle"};
#    print "Found rotamer_$angle\n";
  }

  # read only the relevant timepoints from the file.
  my ($line);
  my @data;
  my $i_time = $headers{"time"};
  my $i_major = $headers{"rotamer_major"};
  my ($angle, $index);
  my $last_time;
  my $time_limit;
  if ( not $intervals->{'include_all'} ) { $time_limit = $intervals->{'last_time'} }

  #print Dumper ( $intervals );
  #exit;

  #print "WORKING WITH TIME LIMIT $time_limit\n";
  
  while ( defined ( $line = <$IN> ) ) {
    
    # skip comments and empty lines.
    next if ( ( $line = Clean_Line ( $line ) ) eq '' );
    
    my @line = split /\s+/, $line;
    my $time = $line[$i_time];
    if ( not defined $last_time ) { $last_time = $time }
    else { 
      if ( $last_time >= $time ) { 
        #print "Rotamer file '$file' backtracks from $last_time to $time\n";
        next;
      }
      $last_time = $time;
    }

    last if defined $time_limit and $time_limit < $time;
    #last if   Past_Range ( $intervals, $line[$i_time] ); # SLOWER WAY
    next if not In_Range ( $intervals, $time );

    foreach my $i ( 0..$#angle_names ) {
      $index = $angle_indecies[$i];
      return if not $line[$index] =~ /^[\d\.-]+$/;
      #return if not Is_Number ($line[$index]);
      $line[$index] = (sprintf "%.0f", $line[$index]) % 360;
    }

    $line[$i_major]--;
    if ( $line[$i_major] < 0 ) {
      die "rotamer [read_rotamer_file]: read file with incorrectly indexed major rotamers for $res_name in $res_dir.\n";
    }
    push @data, \@line;

    
  }

  close $IN or die "Couldn't close rotamer data file, '$file'. $!";
  
  if ( not @data ) { 
    return;
    die "No data in Read_Rotamer_File for this protein.  $res_name, $res_dir\n";
  }

  my $RES = 
    {
      res_name        => $res_name,
      timepoints      => scalar @data,
      headers         => \%headers,
      angle_names     => $res->{'angle_names'},
      rotamer_names   => $res->{'rotamer_names'},
      data            => \@data,
      files           => [$file],
      residue         => $res,
      for             => 'residue',
      total_rotamers  => $res->{'total_rotamers'},
      populated       => 1,
    };
  
  $self->Generate_Histograms ( $RES );
  
#  print Dumper ($RES);
#  die;
  
  return ( $RES );
}



############################################## IMPORT DIHEDRAL FILES
# Reads raw dihedral data and bins it into rotamers.
#
# I'm going to make a bunch of assumptions about the files,
# and I'm not going to do a lot of data checking, that way this is
# faster, cuz this is the slow part.
# Assumptions:
#   1) the join column contains exactly the same data in each file
#   2) the join column data is in the same order in each file.
#   3) the columns are delimited by something matching /\s+/.
#   4) there is only one data column per file.

sub Read_Dihedral_Files {

  my $self        = $_[0];
  my $res_name    = $_[1];
  my $res_dir     = $_[2];
  my $overwrite   = $_[3] || $self->{'overwrite_dihedrals'} || 0;

  #========== Look for the aggregate file
  my $file = "$res_dir/$self->{'files'}{'cache'}";
  my $file_old = "$res_dir/$self->{'files'}{'dihedral_time'}";
  my @files_raw = glob ( "$res_dir/$self->{'files'}{'dihedral_examples'}" );
  if ( (not -e $file) or $overwrite) {

    # Possibly have compressed the file to save space
    if ( -e "$file.gz" ) { 
      $file = "$file.gz";

    # Check for old summary file or zipped version
    } elsif ( (-e $file_old) and (not $overwrite) ) { 
      $file = $file_old; 
    } elsif ( ( -e "$file_old.gz" ) and ( not $overwrite ) ) { 
      $file = "$file_old.gz";

    # Better test if the data even exists
    } elsif ( not scalar @files_raw ) {
      die "There is not dihedral data in the dihedral directory, '$res_dir/$self->{'files'}{'dihedral_examples'}'\n"
        . "while looking for file, $file\n";
    # And then run the aggregator over the raw data
    } else {
      $self->Run_Aggregator($res_dir)
        or die "Could not aggregate dihedral data to summary file in directory, '$res_dir'.\n";
      $file = $file_old;
    }
  }

  # If it's still not there, something is really broken.
  if ( not -e $file ) { die "Rotamer [Read_Dihedral_File]: Could not find/create dihedral data file in dir: '$res_dir'.\n"; }
  


  #========== Convert the dihedral data to rotamer bins
  
  my $IN;
  open $IN, "<gzip(autopop)", $file or die "Couldn't open file, '$file'. $!";

  #my $IN = Open_File ( $file );
  #open IN, $file or die "Couldn't open dihedral data file, '$file'. $!";

  # for looking up which column is which
  my $headers = lc <$IN>;
  my @headers = split /\s+/, $headers;
  my %headers = ();
  my $i = 0;
  foreach my $h (@headers) {
    $h =~ s/#//g;
    $headers{$h} = $i++;
  }
  my $i_time = $headers{"time"};

  # Figure out which angles are important
  my $AA = $self->{'AA'};
  my $res               = $AA->{$res_name};
  my @angle_names       = grep { $headers{$_} } @{$res->{'angle_names'}};
  my @rotamer_names     = @{$res->{'rotamer_names'}};
  my @angle_indecies    = map { $headers{$_} } @angle_names;
  my @rotamer_indecies  = map { $headers{$_} } @rotamer_names;
  my %confs             = %{$res->{'conformations'}};

  # Make sure we're not missing any columns of dihedral data.
  my @missing;
  if ( @missing = grep { not exists $headers{$_} } @rotamer_names ) {
    die "Dihedral data is incomplete for directory, '$res_dir'.  Missing angles " 
      . (join ', ', @missing) . ".\n";
  }

  # Add the minor rotamer indecies to the header lookup hash
  foreach ( @rotamer_names, 'major' ) { 
    if ( not exists $headers{"rotamer_$_"} ) {
      $headers{"rotamer_$_"} = $i++; 
    }
  }

  

  
  # Cache up the conversion data
  my @quickstarts     = @{$res->{'quickstarts' }};
  my @quicknumbers    = @{$res->{'quicknumbers'}};
  my @quickmults      = @{$res->{'quickmults'  }};
  
 
  my ( @data, $datum, $major, $angle, $index );
  my ( $line, $q, $minor, $k );
  my $last_time;

  while ( defined ( $line = <$IN> ) ) {
    
    next if ( ( $line = Clean_Line ( $line ) ) eq '' );

    $major = 0;
    my @line = split /\s+/, $line;

    my $time = $line[$i_time];
    if ( not defined $last_time ) { $last_time = $time }
    else { 
      if ( $last_time >= $time ) { 
        #print "Rotamer file '$file' backtracks from $last_time to $time.  Skipping line.\n";
        next;
      }
      $last_time = $time;
    }


    foreach my $i ( 0..$#angle_names ) {
      $angle = $angle_names[$i];
      $index = $angle_indecies[$i];
      return if not Is_Number ($line[$index]);
      $line[$index] = (sprintf "%.0f", $line[$index]) % 360;
    }

    foreach my $i ( 0..$#rotamer_names ) {

      $angle = $rotamer_names[$i];
      $index = $rotamer_indecies[$i];
      $datum = $line[$index];
    
      $q = $quickstarts[$i];
      $minor = $#$q;
      foreach $k ( 0..$minor ) {
        if ( $datum < $q->[$k] ) {
          $minor = $k - 1;
          last;
        }
      }
      $line[$headers{"rotamer_$angle"}] = $confs{$angle}[$minor];
      $major += ( $quicknumbers[$i][$minor] - 1 ) * $quickmults[$i];
    }
    if ( $major < 0 ) {
      die "rotamer [read_rotamer_file]: read file with incorrectly indexed major rotamers for $res_name in $res_dir.\n";
    }
    $line[$headers{"rotamer_major"}] = $major;
    push @data, \@line;
  }

  # Close filehandles
  close $IN or die "Couldn't close dihedral data file, '$file'. $!";

  my $RES = 
    {
      res_name        => $res_name,
      headers         => \%headers,
      timepoints      => scalar @data,
      rotamer_names   => $res->{'rotamer_names'},
      angle_names     => $res->{'angle_names'},
      data            => \@data,
      files           => [$file],
      residue         => $res,
      total_rotamers  => $res->{'total_rotamers'},
      for             => 'residue',
      populated       => 1,
    };
  
    #$self->Generate_Histograms ( $RES );
    #
  
  return ( $RES );

}

####################################### CALCULATE SINGLE ROTAMER
# Simplified subroutine that just calculates one rotamer at a time.
#

#sub Calculate_Single_Rotamer { 
# 
#  my $self = $_[0]; 
#  my $resn = $_[1];
#  my $chis = $_[2];
#
#  my $res  = $self->{'AA'}{$resn};
#  
#  # Cache up the conversion data
#  my @quickstarts     = @{$res->{quickstarts }};
#  my @quicknumbers    = @{$res->{quicknumbers}};
#  my @quickmults      = @{$res->{quickmults  }};
#  
#  my @angle_names     = @{$res->{angle_names}};
#  my @rotamer_names   = @{$res->{rotamer_names}};
#  
#  my @minors = ();
#  my $major = 0;
#
#  foreach my $i ( 0..$#rotamer_names ) { 
#    $angle = $rotamer_names[$i];
#    $index = $rotamer_indecies[$i];
#  }
#
#}

####################################### GENERATE DIHEDRAL HISTOGRAMS
#TODO: 

sub Generate_Dihedral_Histograms {
}

######################################## GENERATE ROTAMER HISTOGRAMS
# TODO: generate histograms for individual angles
#
sub Generate_Rotamer_Histograms {

  my $self    = $_[0];
  my $RES     = $_[1];

  my $i_major = $RES->{'headers'}{'rotamer_major'};
  my $data    = $RES->{'data'};

  #----- Major Rotamer Historgram
  my @majors  = ();
  my $major   = 0;
  for my $line ( @$data ) {
    $major = $line->[$i_major];
    if ( not defined $majors[$major] ) { $majors[$major] = 1; }
    else                               { $majors[$major]++;   }
  }
  
  #----- Initialize all the empty ones to 0
  foreach my $i ( 0..$#majors ) { 
    if ( not defined $majors[$i] ) { $majors[$i] = 0; }
  }
  
  #----- Save it in the residue
  $RES->{'histograms'}{'rotamer_major'} = \@majors;
}



####################################################################
#                                                       OUTPUT CACHE
####################################################################


################################################# WRITE ROTAMER FILE

sub Write_Rotamer_Cache_Files {

  my $self  = $_[0];
  my $res   = $_[1];
  my $dir   = $_[2];
  
  
  $self->Write_Rotamer_Cache_Time_File ( $res, $dir );
  return;
  $self->Write_Rotamer_Cache_Histogram_Files ( $res, $dir );
}

sub Write_Rotamer_Cache_Time_File {
  
  my $self  = $_[0];
  my $res   = $_[1];
  my $dir   = $_[2];

  my $file = "$dir/$self->{'files'}{'cache'}.gz";
  #my $file = "$dir/$self->{'files'}{'dihedral_time'}";
  my $file_temp = "$file.temp.gz";
  #print "$file\n$file_temp\n";

  my $res_name        =   $res->{'res_name'};
  my %headers         = %{$res->{'headers'  }};
  my @rotamer_names   = @{$res->{'rotamer_names' }};
  my @angle_names     = @{$res->{'angle_names'}};
  my $data            =   $res->{'data'     };
  my $AA              =   $self->{'AA'};

  open OUT, ">gzip", $file_temp or die "Couldn't open rotamer time file '$file_temp' for writing. $!";

  my $i_time      = $headers{'time'};
  my $i_major     = $headers{'rotamer_major'};
  my @i_rotamers  = map { $headers{"rotamer_$_"} } @rotamer_names;
  my @i_angles    = map { $headers{$_} } @angle_names;
  
  my %conformations = %{$res->{'residue'}{'conformations'}};

  print OUT "#" . ( join "\t", "time", @angle_names, "rotamer_major", map { "rotamer_$_" } @rotamer_names ) . "\n";
  
  for my $line ( @$data ) {
    print OUT join "\t", 
                  $line->[$i_time],
                  ( map { $line->[$_] } @i_angles ),
                  $line->[$i_major]+1,
                  (map { $line->[$headers{"rotamer_$_"}] } @rotamer_names),
                  ;
    print OUT "\n";
  }

  close OUT or die "Couldn't close rotamer time file '$file_temp' after writing. $!";

  `mv $file_temp $file`;
  #`gzip $file`;
  #print "cp $file_temp $file\n";

}




1;
__END__


