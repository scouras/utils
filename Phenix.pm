package Phenix;

BEGIN { 
  push @INC, "$ENV{'CODE_BASE'}/scouras";
  push @INC, "$ENV{'CODE_BASE'}/ParsePDB";
}
use warnings;
use strict;
use Carp;
use Utility;
use Cwd;
use Storable;
use Data::Dumper;
use IPC::Open2;
use IPC::Open3;
use Time::HiRes;
use ParsePDB;
use POSIX ":sys_wait_h";
$Data::Dumper::Indent = 1;
$Data::Dumper::Terse  = 1;
$Data::Dumper::Maxdepth = 6;

#my $MODE = 'RNASEA';
my $MODE = 'LIPOXYGENASE';


###############################################################################
# Phenix.pm
# Alexander Scouras
# University of California, Berkeley
# scouras@berkeley.edu
# Modified 14 April 2013
#------------------------------------------------------------------------------
# TODO:
#   Parse table 1?
#   Change Output_Analysis_Matrix to include altconfs
#   Consider altconfs in Output_Analysis_Regroup
#     Think cannot do so, since not all proteins have same altconfs.
#     But instead of writing main altconf, could write averages over all
#     altconfs, which might be smarter.
#     Just change _ to AA?
#       Need to make sure AA is calculated correctly and reliably.
#
#   Make all analyses detect if input is newer than output.
#   Output Dihedral Angle data?
#     Ah, there it is in DATA/*protein.angle.txt
#------------------------------------------------------------------------------
# NOTES
# -----
# For this module, the special altconf '_' refers to the primary altconf, which
# is the one with the highest occupancy.  There is also a special altconf '*', 
# which is the summary over all altconfs for that atom. The summary function
# is specific to each property, as determined when Collapse_Atom_Properties
# is called.  For instance, B-Factors can be averaged, whereas peak counts 
# would instead be summed, and PISA contact area might be summed.
###############################################################################



##### READ ENVIRONMENT VARIABLES AND OTHER DATA
our $VERBOSE  = 0;
our $START     = time;
our $PWD       = Cwd::getcwd;
our $CODE_BASE = $ENV{'CODE_BASE'};
our $ROSETTA   = $ENV{'ROSETTA3_BIN'};

##### EXPORTS
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw(
        &Can_Run  
        &Create_Structure_Configuration 
        &Is_Structure_Config  

        &Fetch_PDB_Redo  
        &Parse_PDB_Meta 
        &Fetch_PDB_General  
        &Fetch_PDB_and_CIF  

        &XTriage 
        &XTriage_to_String 
        &Print_XTriage  

        &MTZ_Dump  
        &Print_MTZ_Dump 

        &IOBS_to_FOBS  
        &Generate_R_Free  
        &Maps  

        &Rosetta_Holes  
        &Parse_Holes_Data  
        &DSSP  
        &Pisa  
        &Parse_Pisa_Data  
        &Ringer  
        &Parse_Ringer_Data  
        &Aggregate_Structure_Data 


        &diffdump  
        &Spotfinder  
        &Reap_Spotfinder_Job  
        &Parse_Spotfinder  
        &Cleanup 
) ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw();


##############################################################################
#                             KNOWN PROTEINS
#-----------------------------------------------------------------------------
# Just some convenient properties of known proteins for expediting various
# analyses. Should not be relied upon for serious analysis, more intended
# for programs like mtz_summarize_analysis_files.pl which just does a 
# quick parse of an mtz file. 
##############################################################################

our %KNOWN_PROTEINS = (
  RNASEA => { 
    length => 124,
  },

  LIPOXYGENASE => {
    length => 839,
  },
);

##### Adding prefixes to make heuristic lookup of these simpler
$KNOWN_PROTEINS{lipo} = $KNOWN_PROTEINS{LIPOXYGENASE};
$KNOWN_PROTEINS{slo } = $KNOWN_PROTEINS{LIPOXYGENASE};

##############################################################################
#                            OTHER PARAMETERS
##############################################################################

our $ALTCONF_DISTINCTION_DISTANCE = 1.0; # Angstroms


##############################################################################
#                         FORMATS OF ANALYSIS OUTPUT
#-----------------------------------------------------------------------------
# These are not directly the formats used, but general formats that the 
# analysis parsing subroutines can query based on what results they actually
# find. For instance, PISA may have results for specific ligands that cannot
# be predicted, but is always reported similary in absolute and relative
# surface areas.
##############################################################################

# Convenience Table for more condensed data entry
# Column Order:
#   F     | Format for sprintf
#   D     | Default Value
#   PMIN  | Minimum Plot/Histogram Value
#   PMAX  | Maximum Plot/Histogram Value
#   PZS   | Plot ZScale 
#         | lin or log
#   PCS   | Plot Colorscheme
#         | bow, rb, rwb, cyc
#   PZSC  | Plot ZScale Colors
#   PR    | Plot Colorscheme Reversal (boolean)
#   PLOT  | Boolean whether to even plot it
#   HBW   | Histogram Bin Width
#   DESC  | Description
#
our %PLOT_COLORSCHEMES = (
  bow => 'rainbow',
  cyc => 'rainbow-cycle',
  rb  => 'red-blue',
  rwb => 'red-white-blue',
  gry => 'grayscale',
);

our %PLOT_ZSCALES = (
  lin => 'linear', 
  log => 'log',
);

our %META_PROPERTIES_TABLE_RNASEA = (
#                      F          D   PMIN   PMAX   DMIN  DMAX    HBW   ALTS
  PDB     => {
    X             => [ "%.2f",  -10, -35.0,  35.0,-100.0, 100.0,  1.00, '*_'],
    Y             => [ "%.2f",  -10, -35.0,  35.0,-100.0, 100.0,  1.00, '*_'],
    Z             => [ "%.2f",  -10, -35.0,  35.0,-100.0, 100.0,  1.00, '*_'],
    B             => [ "%.2f",  -10,   0.0,  30.0,   0.0, 100.0,  1.00, '*_'],
    Q             => [ "%.2f",   -1,   0.0,   1.0,   0.0,   1.0,  0.10, '_' ],
    DIHEDRALS     => [ "%.2f", -999,  -180,   180,  -180, 180.0,  5.00, '_' ],
    ALT_COUNTS    => [ "%d",     -1,     0,     4,     0,    10,  1.00, '*' ],
    DISTINCT_ALTS => [ "%d",     -1,     0,     4,     0,    10,  1.00, '*' ],
  },
  HOLES   => {
    SCORE         => [ "%.2f",   -1,     0,   1.0,   0.0,   1.0,  0.05, '*' ],
    DEPTH         => [ "%.2f",   -1,     0,   5.0,   0.0,   5.0,  0.50, '*' ],
  },
  PISA    => {
    REL           => [ "%d",   -100,     0,   100,   0.0, 100.0,  5.00, '_'],
    ABS           => [ "%d",   -100,     0,   200,   0.0, 500.0,  5.00, '_'],
    CLASS         => [ "%s",   'NA', undef, undef, undef, undef, undef, '_'],
  },
  RINGER  => { 
    PEAK_COUNTS   => [ "%d",     -1,     0,     5,     0,     5,  1.00, '_'],
    DENSITY_MEAN  => [ "%.2f",  -10,   0.0,   1.0,   0.0,   1.0,  0.10, '_'],
    DENSITY_MIN   => [ "%.2f",  -10,  -1.0,   1.0,  -1.0,   1.0,  0.10, '_'],
    DENSITY_MAX   => [ "%.2f",  -10,   0.0,  10.0,   0.0,  10.0,  1.00, '_'],
  },
  COOT  => {
    ISOr          => [ "%.2f", -100,   0.0,  10.0, -10.0, 100.0,  0.50, '_'],
    DIFr          => [ "%.2f", -100,  -5.0,   5.0, -10.0, 100.0,  0.25, '_'],
    VOLr          => [ "%.2f", -100,   0.0,  10.0, -10.0, 100.0,  0.50, '_'],
    PEAKr         => [ "%.2f", -100,  -5.0,   5.0, -10.0, 100.0,  0.50, '_'],

    ISOa          => [ "%.2f", -100,   0.0,   5.0, -10.0, 100.0,  1.00, '_'],
    DIFa          => [ "%.2f", -100,  -1.0,   1.0, -10.0, 100.0,  0.10, '_'],
    VOLa          => [ "%.2f", -100,   0.0,   5.0, -10.0, 100.0,  1.00, '_'],
    PEAKa         => [ "%.2f", -100,   0.0,   5.0, -10.0, 100.0,  1.00, '_'],

    ANISO         => [ "%.2f",   -1,   0.0,   1.0,   0.0,   1.0,  0.05, '*'],
    ROTAMER       => [ "%s",   'NA', undef, undef, undef, undef, undef, '_'],
  },
);

our %META_PROPERTIES_TABLE_LIPOXYGENASE = (
#                      F          D   PMIN   PMAX   DMIN  DMAX    HBW   ALTS
  PDB     => {
    X             => [ "%.2f",  -10,-100.0, 100.0,-100.0, 100.0,  1.00, '*_'],
    Y             => [ "%.2f",  -10,-100.0, 100.0,-100.0, 100.0,  1.00, '*_'],
    Z             => [ "%.2f",  -10,-100.0, 100.0,-100.0, 100.0,  1.00, '*_'],
    B             => [ "%.2f",  -10,   0.0,  75.0,   0.0, 100.0,  1.00, '*_'],
    Q             => [ "%.2f",   -1,   0.0,   1.0,   0.0,   1.0,  0.05, '_' ],
    DIHEDRALS     => [ "%.2f", -999,  -180,   180,  -180, 180.0,  5.00, '_' ],
    ALT_COUNTS    => [ "%d",     -1,     0,     4,     0,    10,  1.00, '*' ],
    DISTINCT_ALTS => [ "%d",     -1,     0,     4,     0,    10,  1.00, '*' ],
  },
  HOLES   => {
    SCORE         => [ "%.2f",   -1,     0,   1.0,   0.0,   1.0,  0.05, '*' ],
    DEPTH         => [ "%.2f",   -1,     0,   5.0,   0.0,   5.0,  0.50, '*' ],
  },
  PISA    => {
    REL           => [ "%d",   -100,     0,   100,   0.0, 100.0,  5.00, '_'],
    ABS           => [ "%d",   -100,     0,   200,   0.0, 500.0,  5.00, '_'],
    CLASS         => [ "%s",   'NA', undef, undef, undef, undef, undef, '_'],
  },
  RINGER  => { 
    PEAK_COUNTS   => [ "%d",     -1,     0,     5,     0,     5,  1.00, '_'],
    DENSITY_MEAN  => [ "%.2f",  -10,   0.0,   1.0,   0.0,   1.0,  0.10, '_'],
    DENSITY_MIN   => [ "%.2f",  -10,  -1.0,   1.0,  -1.0,   1.0,  0.10, '_'],
    DENSITY_MAX   => [ "%.2f",  -10,   0.0,  10.0,   0.0,  10.0,  1.00, '_'],
  },
  COOT  => {
    ISOr          => [ "%.2f", -100,   0.0,  10.0, -10.0, 100.0,  0.50, '_'],
    DIFr          => [ "%.2f", -100,  -5.0,   5.0, -10.0, 100.0,  0.25, '_'],
    VOLr          => [ "%.2f", -100,   0.0,  10.0, -10.0, 100.0,  0.50, '_'],
    PEAKr         => [ "%.2f", -100,  -5.0,   5.0, -10.0, 100.0,  0.50, '_'],

    ISOa          => [ "%.2f", -100,   0.0,   5.0, -10.0, 100.0,  1.00, '_'],
    DIFa          => [ "%.2f", -100,  -1.0,   1.0, -10.0, 100.0,  0.10, '_'],
    VOLa          => [ "%.2f", -100,   0.0,   5.0, -10.0, 100.0,  1.00, '_'],
    PEAKa         => [ "%.2f", -100,   0.0,   5.0, -10.0, 100.0,  1.00, '_'],

    ANISO         => [ "%.2f",   -1,   0.0,   1.0,   0.0,   1.0,  0.05, '*'],
    ROTAMER       => [ "%s",   'NA', undef, undef, undef, undef, undef, '_'],
  },
);

our %META_PLOT_TABLE = (
#                      O   P   PZS    PCS  PR  PZSC  DESC          
  PDB     => {
    #B             => [ 1, 1, 'lin', 'bow', 1,    4, 'B-Factors',           ],
    #Q             => [ 1, 1, 'lin', 'bow', 0,    4, 'Occupancies',         ],
    X             => [ 1, 0, 'lin', 'rwb', 1,    4, 'X Coordinate',        ],
    Y             => [ 1, 0, 'lin', 'rwb', 1,    4, 'Y Coordinate',        ],
    Z             => [ 1, 0, 'lin', 'rwb', 1,    4, 'Z Coordinate',        ],
    B             => [ 1, 1, 'lin', 'rwb', 1,    4, 'B-Factors',           ],
    Q             => [ 1, 1, 'lin', 'rwb', 0,    4, 'Occupancies',         ],
    DIHEDRALS     => [ 1, 0, 'lin', 'cyc', 0,    0, 'Dihedral Angles',     ],
    ALT_COUNTS    => [ 1, 1, 'lin', 'rwb', 0,    1, 'Alt Confs',           ],
    DISTINCT_ALTS => [ 1, 1, 'lin', 'rwb', 0,    1, 'Distinct AltConfs',   ],
  },
  HOLES   => {
    SCORE         => [ 0, 0, 'lin', 'rwb', 0,    0, 'Holes Score',         ],
    DEPTH         => [ 0, 0, 'lin', 'rwb', 0,    0, 'Holes Depth',         ],
  },
  PISA    => {
    REL           => [ 1, 1, 'lin', 'rwb', 1,    0, undef,                 ],
    ABS           => [ 0, 0, 'lin', 'rwb', 1,    0, undef,                 ],
    CLASS         => [ 1, 0, 'lin', 'rwb', 1,    0, undef,                 ],
  },
  RINGER  => { 
    PEAK_COUNTS   => [ 1, 1, 'lin', 'rwb', 1,    0, 'Ringer Peak Counts',  ],
    DENSITY_MEAN  => [ 1, 1, 'lin', 'rwb', 0,    0, 'Mean Angle Density',  ],
    DENSITY_MIN   => [ 1, 1, 'lin', 'rwb', 0,    0, 'Min Angle Density',   ],
    DENSITY_MAX   => [ 1, 1, 'lin', 'rwb', 0,    0, 'Max Angle Density',   ],
  },
  COOT  => {
    ISOr          => [ 1, 1, 'lin', 'rwb', 0,    4, '2Fo-Fc RMSD',         ],
    DIFr          => [ 1, 1, 'lin', 'rwb', 0,    2, 'Fo-Fc RMSD',          ],
    VOLr          => [ 1, 1, 'lin', 'rwb', 0,    4, 'Volume Density RMSD', ],
    PEAKr         => [ 1, 1, 'lin', 'rwb', 0,    4, 'Peak Density RMSD',   ],
    
    ISOa          => [ 0, 0, 'lin', 'rwb', 0,    0, '2Fo-Fc Absolute',     ],
    DIFa          => [ 0, 0, 'lin', 'rwb', 0,    0, 'Fo-Fc Absolute',      ],
    VOLa          => [ 0, 0, 'lin', 'rwb', 0,    0, 'Volume Density Absolute',     ],
    PEAKa         => [ 0, 0, 'lin', 'rwb', 0,    0, 'Peak Density Absolute',     ],
    
    ANISO         => [ 1, 1, 'lin', 'rwb', 0,    0, 'Anisotropy',          ],
    ROTAMER       => [ 1, 0, undef, undef, 0,    0, 'Rotamers',            ],
  },
);

my %META_PROPERTIES_TABLE = %META_PROPERTIES_TABLE_RNASEA;
if ( $MODE eq 'LIPOXYGENASE' ) { 
  %META_PROPERTIES_TABLE = %META_PROPERTIES_TABLE_LIPOXYGENASE;
}

our %META_PROPERTIES = ();
foreach my $anal ( keys %META_PROPERTIES_TABLE ) { 
  foreach my $prop ( keys %{$META_PROPERTIES_TABLE{$anal}} ) { 
    my $t_prop = $META_PROPERTIES_TABLE{$anal}{$prop};
    my $t_plot = $META_PLOT_TABLE{$anal}{$prop};

    my ( $format, $default, $pmin, $pmax, $dmin, $dmax, $hist_bin_width, $alt ) = @$t_prop;
    my ( $do_output, $do_plot, $zscale, $cs, $cs_rev, $zcolors, $desc ) = @$t_plot;

    $META_PROPERTIES{$anal}{$prop} = {
      F     => $format,
      D     => $default,
      PMIN  => (sprintf "%.4f", ($pmin || 0)),
      PMAX  => (sprintf "%.4f", ($pmax || 0)),
      DMIN  => (sprintf "%.4f", ($dmin || 0)),
      DMAX  => (sprintf "%.4f", ($dmax || 0)),
      PZS   => ($zscale ? $PLOT_ZSCALES     {$zscale} : 'linear' ),
      PCS   => ($cs     ? $PLOT_COLORSCHEMES{$cs    } : 'rainbow'), 
      PR    => $cs_rev,
      HBW   => $hist_bin_width,
      PZSC  => $zcolors,
      PLOT  => $do_plot,
      OUT   => $do_output,
      DESC  => $desc,
      ALT   => $alt,
    };
  }
}

$META_PROPERTIES{PISA}{CLASS  }{TEXT} = 1;
$META_PROPERTIES{COOT}{ROTAMER}{TEXT} = 1;

##### Customize some PISA properties
#IPA BTB EPE FE
for my $pisa_text ( qw(SASA Biomer Crystal Ligand Metal Anion Cation TASA Buried ASym Sym), keys %Utility::KNOWN_RESIDUES ) { 
  my $pisa_key = uc $pisa_text;
  for my $prefix ( qw(ABS REL) ) { 
    my $desc = "Pisa $pisa_text";
    if ($prefix eq 'ABS') { $desc .= " Absolute" }
    $META_PROPERTIES{PISA}{"$prefix\_$pisa_key"} = { %{$META_PROPERTIES{PISA}{$prefix}} };
    $META_PROPERTIES{PISA}{"$prefix\_$pisa_key"}{DESC} = $desc;
  }
}
for my $pisa ( qw(CRYSTAL LIGAND EPE IPA BTB ACT SYM ASYM) ) { 
  $META_PROPERTIES{PISA}{"REL_$pisa"}{PMAX} = 25.0;
}

foreach my $pisa ( qw(BIOMER ASYM SYM TASA BTB ) ) { 
  $META_PROPERTIES{PISA}{"REL_$pisa"}{PLOT} = 0;
}

##### Customize Ringer properties
$META_PROPERTIES{RINGER}{PEAK_COUNTS_CHI1         } = { %{$META_PROPERTIES{RINGER}{PEAK_COUNTS}}};
$META_PROPERTIES{RINGER}{PEAK_COUNTS_DISCRETE     } = { %{$META_PROPERTIES{RINGER}{PEAK_COUNTS}}};
$META_PROPERTIES{RINGER}{PEAK_COUNTS_DISCRETE_CHI1} = { %{$META_PROPERTIES{RINGER}{PEAK_COUNTS}}};

$META_PROPERTIES{RINGER}{PEAK_COUNTS_CHI1         }{DESC} = "Peak Counts Chi1";
$META_PROPERTIES{RINGER}{PEAK_COUNTS_DIESCRETE    }{DESC} = "Peak Counts Discrete";
$META_PROPERTIES{RINGER}{PEAK_COUNTS_DISCRETE_CHI1}{DESC} = "Peak Counts Discrete Chi1";


################################################################ DEFAULT PATHS
our %DEFAULT_PATHS = (
  DIRS => {
    PROPERTIES  => "PROPERTIES",
    HOLES       => "HOLES",
    PISA        => "PISA",
    MAPS        => "MAPS",
    DATA        => "DATA",
  },
);
    
# Occupancies of alternate conformations should sum to 1.0, to within this tolerance
our $OCCUPANCY_DIFFERENCE_THRESHOLD = 0.05; 

##############################################################################
#                          FLAGS TO UTILITY FUNCTIONS
##############################################################################

our %FLAGS_TRY_COMMAND = (
  verbose => 1,
);

our %FLAGS_REPORT_ERROR = (
  write_file  => 1,
  fatal       => 0,
);


#----------------- Pre and post-verify analysis files
our %FLAGS_VERIFY_PRE = (
  verify      => 1,
  fatal_in    => 1,
  fatal_out   => 0,
);

our %FLAGS_VERIFY_POST = (
  verify      => 1,
  fatal_in    => 0,
  fatal_out   => 0,
);

our $DEFAULT_IOBS_LABEL = "IOBS";
our $DEFAULT_FOBS_LABEL = "FOBS";
our $DEFAULT_FREE_LABEL = "R-free-flags";


##############################################################################
#                        MISCELLANEOUS OTHER VARIABLES
##############################################################################
my $DISCRETE_PEAK_CUTOFF = 15;   # Minimum angle between discrete peaks
my $PEAK_DENSITY_MINIMUM =  0.3; # Minimum density for a peak




##############################################################################
#                              PROGRAM DETECTION
##############################################################################

our %PROGRAMS = ();
$PROGRAMS{FETCH_PDB     } = Utility::Get_Exe("phenix.fetch_pdb");
$PROGRAMS{CIF_AS_MTZ    } = Utility::Get_Exe("phenix.cif_as_mtz");
$PROGRAMS{MTZ_DUMP      } = Utility::Get_Exe("phenix.mtz.dump");
$PROGRAMS{READY_SET     } = Utility::Get_Exe("phenix.ready_set");
$PROGRAMS{MAPS          } = Utility::Get_Exe("phenix.reflection_file_converter");
$PROGRAMS{RINGER        } = Utility::Get_Exe("$CODE_BASE/ringer/ringer");
$PROGRAMS{HOLES         } = Utility::Get_Exe("$CODE_BASE/scouras/holes.pl");
$PROGRAMS{PISA          } = Utility::Get_Exe("$CODE_BASE/scouras/parse_pisa.pl");
$PROGRAMS{XTRIAGE       } = Utility::Get_Exe("phenix.xtriage");
$PROGRAMS{DIFFDUMP      } = Utility::Get_Exe("diffdump");
$PROGRAMS{SPOT_SERVER   } = Utility::Get_Exe("distl.mp_spotfinder_server_read_file");
$PROGRAMS{SPOT_CLIENT   } = Utility::Get_Exe("distl.thin_client");
$PROGRAMS{COOT          } = Utility::Get_Exe("/Library/Coot/bin/coot-real");
if ( not defined $PROGRAMS{COOT} ) { 
  $PROGRAMS{COOT        } = Utility::Get_Exe("/programs/x86_64-linux/coot/0.7-pre-1-r4135/bin/coot");
}
$PROGRAMS{COOT_PROPS    } = Utility::Get_Exe("$CODE_BASE/coot/coot_atom_density.py");
$PROGRAMS{PYMOL         } = Utility::Get_Exe("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL");
$PROGRAMS{PYMOL_UTILITY } = Utility::Get_Exe("$CODE_BASE/pymol/PyMOL.py");
$PROGRAMS{PYMOL_WATERS  } = Utility::Get_Exe("$CODE_BASE/pymol/pymol_waters.py");

if ( not defined $PROGRAMS{COOT_PROPS} ) { 
  confess "No COOT_PROPS";
}


sub Can_Run { 
  my $program = $_[0];
  my $verbose = $_[1];
  my $fatal   = $_[2] || 0;


  my $prog_uc = uc $program;
  my $error   = "";
  if ( not $PROGRAMS{$prog_uc} ) { 
    $error = "Cannot run program $program. Cannot find executable.\n";
    confess $error if $fatal;
    print $error;
    return 0;
  }


  if ( $prog_uc eq 'HOLES' ) { 
    if ( not defined $ROSETTA ) { 
      $error = "Cannot run $program. Requires an installation of rosetta.\n";
      confess $error if $fatal;
      print $error;
      return 0;
    }

    my @holes = glob("$ROSETTA/holes.*");
    if ( scalar @holes == 0 ) {
      $error = "Cannot run $program. Requires an installation of rosetta holes.\n";
      confess $error if $fatal;
      print $error;
      return 0;
    }
  }


  return 1;

}




##############################################################################
##############################################################################
#                         CONFIGURATION OBJECT
##############################################################################
##############################################################################

# Configuration Object Structure
# TITLE - Main title for this PDB4 or novel crystal structure
# JOB   - Suffix for analysis directories for custom parameters 
#       - or stages in refinement, etc.
# BASE  - Base filename for intermediate data files created, such as ccp4 maps
#
# FILES - Paths for standard files
#   PDB 
#   MTZ 
#   MAP     - Map Coefficients for Ringer
#   DUMP    - MTZ dump file
#   XTRIAGE - XTriage log file
#   SCALA
#   POINTLESS
#   TRUNCATE
#   TABLE1
# 
# --- Analyses below are only initialized as they are run/read
# MTZ     - MTZ Properties, mostly from mtz.dump
# XTRIAGE - XTriage Analysis
# META    - Meta data from the PDB file
# PISA    - Pisa contact surface area data
# HOLES   - 
# RINGER  - 
# ANGLES  - Omega, phi, psi, and chi angles and rotamers 
# 
# 
# DO      - Which analyses to perform this round
# FLAGS   - Flags to configure error reporting, etc.
#   pdb_data    - only parse pdb coordinates, not metadata
#   pdb_verbose - checkpoints in pdb parsing
#
#------------------------------------------------------------------------------
# TODO:
#   Create analysis directory
#     Initialize ringer, etc. analysis dirs and titles
#     Sort out what belongs in Dump vs. Files vs. PDB vs. Main
#      

sub Create_Structure_Configuration {

  my %args = @_;

  # Files
  my $dir           = $args{dir         } || Cwd::cwd(); # Output for all data
  my $pdb           = $args{pdb         };
  my $mtz           = $args{mtz         };
  my $subdirs       = $args{subdirs     } || {} ; # Optionally specify analysis subdirectories

  # Names and Such
  my $name          = $args{name        }; # A human meaningful name for the structure
  my $pdb4          = $args{pdb4        }; # A pdb4 code, though this will be autodetected in title
  my $desc          = $args{desc        }; # A description of the structure
  my $title         = $args{title       }; # A full title for all output with this structure
  my $base          = $args{base        }; # The title converted to filesytem friendly text

  # Attributes of this Structure
  my $crystal       = $args{crystal     }; # Specifies the crystal and dataset
  my $dfxntemp      = $args{dfxntemp    }; # Diffraction Temperature
  my $job           = $args{job         }; # Specifies the refinement state, typically 
    #TODO: delete jobs? Maybe useful for sorting analysis output, but tricky.
  my $annotations   = $args{annotations }; # Generalized annotations
  my $anno_list     = $args{anno_list   }; # List and order of annotation output

  # Flags
  my $params        = $args{params      };
  my $do            = $args{do          };
  my $flags         = $args{flags       };


  #============================================== Validate Arguments

  #---------------- Directory
  if ( not -e $dir ) { confess "Cannot find directory, '$dir'" }
  if ( not -d $dir ) { confess "Directory '$dir' is not a directory." }
  chdir $dir or confess "Couldn't chdir to '$dir'. $!";
  
  #---------------- PDB File
  my ( $pdb_path, $pdb_dir, $pdb_file );
  if ( defined $pdb ) { 
    if ( (not -e $pdb) and ( -e "$dir/$pdb" ) ) { 
      $pdb = "$dir/$pdb";
    }
    ( $pdb_path, $pdb_dir, $pdb_file ) = Utility::Split_Path($pdb);
  }

  #---------------- MTZ File
  my ( $mtz_path, $mtz_dir, $mtz_file );
  if ( defined $mtz ) { 
    if ( (not -e $mtz) and ( -e "$dir/$mtz" ) ) { 
      $mtz = "$dir/$mtz";
    }
    ( $mtz_path, $mtz_dir, $mtz_file ) = Utility::Split_Path($mtz);
  }

  #---------------- Title
  if ( not defined $title ) { 
    if ( defined $pdb_file ) { 
      $title = $pdb_file;
      $title =~ s/\.pdb//;
    } elsif ( defined $mtz_file ) { 
      $title = $mtz_file;
      $title =~ s/\.mtz//;
    } else { 
      $title = "unknown";
    }
  }

  #---------------- Annotations
  $annotations = {} if not defined $annotations;
  $anno_list   = [] if not defined $anno_list;

  #---------------- Clean Base Filename for Analyses
  if ( not defined $base ) { 
    if ( defined $title ) { 
      $base = Utility::Clean_Filename($title);
    } else { 
      $base = "";
      confess "We really need a base filename.\n";
    }
  }

  #---------------- TODO: Set default flags


  #============================================== Construct CFG Object
  my %CFG = (

    DIR       => $dir,
    TITLE     => $title,
    PDB4      => $pdb4,
    NAME      => $name,
    DESC      => $desc,

    CRYSTAL   => $crystal,
    DFXNTEMP  => $dfxntemp,
    JOB       => $job,

    FILES => {},
    PATHS => {},
    DIRS  => {},

    ANNOTATIONS => $annotations,
    ANNO_LIST   => $anno_list,

    PROPS => { 
      LIST      => [],
      ATOMS     => [],
      ANGLES    => [],
      RESIDUES  => [],
      PROTEIN   => [],
      LISGAND   => [],
      WATER     => [],
      LOOKUP    => {},
    },

    DO        => $do,
    FLAGS     => $flags,
    PARAMS    => $params,

    STRUCTURE_CONFIG => 1,
  );

  
  # Add paths
  my $base_path = "$dir/$base";
  if ( defined $pdb_path  ) { Utility::Add_Path_to_Config(\%CFG, 'PDB',  $pdb_path) }
  if ( defined $mtz_path  ) { Utility::Add_Path_to_Config(\%CFG, 'MTZ',  $mtz_path) }
  if ( defined $base_path ) { Utility::Add_Path_to_Config(\%CFG, 'BASE', $base_path) }

    


  # Autodetect pdb codes
  if ( (not $CFG{PDB4}) and $title =~ /^\d\w\w\w$/ ) { $CFG{PDB4} = $title }

  # Add optional directories
  foreach my $anal ( keys %$subdirs ) { 
    my $anal_subdir = $subdirs->{$anal};
    my $anal_dir    = "$dir/$subdirs->{$anal}";
    if ( not -e $anal_dir ) { 
      mkdir $anal_dir or confess "Couldn't mkdir anal dir '$anal_dir'. $!";
    }
    $CFG{SUBDIRS}{$anal} = $anal_subdir;
    $CFG{DIRS   }{$anal} = Cwd::abs_path($anal_dir);

  }

  # Check directory for previous error flags
  if ( -e "$dir/NO_PDB_REDO"          ) { $CFG{ERRORS}{NO_REDO} = 1 }
  if ( -e "$dir/NO_STRUCTURE_FACTORS" ) { $CFG{ERRORS}{NO_SF  } = 1 } 
 
  if (wantarray) { return  %CFG }
  else           { return \%CFG }

}

#=========================================== Quick Test if Argument is a Config
sub Is_Structure_Config { 
  my $CFG = $_[0];
  return 0 if not (ref $CFG eq 'HASH');
  return 0 if not exists $CFG->{STRUCTURE_CONFIG};
  return 0 if not $CFG->{STRUCTURE_CONFIG};
  return 1;
}


###############################################################################
#            Register Atom, Angle, and Residue Analyses
#------------------------------------------------------------------------------
# All analysis data must be registered in a standardized format after it is
# added to the CFG object.  Each prop/anal pair is stored in the CFG as:
#   $CFG->{ANAL}{PROP}{CHAIN}[RESI]{ATOM}{ALT} = $VALUE
#
# Registration indicates what resolution level the property applies to
# (atom, angle, residue) which changes what "atom names" are used to
# access the data.  The anal/prop names are given. Custom formats may be 
# specified for each, but this is being deprecated in favor of centralizing
# all the formats on the global %META_PROPERTIES hash,
# Finally, registration needs to specify what types of residues the analysis
# might apply to.  Possible types are taken from Utility.pm but are listed below.
#
#  PROTEIN   => 'P',
#  WATER     => 'W',
#  SOLVENT   => 'S', # Not using currently
#  ANION     => 'A', # Not using
#  CATION    => 'C', # Not using
#  METAL     => 'M', # Not using
#  LIGAND    => 'L',
#  UNKOWN    => 'U', 
###############################################################################


our %ANALYSIS_LEVELS = (
  ATOM => 1,
  ANGLE => 1,
  RESIDUE => 1,
);

our %ANALYSIS_TYPES = (
  PROTEIN => 1, P => 1,
  WATER   => 1, W => 1,
  LIGAND  => 1, L => 1,
);


sub Register_Analyses {
  my $CFG     = $_[0];
  my $levels  = $_[1];
  my $anal    = $_[2];
  my $props   = $_[3];
  my $types   = $_[4] || ['P'];
  my $formats = $_[5] || {};


   
  foreach my $prop ( @$props ) { 
    my $prop_hash = {
      anal    => $anal,
      prop    => $prop,
      types   => { map  { $_=>1, $Utility::RESIDUE_TYPES{$_}=>1 } @$types },
      levels  => { (map { $_=>1 } @$levels) },
      #props   => { (map { $_=>1 } @$props) },
    };
    push @{$CFG->{PROPS}{LIST}}, $prop_hash;
    $CFG->{PROPS}{LOOKUP}{$anal}{$prop} = $prop_hash;

    foreach my $level ( @$levels ) { 
      if ( not exists $ANALYSIS_LEVELS{$level} ) { 
        confess "Unknown analysis level '$level'.\n";
      }
      push @{$CFG->{PROPS}{$level}}, $prop_hash;
    }
  }
  foreach my $key ( keys %$formats ) { 
    $CFG->{FORMATS}{$anal}{$key} = $formats->{$key};
  }
}


###############################################################################
#                       CONSENSUS SEQUENCE AND PROPERTY SET
#------------------------------------------------------------------------------
# In order to have standardized outputs when several structures are analysed,
# thus entailing several separate $CFGs, this function registers all the
# analyses from a $CFG with a passed hash and list.  These "global" structs
# can then be passed to output functions (particularly 
# Output_Analysis_Regroup()) or used in other ways.  
#
# Collect Properties() and Unify_Output_Properties() perform similar, but 
# distinct functions.
#
# $REGROUP{LEVEL}{BASE}{ANAL}{PROP} = ANAL_PROP_FROM_BASE
# $LIST{LEVEL}{ANAL}{PROP} = FIRST_ANAL_PROP
#
# REGROUP actually contains the data from every structure/analysis/property.
#
###############################################################################

#==================================================== Generate Consensus Config
# The consensus object keeps a list of the consensus sequence in a standard
# property structure. Thus it needs to track:
#   PDB
#     CHAINS
#     RESI
#     RESN
#     TYPE
#     ALTS
#     ATOMS
# It gathers these properties from each CFG, inserting new atoms when they're
# undefined, and checking that the atoms are the same otherwise.

sub Generate_Consensus_Config {

  my @CFGS = @_;

  my %PDB = (
    CHAINS  => [],
    RESI    => {},
    RESN    => {},
    TYPE    => {},
    ALTS    => {},
    ATOMS   => {},
  );

  my $discrepancies = 0;

  #================================================================== Loop CFGS
  foreach my $CFG ( @CFGS ) {

    #============================================================== Loop Chains
    foreach my $chain ( @{$CFG->{PDB}{CHAINS}} ) {

      #----------- Add missing chains
      if ( not grep { $_ eq $chain } @{$PDB{CHAINS}} ) { 
        push @{$PDB{CHAINS}}, $chain;
        $PDB{RESI }{$chain} = [];
        $PDB{RESN }{$chain} = [];
        $PDB{TYPE }{$chain} = [];
        $PDB{ALTS }{$chain} = [];
        $PDB{ATOMS}{$chain} = [];
      }

      #========================================================== Loop Residues
      foreach my $resi ( 0..$#{$CFG->{PDB}{RESI}{$chain}} ) {
        
        my $resn  = $CFG->{PDB}{RESN }{$chain}[$resi];
        my $type  = $CFG->{PDB}{TYPE }{$chain}[$resi];
        my $alts  = $CFG->{PDB}{ALTS }{$chain}[$resi];
        my $atoms = $CFG->{PDB}{ATOMS}{$chain}[$resi];
        next if not defined $resn;

        #--------- Check/Add RESI/RESN/TYPE
        my $con_resi  = $PDB{RESI }{$chain}[$resi];
        my $con_resn  = $PDB{RESN }{$chain}[$resi];
        my $con_type  = $PDB{TYPE }{$chain}[$resi];
        my $con_alts  = $PDB{ALTS }{$chain}[$resi];
        my $con_atoms = $PDB{ATOMS}{$chain}[$resi];

        #--------- Newly found residue
        if ( not defined $con_resn ) { 
          $PDB{RESI }{$chain}[$resi] = $resi;
          $PDB{RESN }{$chain}[$resi] = $resn;
          $PDB{TYPE }{$chain}[$resi] = $type;
          $PDB{ALTS }{$chain}[$resi] = $alts;
          $PDB{ATOMS}{$chain}[$resi] = $atoms;
       
        #--------- Compare to existing residue 
        } else { 
          if ( $con_resi ne $resi ) { 
            print "Sequence discrepancy in $CFG->{TITLE} resi $resi $resn: resi $resi vs. consensus $con_resi\n";
            #$discrepancies++;
          } elsif ( $con_resn ne $resn ) { 
            print "Sequence discrepancy in $CFG->{TITLE} resi $resi $resn: resn $resn vs. consensus $con_resn\n";
            #$discrepancies++;
          } elsif ( $con_type ne $type ) { 
            print "Sequence discrepancy in $CFG->{TITLE} resi $resi $resn: type $type vs. consensus $con_type\n";
            #$discrepancies++;
          }
          my @atoms     = sort ( grep { not /H/ } @$atoms );
          my @con_atoms = sort ( grep { not /H/ } @$con_atoms );
          if ( not Utility::Identical_Arrays( \@atoms, \@con_atoms ) ) { 
            print "Sequence discrepency in $CFG->{TITLE} resi $resi $resn: atoms @$atoms vs. consensus @$con_atoms\n";
            #$discrepancies++;
          }
          
          
          #------- Alts do not need to be identical. Add all of them together.
          $PDB{ALTS}{$chain}[$resi] = Utility::Unique( $con_alts, $alts );

        }
      }
    }
  }

  if ( $discrepancies > 1 ) { 
    die "$discrepancies discrepancies found in protein sequences.\n";
  }

  my ( $DATA, $PROPS ) = Phenix::Unify_Output_Properties( @CFGS );
  my %CON = (
    PDB   => \%PDB,
    DATA  => $DATA,
    PROPS => $PROPS,
    CFGS  => \@CFGS,
    TITLES  => [ map { $_->{TITLE} } @CFGS ],
  );

  if ( wantarray ) { return  %CON }
  else             { return \%CON }

}


#============================================================= Regroup Analyses

sub Regroup_Analyses { 

  my $CFG     = $_[0];
  my $regroup = $_[1];
  my $list    = $_[2];

  my $base = $CFG->{FILES}{BASE};

  foreach my $level ( keys %ANALYSIS_LEVELS ) { 
    foreach my $analprop ( @{$CFG->{"PROP_$level"}} ) { 
      my $anal = $analprop->{anal};
      my $prop = $analprop->{prop};
      $regroup->{$level}{$base}{$anal}{$prop} = $CFG->{$anal}{$prop};
      if ( not exists $list->{$level}{$anal}{$prop} ) { 
        $list->{$level}{$anal}{$prop} = $analprop;
      }
      $list->{$level}{$anal}{$prop}{count}++;
    }
  }
}


#============================================== Collect Properties from CFG Set
# Collects all properties across all CFGs, then copies them to an OUT_PROP_$LEVEL
# array for use with Output_Analysis_Matrix subroutine.  This ensures all
# structures output consistent properties for import into excel or whatever.
#
sub Collect_Properties { 

  my @CFGs = @_;
  if ( ref($CFGs[0]) eq 'ARRAY' ) { @CFGs = @{$CFGs[0]} }

  my %DATA    = ();
  my %PROPS   = ();

  foreach my $cfg ( @CFGs ) { 
  
    my $structure = $cfg->{TITLE};
    my @analprops = @{$cfg->{PROPS}{LIST}};

    foreach my $analprop ( @analprops ) { 
      my $anal = $analprop->{anal};
      my $prop = $analprop->{prop};

      # Save the analprop to the search list and lookup hash
      if ( not exists $PROPS{LOOKUP}{$anal}{$prop} ) { 
        push @{$PROPS{LIST}}, $analprop;
        foreach my $level ( keys %{$analprop->{levels}} ) { 
          push @{$PROPS{$level}}, $analprop;
        }
        $PROPS{LOOKUP}{$anal}{$prop} = $analprop;
        $analprop->{structures} = [];
        $analprop->{count} = 0;
        #print "Added $anal $prop $structure\n";
      }

      # Add this structure to the list of structures
      push @{$analprop->{structures}}, $structure;
      $analprop->{count}++;

      # Save the data for this structure
      $DATA{$structure}{$anal}{$prop} = $cfg->{$anal}{$prop};
    }
  }

  return \%DATA, \%PROPS;
}


sub Unify_Output_Properties { 
  my @CFGs = @_;
  if ( ref($CFGs[0]) eq 'ARRAY' ) { @CFGs = @{$CFGs[0]} }

  my ($DATA, $PROPS) = Collect_Properties(\@CFGs);

  foreach my $cfg ( @CFGs ) { 
    $cfg->{PROPS} = $PROPS;
    #$cfg->{PROPS_OUT} = $PROPS;
  }
  return $DATA, $PROPS;
}




#================================================== Get/Verify Appropriate Protperties
sub Get_Appropriate_Properties {
  my $props = $_[0];
  my $level = $_[1];
  my $type  = $_[2];

  # Verify analysis level exists
  if ( not exists $ANALYSIS_LEVELS{$level} ) { 
    confess "Unknown analysis resolution level '$level'.\n";
  }

  # Verify analysis type exists
  if ( $type ne '*' and not exists $ANALYSIS_TYPES{$type} ) { 
    confess "Unknown analysis atom type '$type'.\n";
  }

  my @props = ();
  foreach my $analprop ( @{$props->{LIST}} ) { 
    next if not $analprop->{levels}{$level};
    next if not $analprop->{types }{$type };
    push @props, $analprop;
  }
  @props = sort { ($a->{anal} cmp $b->{anal}) or ($a->{prop} cmp $b->{prop}) } @props;

  #print "Appropriate Properties for $level/$type\n" . Dumper(\@props);
  if ( wantarray ) { return  @props }
  else             { return \@props }
}

#=============================================== Verify Appropriate Properties
# For a single analysis, verifies that the properties are appropriate for
# a list of resolution levels and atom types.
#
sub Verify_Appropriate_Properties { 
  my $master  = $_[0];
  my $anal    = $_[1];
  my $props   = $_[2];
  my $levels  = $_[3];
  my $types   = $_[4];

  my @props  = ref $props  eq 'ARRAY' ? @$props  : ($props );
  my @levels = ref $levels eq 'ARRAY' ? @$levels : ($levels);
  my @types  = ref $types  eq 'ARRAY' ? @$types  : ($types );

  #---- Verify known property for analysis
  foreach my $prop ( @props ) { 
    if ( not exists $master->{LOOKUP}{$anal}{$prop} ) { 
      confess "Unknown analysis $anal $prop.";
    }
  }

  #---- Verify known resolution level
  foreach my $level ( @levels ) { 
    if ( not exists $ANALYSIS_LEVELS{$level} ) { 
      confess "Unknown analysis level '$level'.";
    }
  }
 
  #---- Verify atom types and levels are appropriate to properties
  foreach my $prop ( @props ) {     
    foreach my $level ( @levels ) { 
      if ( not exists $master->{LOOKUP}{$anal}{$prop}{levels}{$level} ) { 
        confess "Analysis $anal $prop is not appropriate for resolution level $level\n";
      }
    } 

    foreach my $type ( @types ) { 
      if ( not exists $master->{LOOKUP}{$anal}{$prop}{types}{$type} ) { 
        confess "Analysis $anal $prop is not appropriate for atom type $type\n";
      }
    }
  }
  return 1;
}






###############################################################################
#                        Assign Default Altconfs
#------------------------------------------------------------------------------
# For analprops where it is inconvenient or not sensical to Collapse the 
# atom properties, this subroutine takes care of the necessary task of
# giving assignments to '_' and '*', that is making them the same at the main
# altconf.  If not alts are passed, assume that '_' is assigned and copy that
# to '*'.
###############################################################################

sub Assign_Default_Altconfs { 
  my $data = $_[0];
  my $alts = $_[1];

  foreach my $chain ( keys %$data ) { 
    next if not defined $data->{$chain};
    foreach my $resi ( 0..$#{$data->{$chain}} ) { 
      next if not defined $data->{$chain}[$resi];
      foreach my $atom ( keys %{$data->{$chain}[$resi]} ) { 
        next if not defined $data->{$chain}[$resi]{$atom};
        my $a = $data->{$chain}[$resi]{$atom};
        my @alts = ('_');
        if ( defined $alts and defined $alts->{$chain} and defined $alts->{$chain}[$resi] ) 
            { @alts = @{$alts->{$chain}[$resi]} }
        #else  { carp ("Assign Default Altconfs !! Missing alt list for Chain $chain, resi $resi, atom $atom\n"); }
        my $main = $alts[0];
        if ( $main eq '' or $main eq ' ' ) { $main = '_' }
        if ( not exists $a->{'_'} ) { $a->{'_'} = $a->{$main} }
        if ( not exists $a->{'*'} ) { $a->{'*'} = $a->{$main} }
      }
    }
  }
}



###############################################################################
#                        Collapse_Atom_Properties
#------------------------------------------------------------------------------
# TODO: Document this crazy thing. 
###############################################################################

sub Collapse_Atom_Properties { 

  my $data = $_[0];
  my $resn = $_[1];
  my $alts = $_[2];
  my $prop = $_[3];
  my $func = $_[4];
  my $type = $_[5] || 'PROTEIN';

  my @groups = ();
  if    ( $type eq 'PROTEIN' ) { @groups = ( @Utility::ATOM_GROUPS, @Utility::ANGLE_GROUPS ) }
  elsif ( $type eq 'WATER'   ) { @groups = ( @Utility::SOLVENT_GROUPS ) }
  elsif ( $type eq 'SOLVENT' ) { @groups = ( @Utility::SOLVENT_GROUPS ) }
  elsif ( $type eq 'LIGAND'  ) { @groups = ( @Utility::SOLVENT_GROUPS ) }
  #my @groups = ( @Utility::ATOM_GROUPS, @Utility::ATOM_SETS, @Utility::ATOM_SUPER );

  my %functions = (
    fst => 1, 
    avg => 1,
    sum => 1,
    min => 1,
    max => 1,
    nan => 1, 
  );

  #---- Test for known functions
  if ( not exists $functions{$func} ) { 
    confess "Phenix::Collapse_Atom_Properties: Unknown function '$func'."
          . "Valid functions are " . (join ',', keys %functions) . "\n";
  }

  my @chains = keys %{$data->{$prop}};
  foreach my $chain (@chains) { #---------------------------------- Loop Chains

    my $atm = $data->{$prop}{$chain};
    foreach my $resi ( 0..$#{$data->{$prop}{$chain}} ) { #------- Loop Residues
      #my $resi = $resi->{$chain}[$resi];
      my $resn = $resn->{$chain}[$resi];
      next if not $resn; # skip missing residues

      #---- Get amino acid group definitions
      my $AA = $Utility::KNOWN_RESIDUES{$resn} || {};
      my @all_atoms = keys %{$atm->[$resi]};
      my @all_sets  = @all_atoms;
      #if ( not exists $AA->{aa_atoms} ) { @{$AA->{aa_atoms}} = [keys %{$atm->[$resi]}] };

      my @alts = @{$alts->{$chain}[$resi]};
      my $alt_main = $alts[0];
      my $alt_back = $alts[1];

      #print "Alts: " . (join ',', @alts) . "\n";


      foreach my $alt ( @alts ) { #---------------------------------- Loop Alts
        #print "  Alt '$alt'\n";

        foreach my $group ( @groups ) { #-------------------------- Loop Groups
          #print "    Group '$group'\n";
          my @atoms = Utility::Get_Atom_Group_Components($group, $resn, \@all_atoms);

          #print "$chain\t$resi\t$resn\t$alt\t$group\t" 
          #    . (join ',', @atoms) . "\t"
          #    . (join ',', @all_atoms) . "\n"
          #    ;
          #exit;
          next if scalar @atoms == 0;
          my @has_atoms = Utility::Intersection(\@atoms, \@all_sets);
          next if scalar @has_atoms == 0;
          push @all_sets, $group;


          my @values    = ();
          my $collapsed = 0;
          my $trash     = 0;

          # Collect all atoms in each group
          foreach my $a ( @atoms ) { #----------------------------- Loop Atoms
            #print "      Atom $a\n";
            next if $a=~/^H/ or $a=~/^\d/; #skip hydrogens
            my $value;
            
            #---- Skip missing atoms
            if ( not defined $atm->[$resi]{$a} ) {
              #print "Missing atom: $chain $resi $resn '$alt' $group $a.\n";
              next;
            }

            #==== Assign default alt to the primary alt if none is assigned
            # Primary alts are on a per atom basis, so might have to search
            # some alts (usually no more than ' ' and 'A') to find the highest
            # occupancy.  So just go through the list until we find one
            
           
             
            # Test primary Alt
            #if ( not exists $atm->[$resi]{$a}{_} ) { 
            #  $atm->[$resi]{$a}{_} = $atm->[$resi]{$a}{$alt_main};
            #  #if ( ($default eq ' ') or ($default eq '') or ($default eq $alts[0]) ) { 
            #  #  $atm->[$resi]{$a}{_} = $atm->[$resi]{$a}{$default};
            #  #}
            #}
            
            # Test backup Alt
            #if ( not defined $atm->[$resi]{$a}{_} ) { 
            #  $atm->[$resi]{$a}{_} = $atm->[$resi]{$a}{$alt_back};
            #}
            
            # Test every alt we can find
            if ( not defined $atm->[$resi]{$a}{_} ) { 
              foreach my $key ( @alts ) { 
              #foreach my $key ( keys %{$atm->[$resi]{$a}} ) { 
                if ( defined $atm->[$resi]{$a}{$key} ) { 
                  $alt = $key;
                  $atm->[$resi]{$a}{_} = $atm->[$resi]{$a}{$key};
                  last;
                }
              }
            }

            #---- If it's still undefined, we have a missing atom
            if ( not defined $atm->[$resi]{$a}{_} ) { 
              confess "There is not data for this atom: Prop $prop  \n"
                  . "\tChain $chain Resi $resi $resn  Alt '$alt' $group $a.\n"
                  . "\tAlts Attempted: " . (join ',', @alts) . "\n"
                  . Dumper($atm->[$resi])
                  ;
            }


            
            #---- Get data from the alt, if it exists, else take from main conf
            if ( exists $atm->[$resi]{$a}{$alt} ) { $value = $atm->[$resi]{$a}{$alt} } 
            else                                  { $value = $atm->[$resi]{$a}{_}    } 
           
            #---- Check and save value 
            if ( not defined $value ) { 
              confess "Property value isn't defined for prop $prop.  Chain $chain $resi $resn '$alt' $group $a.\n"
                  . "Alts Attempted: " . (join ',', @alts) . "\n"
                  . Dumper($atm->[$resi])
                  ;
            }
            push @values, $value;
            #print "Property value: $resi $resi $resn '$alt' $group $a $value.\n"
          }
          
          #---- Perform the collapse function on each group
          if    ( $func eq 'avg' ) { ($collapsed, $trash) = Utility::Mean_and_Stddev(\@values) }
          elsif ( $func eq 'sum' ) { ($collapsed, $trash) = Utility::Sum(\@values) || 0 }
          elsif ( $func eq 'min' ) { ($collapsed, $trash) = Utility::Min(\@values) || 0 }
          elsif ( $func eq 'max' ) { ($collapsed, $trash) = Utility::Max(\@values) || 0 }
          elsif ( $func eq 'nan' ) { ($collapsed, $trash) = (0, 0) }
          elsif ( $func eq 'fst' ) { ($collapsed, $trash) = $atm->[$resi]{$group}{$alt_main} }
          $atm->[$resi]{$group}{$alt} = $collapsed;
          #$atm->[$resi]{$group}{$alt} = sprintf "%.3f", $collapsed;
          
          
          #---- Assign the value for the main alt
          if ( ($alt eq $alt_main) or (not defined $atm->[$resi]{$group}{_})) { 
          #if ( $alt eq ' ' or $alt eq 'A' ) { 
            $atm->[$resi]{$group}{_} = $atm->[$resi]{$group}{$alt};
          }
        } #/Groups
      } #/Alts

      #=========== Use Angle Group XX for Angle Properties
      if ( not exists $atm->[$resi]{AA} ) {
        if ( exists $atm->[$resi]{XX} ) { 
          #print "Subbing XX for AA in $chain, $resi, $resn\n";
          $atm->[$resi]{AA} = $atm->[$resi]{XX};
        } else {  
          #print "No residue summary for $chain, $resi, $resn\n";
        }
      }

      #=========== Compute Summary Altconf * from All Other Altconfs
      @all_atoms = keys %{$atm->[$resi]}; # Now has all groups as well
      foreach my $atom ( @all_atoms ) { 
        my @values = ();
        
        #print "Summarizing * for Prop: $prop   Res: $resi $resn $type  Atom: $atom\n";
        #---- If no alts, freak out
        if ( scalar @alts == 0 ) { 
          confess "No altconfs * for Prop: $prop   Res: $resi $resn $type  Atom: $atom\n";
        }

        #---- If only one alt, just use that value
        my @these_alts = keys %{$atm->[$resi]{$atom}};
        if ( scalar @alts == 1 ) { 
          if ( defined $atm->[$resi]{$atom}{$alt_main} ) { 
            $atm->[$resi]{$atom}{'*'} = $atm->[$resi]{$atom}{$alt_main};
          } elsif ( defined $atm->[$resi]{$atom}{_} ) { 
            $atm->[$resi]{$atom}{'*'} = $atm->[$resi]{$atom}{_};
          } elsif ( scalar @these_alts > 0 ) { 
            foreach my $key ( @these_alts ) { 
              if ( defined $atm->[$resi]{$atom}{$key} ) { 
                $atm->[$resi]{$atom}{'*'} = $atm->[$resi]{$atom}{$key};
              }
            }
          } else { 
            confess "Couldn't find a single alt * for Prop: $prop   Res: $resi $resn $type  Atom: $atom\n";
          }
          next;
        }

        #---- Retrieve the values for each altconf of an atom/group
        foreach my $alt ( @alts ) { 
          my $value = $atm->[$resi]{$atom}{$alt};
          push @values, $value if defined $value;
        }

        #---- Perform the collapse function on each group
        my ( $collapsed, $trash );
        if    ( $func eq 'avg' ) { ($collapsed, $trash) = Utility::Mean_and_Stddev(\@values) }
        elsif ( $func eq 'sum' ) { ($collapsed, $trash) = Utility::Sum(\@values) || 0 }
        elsif ( $func eq 'min' ) { ($collapsed, $trash) = Utility::Min(\@values) || 0 }
        elsif ( $func eq 'max' ) { ($collapsed, $trash) = Utility::Max(\@values) || 0 }
        elsif ( $func eq 'nan' ) { ($collapsed, $trash) = (0, 0) }
        elsif ( $func eq 'fst' ) { ($collapsed, $trash) = $atm->[$resi]{$atom}{$alt_main} }
        $atm->[$resi]{$atom}{'*'} = $collapsed;
      }
    } #/Residues
  } #/Chains
}



###############################################################################
#                         Output Atom Analysis Matrix
#------------------------------------------------------------------------------
# Outputs a complete set of analysis for a given structure/level/residue type
# combination. Detects all appropriate properties to include in the file.
###############################################################################

sub Output_Analysis_Matrix {

  my $CFG       = $_[0];
  my $level     = $_[1];
  my $file      = $_[2];
  my $group     = $_[3];
  my $out_type  = $_[4] || 'P';
  my $flags     = $_[5] || {};

  # Map type to the short name.
  if (length($out_type) > 1) { $out_type = $Utility::RESIDUE_TYPES{$out_type} }
  my @props = Get_Appropriate_Properties($CFG->{PROPS}, $level, $out_type);
  #my @props = Get_Appropriate_Properties($CFG->{PROPS_OUT}, $level, $out_type);


  #------------------------------------ Output Header
  my $output = "#" 
             . ( join "", map { "$_\t" } @{$CFG->{ANNO_LIST}} )
             . "CHAIN\tRESI\tRESN\tATOM\tALT";
  foreach my $prop ( @props ) { 
    $output .= "\t$prop->{anal}_$prop->{prop}";
  }
  $output .= "\n";

  foreach my $chain ( @{$CFG->{PDB}{CHAINS}} ) { #----------------- Loop Chains
    foreach my $resi ( 0..$#{$CFG->{PDB}{RESI}{$chain}} ) { #---- Loop Residues
      
      my $resn      = $CFG->{PDB}{RESN}{$chain}[$resi];
      my $type      = $CFG->{PDB}{TYPE}{$chain}[$resi];
      my $X         = $CFG->{PDB}{X}{$chain}[$resi];
      my @all_atoms = ();
      my @set       = ();

      # Fill in missing residues, for my IPAs that are non-consequtively numbered
      if ( $flags->{fill_missing_residues} ) { 
        if ( not defined $resn ) { 
          $resn       = 'XXX';
          $type       = $out_type;
          @all_atoms  = ('X');
          @set        = ('X');
        } 
      } else { 
        next if not defined $resn;
        next if not defined $out_type;
        next if $out_type ne '*' and $out_type ne (uc $type);
        @all_atoms = keys %$X;
        @set = Utility::Get_Atom_Set_Elements($group, $resn, \@all_atoms);
      }


      #print "$chain\t$resi\t$resn\t$type\t$group\n" . Dumper(\@all_atoms) . Dumper(\@set);

      my @alts = @{$CFG->{PDB}{ALTS}{$chain}[$resi] || ['_']};
      if ( $level eq 'RESIDUE' ) { @alts = ("_") }

      foreach my $alt ( @alts ) { #---------------------------------- Loop Alts
        my $a = $alt; $a = '_' if $a eq ' ';

        foreach my $atom ( @set ) { #------------------------------- Loop Atoms

          #print "$chain $resi $resn $type $alt $atom\n";
          next if Utility::Is_Hydrogen_Atom($atom);
          if ( not defined $X->{$atom} ) { 
            next;
            confess "Unknown atom/set/group type '$group' - '$atom'.\n"
                  . "Available sets:   " . (join ',', @Utility::ATOM_SETS  ) . "\n"
                  . "Available groups: " . (join ',', @Utility::ATOM_GROUPS) . "\n"
                  #. Dumper(\%Utility::ATOM_GROUPS)
                  . Dumper($X)
                  ;
          }
          
          if ( not exists $X->{$atom}{$a} ) { 
            print "Missing X-Coordinate for alt '$a' in $chain $resi $resn $alt $atom\n"
                . Dumper($X->{$atom});
            next;
          }

          # Output Crystal Structure Annotations
          $output .= (join "", map { $CFG->{ANNOTATIONS}{$_} . "\t" } @{$CFG->{ANNO_LIST}} );
          
          # Output Atom Spec
          $output .= sprintf "%s\t%d\t%s\t%s\t%s",
                              $chain,
                              $resi,
                              $resn,
                              $atom,
                              $a,
                              ;
         
          foreach my $analprop ( @props ) { 
            my $anal = $analprop->{anal};
            my $prop = $analprop->{prop};
            my $value = $CFG->{$anal}{$prop}{$chain}[$resi]{$atom}{$a};
            if ( not defined $value ) {
              $value = $CFG->{$anal}{$prop}{$chain}[$resi]{$atom}{_};
            }  
            if ( not defined $value ) {  
              if ( defined $META_PROPERTIES{$anal}{$prop}{D} ) { 
                $value = $META_PROPERTIES{$anal}{$prop}{D}; 
              #} elsif ( defined $default ) { 
              #  $value = $default;
              } else { 
                confess "Cannot find a value nor a default value to use for:\n"
                      . " $CFG->{TITLE} $anal $prop $chain $resi $resn $atom $a\n"
                      . Dumper($CFG->{$anal}{$prop}{$chain}[$resi]);
              }
            }
            #print "\t$anal\t$prop\t$value\n";
            #$output .= "\t$anal\t$prop\t$value";
            if ( not exists $META_PROPERTIES{$anal}{$prop}{F} ) { 
              print "Format $anal $prop '$META_PROPERTIES{$anal}{$prop}{F}'\n";
            }
            $output .= sprintf "\t$META_PROPERTIES{$anal}{$prop}{F}", $value;
            #$output .= sprintf "\t$CFG->{FORMATS}{$anal}{$prop}", $value; # Why am I not using global formats?  Guess we'll find out the hard way.
          }
          $output .= "\n";
        }
      }
    }
  }


  if ( not $file ) { 
    print $output;
  } else { 
    Utility::Write_File($file, $output, 1);
  }
}



###############################################################################
#                          Output Analysis Regroup
#------------------------------------------------------------------------------
# Output only information for primary alternate conformations for selected
# analysis/properties for each pdb. Structure information is taken from a 
# passed structure config, typically the first pdb.
###############################################################################

sub Output_Analysis_Regroup { 

  my $CON         = $_[0] or confess "Need consensus configuration.";
  #my $CFG         = $_[0] or confess "Need CFG.";
  #my $data        = $_[1] or confess "Need data.";
  #my $pdbs        = $_[2] or confess "Need bases list.";
  my $analprop    = $_[1] or confess "Need analprop hash";
  my $output_file = $_[2] or confess "Need output file name";
  my $flags       = $_[3] || {};

  my $verbose = $flags->{verbose} || 0;
  #print Dumper($flags);


  ############################################################## INITIALIZATION

  #===== Retrieve or Set Defaults for Options Flags
  my $use_altconf       = $flags->{use_altconf      } || '_';
  my $structure_level   = $flags->{structure_level  } || 'ATOM';
  my $structure_group   = $flags->{structure_group  } || 'RES';
  my $file_format       = $flags->{file_format      } || 'ALL';
  my $filter_type       = $flags->{filter_type      } || 'P';
  my $filter_resn       = $flags->{filter_resn      } || '';
  my $filter_chain      = $flags->{filter_chain     } || '';
  if ( length($filter_type) > 1 ) { $filter_type = $Utility::RESIDUE_TYPES{$filter_type}}

  #===== Verify Settings
  my $anal = $analprop->{anal};
  my $prop = $analprop->{prop};
  Verify_Appropriate_Properties ( $CON->{PROPS}, $anal, $prop, $structure_level, $filter_type );

  #===== Validate and Aggregate Data for Easier Lookup
  my %data = ();
  my $data_format = $META_PROPERTIES{$anal}{$prop}{F};
  foreach my $pdb ( @{$CON->{TITLES}}) { 
  #foreach my $pdb ( @$pdbs ) { 
    if ( not exists $CON->{DATA}{$pdb} ) { 
      confess "PDB $pdb does not have regroup level entries for $structure_level.\n";
    } 
    if ( not $flags->{fill_missing_pdbs} ) { 
      if ( not exists $CON->{DATA}{$pdb}{$anal} ) { 
        confess "PDB $pdb does not have regroup analysis entries for $structure_level $anal.\n";
      }
      if ( not exists $CON->{DATA}{$pdb}{$anal}{$prop} ) { 
        confess "PDB $pdb does not have regroup property entries for $structure_level $anal $prop.\n";
      }
    }
    $data{$pdb} = $CON->{DATA}{$pdb}{$anal}{$prop};
    #$data{$pdb} = $data->{$pdb}{$anal}{$prop};
  }

   
   
  ################################################################ OUTPUT DATA
   
  #===== Output Headers
  my $output = "#";
  if    ( $file_format eq 'RESI'    ) { $output .= "RESI" }
  elsif ( $file_format eq 'COMPACT' ) { $output .= join "_",  qw(RESI RESN CHAIN ATOM) }
  elsif ( $file_format eq 'PLOT'    ) { $output .= "RESI" }
  elsif ( $file_format eq 'ALL'     ) { $output .= join "\t", qw(RESI RESN CHAIN ATOM) }
  else { confess "Unknown row label $file_format.\n" }
  
  foreach my $pdb ( @{$CON->{TITLES}} ) { $output .= "\t$pdb" }
  #foreach my $pdb ( @$pdbs ) { $output .= "\t$pdb" }
  $output .= "\n";
  
  my @chains = @{$CON->{PDB}{CHAINS}}; 
  foreach my $chain ( @chains ) { #===================================== CHAIN
    next if $filter_chain and $chain ne $filter_chain;
    print "  Chain $chain\n" if $verbose >= 2;
    foreach my $resi ( 0..$#{$CON->{PDB}{TYPE}{$chain}} ) { #============ RESI

      my $resn      = $CON->{PDB}{RESN }{$chain}[$resi];
      my $type      = $CON->{PDB}{TYPE }{$chain}[$resi];
      my $atoms     = $CON->{PDB}{ATOMS}{$chain}[$resi];
      #my $X         = $CFG->{PDB}{X}{$chain}[$resi];
      #my @all_atoms = ();
      my @set       = ();
      my @atoms     = ();

      # Fill in missing residues, for my IPAs that are non-consequtively numbered
      if ( $flags->{fill_missing_residues} ) { 
        if ( not defined $resn ) { 
          $resn       = 'XXX';
          $type       = $filter_type;
          $atoms      = ['X'];
          #@all_atoms  = ('X');
          @set        = ('X');
        } 
      } else { 
        next if not defined $resn;
        next if not defined $type;
        next if $filter_type ne (uc $type);
        next if $filter_resn and ($filter_resn ne $resn);
        #@all_atoms = keys %$X;
        @set = Utility::Get_Atom_Set_Elements($structure_group, $resn, $atoms);
        #@set = Utility::Get_Atom_Set_Elements($structure_group, $resn, \@all_atoms);
      }
      print "    RES $resi $resn $type\n" if $verbose>=2;

      my $set_size = scalar @set;
      foreach my $s ( 0..$#set ) { #====================================== SET
        my $atom = $set[$s];
        next if Utility::Is_Hydrogen_Atom($atom);
        print "      Atom $atom\n" if $verbose >= 2;
        #if ( not grep { $_ eq $atom } @atoms ) { 
        #if ( not exists $X->{$atom} ) { 
        #  next;
        #  confess "Unknown atom/set/group type $anal $prop $chain $resi $resn $atom $structure_group.\n"
        #        . "Available sets:   " . (join ',', @Utility::ATOM_SETS  ) . "\n"
        #        . "Available groups: " . (join ',', @Utility::ATOM_GROUPS) . "\n"
        #        #. Dumper(\%Utility::ATOM_GROUPS)
        #        . Dumper($CON)
        #        ;
        #}

        #if ( not exists $X->{$atom}{$use_altconf} ) { 
        #  confess "Atom is missing altconf '$use_altconf': $anal $prop $chain $resi $resn $atom.\n";
        #}
        #if ( not defined $X->{$atom}{$use_altconf} ) { 
        #  confess "Atom not defined for altconf '$use_altconf': $anal $prop $chain $resi $resn $atom.\n";
        #}

        #---- Output Residue/Atom Identifiers
        if    ( $file_format eq 'RESI'       ) { $output .= $resi }
        elsif ( $file_format eq 'COMPACT'    ) { $output .= join "_",  $resi, $resn, $chain, $atom }
        elsif ( $file_format eq 'ALL'        ) { $output .= join "\t", $resi, $resn, $chain, $atom }
        elsif ( $file_format eq 'PLOT'       ) { 
          if ( $structure_level eq 'RESIDUE' ) { $output .= sprintf "%.0f", $resi + ($s/$set_size) }
          else                                 { $output .= sprintf "%.2f", $resi + ($s/$set_size) }
        }
        
        #print "        HERE: $chain, $resi, $resn, $atom, $use_altconf\n";

        foreach my $pdb ( @{$CON->{TITLES}} ) { #========================== PDB
          
          #--- Get the value
          my $value;
          $value = $data{$pdb}{$chain}[$resi]{$atom}{$use_altconf};

          #--- Try to find a default value if it's missing
          if ( not defined $value ) {  
            #if ( defined $data{$pdb} ) { $value = 0 }  # Think this is just for faking it.
            if ( defined $META_PROPERTIES{$anal}{$prop}{D} ) { 
              $value = $META_PROPERTIES{$anal}{$prop}{D};
            } else { 
              confess "PDB $pdb is missing data for $anal $prop $chain Residue: $resi $resn  Atom: $atom   Alt '$use_altconf'.\n"
                    . Dumper($data{$pdb}{$chain}[$resi])
                    ;
            }
          }

          $output .= sprintf "\t$data_format", $value;
        }
        $output .= "\n";
      }
    }
  }

  
  #---- Write File or Print to STDOUT
  if ( not $output_file ) { 
    print $output;
  } else { 
    Utility::Write_File($output_file, $output, 1);
  }
}



        














###############################################################################
###############################################################################
#                   ANALYSIS SETS AND BATCH OPERATIONS
###############################################################################
###############################################################################


# TODO: 
# Convert prepare_pdbs.pl to initializing with CFGs
# Write batch analysis blocks for 
#   structure analyses: 
#     PISA
#     HOLES
#     ANGLES
#     ALT CONFS
#   density analyses:
#     RINGER
#
#
# Convert rest of phenix subroutines to using CFGs
# Rename Phenix.pm to something sensible
# Batch MTZ cleanup
#   CIF to MTZ
#   FOBS
#   RFree
#   Cut anisotropic data
#
#
#
#
# Usage Cases:
#
# prepare_pdbs:
#   download pdb/mtz
#   clean shit up
#   get meta data: dump, xtriage, meta
#   analyze structure
#   analyze density
#
#
##############################################################################





##############################################################################
##############################################################################
#                        PDB PARSING AND OTHER FUNCTIONS
##############################################################################
##############################################################################



#################################################################### Parse PDB
# NOTE: Assuming for now that there is only one model of interest per PDB.
#       Returning only information for the first model.

sub Parse_PDB { 

  my $CFG = $_[0];
  my $pdb;
  my $flags = $_[1];

  if ( Is_Structure_Config($CFG) ) { 
    $pdb   = $CFG->{PATHS}{PDB};
    $flags = $CFG->{FLAGS};
  } else { 
    $pdb   = $CFG;
  }

  my $flag_do_meta = 1;
  my $flag_verbose = $flags->{pdb_verbose} || 0;
  if ( $flags->{pdb_data} ) { 
    $flag_do_meta = 0;
  }

  my %PDB = ();

  $ParsePDB::NoANISIGDefault = 1;
  my $PDB = ParsePDB->new(FileName => $pdb);

  #----------- Validation
  $PDB{MODELS} = [$PDB->IdentifyModels()];
  if ( scalar @{$PDB{MODELS}} > 1 ) { 
    confess "I haven't decided how to deal with multiple models yet.\n";
  }

  if ( not $PDB->ChainLabelsValid() ) { 
    #print "Invalid chain lables for pdb '$pdb'.\n";
    #exit;
  } 

  
  #----------- Precompute dihedral angles for phi/psi
  $PDB{CHAIN_IDS} = [$PDB->IdentifyChainLabels()];
  $PDB{CHAINS}    = [sort {$a cmp $b} Utility::Unique($PDB->IdentifyChainLabels())];
  #$PDB->GetAngles();
  $PDB->GetElement;

  #----------- Headers and Footers
  $PDB{HEADER_LINES} = [$PDB->GetHeader()];
  $PDB{FOOTER_LINES} = [$PDB->GetFooter()];

  #----------- Metadata
  my $META;
  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{PDB     } = \%PDB;
    $CFG->{ParsePDB} =  $PDB;
    $META = Parse_PDB_Meta($CFG) if $flag_do_meta;
  } else { 
    $META = Parse_PDB_Meta($PDB{HEADER_LINES}) if $flag_do_meta;
  }
  
  # Mark protein residues from the sequence (handles 1DY5 67 IAS)
  if ( exists ( $META->{dbref} ) ) { 
    foreach my $dbref ( @{$META->{dbref}} ) { 
      foreach my $i ( $dbref->{seq_begin}..$dbref->{seq_end} ) { 
        $PDB{TYPE}{$dbref->{chain}}[$i] = 'P';
      }
    }
  }


  #----------- Get Residue Sequence and Numbering
  # TODO: May need to add AtomLocations=>First|All|None|'A'|...
  my $residue_index = 0;
  my $atom_index = 0;
    


  foreach my $c ( 0..$#{$PDB{CHAIN_IDS}} ) { 
    my $chain = $PDB{CHAIN_IDS}[$c];
    print "Extracting residues for chain $c ($chain).\n" if $flag_verbose >= 1;
    my @res = $PDB->Get( Chain=>$c, ResidueIndex=>1 );
    foreach my $res ( @res ) { 
      my $resi = $res->{ResidueNumber};
      my $resn = $res->{ResidueLabel};
      #$PDB{CHAIN  }{$chain}[$resi] = $res;
      $PDB{RESI   }{$chain}[$resi] = $resi;
      $PDB{RESN   }{$chain}[$resi] = $resn;
      $PDB{RESIDUE}{$chain}[$resi] = $res;
      $PDB{ATOMS  }{$chain}[$resi] = Utility::Unique( map { $_->{AtomType} } @{$res->{Atoms}} );
    }

    
    #--- Break out protein vs. other residues
    $PDB{TYPE}{$chain} = [];
    Utility::Assign_Residue_Types($PDB{RESN}{$chain}, $PDB{TYPE}{$chain});
    foreach my $resi ( 0..$#{$PDB{TYPE}{$chain}} ) { 
      my $type = $PDB{TYPE}{$chain}[$resi];
      next if not $type;
      push @{$PDB{$Utility::RESIDUE_TYPES{$type}}{$chain}}, $resi;
    }

    #=== Extract residue properties
    

    #foreach my $resi ( @{$PDB{PROTEIN}{$chain}} ) { 
    foreach my $resi ( 0..$#{$PDB{TYPE}{$chain}} ) { 
      my %alts = ();
      next if not defined $PDB{RESIDUE}{$chain}[$resi];
      my $res  = $PDB{RESIDUE}{$chain}[$resi];
      my $resn = $PDB{RESN}{$chain}[$resi];
      my $type = $PDB{TYPE}{$chain}[$resi];
      #my $resi = $PDB{RESI}{$chain}[$resi];
      my $incd = $res->{InsResidue};
      if ( ($incd ne '') and ($incd ne ' ') ) { 
        confess ("This protein has insertion codes (Chain $chain, Res $resi $resn, ICode '$incd').\n");
      }

      $PDB{TYPE_COUNTS}{$type}{$resn}++;

      #=========== Create Residue Object
      my $obj_residue = {
        chain => $chain,
        resi  => $resi,
        resn  => $resn,
        residue_index => $residue_index,
        label => "$chain\_$resi",
        atoms => [],
      };
      push @{$PDB{RES_MAP}}, $obj_residue; # Residue Map

      my $atom_subindex = 0;
      foreach my $a ( 0..$#{$res->{Atoms}} ) { 
        my $atom = $res->{Atoms}[$a];
        my $name = $atom->{AtomType};
        my $alt  = $atom->{AltLoc  };
        next if Utility::Is_Hydrogen_Atom($name);

        print "Chain $chain\tRESI $resi\tRESN $resn\tAtom $name\tAlt '$alt'\t$atom->{Temp}\t$atom->{Occupancy}\t$atom->{x}\t$atom->{y}\t$atom->{z}\n"
          if $flag_verbose >= 5;

       
        my $old_occ = $PDB{Q}{$chain}[$resi]{$name}{_} || 0;
        if ( $flags->{no_alts} ) { 
          $alt = "_";
          next if $atom->{Occupancy} < $old_occ;
        } elsif ( $flags->{min_occ} ) { 
          if ( $atom->{Occupancy} < $flags->{min_occ} ) { 
            next;
          }
        }
        
        $alts{$alt} = $atom->{Occupancy};
        $PDB{B}{$chain}[$resi]{$name}{$alt} = $atom->{Temp};
        $PDB{Q}{$chain}[$resi]{$name}{$alt} = $atom->{Occupancy};
        $PDB{X}{$chain}[$resi]{$name}{$alt} = $atom->{x};
        $PDB{Y}{$chain}[$resi]{$name}{$alt} = $atom->{y};
        $PDB{Z}{$chain}[$resi]{$name}{$alt} = $atom->{z};
        
        $PDB{COORDS}{$chain}[$resi]{$name}{$alt} = [$atom->{x}, $atom->{y}, $atom->{z}];


        #========= Create Atom Object        
        my $obj_atom = {
          chain   => $chain,
          resi    => $resi,
          resn    => $resn,
          atom    => $name,
          alt     => $alt,

          x             => $atom->{x},
          y             => $atom->{y},
          z             => $atom->{z},
          b             => $atom->{Temp},
          q             => $atom->{Occupancy},
          element       => $atom->{Element},
          number        => $atom->{AtomNumber},

          residue_index => $residue_index,
          atom_index    => $atom_index,
          atom_subindex => $atom_subindex,
          label         => "$chain\_$resi\_$atom\_$alt",
          residue       => $obj_residue,
        };
        push @{$PDB{ATOM_MAP}}, $obj_atom;                   # Atom Map
        push @{$obj_residue->{atoms}}, $obj_atom;               # Residue Map
        push @{$PDB{ATOM_MATRIX}[$residue_index]}, $obj_atom;  # Atom Matrix
        $atom_index++;
        $atom_subindex++;
      } #/END ATOM
      $residue_index++;

      #=== Count Alternate Conformations and Sort By Occupancy, Rather than Alphabetically
      my @alts = sort { - ($alts{$a} <=> $alts{$b}) } keys %alts;
      $PDB{ALTS}{$chain}[$resi] = \@alts;
      $PDB{ALT_COUNTS}{$chain}[$resi]{AA}{'_'} = scalar @alts;
      $PDB{ALT_COUNTS}{$chain}[$resi]{AA}{'*'} = scalar @alts;

      foreach my $name ( keys %{$PDB{COORDS}{$chain}[$resi]} ) { 
        $PDB{COORDS}{$chain}[$resi]{$name}{_} = $PDB{COORDS}{$chain}[$resi]{$name}{$alts[0]};
      }
     
      #=== Assign Default Alt 
      #foreach my $resi ( @{$PDB{PROTEIN}{$chain}} ) { 
      #  foreach my $prop ( qw(B Q X Y Z COORDS) ) { 
      #    $PDB{$prop}{$chain}[$resi]{$name}{_} = $PDB{$prop}{$chain}[$resi]{$name}{$alt};
      #  }
      #}

      #=== Validate Occupancies and Names for Alts for Protein Residues (not ligands, waters)
      if ( $type eq 'P' ) { 
        if ( scalar @alts > 1 ) { 

          my $last_occ;
          my $last_alt;
          if ( not exists $alts{A} ) { 
            warn "Residue $chain $resi $resn is missing alt conf A.\n"
              if $flag_verbose >= 3;
          }
          foreach my $alt ( sort @alts ) { 
            if ( $alt ne uc $alt ) { 
              warn "Residue: $chain $resi $resn alt $alt is not a capital letter.\n"
                if $flag_verbose >= 3;
            }

            if ( $last_alt and ( not $alt eq (chr(ord($last_alt)+1)))) { 
              warn "Residue $chain $resi $resn: missing alt between $last_alt and $alt\n"
                . "\t" . (join " ", @alts) . "\n"
                if $flag_verbose >= 3;
            }
            if ( $last_occ and ( $alts{$alt} > $last_occ ) ) { 
              warn "Residue $chain $resi $resn: Alt $alt > $last_alt ($alts{$alt} vs. $last_occ.)\n"
                if $flag_verbose >= 3;
            } 

            $last_alt = $alt;
            $last_occ = $alts{$alt};
          }

        } elsif ( $alts[0] !~ /^\s*$/ ) { 
          warn "Residue has only one alt, but it is not ' ': $chain $resi $resn\n"
            if $flag_verbose >= 3;
        }

        my $occ_sum = Utility::Sum(values %alts);
        if ( abs($occ_sum - 1.0) > $OCCUPANCY_DIFFERENCE_THRESHOLD ) { 
          warn "Residue $chain $resi $resn occupancy does not sum to 1.0 ($occ_sum)\n"
             . Dumper(%alts)
             if $flag_verbose >= 3;
        }
      } 


      #---- Little Output 
      my $atom = 'CA';
      if ( not grep { /$atom/ } @{$PDB{ATOMS}{$chain}[$resi]} ) { 
        $atom = $PDB{ATOMS}{$chain}[$resi][0];
      }
      print "$PDB{TYPE}{$chain}[$resi]\t"
          . "$PDB{RESI}{$chain}[$resi]\t"
          . "$PDB{RESN}{$chain}[$resi]\t"
          . "$PDB{B   }{$chain}[$resi]{$atom}{$alts[0]}\t"
          . "$PDB{Q   }{$chain}[$resi]{$atom}{$alts[0]}\t"
          . "$PDB{X   }{$chain}[$resi]{$atom}{$alts[0]}\t"
          . "$PDB{Y   }{$chain}[$resi]{$atom}{$alts[0]}\t"
          . "$PDB{Z   }{$chain}[$resi]{$atom}{$alts[0]}\t"
          . "\n"
          if $flags->{pdb_verbose} >= 5
          ;

    }

  }
      


  ######################################### Count Residues and Types
  

#  my %TYPES = ();
#  foreach my $chain ( keys %{$PDB{type}} {
#    $TYPES{$RESN} = 0;
#    freach my $resi ( @{$PDB[$resi]{$chain]}} ) { 
#      my $type = $PDB{TYPE}{$chain}[$resi];
#      my $resn = $PDB{RESN}{$chakin}[$esi];
#      $TTYPES{$types}{$resn} ++;
#    }
#  }


  
  my $meta_path = $CFG->{PATHS}{META};
  if ( defined $meta_path ) { 

    open META, ">>$meta_path" or die "Couldn't open pdb meta data file, '$meta_path'. $!";
    foreach my $type ( sort keys %{$PDB{TYPE_COUNTS}} ) { 
      foreach my $resn ( sort keys %{$PDB{TYPE_COUNTS}{$type}} ) { 
        print META "TYPECOUNT\t$type\t$resn\t$PDB{TYPE_COUNTS}{$type}{$resn}\n";
      }
    }
    close META or die "Couldn't close pdb meta data file, '$meta_path'. $!";

    #my $meta_path = $CFG->{PATHS}{META}
    #Utility::Write_File($meta_type_counts_path, Dumper(\%TYPE_COUNTS), { overwrite => 1 });
  }



  ######################################## Delete unwantee redsduetypes

  # Delete chains that don't actually contain the residue type
  foreach my $type ( @Utility::RESIDUE_TYPES ) { 
    foreach my $chain ( keys %{$PDB{$type}} ) { 
      next if not exists $PDB{$type}{$chain};
      my $count = scalar @{$PDB{$type}{$chain}};
      #print "$type\t$chain\t$count\t" . (join ',', @{$PDB{$type}{$chain}}) . "\n";
      if ( $count == 0 ) { 
        #print "  Deleting $type $chain\n";
        delete $PDB{$type}{$chain};
      }
    }
  }




  Collapse_Atom_Properties(\%PDB, $PDB{RESN}, $PDB{ALTS}, 'B', 'avg');
  Collapse_Atom_Properties(\%PDB, $PDB{RESN}, $PDB{ALTS}, 'Q', 'avg');
  Collapse_Atom_Properties(\%PDB, $PDB{RESN}, $PDB{ALTS}, 'X', 'avg');
  Collapse_Atom_Properties(\%PDB, $PDB{RESN}, $PDB{ALTS}, 'Y', 'avg');
  Collapse_Atom_Properties(\%PDB, $PDB{RESN}, $PDB{ALTS}, 'Z', 'avg');

  Assign_Default_Altconfs($PDB{COORDS}, $PDB{ALTS});

  if ( Is_Structure_Config($CFG) ) { 
    Register_Analyses(
      $CFG, 
      [qw(ATOM ANGLE RESIDUE)], 
      'PDB', 
      [ qw(B Q X Y Z) ],
      ['PROTEIN', 'WATER', 'LIGAND'],
    );
    Register_Analyses ( 
      $CFG,
      ['RESIDUE'],
      'PDB',
      ['ALT_COUNTS'],
      ['PROTEIN', 'WATER', 'LIGAND'],
    );
  }

  #==== Specialty Analyses
  Calculate_Dihedral_Angles($CFG);      # Dihedral Angles
  Detect_Distinct_Altconfs($CFG);       # Distinct Altconfs

  return $PDB, \%PDB;
}







############################################################### Parse PDB Meta
# Parses several pieces of metadata out of the pdb headers:
#   pdb         - name of the pdb file
#   title       - Description of the protein from the HEADER
#   date        - Date of deposition from the HEADER (also day, month, year)
#   exp_type    - XRAY, NMR, etc.
#   r_work
#   r_free
#   monomers
#
#   spacegroup      - name of the space group
#   edges           - array of unit cell edge lengths
#   edges_txt       ^ as text
#   angles          - 
#   angles_txt      ^ as text
#   protein_atoms   - Atom count of each type (protein, nucleic, etc)
#   nucleic_atoms   ^
#   het_atoms       ^
#   solvent_atoms   ^
#
#   chains          - chains listed in the PDB (from DBREFs)
#   dbref           - alignment of this structure to protein database sequences
#   sequence        - protein sequence for each chain
#   compound        - molecule, chain, synonym, and ec from COMPND
#
#   TODO: *'d ones are important
#   formula         * FORMUL data for each hetatm type
#   source          * animal SOURCE for the proteins
#   keywords        - KEYWDS
#   authors         - AUTHOR
#   journals        - JRNL
#   remark 4        - Compliance with PDB format versions
#   remark 200      - experimental details (diffraction)
#   remark 280      * crystal conditions
#   remark 900      - related entries (pdbs from the same paper, it looks like)
#   het and hetname - hetatm stuff?
#   helix           * 
#   sheet           * 
#   ssbond          - 
#   site            - potentially interesting active site info?
#   Anything else I already parse the hard way from pdb_repository
#-----------------------------------------------------------------------------
# Data Structure
# --------------
#
# PDB::HEADER_LINES   - Lines to parse
# 
#
#
#
#-----------------------------------------------------------------------------
# Remark types (http://www.wwpdb.org/documentation/format33/remarks1.html):
#         0 - Rerefinement of old data
#         1 - Publications
#         2 - Resolution
#         3 - Refinement programs and statistics
#         4 * PDB file format version
#         5 - Obsolete statement TODO: Detect this and crash for now
#   6 -  99 - Free text annotations
#       100 - Processing site (RCSB, PDBe, PDBj, BNL)
#
#       200 * Experimental details for X-Ray Diffraction
#       210 ^                          NMR
#       215 ^                          Solution NMR
#       217 ^                          Solid State NMR
#       230 ^                          Neutron Diffraction
#       240 ^                          Electron Crystallography
#       245 ^                          Electron Microscopy
#       247 ^                          Electron Microscopy Again?
#       250 ^                          Other
#       265 ^                          Solution Scattering
#
#       280 * Crystallographic Details 
#       285 - Free text description of CRYST1
#       290 - Crystallographic Symmetry
#
#       300 - Biological Function (Free Text)
#       350 - Coordinate transformations to generate biomolecule
#       375 - Atoms on symmetry boundary
#       400 - Compound information
#       450 - Biological source
#       465 - Missing residues
#       470 - Missing atoms
#       475 - Residues with zero occupancy
#       480 - Atoms with zero occupancy
#       500 - Geometry and sterochemistry
#       525 - Distant solvent atoms
#       600 - Heterogens
#       610 - HETATM missing atoms
#       615 - HETATM atoms with zero occupancy
#       620 - Metal coordination
#       630 - Inhibitor description
#       650 - Helical SS
#       700 - Sheet SS
#       800 - Description of SITE data
#       900 - Related entries (from same experiment set, it looks like)
#       999 - Sequence
#
#



sub Parse_PDB_Meta {

  my $CFG = $_[0];
  my ($pdb_path, $pdb_file);
  my @lines = ();
  my $flags = $_[1] || {};

  if ( Is_Structure_Config($CFG) ) { 
    $pdb_path = $CFG->{PATHS}{PDB};
    $pdb_file = $CFG->{PATHS}{PDB};
    if ( defined $CFG->{PDB}{HEADER_LINES} ) { 
      @lines = @{$CFG->{PDB}{HEADER_LINES}};
    } else { 
      @lines = Utility::Read_File($pdb_path);
    }
    $flags = $CFG->{FLAGS};
  } elsif ( ref $CFG eq 'ARRAY' ) { 
    @lines = @$CFG;
  } else { 
    $pdb_path = $CFG;
    @lines = Utility::Read_File($pdb_path);
  }

  $flags->{meta_verbose} = 0 if not exists $flags->{meta_verbose};
  $flags->{pdb_meta_output  } = 0 if not exists $flags->{pdb_meta_output  };
  
  my @props = qw(
    title
    date day month year
    exp_type resolution
    r_work r_free
    spacegroup monomers edges_txt angles_txt
    protein_atoms nucleic_atoms het_atoms solvent_atoms
    pdb_version pdb_version_date
  
    refinement_program
    dfxn_temperature dfxn_temperature_unit ph
    completeness_all completeness_hr
  );


  my %meta = ( 
    chains    => [],
    dbref     => [],
    sequence  => {},
    compounds => [],
    formulas  => [],
  );


  my ($line, $record, $record_info, $remark, $remark_info);
  my $last = { id => -1 };
  foreach my $line ( @lines ) { 
    chomp $line;
    #print "Line: $line\n";

    $record      = Utility::Trim(substr($line, 0, 6));
    $record_info = substr($line, 7);
    $remark      = Utility::Trim(substr($line, 7, 3)); # only meaningful for remark lines!
    $remark_info = substr($line, 11);

    #========= Header
    if ( $record eq 'HEADER' ) { 
      $line =~ /^HEADER\s+(.*)(\d\d)-(\w\w\w)-(\d\d)\s+(\d\w\w\w)/
        or confess("Couldn't parse $record line for pdb '$pdb_path':\n'$line'\n");
      my ($title, $day, $month, $year, $header_pdb) 
        = map { Utility::Trim($_) } ($1, $2, $3, $4, $5);
      $meta{'title'       } = $title;
      $meta{'date'        } = "$day-$month-$year";
      $meta{'day'         } = $day;
      $meta{'month'       } = $month;
      $meta{'year'        } = $year;
      $meta{'line_header' } = $line;
      next;
    }

    #========= Experiment Type
    if ( $line =~ /^EXPDTA\s+(.*)$/ ) { 
      my $exp = Utility::Trim($1);
      if    ( $exp =~ /NMR/     ) { $meta{'exp_type'} = 'NMR'  }
      elsif ( $exp =~ /X-RAY/   ) { $meta{'exp_type'} = 'XRAY' }
      elsif ( $exp =~ /NEUTRON/ ) { $meta{'exp_type'} = 'XRAY' }
      else                        { $meta{'exp_type'} = $exp   }
      next;
    }



    #========= All Remarks
    if ( $record eq 'REMARK' ) { 

      next if not Utility::Is_Number($remark);

      #======= Remark 2
      #------- Resolution (If XRay)
      if ( $remark == 2 ) {
        #----- From PDB
        if ( $line =~ /RESOLUTION\.\s+([\d\.])+ ANGSTROMS./ ) { 
          $meta{resolution} = $1;
        #----- From PHENIX
        } elsif ( $line =~ /RESOLUTION RANGE HIGH (ANGSTROMS) : ([\d\.]+)/ ) { 
          $meta{resolution} = $1;
        }
      }

      #======= Remark 3: Structure Refinement Details
      if ( $remark == 3 ) { 

        #----- Refinement Program
        if ( $line =~ /PROGRAM\s+:\s+(.+?)\s*$/ ) { 
          $meta{refinement_program} = $1;
        }

        #----- Atom Counts
        if ( $line =~ /PROTEIN ATOMS\s+:\s+(\d+)\s+$/ ) { 
          $meta{protein_atoms} = $1;
          next;
        }
          
        if ( $line =~ /NUCLEIC ACID ATOMS\s+:\s+(\d+)\s+$/ ) { 
          $meta{nucleic_atoms} = $1;
          next;
        }
          
        if ( $line =~ /HETEROGEN ATOMS\s+:\s+(\d+)\s+$/ ) { 
          $meta{het_atoms} = $1;
          next;
        }
          
        if ( $line =~ /SOLVENT ATOMS\s+:\s+(\d+)\s+$/ ) { 
          $meta{solvent_atoms} = $1;
          next;
        }

        #----- R and R FREE
        if ( $line =~ /  R VALUE\s+(.*?)\s+: ([\.\d]+)\s+$/ ) { 
          my ($r_set, $r) = ($1, $2);
          #print "R Work: $r_set / $r\n";
          if ( exists $meta{'r_work'} ) {
            next if $r_set =~ /TEST SET/;
            next if $r_set !~ /WORKING SET/;
          }
          $meta{'r_work'} = $r;
          next;
        }

        if ( $line =~ /FREE R VALUE\s+: ([\.\d]+)\s+$/ ) { 
          $meta{'r_free'} = $1;
          next;
        }
      }

      #======= Remark 4: PDB Version
      if ( $remark == 4 ) { 
        if ( $line =~ /COMPLIES WITH FORMAT V\. (\d)\.(\d\d), (\d\d)-(\w\w\w)-(\d\d)/ ) { 
          $meta{pdb_version      } = "$1.$2";
          $meta{pdb_version_major} = $1;
          $meta{pdb_version_minor} = $2;
          $meta{pdb_version_date } = "$3-$4-$5";
          $meta{pdb_version_day  } = $3;
          $meta{pdb_version_month} = $4;
          $meta{pdb_version_year } = $5;
        }
      }

      #======= Remark 200: Experimental Details
      # TODO: rmerge, rsym, i/sigmai
      if ( $remark == 200 ) { 

        #----- Diffraction Temperature
        if ( $remark_info =~ /^ TEMPERATURE\s+(\(\w+\))? : (\S+)\s*$/ ) { 
          $meta{dfxn_temperature_unit} = $1;
          $meta{dfxn_temperature_unit} =~ s/\(|\)//g;
          $meta{dfxn_temperature     } = $2;
        }

        #----- pH
        if ( $remark_info =~ /^ PH\s+: (\S+)\s*$/ ) { 
          $meta{ph} = $1;
        }

        #----- Overall Completeness
        if ( $remark_info =~ /^ COMPLETENESS FOR RANGE     \(%\) : ([\d\.]+)\s*$/ ) { 
          $meta{completeness_all} = $1;
        }

        #----- Completeness in HR Shell
        if ( $remark_info =~ /^ COMPLETENESS FOR SHELL     \(%\) : ([\d\.]+)\s*$/ ) { 
          $meta{completeness_hr} = $1;
        }


      }


    } #======= END REMARKS


    #========= Space Group and Unit Cell
    my $num = '([\d\.]+)\s+';
    if ( $line =~ /^CRYST1\s+$num$num$num$num$num$num(.*)\s+(\d+)\s*$/ ) { 
      $meta{edges     } = [$1, $2, $3];
      $meta{angles    } = [$4, $5, $6];
      $meta{edges_txt } = sprintf "%.2f,%.2f,%.2f", @{$meta{'edges'}};
      $meta{angles_txt} = sprintf "%.2f,%.2f,%.2f", @{$meta{'angles'}};
      $meta{spacegroup} =  $7;
      $meta{monomers  } =  $8;
    }

    #========= DB Ref
    if ( $line =~ /^DBREF/ ) { 
      my @dbref = split /\s+/, $line;

      my %dbref = ();
      ( $dbref{record   }, $dbref{idcode       }, $dbref{chain}, 
        $dbref{seq_begin}, $dbref{seq_ins_begin}, 
        $dbref{seq_end  }, $dbref{insert_end   },
        $dbref{database }, $dbref{accession    }, $dbref{dbid},
        $dbref{db_begin }, $dbref{db_ins_begin },
        $dbref{db_end   }, $dbref{db_end_ins},
      ) = map { Utility::Trim($_) } 
            unpack "A7A5A2 A4A2 A4A1xA6xA8xA12xA5A1xA5A1", $line;
      
      $dbref{seq_length} = $dbref{seq_end}-$dbref{seq_begin}+1;
      $dbref{db_length } = $dbref{db_end }-$dbref{db_begin }+1;


      if ( not defined $meta{dbref} ) { $meta{dbref}=[] }
      push @{$meta{dbref}}, \%dbref;

      push @{$meta{chains}}, $dbref{chain};
    }

    #========= Ligand Formulas
    # TODO:
    if ( $line =~ /FORMUL \s+(\d+)\s+(\w+)\s+(.*)$/ ) { 

      my $comp    = $1;
      my $id      = $2;
      my $formula = $3;


    }

    #========= Animal Source for Protein
    # TODO:

    #========= SeqRes
    if ( $line =~ /^SEQRES/ ) { 
 
      my ( $record, $N, $chain, $length, @sequence ) = 
        map { Utility::Trim($_) } 
          unpack "A6xA3xA1xA4xx(A4)*", $line;

      @sequence = grep {$_} @sequence;

      if ( not exists $meta{sequence}{$chain} ) { 
        $meta{sequence}{$chain} = [];
      }
      #print "LENGTH: $chain $length\n";
      push @{$meta{sequence}{$chain}}, @sequence;

    }




    #========= Compound
    if ( $line =~ /^COMPND/ ) { 
      if ( $line =~ /MOL_ID: (\d+)/ ) {  
        my $id = $1;
        if ( (not defined $last->{id}) or ($id != $last->{id} ) ) { 
          push @{$meta{compounds}}, {};
          $last = $meta{compounds}[-1];
        }
                                           $last->{mol_id  }  = $1 }
      elsif($line =~ /MOLECULE: (.*);/ ) { $last->{molecule} .= $1 }
      elsif($line =~ /CHAIN: (.*);/    ) { $last->{chain   }  = [split /\s+,\s+/, $1] }
      elsif($line =~ /SYNONYM: (.*);/  ) { $last->{synonym }  = $1 }
      elsif($line =~ /EC: ([\d\.]+)/   ) { $last->{ec      }  = $1
      }
    }
  }



  #=========================== Parse Residues to get it right


  


  #========================== Write summary file to .meta
  

  #-------------------------- Create Text Summary
  my $meta_data = "";
  foreach my $prop (@props) { 
    if ( not exists $meta{$prop} ) { 
      warn "Didn't find or couldn't parse data for $prop.\n" if $flags->{meta_verbose} > 3;
    }
    $meta_data .= sprintf "%-30s%s\n", uc("$prop"), ($meta{$prop}||"NULL");
  }

  #------------------------- Write File
  if ( $flags->{pdb_meta_write} ) { 
    my $meta_path = $pdb_path;
    $meta_path =~ s/\.pdb/\.meta/;
    Utility::Add_Path_to_Config($CFG, 'META', $meta_path);

    Utility::Write_File($meta_path, $meta_data, 1);
    #---- Raw version of the meta data structure
    if ( $flags->{pdb_meta_write_raw} ) { 
      my $meta_raw_path = "$meta_path.raw";
      Utility::Add_Path_to_Config($CFG, 'META_RAW', $meta_raw_path);
      Utility::Write_File($meta_raw_path, Dumper(\%meta), 1);
    }
  }

  #---- And print it 
  print $meta_data if $flags->{pdb_meta_output};


  #========================== Add it to the CFG
  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{META} = \%meta;
  }

  #========================== Return the data
  if ( wantarray ) { return  %meta }
  else             { return \%meta }
}








####################################################### FIND PDB AND MTZ FILES

sub Find_PDB_and_MTZ_Files {

  my $path_pdb = $_[0];
  my $path_mtz = $_[1];
  
  my ($dir_pdb, $file_pdb);
  my ($dir_mtz, $file_mtz);

  #=============== Nothing Passed. Look in Current Directory
  if ( not defined $path_pdb ) { 
    #confess "Need a PDB file or a refinement directory to start from.";
    $path_pdb = "."
  }

  #=============== Passed a PDB File
  if ( -f $path_pdb ) { 

    #------------- Check extension
    if ( not $path_pdb =~ /\.pdb$/ ) { 
      confess "Was passed a file, but it's not a .pdb ($path_pdb)";
    }
    #------------- Generate path/file/dir
    # NOTE: For now, just leaving as (potentially) relative path.
    #($path_pdb, $dir_pdb, $file_pdb) = Utility::Split_Path($path_pdb);


    #============= Now find matching MTZ file
    
    #------------- Easiest: Look for one matching the PDB name
    if ( not defined $path_mtz ) { 
      $path_mtz = $path_pdb;
      $path_mtz =~ s/\.pdb/\.mtz/;
    }

    #------------- Edited PDB in Coot after Phenix Refinement
    if ( not -e $path_mtz ) { 
      $path_pdb =~ /refine_(\d+)-/;
      my $refine = $1;
      $path_mtz = "Refine_$refine/$path_pdb";
      $path_mtz =~ s/-coot-\d*.pdb/.mtz/;
    }

    #------------- No luck!
    if ( not -e $path_mtz ) { 
      confess "Couldn't find a MTZ file '$path_mtz'\n";
    }
  }


  #=============== Passed a Directory. Find common PDB/MTZ pairs
  if ( -d $path_pdb ) { 
    $dir_pdb = $path_pdb;
    #------------- Phenix Refinement Directory
    if ( $path_pdb =~ /Refine_(\d+)/ ) { 
      my $refine = $1;
      $path_pdb = (glob("$dir_pdb/*refine_$refine.pdb"))[0];
      $path_mtz = (glob("$dir_pdb/*refine_$refine.mtz"))[0];

      if ( not -e $path_pdb ) { die "Couldn't find pdb in refine dir '$dir_pdb'" }
      if ( not -e $path_mtz ) { die "Couldn't find mtz in refine dir '$dir_pdb'" }

    #------------- Find ANY PDB/MTZ Pair
    } else { 
      my @pdbs = glob("$dir_pdb/*.pdb");
      foreach my $pdb ( @pdbs ) { 
        my $mtz = $pdb;
        $mtz =~ s/\.pdb/\.mtz/;
        if ( -e $mtz ) { 
          $path_pdb = $pdb;
          $path_mtz = $mtz;
          last;
        }
      }
      if ( not $path_mtz ) { 
        confess "Couldn't find a matching pdb/mtz pair with '$path_pdb'";
      }
    }
  }


  #=============== Return matched PDB/MTZ pair
  return ($path_pdb, $path_mtz);

}




########################################################## DETECT TRUE ALTCONFS
# 
# The algorithm here is basically wrong now, but it helped me think it out.
# Loop is now AltT - Set - AltR - Atom
#
# Should work on all residue types at all resolution levels.
#------------------------------------------------------------------------------
# True altconfs are defined here as a conformation with at least one atom that
# is at least $DIST_THRESHOLD (0.5A) from the corresponding atom in each 
# of the stronger altconfs.  Thus the first altconf is always "real" (no 
# stronger altconfs).
#
# This specific phrasing above is important, as it justifies the odd nested
# loop structure below.  AltT - Atom - AltR.  AltT is the alt being tested,
# so AltR is always a stronger Reference Altconf (by occupancy).  
# If any single atom in AltT is distinct from all AltR copies, then all 
# of AltT is distinct, or "real".  So for the tested AltT, we take one 
# atom at a time and compare it against all previous versions.  At the end 
# of that comparison, if it is distinct from all, then we could terminate 
# the whole AltT loop succesffully. 
#
# However, we do not.  Instead this is an opportunity to also measure "how" 
# distinct an altconf is, so we'll sum the number of atoms which are 
# distinct for the Collapse_Atom_Properties call. Then we'll reexamine the
# structure and, instead of counting distinct atoms in an altconf, we'll
# count distinct altconfs for an atom. By copying the tree before collapsing,
# for each group/altconf, we take the max instead of the sum.  Then for 
# altconf '_', we take the sum.
#
#
# To summarize:
#
# DISTINCT_ATOMS    : For a group/altconf, count the number of distinct atoms
# DISTINCT_ALTCONFS : For a group/(alt '_'), count distinct altconfs
##############################################################################

sub Detect_Distinct_Altconfs { 

  my $CFG = $_[0];

  my @all_sets = @Utility::ATOM_GROUPS;

  foreach my $chain ( @{$CFG->{PDB}{CHAINS}} ) { #====================== Chains
    foreach my $resi ( 0..$#{$CFG->{PDB}{RESI}{$chain}} ) { #========= Residues
      my $resn = $CFG->{PDB}{RESN}{$chain}[$resi];
      #my $type = $CFG->{PDB}{TYPE}{$chain}[$resi];
      next if not defined $resn;

      #next if $resi != 39; 

      #if ( not defined $type ) { die "$chain  $resi  $resn\n" }
      #next if not uc $type eq "P";
      my $XYZ = $CFG->{PDB}{COORDS}{$chain}[$resi];
      my @alts = @{$CFG->{PDB}{ALTS}{$chain}[$resi]};
      my $AA = $Utility::AA{$resn};

      #print "\nWorking on residue $resi $resn with alts @alts\n";
      my $alt_main = $alts[0];
      $alt_main = '_' if not $alt_main or $alt_main eq ' ';
      
      #===== Remove ' ' or '_' alts, as they're not meant to be distinct
      @alts = grep { not ($_ eq ' ' or $_ eq '_') } @alts;

      
      my @all_atoms = keys %$XYZ;
      my @sets = (@all_atoms, @all_sets);

      #====== Create a cache on {atom}[altT][altR] to save distance calculations
      my %DISTINCT_CACHE = ();

      #====== Set first alt as distinct
      foreach my $set ( @sets ) { 
        #my @atoms = Utility::Get_Atom_Set_Elements($set, $resn, \@all_atoms);
        my @atoms = Utility::Get_Atom_Group_Components($set, $resn, \@all_atoms);
        my $gtype = $Utility::ATOM_GROUPS{$set} || 'unknown';
        next if not @atoms;
        #print "C: $chain Res: $resi $resn Set: $set Atoms: @atoms GType: $gtype\n" if $chain eq 'B';
        foreach my $atom ( @atoms, $set ) { 
          $CFG->{PDB}{DISTINCT_ALTS}{$chain}[$resi]{$atom}{$alt_main} = 1;
          $CFG->{PDB}{DISTINCT_ALTS}{$chain}[$resi]{$atom}{'*'      } = 1;
          $CFG->{PDB}{DISTINCT_ALTS}{$chain}[$resi]{$atom}{'_'      } = 1;
          #print "Distinct Alts: $chain  $resi $resn  $set $atom\n";
        }
      }
   
      #====== Test Subsequent AltConfs 
      foreach my $T ( 1..$#alts ) { #========================== Alt Testing
        my $altT = $alts[$T];
            
        foreach my $set ( @sets ) { #==================================== Sets
          my @atoms = Utility::Get_Atom_Group_Components($set, $resn, \@all_atoms);
          next if not scalar @atoms; # Set doesn't apply to this residue probably.
          #my @atoms = Utility::Get_Atom_Set_Elements($set, $resn, \@all_atoms);
          my $distinct_set = 1;
        
          foreach my $R ( 0..($T-1) ) { #================= Alt Reference
            my $altR = $alts[$R];
            my $distinct_alt = 0;
            
            foreach my $atom ( @atoms ) { #============================= Atoms

              #=== Skip Hydrogen Atoms
              next if Utility::Is_Hydrogen_Atom($atom);
           
              #=== Lookup distinction in cache 
              my $distinct_atom = $DISTINCT_CACHE{$atom}[$T][$R];
             
              #print "Comparing distinct alts for $resi $resn  target $T '$altT'  reference $R '$altR'   set $set  atom $atom\n";
              #=== If not there, check distinction
              if ( not defined $distinct_atom ) { 
                #--- Get Coordinates
                
                if ( not exists $XYZ->{$atom}{$altT} ) { next }
                if ( not exists $XYZ->{$atom}{$altR} ) { next }
                my $xyzT = $XYZ->{$atom}{$altT};
                my $xyzR = $XYZ->{$atom}{$altR};
                if ( not $xyzT ) { die "Missing xyzT for Chain $chain   Res $resi $resn   Atom Set $set $atom   Alt $altT." }
                if ( not $xyzR ) { die "Missing xyzR for Chain $chain   Res $resi $resn   Atom Set $set $atom   Alt $altR." }
                #--- Get Distance and Test
                my $dist = Utility::Distance( $xyzT, $xyzR );
                $distinct_atom = ( $dist > $ALTCONF_DISTINCTION_DISTANCE ) ? 1 : 0;
                $DISTINCT_CACHE{$atom}[$T][$R] = $distinct_atom;

              }
              #print "   Distinct? $distinct_atom\n";
  
              
              #=== Distinct atom implies distinct alt
              $distinct_alt ||= $distinct_atom;
            } #/END ATOM

            #=== A setmust be distinct from all altconfs
            $distinct_set &&= $distinct_alt;
            #print "  distinct set: $distinct_set\n";
          } #/END ALT REF
        
          #=== Record Final Observation on Atom/Set Distinction
          $CFG->{PDB}{DISTINCT_ALTS}{$chain}[$resi]{$set}{$altT} = $distinct_set;
          #print "DISTINCT:  $chain   $resi $resn  $set  '$altT'  '$distinct_set'\n";
          
          #=== Increment Count for Summary Alt
          $CFG->{PDB}{DISTINCT_ALTS}{$chain}[$resi]{$set}{'*'} += $distinct_set;

        } #/END SET
      } #/END GROUP
    } #/END RESI
  } #/END CHAIN


  #===== Properties are already collapsed, so just report it
  #Collapse_Atom_Properties($CFG->{PDB}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DISTINCT', 'sum');

  Register_Analyses(
    $CFG,
    [qw(ATOM ANGLE RESIDUE)],
    'PDB',
    [qw(DISTINCT_ALTS)],
    [qw(PROTEIN WATER LIGAND)],
  );

  #print "Distinct Alts: \n" . Dumper($CFG->{PDB}{DISTINCT_ALTS});
  #exit;
}
      


#################################################### Calculate Dihedral Angles
# For all defined residue dihedral angles (phi, psi, chi1, etc), calculate
# the angles for every altconf of every residue.
##############################################################################

sub Calculate_Dihedral_Angles { 

  my $CFG = $_[0];
  my %DH = %Utility::DIHEDRAL_ANGLES;

  foreach my $chain ( @{$CFG->{PDB}{CHAINS}} ) { #----------------- Loop Chains
    #foreach my $resi ( 1, 2, 121, 122, 123, 124 ) { 
    foreach my $resi ( 0..$#{$CFG->{PDB}{RESI}{$chain}} ) { #---- Loop Residues
      my $resn = $CFG->{PDB}{RESN}{$chain}[$resi];
      my $type = $CFG->{PDB}{TYPE}{$chain}[$resi];
      next if not defined $resn;
      if ( not defined $type ) { die "$chain  $resi  $resn\n" }
      next if not uc $type eq "P";
      my $X = $CFG->{PDB}{X}{$chain}[$resi];
      my @alts = @{$CFG->{PDB}{ALTS}{$chain}[$resi]};
      my $AA = $Utility::AA{$resn};
      foreach my $alt ( @alts ) { #---------------------------------- Loop Alts
        my $a = $alt; $a = '_' if $a eq ' ';
    
        foreach my $dh ( sort keys %DH ) { 
          if ( exists ( $DH{$dh}{uc $resn} ) ) { 
            my @atom_names = @{$DH{$dh}{uc $resn}};
            my @atom_xyzs = ();
            foreach my $nameoff ( @atom_names ) { 
              my $offset = 0;
              if ( $nameoff =~ /\-/ ) { $offset = -1 }
              if ( $nameoff =~ /\+/ ) { $offset =  1 }
              my $name = $nameoff;
              $name =~ s/\+//;
              $name =~ s/\-//;


              # Skip N- and C-terminal residues
              last if ($resi+$offset) < 0;
              last if not defined $CFG->{PDB}{RESN}{$chain}[$resi+$offset];
             
              # Get the coordinates 
              my $coords = $CFG->{PDB}{COORDS}{$chain}[$resi+$offset]{$name}{$alt};
              # If not defined, see if there is a major alt we can use.
              if ( not defined $coords ) { 
                $coords = $CFG->{PDB}{COORDS}{$chain}[$resi+$offset]{$name}{_};
              }
              # Skip residues that are essentially missing an atom.
              #print "Checking coords defined\n";
              last if not defined $coords; # Skip undefined coordinates
              #print "Checking coords array\n";
              last if not ref $coords eq 'ARRAY';
              #print "Checking coords is 3\n";
              last if scalar @$coords != 3;
              push @atom_xyzs, $coords;
            }
            # Skip dihedrals that are missing an atom.
            #print "Check dihedral Angle of $resi, $resn, $type, $alt, $dh, @atom_names\n";
            #print "Check dihedral Angle of $resi, $resn, $type, $alt, $dh, @atom_names " . Dumper(\@atom_xyzs) . "\n";
            next if ((scalar @atom_xyzs) != 4);
            #print "Dihedral Angle of " . Dumper(\@atom_xyzs) . "\n";
            my $angle = Utility::Dihedral_Angle(@atom_xyzs);

            # Account for symmetric Chi Angles.
            if ( exists $AA->{sym}{lc $dh} ) { 
              #print "Chopping down an angle: $resn, $dh, $angle\n";
              #exit;
              $angle %= (360.0 / $AA->{sym}{lc $dh});
            }

            my $dh_map = $Utility::ATOM_SET_MAP{$dh};
            #printf "Angle of $resi, $resn, $type, $alt, $dh, @atom_names: %8.2f\n", $angle;
            if ( not defined $dh_map ) { 
              confess "Couldn't get atom super set definition for dihedral angle $dh.\n";
            }
            #$CFG->{PDB}{DIHEDRALS}{$chain}[$resi]{$dh}{$a} = $angle;
            $CFG->{PDB}{DIHEDRALS}{$chain}[$resi]{$dh_map}{$a} = $angle;
            if ( not defined $angle ) { 
              die "This isn't an angle.\n";
            }
          }
        }
      }
    }
  }
  Collapse_Atom_Properties($CFG->{PDB}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DIHEDRALS', 'nan');
  Register_Analyses(
    $CFG,
    [qw(ANGLE)],
    'PDB',
    [ qw(DIHEDRALS) ],
    ['PROTEIN'],
  );

}
      







##############################################################################
#                                                                  CCP4 ANGLES
##############################################################################

sub CCP4_Angles { 

  die "Don't use this.  Dihedral angles is better. Just here for nostalgia or something.\n";

#  print "Doing Angles\n";
#  my $angles_target = "$title.angles.dat";
#  my $angplt_target = "$title.angles.plt";
#  my $angps_target  = "$title.angles.ps";
#  my $angpng_target = "$title.angles.png";
#
#  my $dir_angles = $CFG->{DIRS}{ANGLES};
#  Utility::Make_Directories($dir_angles);
#  my $cwd = Cwd::getcwd();
#  chdir $dir_angles or croak "Chdir $dir_angles $!";
#  if ( not -e $angles_target ) { 
#    print "  considering angles\n";
#
#
#    my $cmd = "angles XYZIN ../$pdb ANGOUT $angles_target PLOT $angplt_target << EOF\nDIHEDRAL\nRESIDUE 1 124\nTITLE $title\nCHAIN A\nEOF\n";
#
#    print "$cmd\n";
#    `$cmd`;
# 
#    my $plot = "pltdev -dev ps -xp 0 -yp 0 -pnt 2.0 -lnw 1.0 -lan -i $angplt_target -o $angps_target"; 
#    `$plot`;
#
#    my $cnvt = "convert $angps_target -rotate 90 $angpng_target";
#    `$cnvt`;
#    
#  }
#  chdir $cwd or croak "Couldn't chdir '$cwd'. $!";
}








##############################################################################
##############################################################################
#                              PHENIX INTERFACE
##############################################################################
##############################################################################




##############################################################################
#                                                            FETCH PDB AND MTZ
##############################################################################


# Generic download subroutine for both PDB and MTZs.

sub Fetch_PDB_General { 

  my $PDB   = $_[0];
  my $type  = $_[1] || 'pdb';
  my $file  = $_[2] || "$PDB.$type";
  my $flags = $_[3];

  my $do_sf = ($type eq 'mtz' or $type eq 'cif') ? "-x" : "";

  # Check if it already exists
  if ( -e $file and not $flags->{'overwrite'} ) { 
    #if ( $flags->{'verbose'} ) { 
      print "File $file already downloaded.\n";
    #}
    return;
  }


  # Download it
  my $command = "$PROGRAMS{'FETCH_PDB'} $do_sf $PDB";
  my %result = Utility::Try_Command ( 
                  $command,
                  "Downloading pdb $PDB",
                  $flags,
                );

  # Check standard structure factors not existing error
  if ( $result{'err'} =~ /Couldn't download structure factors/m ) { 
    print "Structure factors do not exist for pdb $PDB.\n";
    return;
  }

  # Check generally that file could be downloaded
  if ( $result{'value'} ) { 
    Utility::Report_Error ( 
      "Couldn't download $type file for pdb '$PDB', '$file'.",
      \%result,
      $flags,
    );
    return;
  }

  # Double check file is downloaded
  my $out_file;
  if ( $result{'out'} =~ /saved to (.*)$/ ) { 
    $out_file = (Utility::Split_Path($1))[2];
  }
  if ( (not $out_file) or (not -e $out_file) ) { 
    Utility::Report_Error ( 
      "Could not download $type file for pdb '$PDB', '$file' (invisibly).",
      "",
      $flags,
    );
    return;
  }

  # Rename file if indicated
  if ( $file ne $out_file ) {  
    rename $out_file, $file
      or confess "Couldn't rename file $out_file to $file. $!";
  }

  return;
}



# Fetches the PDB and then decides whether to fetch the MTZ.
sub Fetch_PDB_and_CIF { 

  my $pdb   = $_[0];
  my $f_pdb = $_[1];
  my $f_cif = $_[2];
  my $flags = $_[3];

  my $flag_no_sf = "NO_STRUCTURE_FACTORS_FOR_$f_pdb";

  Fetch_PDB_General ( $pdb, 'pdb', $f_pdb, $flags );
  my $meta = Parse_PDB_Meta($f_pdb);
  if ( $meta->{'exp_type'} eq 'XRAY' or $meta->{exp_type} eq 'NEUTRON') { 
    if ( not -e $flag_no_sf ) { 
      Fetch_PDB_General ( $pdb, 'cif', $f_cif, $flags );
    } else { 
      print "PDB $pdb was previously found to have no structure factors.\n";
    }
  } else { 
    print "PDB $pdb is not a crystal structure, so has no structure factors.\n";
  }
  $meta->{'SF'} = (-e $f_cif );
  if ( not $meta->{'SF'} ) { 
    `touch $flag_no_sf`;
  }
  if ( wantarray ) { return %$meta }
  else             { return  $meta }

}


##############################################################################
#                         FETCH PDB REDO FILES
##############################################################################

sub Fetch_PDB_Redo { 

  my $CFG = $_[0];
  
  
  my $pdb4;
  my $job = 'redo';

  if ( Is_Structure_Config($CFG) ) { 
    $pdb4 = $CFG->{PDB4};
    $job  = $CFG->{JOB} || $job;
  } else { 
    $pdb4 = $CFG;
    $job  = $_[1] || $job;
  }

  if ( not defined $pdb4     ) { confess("Maybe this isn't a PDB4 file $CFG->{TITLE}?"); }
  if ( $pdb4 !~ /^\d\w\w\w$/ ) { confess("Probably not a pdb4 file. $pdb4, $CFG->{TITLE}"); }

  if ( -e "NO_PDB_REDO" ) { 
    return;
  }
  
  my $p = lc $pdb4;
  my $q = substr($p,1,2);
  my $p_pdb_final   = "$p\_final.pdb";
  my $p_mtz_final   = "$p\_final.mtz";
  my $p_pdb_besttls = "$p\_besttls.pdb";
  my $p_mtz_besttls = "$p\_besttls.mtz";

  my $f_pdb_redo = "$pdb4.$job.pdb";
  my $f_mtz_redo = "$pdb4.$job.mtz";

  my $result;
  if ( not -e $f_pdb_redo ) { 
    $result = Utility::Try_Command(
      "rsync -avz rsync://rsync.cmbi.ru.nl/pdb_redo/$q/$p/$p_pdb_final $f_pdb_redo",
      "Downloading PDB Redo $p_pdb_final",
    );
  }
  if ( not -e $f_mtz_redo ) { 
    $result = Utility::Try_Command(
      "rsync -avz rsync://rsync.cmbi.ru.nl/pdb_redo/$q/$p/$p_mtz_final $f_mtz_redo",
      "Downloading PDB Redo $p_mtz_final",
    );
  }
  if ( not -e $f_pdb_redo ) { 
    $result = Utility::Try_Command(
      "rsync -avz rsync://rsync.cmbi.ru.nl/pdb_redo/$q/$p/$p_pdb_besttls $f_pdb_redo",
      "Downloading PDB Redo $p_pdb_besttls",
    );
  }
  if ( not -e $f_mtz_redo ) { 
    $result = Utility::Try_Command(
      "rsync -avz rsync://rsync.cmbi.ru.nl/pdb_redo/$q/$p/$p_mtz_besttls $f_mtz_redo",
      "Downloading PDB Redo $p_mtz_besttls",
    );
  }
  
  if ( not -e $f_pdb_redo ) { 
   `touch NO_PDB_REDO`;
   Utility::Report_Error("Could not download PDB REDO $f_pdb_redo", "No file found");
   return;
  }


 
  my %REDO = Create_Structure_Configuration ( 
    $CFG->{DIR},
    $CFG->{TITLE},
    $job,
    $f_pdb_redo,
    $f_mtz_redo,
    $CFG->{FLAGS},
    $CFG->{DO},
  );

  if ( wantarray ) { return  %REDO } 
  else             { return \%REDO }

}


###################################################################### XTriage
##############################################################################

sub XTriage {
  
  print "----- XTriage\n";

  my $CFG = $_[0];

  my ($mtz_path, $xtriage_path, $mtz_arrays);
  
  if ( Is_Structure_Config($CFG) ) { 
    $mtz_path     = $CFG->{PATHS}{MTZ};
    $xtriage_path = $CFG->{PATHS}{XTRIAGE};
    $mtz_arrays   = $CFG->{DUMP}{iobs} || $CFG->{DUMP}{fobs};
  } elsif ( $CFG =~ /\.mtz/ ) { 
    $mtz_path = $CFG;
  } elsif ( $CFG =~ /\.xtriage\.log/ ) { 
    $xtriage_path = $CFG;
  } else { 
    confess("Phenix::XTriage must receive a CFG, an mtz, or an xtriage.log file.");
  }

  if ( not defined $xtriage_path) { 
    if ( not defined $mtz_path ) {
      carp("No MTZ file defined for xtriage analysis.\n");
      return;
    }
    $xtriage_path = $mtz_path;
    $xtriage_path =~ s/\.mtz/\.xtriage\.log/;
    if ( Is_Structure_Config($CFG) ) { 
      Utility::Add_Path_to_Config($CFG, 'XTRIAGE', $xtriage_path);
    }
  } 

  if ( -e $xtriage_path ) { 
    if ( -s $xtriage_path == 0 ) { 
      print "XTriage file is empty. Regenerating.\n";
      unlink $xtriage_path or confess "Couldn't delete mtz dump, '$xtriage_path'. $!";
    } elsif ( -M $xtriage_path > -M $mtz_path ) { 
      print "XTriage file is older than mtz. creating new replacement.";
      unlink $xtriage_path or confess "Couldn't delete mtz dump, '$xtriage_path'. $!";
    }
  }

  # Parse protein length from the CFG file, or guess.
  if ( not exists $CFG->{TOTAL_RESIDUES} ) { 
    print " Guessing protein name from path $CFG->{PATHS}{MTZ}\n";
    foreach my $p ( keys %KNOWN_PROTEINS ) { 
      if ( $CFG->{PATHS}{MTZ} =~ /$p/i ) { 
        print "  Found protein $p\n\n";
        $CFG->{TOTAL_RESIDUES} = $KNOWN_PROTEINS{$p}{length};
        last;
      }
    }
  }


  # No existing xtriage analysis
  if ( not -e $xtriage_path ) {  

    if ( not -e $mtz_path ) {
      carp("MTZ file $mtz_path and XTriage file $xtriage_path not found.\n");
      return;
    }
    
    # Run XTriage if necessary 
    my $labels = $CFG->{DUMP}{iobs} || $CFG->{DUMP}{root};
    my $cmd_xtriage = "$PROGRAMS{XTRIAGE} "
                      . "'$CFG->{PATHS}{MTZ}' "
                      . "log='$CFG->{PATHS}{XTRIAGE}' "
                      . "obs_labels=$labels "
                      . (exists $CFG->{TOTAL_RESIDUES} ? "n_residues=$CFG->{TOTAL_RESIDUES} " : "")
                      ;

    my $result = Utility::Try_Command($cmd_xtriage, "XTriage");
    if ( $result ) { Utility::Report_Error("XTriage Failed.", 
                     $result, 
                     $CFG->{report_error_flags}) }
  }
   
  
  
  my %xtriage = (
    notes   => [],
    issues  => [],
  );
  my $contents = Utility::Read_File($CFG->{PATHS}{XTRIAGE});

  ##### Extract Date
  $contents =~ /Date .* \(([\d\.]+) s\)/ 
    or carp("XTriage couldn't parse the date in $CFG->{PATHS}{XTRIAGE}.");
  $xtriage{timestamp} = $1;
  $xtriage{localtime} = localtime($1);
  
  ##### Indices
  (($xtriage{mtz_path}, $xtriage{mtz_array})
    = ($contents =~ /Miller array info: (.*):(.*)/))
    or carp("XTriage couldn't parse the miller indices in $CFG->{PATHS}{XTRIAGE}.");

  (($xtriage{indices}) = ($contents =~ /Miller indices: (\d+)/))
    or carp("XTriage couldn't parse miller indices in $CFG->{PATHS}{XTRIAGE}.");

  $contents =~ /Unit cell: \(([\d\., ]+)\)/
    or carp("XTriage couldn't parse the unit cell in $CFG->{PATHS}{XTRIAGE}.");
  my $uc = $1;
  $uc =~ s/,//;
  $xtriage{unitcell} = [split /\s+/, $uc];


  (($xtriage{spacegroup}) = ($contents =~ /Space group: (.*) \(/))
    or carp("XTriage couldn't parse the space group in $CFG->{PATHS}{XTRIAGE}.");

  (($xtriage{pointgroup}) = ($contents =~ /likely point group of this data is: (.*)/))
    or $xtriage{pointgroup} = "N/D";
    #or warn("XTriage couldn't parse the point group in $CFG->{PATHS}{XTRIAGE}.");

  (($xtriage{rsln_low}, $xtriage{rsln_high})
    = ($contents =~ /Resolution range: ([\d\.]+) ([\d\.]+)/))
    or carp("XTriage couldn't parse the resolution in $CFG->{PATHS}{XTRIAGE}.");


  ################ COMPLETENESS

  #=============== Overall Completeness
  (($xtriage{comp_all})
    = ($contents =~ /Completeness in resolution range: ([\d\.]+)/))
    or carp("XTriage couldn't parse the resolution in $CFG->{PATHS}{XTRIAGE}.");



  #=============== Full CCP4 Completeness Table
#print "Header: $comp_header\n";
#print "Data:\n$comp_data\n";
  $xtriage{COMP_CUTOFF} = 0.90;
  if ( $contents =~ /TABLE: Z scores and comp.*?\$\$\n(.*?)\s+\$\$\n(.*?)\n\$\$/s ) { 
    my ($comp_header, $comp_data) = ($1, $2);
    my @comp_data = Utility::Text_to_Matrix($comp_data, { delimiter=>'\s+' });
    map { $_->[3] = 1/sqrt($_->[0]) } @comp_data;
    $xtriage{comp_lr_rsln  } = $comp_data[ 0][3];
    $xtriage{comp_lr_value } = $comp_data[ 0][2];
    $xtriage{comp_hr_rsln  } = $comp_data[-1][3];
    $xtriage{comp_hr_value } = $comp_data[-1][2];
    $xtriage{comp_co_rsln  } = $comp_data[-1][3];
    $xtriage{comp_co_value } = $comp_data[-1][2];
    if ( $xtriage{comp_hr_value} < $xtriage{COMP_CUTOFF} ) { 
      foreach my $i (reverse 0..$#comp_data) { 
        if ( $comp_data[$i][2] >= $xtriage{COMP_CUTOFF} ) { 
          $xtriage{comp_co_rsln  } = $comp_data[$i][3];
          $xtriage{comp_co_value } = $comp_data[$i][2];
          last;
        }
      }
    }
  } else { 
    push @{$xtriage{notes}}, "Full completeness data not included in file.";
  } 

  #=============== Human Readable Completeness Table for HR Shell
  if ( $contents =~ /(\|.*\|)\n(\-+)\n(\s+The completeness of data for which)/ ) { 
    #exit;
    my $line = $1;
    print "Found line: \n$line\n";
    my ( $trash, $res, @comps ) = split /\s*\|\s*/, $line;
    map { s/%// } @comps;
    my ( $res_low, $res_high ) = split /\s*-\s*/, $res;
    my ( $comp_1s, $comp_2s, $comp_3s, $comp_5s, $comp_10s, $comp_15s ) = @comps;

    $xtriage{comp_hr_rsln  } = $res_high;
    $xtriage{comp_hr_value } = $comp_1s  / 100.0;
    $xtriage{comp_hr_1sig  } = $comp_1s  / 100.0;
    $xtriage{comp_hr_2sig  } = $comp_2s  / 100.0;
    $xtriage{comp_hr_3sig  } = $comp_3s  / 100.0;
    $xtriage{comp_hr_5sig  } = $comp_5s  / 100.0;
    $xtriage{comp_hr_10sig } = $comp_10s / 100.0;
    $xtriage{comp_hr_15sig } = $comp_15s / 100.0;
  }



  ################ I/SigI
  $xtriage{ISIGI_CUTOFF} = 1.5;
  if ( $contents =~ /TABLE: <I\/sigma_I>.*?\$\$\n(.*?)\s+\$\$\n(.*?)\n\$\$/s ) { 
    my ($isigi_header, $isigi_data) = ($1, $2);
    my @isigi_data = Utility::Text_to_Matrix($isigi_data, { delimiter=>'\s+' });
    map { $_->[2] = 1/sqrt($_->[0]) } @isigi_data;
    $xtriage{isigi_hr_rsln } = $isigi_data[-1][2];
    $xtriage{isigi_hr_value} = $isigi_data[-1][1];
    $xtriage{isigi_co_rsln } = $isigi_data[-1][2];
    $xtriage{isigi_co_value} = $isigi_data[-1][1];
    
    if ( $xtriage{isigi_co_value} < $xtriage{ISIGI_CUTOFF} ) { 
      foreach my $i (reverse 0..$#isigi_data) { 
        if ( $isigi_data[$i][1] >= $xtriage{ISIGI_CUTOFF} ) { 
          $xtriage{isigi_co_rsln } = $isigi_data[-1][2];
          $xtriage{isigi_co_value} = $isigi_data[-1][1];
          last;
        }
      }
    }
  } else { 
    push @{$xtriage{notes}}, "Full I/SigI data not included in file.";
  } 


    
  ###### Intensities
  if ( $contents =~ /TABLE: Intensity plots:.*?\$\$\n(.*?)\s+\$\$\n(.*?)\n\$\$/s ) { 
    my ($int_header, $int_data) = ($1, $2);
    my @int_data  = Utility::Text_to_Matrix($int_data, { delimiter=>'\s+' });
    map { $_->[3] = 1/sqrt($_->[0]) } @int_data;
    $xtriage{int_hr_rsln  } = $int_data[-1][3];
    $xtriage{int_hr_value } = $int_data[-1][2];
    $xtriage{int_min_rsln } = $int_data[-1][3];
    $xtriage{int_min_value} = $int_data[-1][2];

    foreach my $i (0..$#int_data) { 
      my $v = $int_data[$i][2];
      next if $v == 0;
      if ( $v < $xtriage{int_min_value} ) { 
        $xtriage{int_min_value} = $v;
        $xtriage{int_min_rsln } = $int_data[$i][3];
      }
    }
  } else { 
    push @{$xtriage{notes}}, "Full intensity data not included in file.";
  }


  ##### Matthews Solvent Content and Copies
  $xtriage{matthews_type} = "undetermined";
  if ( $contents =~ /Matthews coefficient(.*)Best guess :\s+(\d+)\s+(\w+)\s/is ) {
    #print "*** Fount matthews stuff ***\n";
    my $table   = $1;
    my $number  = $2;
    my $type    = $3;
    #print "Number '$number'  Type '$type'.\n";
    if ( $type eq 'residues' ) { 
      $xtriage{matthews_solvent} = .50;
      $xtriage{matthews_type   } = "residues";
      $xtriage{matthews_number } = $number;
    } else { 
      my @table = split "\n", $table;
      $xtriage{matthews_type  } = "monomers";
      $xtriage{matthews_number} = $number;
      foreach my $line ( @table ) { 
        #print "  Working line $line\n";
        next if not $line =~ /^\|\s+\d/;
        #print "    passed\n";
        my $sep = '\s*\|\s*';
        my $regnum = '[\d\.-]+';
        $line =~ /^$sep($regnum)$sep($regnum)$sep($regnum)$sep($regnum)$sep$/
          or die "Couldn't parse line $line.";
        my $copies  = $1;
        my $content = $2;
        my $coef    = $3;
        my $prob    = $4;
        if ( (not exists $xtriage{matthews_prob}) or ($xtriage{matthews_prob} < $prob) ) { 
          $xtriage{matthews_copies  } = $copies;
          $xtriage{matthews_content } = $content;
          $xtriage{matthews_coef    } = $coef;
          $xtriage{matthews_prob    } = $prob;
        }
      }
    }
  }



  
  
  ##### Quality of File 
  $xtriage{quality  } = 'Good';
  $xtriage{issues   } = [];
  $xtriage{rsln_sug } = $xtriage{rsln_high};

  # Completeness issues
  if ( $xtriage{comp_all} < $xtriage{COMP_CUTOFF} ) 
    { push @{$xtriage{issues}}, "Poor completeness overall." }

  if ( $xtriage{comp_hr_value} < $xtriage{COMP_CUTOFF} ) { 
    push @{$xtriage{issues}}, "Poor completeness in HR shell.";
    if ( $xtriage{comp_co_rsln} > $xtriage{rsln_sug} ) { 
      $xtriage{rsln_sug} = $xtriage{comp_co_rsln};
    }
  }

  if ( $xtriage{comp_hr_rsln} > $xtriage{rsln_high} ) { 
    push @{$xtriage{notes}}, "Completeness not measured to full resolution of file.";
  }

  # I/SigI Issues
  if ( $xtriage{isigi_hr_value} < $xtriage{ISIGI_CUTOFF} ) { 
    push @{$xtriage{issues}}, "Poor I/SigI in HR shell.";
    if ( $xtriage{isigi_co_rsln} > $xtriage{rsln_sug} ) { 
      $xtriage{rsln_sug} = $xtriage{isigi_co_rsln};
    }
  }

  if ( $xtriage{isigi_hr_rsln} > $xtriage{rsln_high} ) { 
    push @{$xtriage{notes}}, "I/SigI not measured to full resolution of file.";
  }

  if ( scalar @{$xtriage{issues}} > 0 ) { $xtriage{quality} = "Bad" }
  push @{$xtriage{notes}}, "The quality of this structure is $xtriage{quality}.";

  if ( $xtriage{rsln_sug} > $xtriage{rsln_high} ) { 
    push @{$xtriage{notes}}, "Resolution should be truncated to $xtriage{rsln_sug}.";
  }
  
  $CFG->{XTRIAGE} = \%xtriage;
  return %xtriage;

}

#===================================================================== PRINT XTRIAGE

sub XTriage_to_String {

  my $CFG = $_[0];

  my $xtriage;

  if ( Is_Structure_Config($CFG) ) { 
    if ( $CFG->{XTRIAGE} ) { 
      $xtriage = $CFG->{XTRIAGE};
    } else { 
      carp("CFG has no xtriage data.\n");
    }
  } elsif ( ref $CFG eq 'HASH' and exists $CFG->{xtriage_data} ) { 
    $xtriage = $CFG;
  } elsif ( $CFG =~ /\.mtz/ or $CFG =~ /\.xtriage\.log/ ) { 
    $xtriage = XTriage($CFG);
  } else { 
    confess "XTriage_to_String was not passed a Structure object.\n";
  }

  my $string = "";


  $string .= sprintf "%-25s: %s\n",        "MTZ Path",           $xtriage->{mtz_path};
  $string .= sprintf "%-25s: %s\n",        "Date",               $xtriage->{localtime};
  $string .= sprintf "%-25s: %s\n",        "Array",              $xtriage->{mtz_array};
  $string .= sprintf "%-25s: %s\n",        "Space Group",        $xtriage->{spacegroup};
  #$string .= sprintf "%-25s: %s\n",        "Point Group",        $xtriage->{pointgroup};
  $string .= sprintf "%-25s: %8d\n",        "Miller Indices",     $xtriage->{indices};
  $string .= sprintf "%-25s: %8.2fA\n",   "Resolution (Low)",   $xtriage->{rsln_low};
  $string .= sprintf "%-25s: %8.2fA\n",   "Resolution (High)",  $xtriage->{rsln_high};
  $string .= sprintf "\n";

  $string .= sprintf "%-25s: %8.0f%%\n",           "Comp Overall",   $xtriage->{comp_all}*100;
  $string .= sprintf "%-25s: %8.0f%% @ %5.2fA\n", "Comp HR",      $xtriage->{comp_hr_value}*100, $xtriage->{comp_hr_rsln};
  $string .= sprintf "%-25s: %8.0f%% @ %5.2fA\n", "Comp Cutoff",  $xtriage->{comp_co_value}*100, $xtriage->{comp_co_rsln};
  $string .= sprintf "\n";
 
  $string .= sprintf "%-25s: %8.2f  @ %5.2fA\n", "I/SigI HR",     $xtriage->{isigi_hr_value}, $xtriage->{isigi_hr_rsln};
  $string .= sprintf "%-25s: %8.2f  @ %5.2fA\n", "I/SigI Cutoff", $xtriage->{isigi_co_value}, $xtriage->{isigi_co_rsln};
  $string .= sprintf "\n";

  $string .= sprintf "%-25s: %8.0f  @ %5.2fA\n", "Intensity HR",     $xtriage->{int_hr_value}*100,     $xtriage->{int_hr_rsln};
  $string .= sprintf "%-25s: %8.0f  @ %5.2fA\n", "Intensity Minimum", $xtriage->{int_min_value}*100,     $xtriage->{int_min_rsln};
  $string .= sprintf "%-25s: %8.0f %%\n",  "Intensity HR/Min Ratio",   $xtriage->{int_hr_value}/$xtriage->{int_min_value} * 100;
  $string .= "\n";

  if ( $xtriage->{matthews_type} eq 'residues' ) { 
    $string .= sprintf "%-25s: %8d Residues\n", "ASU Content", $xtriage->{matthews_number};
  } elsif ( $xtriage->{matthews_type} eq 'monomers' ) { 
    $string .= sprintf "%-25s: %8d Monomers\n", "ASU Content", $xtriage->{matthews_number};
    $string .= sprintf "%-25s: %8.2f%%\n", "Solvent %", $xtriage->{matthews_content};

  } else { 
    $string .= "Unknown/undetermined matthews type '$xtriage->{matthews_type}'.\n";
  }
  $string .= "\n";

=pod

  $string .= sprintf "%-25s: %8.2f A\n",   "Comp HR A",          $xtriage->{comp_hr_rsln};
  $string .= sprintf "%-25s: %8.0f %%\n",  "Comp HR Value",      $xtriage->{comp_hr_value}*100;
  $string .= sprintf "%-25s: %8.2f A\n",   "Comp Cutoff A",      $xtriage->{comp_co_rsln};
  $string .= sprintf "%-25s: %8.0f %%\n",  "Comp Cutoff Value",  $xtriage->{comp_co_value}*100;
  $string .= sprintf "\n";
 
  $string .= sprintf "%-25s: %8.2f A\n",   "I/SigI HR A",          $xtriage->{isigi_hr_rsln};
  $string .= sprintf "%-25s: %8.2f \n",    "I/SigI HR Value",      $xtriage->{isigi_hr_value};
  $string .= sprintf "%-25s: %8.2f A\n",   "I/SigI Cutoff A",      $xtriage->{isigi_co_rsln};
  $string .= sprintf "%-25s: %8.2f \n",    "I/SigI Cutoff Value",  $xtriage->{isigi_co_value};
  $string .= sprintf "\n";

  $string .= sprintf "%-25s: %8.2f A\n",   "Intensity HR A",          $xtriage->{int_hr_rsln};
  $string .= sprintf "%-25s: %8.0f\n",     "Intensity HR Value",      $xtriage->{int_hr_value}*100;
  $string .= sprintf "%-25s: %8.2f A\n",   "Intensity Minimum A",      $xtriage->{int_min_rsln};
  $string .= sprintf "%-25s: %8.0f\n",     "Intensity Minimum Value",  $xtriage->{int_min_value}*100;
  $string .= sprintf "%-25s: %8.0f %%\n",  "Intensity HR/Min Ratio",   $xtriage->{int_hr_value}/$xtriage->{int_min_value} * 100;

=cut


  $string .= join "\n", @{$xtriage->{notes}};
  $string .= "\n\n";
  if ( scalar @{$xtriage->{issues}} > 0 ) { 
    $string .= "This structure has the following issues:\n"
        . (join "", map { "  $_\n" } @{$xtriage->{issues}})
        ;
  }
  #$string .= "---------------------------------------------------------------------\n\n";

  return $string;
}

sub Print_XTriage { 
  print XTriage_to_String($_[0]);
}


###################################################################### TABLE 1
# Parse table1.txt files
##############################################################################

sub Table1 { 

  print "----- Table 1\n";
  my $CFG = $_[0];
  my $path_table1;
  if ( Is_Structure_Config($CFG) ) { 
    $path_table1 = $CFG->{PATHS}{TABLE1};
  } else { 
    $path_table1 = $CFG;
  }

  if ( not $path_table1 ) { 
    confess("Phenix::Table1 did not received a CFG or file");
  }

  if ( not -e $path_table1 ) { 
    confess("Phenix::Table1 cannot fine file '$path_table1'.");
  }


  my %TABLE1 = ();
  my $table1 = Utility::Read_File($path_table1);


  
  ($TABLE1{rsln}                  ) = ($table1 =~ /Resolution, A\s+([\d\.]+)/);
  ($TABLE1{rsym}, $TABLE1{rsym_hr}) = ($table1 =~ /Rsym\s+([\d\.]+)\(([\d\.]+)\)/);
  ($TABLE1{comp}, $TABLE1{comp_hr}) = ($table1 =~ /Completeness, %\s+([\d\.]+)\(([\d\.]+)\)/);
  ($TABLE1{mult}, $TABLE1{mult_hr}) = ($table1 =~ /Multiplicity\s+([\d\.]+)\(([\d\.]+)\)/);
  ($TABLE1{isgi}, $TABLE1{isgi_hr}) = ($table1 =~ /I\/SD\s+([\d\.]+)\(([\d\.]+)\)/);

  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{TABLE1} = \%TABLE1;
  }

  return \%TABLE1;


}





##################################################################### MTZ DUMP
# Parse an MTZ file using phenix.mtz.dump
##############################################################################

sub MTZ_Dump { 

  print "----- MTZ Dump\n";

  my $CFG   = $_[0];
  my $flags = $_[1] || {};

  my $mtz_path;
  if ( ref $CFG eq 'HASH' ) { 
    $mtz_path = $CFG->{PATHS}{MTZ};
  } else { 
    $mtz_path = $CFG;
  }
 
  #---- Verify parameters 
  if ( not defined $mtz_path ) { confess "MTZ_Dump !! No MTZ file passed" }
  if ( not -e $mtz_path      ) { confess "MTZ_Dump !! MTZ file not found, '$mtz_path'." }

  #---- Prep for log file
  my $log_path = $mtz_path;
  $log_path =~ s/\.mtz/\.dump\.log/;
  if ( Is_Structure_Config($CFG) ) { 
    Utility::Add_Path_to_Config($CFG, 'MTZ_DUMP', $log_path);
  }
  

  #=============== Check if already complete, or run if necessary
  my ( $errors_in, $errors_out ) = Utility::Verify_Analysis_Files ( 
    "MTZ Dump",
    [ $mtz_path ],
    [ $log_path ],
    \%FLAGS_VERIFY_PRE,
  );

  #---- Yes, complete
  if ( ( not $flags->{overwrite} ) and not @$errors_in and not @$errors_out ) { 
      print "MTZ Dump complete for $mtz_path\n" if $VERBOSE >= 1;


  #--------------- Run MTZ Dump    
  } else { 
    my $cmd_mtz_dump = "'$Phenix::PROGRAMS{MTZ_DUMP}' '$mtz_path' > '$log_path'";
    my %result = Utility::Try_Command($cmd_mtz_dump, "Phenix - MTZ Dump", 
                  { %$flags, validate=>1, inputs=>[$mtz_path], outputs=>[$log_path] } );
    if ( $result{value} ) { return %result }
  }

  

  #=============== Parse Output  
  my (@mtz_dump, $mtz_dump);
  @mtz_dump = Utility::Read_File($log_path);
  $mtz_dump = join "\n", grep { $_ } @mtz_dump;

  #--------------- Basic Properties
  my ($space_group, $space_group_number, $point_group);
  #----- Space Group
  if ( $mtz_dump =~ /Space group symbol from file\:\s*([\w\d]+)/ ) { 
    $space_group = $1;
  }
  if ( $mtz_dump =~ /Space group number from file\:\s*([\w\d]+)/ ) { 
    $space_group_number = $1;
  }
  if ( $mtz_dump =~ /Point group symbol from file\:\s*([\w\d]+)/ ) { 
    $point_group = $1;
  }



  #--------------- Array Data

  my @mtz_data = grep { /(Resolution)|( \w: )|Crystal|Unit cell/ } @mtz_dump;
  my $mtz_data = (join "\n", @mtz_data) . "\n";

  if ( not scalar @mtz_data ) { 
    confess "Couldn't parse mtz_dump:\n\t$mtz_dump\n\nGot:\n" . (join "\n", @mtz_data) . "\n";
  }

  my ($label_fobs, $label_fobs_sig, $label_iobs, $label_iobs_sig);
  my $label_free;
 
  my @errors = ();
  
  foreach my $line ( @mtz_data ) { 

    #print "Parsing line, '$line'\n";

    #------------- Find R-Free
    if ( $line =~ /^\s+(\S*free\S*)/im ) { 
      $label_free = $1;
    }

    #------------- Find Amplitudes
    if ( $line =~ /^\s+(\S+)\s+.*amplitude/i ) { 
      next if defined $label_fobs;
      $label_fobs = $1;
      $label_fobs_sig = "SIG$label_fobs";
      #print "Found amplitude $label_fobs\n";
      if ( not $mtz_dump =~ /$label_fobs_sig/m ) {
        push @errors, "Missing amplitude SD ($label_fobs_sig)";
        undef $label_fobs_sig;
      }
    }

    #------------- Find Intensities
    if ( $line =~ /^\s+(\S+)\s+.*intensity/i 
      or $line =~ /^\s+(\S+)\s+.*I\(\+\) or I\(\-\)/i ) { 
      next if defined $label_iobs;
      $label_iobs = $1;
      $label_iobs =~ s/\(.*\)//;
      $label_iobs_sig = "SIG$label_iobs";
      #print "Found intensity $label_fobs\n";
      if ( not $mtz_dump =~ /$label_iobs_sig/m ) {
        push @errors, "Missing intensity SD ($label_iobs_sig)";
        undef $label_iobs_sig;
      }
    }
  }


  my ($label_root, $label_root_sig);
  if ( $label_fobs and $label_fobs_sig ) { 
    ($label_root, $label_root_sig) = ($label_fobs, $label_fobs_sig);
  } elsif ( $label_iobs and $label_iobs_sig ) { 
    ($label_root, $label_root_sig) = ($label_iobs, $label_iobs_sig);
  } 


  my %mtz_dump = ( fobs => $label_fobs,  fsig => $label_fobs_sig,
                   iobs => $label_iobs,  isig => $label_iobs_sig,
                   root => $label_root,  rsig => $label_root_sig,
                   free => $label_free,
                   dump => $mtz_dump,    errors  => \@errors,
                   data => $mtz_data,
                   sg   => $space_group,
  );

  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{DUMP} = \%mtz_dump;
  }

  if ( wantarray ) { return  %mtz_dump}
  else             { return \%mtz_dump}

}



############################################################### PRINT MTZ DUMP

sub Print_MTZ_Dump {
  my $CFG = $_[0];
  if ( Is_Structure_Config($CFG) ) { print $CFG->{DUMP}{data} }
  else                             { confess("XTriage::Print_MTZ_Dump is confused.") }
}


##############################################################################
#                                                       EDIT REFLECTIONS FILES
##############################################################################

############################################ CONVERT INTENSITIES TO AMPLITUDES
# --massage_intensities is analagous to the French & Wilson method implmented
# in ctruncate, though a bit different. Still should be fine.

sub IOBS_to_FOBS { 

  my $CFG       = $_[0];
  my $mtz_in;
  my $mtz_out   = $_[1];
  my $iobs      = $_[2];
  my $fobs      = $_[3];

  my $overwrite = 0;

  if ( Is_Structure_Config($CFG) ) { 
    $mtz_in   = $CFG->{FILES}{MTZ};
    $iobs     = $CFG->{MTZ}{iobs};
    $fobs     = $DEFAULT_FOBS_LABEL;
  } else { 
    $mtz_in = $CFG;
  }

  # Overwrite existing MTZ file (indirectly, after checking for errors).
  if ( (not defined $mtz_out) or ($mtz_out eq $mtz_in) ) { 
    $overwrite = 1;
    $mtz_out = "$mtz_in.temp";
  }


  my $command = "phenix.reflection_file_converter "
              . "$mtz_in "
              #. "--symmetry=$mtz_in "
              . "--non-anomalous "
              . "--write_mtz_amplitudes "
              . "--massage_intensities "
              . "--mtz-root-label=$fobs "
              . "--label=$iobs "
              . "--mtz $mtz_out "
              ;

  my %result = Utility::Try_Command($command, "Phenix - Convert IOBS to FOBS");

  if ( $result{'value'} ) { return %result }

  if ( not -e $mtz_out ) { 
    Utility::Report_Error ( "Phenix did not generate MTZ file.", "", "" );
    return;
  }

  
  if ( $overwrite ) { 
    rename $mtz_out, $mtz_in 
      or confess("Phenix::IOBS_to_FOBS:: Couldn't rename '$mtz_out' to '$mtz_in'. $!");
    $mtz_out = $mtz_in;
  }

      
  my %mtz_dump = Phenix::MTZ_Dump($mtz_out);
  if ( not $mtz_dump{'fobs'} eq $fobs ) {
    Utility::Report_Error ( 
      "Phenix could not convert iobs to fobs", 
      [$mtz_dump{'dump'}, $mtz_dump{'errors'}],
    );
    return;
  }


  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{MTZ_IOBS} = $CFG->{MTZ};
    $CFG->{MTZ} = \%mtz_dump;
  }

  if ( wantarray ) { return  %mtz_dump}
  else             { return \%mtz_dump}

}


######################################################## GENERATE R FREE FLAGS

sub Generate_R_Free { 
  my $CFG       = $_[0];
  my $mtz_in;
  my $mtz_out   = $_[1];
  my $root      = $_[2];

  my $overwrite = 0;

  if ( Is_Structure_Config($CFG) ) { 
    $root = $CFG->{MTZ}{root};
  } else { 
    $mtz_in = $CFG;
  }

  if ( (not defined $mtz_out) or ($mtz_out eq $mtz_in) ) { 
    $overwrite = 1;
    $mtz_out = "$mtz_in.temp";
  }



  my $command = "phenix.reflection_file_converter "
              . "$mtz_in "
              #. "--symmetry=$pdb "
              . "--label=$root "
              . "--non-anomalous "
              . "--mtz-root-label=$root "
              . "--generate-r_free_flags "
              . "--r-free-flags-fraction=0.05 "
              . "--r-free-flags-max-free=100000 "
              . "--r-free-label=$DEFAULT_FREE_LABEL "
              . "--mtz $mtz_out "
              ;

  my %result = Utility::Try_Command($command, "Phenix - Generate R-Free-Flags");

  if ( $result{'value'} ) { return %result }

  if ( not -e $mtz_out ) { 
    Utility::Report_Error ( "Phenix did not generate MTZ file '$mtz_out'.", "", "" );
    return;
  }

  
  if ( $overwrite ) { 
    rename $mtz_out, $mtz_in 
      or confess("Phenix::IOBS_to_FOBS:: Couldn't rename '$mtz_out' to '$mtz_in'. $!");
    $mtz_out = $mtz_in;
  }

      
  my %mtz_dump = Phenix::MTZ_Dump($mtz_out);
  if ( not $mtz_dump{'root'} eq $root ) {
    Utility::Report_Error ( 
      "Phenix could not convert iobs to fobs", 
      [$mtz_dump{'dump'}, $mtz_dump{'errors'}],
    );
    return;
  }

  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{MTZ_FREE} = $CFG->{MTZ};
    $CFG->{MTZ} = \%mtz_dump;
  }

  if ( wantarray ) { return  %mtz_dump}
  else             { return \%mtz_dump}

}


################################################################## PHENIX.MAPS
# Create electron density maps with intensities/amplitudes and a pdb file.
# If more precision is needed, use the ALL_MAPS file below for reference.
##############################################################################

my $ALL_MAP_PARAMS = <<__ALL_MAPS__
 maps {
   input {
     pdb_file_name = None
     reflection_data {
       file_name = None
       labels = None
       high_resolution = None
       low_resolution = None
       outliers_rejection = True
       french_wilson_scale = True
       french_wilson {
         max_bins = 60
         min_bin_size = 40
       }
       sigma_fobs_rejection_criterion = 0.0
       sigma_iobs_rejection_criterion = 0.0
       r_free_flags {
         file_name = None
         label = None
         test_flag_value = None
         ignore_r_free_flags = False
       }
     }
   }
   output {
     directory = None
     prefix = None
     title = None
     fmodel_data_file_format = mtz cns
   }
   scattering_table = wk1995 it1992 *n_gaussian neutron
   bulk_solvent_correction = True
   #apply_back_trace_of_b_cart = False
   anisotropic_scaling = True
   skip_twin_detection = False
   omit {
     method = *simple
     selection = None
   }

##### The basic F/SigFs
   map_coefficients {
     map_type = mFo
     format = *mtz phs
     mtz_label_amplitudes = F-model
     mtz_label_phases = PHIF-model
     kicked = False
     fill_missing_f_obs = False
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }

##### Coot Readable
   map_coefficients {
     map_type = 2mFo-DFc
     format = *mtz phs
     mtz_label_amplitudes = 2FOFCWT
     mtz_label_phases = PH2FOFCWT
     kicked = False
     fill_missing_f_obs = False
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }
   map_coefficients {
     map_type = mFo-DFc
     format = *mtz phs
     mtz_label_amplitudes = FOFCWT
     mtz_label_phases = PHFOFCWT
     kicked = False
     fill_missing_f_obs = False
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }

##### Default phenix.maps output
   map_coefficients {
     map_type = 2mFo-DFc
     format = *mtz phs
     mtz_label_amplitudes = 2mFoDFc
     mtz_label_phases = P2mFoDFc
     kicked = False
     fill_missing_f_obs = False
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }
   map_coefficients {
     map_type = 2mFo-DFc
     format = *mtz phs
     mtz_label_amplitudes = 2mFoDFc_fill
     mtz_label_phases = P2mFoDFc_fill
     kicked = False
     fill_missing_f_obs = True
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }
   map_coefficients {
     map_type = mFo-DFc
     format = *mtz phs
     mtz_label_amplitudes = mFoDFc
     mtz_label_phases = PmFoDFc
     kicked = False
     fill_missing_f_obs = False
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }
   map_coefficients {
     map_type = anomalous
     format = *mtz phs
     mtz_label_amplitudes = ANOM
     mtz_label_phases = PANOM
     kicked = False
     fill_missing_f_obs = False
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }
   map {
     map_type = 2mFo-DFc
     format = xplor *ccp4
     file_name = None
     kicked = False
     fill_missing_f_obs = False
     grid_resolution_factor = 1/4.
     region = *selection cell
     atom_selection = None
     atom_selection_buffer = 3
     sharpening = False
     sharpening_b_factor = None
     exclude_free_r_reflections = False
     isotropize = True
   }
 }
__ALL_MAPS__
;


#========================================================================= Maps

sub Maps { 

  print "----- Maps\n";

  my $CFG = $_[0];
  my ($pdb_path, $pdb_dir, $pdb_file); 
  my ($mtz_path, $mtz_dir, $mtz_file); 
  my ($map_path, $map_dir, $map_file); 
  my ($mpz_path, $mpz_dir, $mpz_file); # mtz version for coot

     $mtz_file = $_[1];
  my $fobs     = $_[2];
  my $prefix   = $_[3];
  my $flags    = $_[4];
     $map_dir  = $_[5];
  my $free     = $_[6];

  my $r_free   = $DEFAULT_FREE_LABEL;
 
  if ( Is_Structure_Config($CFG) ) { 
    $pdb_file = $CFG->{FILES}{PDB};
    $mtz_file = $CFG->{FILES}{MTZ};
    $pdb_path = $CFG->{PATHS}{PDB};
    $mtz_path = $CFG->{PATHS}{MTZ};
    $prefix   = $CFG->{FILES}{BASE};
    $fobs     = $CFG->{DUMP }{fobs} || $CFG->{DUMP}{iobs};
    $free     = $CFG->{DUMP }{free};
    $flags    = $CFG->{FLAGS};
    $map_dir  = $CFG->{SUBDIRS}{MAPS} || $CFG->{DIRS}{MAPS};
  } else { 
    $pdb_path = $CFG;
    ($pdb_path, $pdb_dir, $pdb_file) = Utility::Split_Path($pdb_path);
    ($mtz_path, $mtz_dir, $mtz_file) = Utility::Split_Path($mtz_path);
  }
  
 
  #----------------------------------------------------------------- Validation
  $map_dir ||= ".";
  if ( not -d $map_dir ) { mkdir $map_dir or confess "Couldn't mkdir holes dir, '$map_dir'. $!" }

  #print "Thinking map dir is $map_dir.\n";
  $map_dir = Cwd::abs_path($map_dir);
  my $cwd = getcwd();
  $map_file = "$prefix\_2mFo-DFc_map.ccp4";
  $mpz_file = "$prefix\_map_coeffs.mtz";
  $map_path = "$map_dir/$map_file";
  $mpz_path = "$map_dir/$mpz_file";
  #print "Making mappath $map_path in dir $map_dir\n";

  #### Check if we're already done
  if ( Is_Structure_Config($CFG) ) { 
    if ( not exists $CFG->{PATHS}{MAP} ) { 
      Utility::Add_Path_to_Config($CFG, 'MAP', $map_path);
      Utility::Add_Path_to_Config($CFG, 'MPZ', $mpz_path);
    } elsif ( $CFG->{PATHS}{MAP} ne $map_path ) { 
      confess "Preexisting MAP Path differs from that expected:\n"
        . "\t$map_path\n"
        . "\t$CFG->{PATHS}{MAP}\n"
        ;
    }
  }

  #--------------------------------------------------------- Holes Already Done
  
  my @input_files  = ( $pdb_path, $mtz_path );
  my @output_files = ( $map_path, $mpz_path );
  my ($errors_in, $errors_out) = Utility::Verify_Analysis_Files ( 
    "Maps", \@input_files, \@output_files, \%FLAGS_VERIFY_PRE );
  
  if ( (not $flags->{'overwrite'})
      and scalar @$errors_in == 0 and scalar @$errors_out == 0 ) { 
    print "Maps complete for $pdb_path\n"
        . "                  $mtz_path.\n"
        . "                  $map_path.\n"
        . "                  $mpz_path.\n"
        if $VERBOSE >= 1;
    return
  }


  #-------------------------------------------------------------- Execute Maps

  #----- Ensure map parameters were parsable
  if ( not defined $fobs ) {
    Utility::Report_Error ( "Maps fobs parameter name is not defined for mtz $mtz_path", "", { fatal => 1 } );
    return;
  } elsif ( not defined $free ) { 
    Utility::Report_Error ( "Maps r-free-flags parameter name is not defined for mtz $mtz_path", "", { fatal => 1 } );
    return;
  }

  print "Running Maps because:\n"
      . (join "", map { "  $_\n" } @$errors_in, @$errors_out)
      . "\n";
 
  

  #----- Write the Maps Parameter file it's missing
  chdir $map_dir or confess "Couldn't chdir to map directory, '$map_dir'. $!";
  my $params = "maps.standard.params";
  if ( not -e $params ) { 
    print "Writing map parameters file, '$params'.\n";
    Utility::Write_File($params, $ALL_MAP_PARAMS);
  }


  #----- Build Command
  my $command = "phenix.maps "
              . "maps.input.pdb_file_name=$pdb_path "
              . "maps.input.reflection_data.labels=$fobs "
              . "maps.input.reflection_data.file_name=$mtz_path "
              . "maps.input.reflection_data.r_free_flags.label=$free "
              . "maps.output.prefix=$prefix "
              #. "maps.output.title=$title "
              . "$params "
              ;

 
  #------ Execute and Verify 
  my %result = Utility::Try_Command (
                  $command,
                  "Creating electron density maps for $pdb_file / $map_file.",
                  { %FLAGS_VERIFY_POST, %$flags, inputs=>\@input_files, outputs=>\@output_files },
               );
 
  chdir $cwd or confess "Couldn't chdir back to original directory, '$cwd'. $!";

  if ( wantarray ) { return  %result }
  else             { return \%result }

}



##############################################################################
##############################################################################
#                               HOLES
##############################################################################
##############################################################################


sub Rosetta_Holes { 

  print "----- Holes\n";

  my $CFG     = $_[0];
  my $pdb;
  my $title   = $_[1];
  my $dir     = $_[2];
  my $flags   = $_[3] || {};


  #---- Detect CFG vs. parameters
  if ( Is_Structure_Config($CFG) ) { 
    $pdb    = $CFG->{PATHS}{PDB};
    $title  = $CFG->{FILES}{BASE};
    $dir    = $CFG->{DIRS }{HOLES};
    $flags  = $CFG->{FLAGS};
  } else { 
    $pdb = $CFG;
  }
  $dir ||= $DEFAULT_PATHS{DIRS}{HOLES};

  #------------------------------------------------------ Files and Directories
  if ( not defined $dir ) { confess "Phenix::Rosetta_Holes !! Undefined holes dir" }
  if ( not defined $title ) { confess "Phenix::Rosetta_Holes !! Undefined file base" }
  Utility::Make_Directories($dir);
  Utility::Add_Path_to_Config($CFG, 'HOLES', $dir );
  Utility::Populate_Paths( 
                $CFG, 
                $dir,
                { 
                  HOLES_PDB => "$title.holes.pdb", # the holes pdb, not the main one
                  HOLES_LOG => "$title.holes.log",
                  HOLES_OUT => "$title.holes.out",
                  HOLES_CA  => "$title.holes.ca.dat",
                  HOLES_RES => "$title.holes.res.dat",
                },
  );
  


  #--------------------------------------------------------- Holes Already Done
  my ($errors_in, $errors_out) = Utility::Verify_Analysis_Files ( 
    "Rosetta Holes",
    [ $pdb ],
    [ map { $CFG->{PATHS}{$_} } qw(HOLES_PDB HOLES_LOG HOLES_OUT HOLES_CA HOLES_RES) ],
    \%FLAGS_VERIFY_PRE,
  );
  
  if ( scalar @$errors_in == 0 and scalar @$errors_out == 0 ) { 
    print "Holes complete for $pdb.\n" if $VERBOSE >= 1;
  }

  #-------------------------------------------------------------- Execute Holes
  else {
    print "Running Rosetta Holes because:\n"
        . (join "", map { "  $_\n" } @$errors_in, @$errors_out)
        . "\n";
    my $cwd = getcwd();
    chdir $dir or confess "Couldn't chdir to $dir. $!";
    my %result = Utility::Try_Command(
                  "$Phenix::PROGRAMS{HOLES} $pdb $title",
                  "Calculating Rosetta Holes for $pdb",
                );
    chdir $cwd or confess "Coudln't chdir to $cwd. $!";

    #------------------------------------------------- Detect and Report Errors
    
    #---- Command failed to execute
    if ( $result{command_error} ) { 
      # Unkown atom
      if ( $result{err} =~ /unknown atom_name:\s+(.*)/ ) { 
        print "Rosetta didn't understand the atom $1.\n";
        return;
      }

      # Other unknown failure
      Utility::Report_Error ( 
        "Rosetta holes failed for $pdb",
        \%result,
        $flags,
      );
      return;
    }
   
    #---- Command seemed to run, but output data is missing 
    if ( not -e $CFG->{PATHS}{HOLES_PDB} ) { 
      Utility::Report_Error (
        "Rosetta holes failed for $pdb (invisibly, missing '$CFG->{PATHS}{HOLES_PDB}')",
        \%result,
        $flags,
      );
      return;
    }

  }

  #--------------------------------------------------------------- Populate CFG
  if ( Is_Structure_Config($CFG) ) { 
    Parse_Holes_Data($CFG);
    return;
  }

  return Parse_Holes_Data($title, $dir);
  
}


#============================================================= Parse Holes Data

sub Parse_Holes_Data { 
  
  my $CFG     = $_[0]; 
  my $title;
  my $dir     = $_[1];
  
  my ($path_log, $path_res);

  if ( Is_Structure_Config($CFG) ) { 
    $dir      = $CFG->{DIRS }{HOLES} || $CFG->{SUBDIRS}{HOLES};
    $title    = $CFG->{FILES}{BASE };
    $path_log = $CFG->{PATHS}{HOLES_LOG};
    $path_res = $CFG->{PATHS}{HOLES_RES};
  }

  $dir = Cwd::abs_path($dir);
  $path_log = "$dir/$title.holes.log"       if not defined $path_log;
  $path_res = "$dir/$title.holes.res.dat"   if not defined $path_res;
  if ( not defined $dir ) { confess "Rosetta Holes !! Dir not defined, '$dir'" }

  my @lines = Utility::Read_File($path_log);
  my ($score_overall, $rsln_overall) = (split /\t/, $lines[1])[2..3];
  my ($data_res, $header_res) = Utility::Read_Hash_Array($path_res);

    
  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{HOLES}{overall_score} = $score_overall;
    $CFG->{HOLES}{overall_rsln } = $rsln_overall;
    foreach my $row ( @$data_res ) { 
      my %row = %$row;
      my ( $chain, $resi, $resn, $atom, $alt, $score, $depth ) 
        = @row{qw(CHAIN RESI RESN ATOM ALT PDB_B PDB_Q)};
      $CFG->{HOLES}{SCORE}{$chain}[$resi]{$atom}{_} = $score;
      $CFG->{HOLES}{DEPTH}{$chain}[$resi]{$atom}{_} = $depth;
    }
   
    Collapse_Atom_Properties($CFG->{HOLES}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'SCORE', 'avg');
    Collapse_Atom_Properties($CFG->{HOLES}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DEPTH', 'avg');


    Register_Analyses(
      $CFG, 
      [qw(ATOM ANGLE RESIDUE)], 
      'HOLES', 
      [ qw(SCORE DEPTH) ],
      ['PROTEIN'],
    );
  }

  return ($score_overall, $rsln_overall, $data_res);
}



##############################################################################
#                          DSSP SECONDARY STRUCTURE
##############################################################################

sub DSSP { 
  my $pdb = $_[0];
  if ( is_structure_cfg($pdb) ) { 
  }
}





##############################################################################
#                              PISA INTERFACES
##############################################################################

sub Pisa { 

  print "----- Pisa\n";

  my $CFG     = $_[0];
  my $pdb;
  my $name    = $_[1];
  my $dir     = $_[2];
  my $chains  = $_[3];
  my $flags   = $_[4];

  #---- Detect CFG vs. parameters
  if ( Is_Structure_Config($CFG) ) { 
    $pdb    = $CFG->{PATHS}{PDB};
    $name   = $CFG->{FILES}{BASE};
    $dir    = $CFG->{DIRS }{PISA} || $CFG->{SUBDIRS}{PISA};
    if ( defined $CFG->{PDB}{PROTEIN} ) { 
      $chains = [keys %{$CFG->{PDB}{PROTEIN}}];
    } else { 
      $chains = $CFG->{PDB  }{CHAINS};
    }
    $flags  = $CFG->{FLAGS};
  } else { 
    $chains ||= ['A'];
  }


  #---- Set Default Parameters
  if ( not -e $pdb ) { confess "Cannot find pdb file, '$pdb'.\n" }
  if ( not defined $name ) { $name = $pdb; $name =~ s/\.pdb//; $name =~ s/.*\///; }
  if ( not defined $dir ) { $dir = "$ENV{REFINE_DIR}/PISA/" };
  if ( not -d $dir ) { mkdir $dir or confess "Couldn't mkdir holes dir, '$dir'. $!" }
  $dir = Cwd::abs_path($dir);
  my $cwd = getcwd();

  #---- Correct for passing 'A' instead of ['A']
  if ( not (ref $chains eq 'ARRAY' ) ) { $chains = [$chains] }

  #---- Loop over chains
  foreach my $chain ( @$chains ) { 

    my $target = "$dir/$name.interfaces.txt";
    
    my ($errors_in, $errors_out) = Utility::Verify_Analysis_Files ( 
        "PISA",
        [ $pdb ],
        [ $target ],
        \%FLAGS_VERIFY_PRE,
      );

    #------------------------------------------------------------- Already Done
    if ( scalar @$errors_in == 0 and scalar @$errors_out == 0 ) { 
      print "PISA complete for $pdb.\n" if $VERBOSE >= 1;

    #-------------------------------------------------------------- Execute Pisa
    } else {
      print "Running PISA because:\n"
          . (join "", map { "  $_\n" } @$errors_in, @$errors_out)
          . "\n";

      chdir $dir or confess "Couldn't chdir to '$dir'. $!";
      my %result = Utility::Try_Command ( 
                      "$PROGRAMS{PISA} $pdb $chain ",
                      "Executed PISA - $pdb",
                      $flags,
                    );
      chdir $cwd or confess "Couldn't chdir to '$dir'. $!";

      #---- Detect Errors
      if ( $result{'value'} ) { 
        Utility::Report_Error ( 
          "PISA failed for $pdb",
          \%result,
          $flags,
        );
        return;
      } elsif ( not -e $target ) { 
        Utility::Report_Error (
          "PISA failed for $pdb (invisibly, no target '$target')",
          \%result,
          $flags,
        );
        return;
      }
    }
  }

  #---- Parse and return results
  my %DATA = Parse_Pisa_Data($CFG, $name, $dir, $chains, $flags);
  if ( wantarray ) { return  %DATA }
  else             { return \%DATA }


}


############################################################## Parse Pisa Data

sub Parse_Pisa_Data { 

  my $CFG    = $_[0];
  my $pdb;
  my $name   = $_[1];
  my $dir    = $_[2];
  my $chains = $_[3];

  if ( Is_Structure_Config($CFG) ) { 
    $pdb    = $CFG->{PATHS}{PDB};
    $name   = $CFG->{FILES}{BASE};    
    if ( defined $CFG->{PDB}{PROTEIN} ) { 
      $chains = [keys %{$CFG->{PDB}{PROTEIN}}];
    } else { 
      $chains = $CFG->{PDB  }{CHAINS};
    }
    $dir    = $CFG->{DIRS }{PISA};
  } else { 
    $chains ||= ['A'];
  }


  my $int = "$dir/$name.interfaces.pisa";
  my $ass = "$dir/$name.assemblies.pisa";
 
  my %PISA = ();
  foreach my $chain (@$chains) {
    my $rel_file = "$dir/$name.$chain.rel.asa";
    my $abs_file = "$dir/$name.$chain.abs.asa";
    my $int_file = "$dir/$name.interfaces.txt";

    my ( $rel_data, $rel_head ) = Utility::Read_Hash_Array($rel_file);
    my ( $abs_data, $abs_head ) = Utility::Read_Hash_Array($abs_file);
    my ( $int_data, $int_head ) = Utility::Read_Hash_Array($int_file);

    foreach my $i ( 0..$#$rel_data ) { 
      my $resi = $rel_data->[$i]{RESI} || $rel_data->[$i]{NUM};
      my $resn = $rel_data->[$i]{RESN} || $rel_data->[$i]{NAME};
      foreach my $h ( 2..$#$rel_head ) { 
        my $head = uc $rel_head->[$h];
        $PISA{"REL_$head"}{$chain}[$resi]{AA}{_} = $rel_data->[$i]{$head};
      }
    }
    
    foreach my $i ( 0..$#$abs_data ) { 
      my $resi = $abs_data->[$i]{RESI} || $abs_data->[$i]{NUM};
      my $resn = $abs_data->[$i]{RESN} || $abs_data->[$i]{NAME};
      foreach my $h ( 2..$#$abs_head ) { 
        my $head = uc $abs_head->[$h];
        $PISA{"ABS_$head"}{$chain}[$resi]{AA}{_} = $abs_data->[$i]{$head};
      }
    }
    


    
    # Move one of the CLASS properties to just CLASS, since they're identical
    $PISA{CLASS} = $PISA{REL_CLASS};
    Assign_Default_Altconfs($PISA{CLASS});
    Register_Analyses(
      $CFG,
      [qw(RESIDUE)],
      'PISA',
      [ 'CLASS' ],
      [ 'PROTEIN' ],
    );

    # Register properties with the data structure.
    if ( Is_Structure_Config($CFG) ) { 
      foreach my $h ( 2..$#$rel_head ) { 
        my $head = uc $rel_head->[$h];
        next if $head eq 'CLASS';
        Assign_Default_Altconfs($PISA{"REL_$head"});
        Assign_Default_Altconfs($PISA{"ABS_$head"});
        Register_Analyses(
          $CFG, 
          [qw(RESIDUE)], 
          'PISA', 
          [ "REL_$head", "ABS_$head" ],
          ['PROTEIN'],
        );

        # Add formats to global format hash (in case of unknown ligands)
        $META_PROPERTIES{PISA}{"REL_$head"}{F} = $META_PROPERTIES{PISA}{REL}{F};
        $META_PROPERTIES{PISA}{"ABS_$head"}{F} = $META_PROPERTIES{PISA}{ABS}{F};
      }
    }
  }

  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{PISA} = \%PISA;
  }

  if ( wantarray ) { return  %PISA }
  else             { return \%PISA }

}

##############################################################################
#                            COOT ATOM PROPERTIES
#-----------------------------------------------------------------------------
# /Library/Coot/bin/coot-real --python 
#                             --pdb IPA_100K_C.pdb  
#                             --auto MAPS/IPA_100K_C_test_map_coeffs.mtz 
#                             --no-guano 
#                             --no-graphics 
#                             -s ~/lab/code/pymol/coot_atom_density.py
##############################################################################

sub Coot_Properties {

  print "----- Coot\n";
  my $CFG     = $_[0];
  my $pdb_path;
  my $mpz_path;
  my $base_path;
  my $flags;

  if ( Is_Structure_Config($CFG) ) { 
    $pdb_path    = $CFG->{PATHS}{PDB};
    $mpz_path    = $CFG->{PATHS}{MPZ};
    $base_path   = $CFG->{PATHS}{BASE};
    $flags       = $CFG->{FLAGS};
  }


  #==== Determine if analysis already run
  my $atm_path = "$base_path.coot.atoms.dat";
  my $res_path = "$base_path.coot.residues.dat";
  my @input_files  = ( $pdb_path, $mpz_path );
  my @output_files = ( $atm_path, $res_path );

  my ($errors_in, $errors_out) = Utility::Verify_Analysis_Files ( 
      "Coot", \@input_files, \@output_files, \%FLAGS_VERIFY_PRE );

  #------------------------------------------------------------- Already Done
  if ( scalar @$errors_in == 0 and scalar @$errors_out == 0 ) { 
    print "Coot complete for $pdb_path.\n" if $VERBOSE >= 1;

  #-------------------------------------------------------------- Execute Pisa
  } else {
    print "Running Coot because:\n"
        . (join "", map { "  $_\n" } @$errors_in, @$errors_out)
        . "\n";


    my $coot_cmd = "$Phenix::PROGRAMS{COOT} "
                 . "--python "
                 . "--pdb $pdb_path "
                 . "--auto $mpz_path "
                 . "--no-guano "
                 . "--no-graphics "
                 . "-s $Phenix::PROGRAMS{COOT_PROPS} "
                 . "2>&1"
                 ;  
  
    my %result = Utility::Try_Command ( 
                    $coot_cmd,
                    "Execute coot for atom properties",
                    { %FLAGS_VERIFY_POST, %$flags, inputs=>\@input_files, outputs=>\@output_files },
                  );
  }
  
  if ( not exists $CFG->{PATHS}{COOT_ATM} ) { 
    Utility::Add_Path_to_Config($CFG, 'COOT_ATM', $atm_path);
    Utility::Add_Path_to_Config($CFG, 'COOT_RES', $res_path);
  }
  Parse_Coot_Properties($CFG);
}


sub Parse_Coot_Properties {

  my $CFG = $_[0];

  my ( $data_atm, $head_atm ) = Utility::Read_Hash_Array($CFG->{PATHS}{COOT_ATM});
  my ( $data_res, $head_res ) = Utility::Read_Hash_Array($CFG->{PATHS}{COOT_RES});

  foreach my $row (@$data_atm) { 
    $CFG->{COOT}{ISOr   }{$row->{CHAIN}}[$row->{RESI}]{$row->{ATOM}}{$row->{ALT}} = $row->{ISOr };
    $CFG->{COOT}{ISOa   }{$row->{CHAIN}}[$row->{RESI}]{$row->{ATOM}}{$row->{ALT}} = $row->{ISOa };
    $CFG->{COOT}{DIFr   }{$row->{CHAIN}}[$row->{RESI}]{$row->{ATOM}}{$row->{ALT}} = $row->{DIFr };
    $CFG->{COOT}{DIFa   }{$row->{CHAIN}}[$row->{RESI}]{$row->{ATOM}}{$row->{ALT}} = $row->{DIFa };
    $CFG->{COOT}{ANISO  }{$row->{CHAIN}}[$row->{RESI}]{$row->{ATOM}}{$row->{ALT}} = $row->{ANISO};
  }
  foreach my $row (@$data_res) { 
    $CFG->{COOT}{ROTAMER}{$row->{CHAIN}}[$row->{RESI}]{AA}{$row->{ALT}} = $row->{ROTAMER};
  }

  Collapse_Atom_Properties($CFG->{COOT}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'ISOr',    'avg');
  Collapse_Atom_Properties($CFG->{COOT}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'ISOa',    'avg');
  Collapse_Atom_Properties($CFG->{COOT}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DIFr',    'avg');
  Collapse_Atom_Properties($CFG->{COOT}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DIFa',    'avg');
  Collapse_Atom_Properties($CFG->{COOT}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'ANISO',   'avg');
  
  Assign_Default_Altconfs($CFG->{COOT}{ROTAMER}, $CFG->{PDB}{ALTS});

  Register_Analyses(
    $CFG, 
    [qw(ATOM ANGLE RESIDUE)], 
    'COOT', 
    [ qw(ISOr DIFr ISOa DIFa ANISO) ],
    ['PROTEIN', 'WATER', 'LIGAND'],
  ); 

  Register_Analyses(
    $CFG, 
    [qw(RESIDUE)], 
    'COOT', 
    [ qw(ROTAMER) ],
    ['PROTEIN'],
  ); 

}


##############################################################################
#                           PyMOL Water Properties
#-----------------------------------------------------------------------------
# pymol -c            # command line mode
#       -q            # quiet launch
#       -r <file>     # script to execute
#       -d <command>  # command to execute
#       <pdb>         # pdb file to load
#       <map>         # map file to load
#
##############################################################################

sub PyMOL_Water_Properties {

  print "----- Pymol Water Properties\n";
  my $CFG     = $_[0];
  my $pdb_path;
  my $base_path;
  my $flags;
  my $water_dir;

  if ( Is_Structure_Config($CFG) ) { 
    $pdb_path   = $CFG->{PATHS}{PDB};
    $base_path  = $CFG->{PATHS}{BASE};
    $water_dir  = $CFG->{DIRS}{WATER} || "WATER";
    $flags      = $CFG->{FLAGS};
  }


  $water_dir = Cwd::abs_path($water_dir);
  my $cwd = getcwd();
  if ( not -e $water_dir ) { 
    mkdir $water_dir or confess "Couldn't mkdir '$water_dir'. $!";
  }
  chdir $water_dir or confess "Couldn't mkdir '$water_dir'. $!";


  #==== Determine if analysis already run
  my $global_path   = "$base_path.water.global.dat";
  my $shells_path   = "$base_path.water.shells.dat";
  my $residues_path = "$base_path.water.residues.dat";
  my @input_files   = ( $pdb_path );
  my @output_files  = ( $global_path, $shells_path, $residues_path );

  my ($errors_in, $errors_out) = Utility::Verify_Analysis_Files ( 
      "PyMOL Waters", \@input_files, \@output_files, \%FLAGS_VERIFY_PRE );

  #------------------------------------------------------------- Already Done
  if ( scalar @$errors_in == 0 and scalar @$errors_out == 0 ) { 
    print "PyMOL Waters complete for $pdb_path.\n" if $VERBOSE >= 1;

  #-------------------------------------------------------------- Execute Pisa
  } else {
    print "Running PyMOL Waters because:\n"
        . (join "", map { "  $_\n" } @$errors_in, @$errors_out)
        . "\n";


    my $pymol_cmd = "$Phenix::PROGRAMS{PYMOL} "
                  . "-c "
                  . "-q "
                  . "-r $Phenix::PROGRAMS{PYMOL_UTILITY} "
                  . "$pdb_path "
                  . "-r $Phenix::PROGRAMS{PYMOL_WATERS} "
                  . "2>&1"
                  ;  
  
    my %result = Utility::Try_Command ( 
                    $pymol_cmd,
                    "Executing PyMOL Water Analysis",
                    { %FLAGS_VERIFY_POST, %$flags, inputs=>\@input_files, outputs=>\@output_files },
                  );
  }
  
  chdir $cwd or die "Couldn't chdir '$cwd'. $!";
  
  #Utility::Add_Path_to_Config($CFG, 'WATER_GLOBAL',   $global_path);
  #Utility::Add_Path_to_Config($CFG, 'WATER_SHELLS',   $shells_path);
  #Utility::Add_Path_to_Config($CFG, 'WATER_RESIDUES', $residues_path);
  #Parse_PyMOL_Water_Properties($CFG);
}


sub Parse_PyMOL_Water_Properties {

  my $CFG = $_[0];

  my ( $data_shl, $head_shl ) = Utility::Read_Hash_Array($CFG->{PATHS}{WATER_SHELLS});
  my ( $data_res, $head_res ) = Utility::Read_Hash_Array($CFG->{PATHS}{WATER_RESIDUES});

  return;

  foreach my $row (@$data_shl) { 
    $CFG->{WATERS}{"SHELL_$row->{SHELL}"} = $row->{$CFG->{TITLE}};
  }

  foreach my $row (@$data_shl) { 
    $CFG->{COOT}{ISOr   }{$row->{CHAIN}}[$row->{RESI}]{$row->{ATOM}}{$row->{ALT}} = $row->{ISOr };
  }
  foreach my $row (@$data_res) { 
    $CFG->{COOT}{ROTAMER}{$row->{CHAIN}}[$row->{RESI}]{AA}{$row->{ALT}} = $row->{ROTAMER};
  }

  Collapse_Atom_Properties($CFG->{COOT}, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'ISOr',    'avg');
  Assign_Default_Altconfs($CFG->{COOT}{ROTAMER}, $CFG->{PDB}{ALTS});

  Register_Analyses(
    $CFG, 
    [qw(ATOM ANGLE RESIDUE)], 
    'COOT', 
    [ qw(ISOr DIFr ISOa DIFa ANISO) ],
    ['PROTEIN', 'WATER', 'LIGAND'],
  ); 

  Register_Analyses(
    $CFG, 
    [qw(RESIDUE)], 
    'COOT', 
    [ qw(ROTAMER) ],
    ['PROTEIN'],
  ); 

}




##############################################################################
#                              RINGER INTERFACE
##############################################################################

my $CHIMERA_HOME = "/programs/x86_64-linux/chimera/1.6.1";

sub Ringer { 

  print "----- Ringer\n";
  my $CFG     = $_[0];
  my $pdb;
  my $map     = $_[1];
  my $dir     = $_[2] || 'ringer';
  my $config  = $_[3] || {};
  my $flags   = $_[4] || {};

  if ( Is_Structure_Config($CFG) ) { 
    $pdb    = $CFG->{PATHS}{PDB};
    $map    = $CFG->{PATHS}{MAP};
    $dir    = $CFG->{SUBDIRS}{RINGER};
    $config = $CFG->{PARAMS}{RINGER};
    $flags  = $CFG->{FLAGS};
    if ( not defined $dir ) { $dir = "ringer_$CFG->{FILES}{BASE}" }
  }

  #--------------------------------------------------------------- Verify Files
  if ( not defined $pdb ) { confess("Phenix::Ringer not passed a pdb file name.\n") }
  if ( not defined $map ) { confess("Phenix::Ringer not passed a map file name.\n") }
  if ( not -e $map      ) { confess "Phenix::Ringer The map ('$map') doesn't exist.\n"; }
  if ( not defined $dir ) { $dir = "ringer" }

  if ( not -d $dir ) { mkdir $dir or confess "Couldn't mkdir holes dir, '$dir'. $!" }
  $dir = Cwd::abs_path($dir);
  my $cwd = getcwd();

  # Prefer .ccp4 maps
  if ( $map !~ /.ccp4/ ) { 
    warn "  This map type probably won't work ($map). Use a .ccp4 map. Will give it a shot though.\n";
  }

  # Works best with Chimera 1.6.1 currently (2012 07 11)
  if ( $ENV{CHIMERA_HOME} =~ /\/programs/ ) { 
    if ( not $ENV{CHIMERA_HOME} =~ /1.6.1/ ) { 
      $ENV{CHIMERA_HOME} = $CHIMERA_HOME;
    }
  }
  
  #------------------------------------------------------- Ringer Configuration
  my %config = (
    #configuration           => 'standard',
    map_type                => 'sigma',
    sigma_map_scaling       => 'mean',
    skip_multi_conf         => 'off',
    atom_sample_type        => 'dynamic',
    chi_sample_degree       => '1',
    lower_sigma_cutoff      => '0.3',
    write_alt_conformers    => 'off',
    write_sigma_plot        => 'on',
    write_peak_list         => 'on',
    write_ring_vectors      => 'off',
    write_plot              => 'on',
    write_rings             => 'on',
    write_chi2chi1          => 'on',
    verbose_outfile_prefix  => 'verbose',
  );
  # Apply any user customization to main config.
  map { $config{$_} = $config->{$_} } keys %$config;



  my $plot_file = "$config{verbose_outfile_prefix}.signal_plot.txt";
  my $plot_path = "$dir/$plot_file";

  my @input_files = ( $pdb, $map );
  my @output_files = ( $plot_path );
  my ($errors_in, $errors_out) = Utility::Verify_Analysis_Files ( 
      "Ringer", \@input_files, \@output_files, \%FLAGS_VERIFY_PRE );
    
  #------------------------------------------------------------- Already Done
  if ( scalar @$errors_in == 0 and scalar @$errors_out == 0 ) { 
    print "Ringer complete for $pdb.\n" if $VERBOSE >= 1;

  #-------------------------------------------------------------- Execute Pisa
  } else {
    print "Running Ringer because:\n"
        . (join "", map { "  $_\n" } @$errors_in, @$errors_out)
        . "\n";


    #------------------------------------------------- Change to Ringer Directory
    if ( not -e $dir ) { 
      mkdir $dir or confess "Couldn't mkdir $dir. $!";
    }
    chdir $dir or confess "Couldn't chdir $dir. $!";

    #--------------------------------------------------- Write Ringer Config File
    my $ringer_in = <<__RINGER_IN__
pdb_name                $pdb
map_name                $map
map_type                $config{'map_type'}
sigma_map_scaling       $config{'sigma_map_scaling'}
skip_multi_conf         $config{'skip_multi_conf'}
atom_sample_type        $config{'atom_sample_type'}
chi_sample_degree       $config{'chi_sample_degree'}
lower_sigma_cutoff      $config{'lower_sigma_cutoff'}
write_alt_conformers    $config{'write_alt_conformers'}
write_sigma_plot        $config{'write_sigma_plot'}
write_peak_list         $config{'write_peak_list'}
write_ring_vectors      $config{'write_ring_vectors'}
write_plot              $config{'write_plot'}
write_rings             $config{'write_rings'}
write_chi2chi1          $config{'write_chi2chi1'}
verbose_outfile_prefix  $config{'verbose_outfile_prefix'}
__RINGER_IN__
;

    my $ringer_file = "ringer.in";
    Utility::Write_File($ringer_file, $ringer_in, 1);
   

    #------------------------------------------------------------- Execute Ringer
    my %result = Utility::Try_Command ( 
                    "$PROGRAMS{RINGER} -i $ringer_file 2>&1",
                    "Execute standard ringer on $pdb / $map",
                    { %FLAGS_VERIFY_POST, %$flags, inputs=>\@input_files, outputs=>\@output_files },
                 );
    chdir $cwd or confess "Couldn't chdir to cwd '$cwd'. $!";


    #----------------------------------------------------------- Check for Errors
    if ( $result{'value'} ) { 
      Utility::Report_Error ( 
        "Ringer failed for $pdb / $map",
        \%result,
        $flags,
      );
      return;
    } elsif ( not -e $plot_path ) { 
      Utility::Report_Error (
        "Ringer failed for $pdb / $map (invisibly, no plot '$plot_file')",
        \%result,
        $flags,
      );
      return;
    }
  }
 
  #========================================================== Parse the Results 
  my %RINGER = ();
  if ( Is_Structure_Config($CFG) ) { 
    $CFG->{DIRS}{RINGER} = $dir;
    $CFG->{FILES}{RINGER_BASE} = $config{'verbose_outfile_prefix'};
    Parse_Ringer_Data($CFG);
    %RINGER = %{$CFG->{RINGER}};
  } else { 
    #TODO: %RINGER = Parse_Ringer_Data($title, $dir, $prefix);
  }


  if ( wantarray ) { return  %RINGER }
  else             { return \%RINGER }

}

#=========================================================== Parse Ringer Data
# Ringer Properties and Files Parsed:
#   FILES and PATHS: filename and full paths to ringer data files
#     dir, out, ala, chi21, chi1, chi2, chi3, chi4, sig
#
#   MULTI: list of residues having multiple conformations in the original pdb
#     { resn, chain, resi }
#
#   UNMATCHED: residues with chi1 angles that do not match a density peak
#     { resn, resi, chain, built (conformation in pdb), found (closest peak) }
#
#   LOW_DENSITY: residues with peak density less than 1.0 sigma
#     { resn, resi, chain, density }
#
#   COUNT_*: counts of various side chains sampled (NOT RESIDUES)
#     total, read, ala, methyl, chi1, chi2, chi3, chi4, ring
# 
#   PEAKS: every peak found for each chi angle
#          read from the chi1-chi4 data files
#     $RINGER{PEAKS}{$chain}[$resi]{"X$chi"} = { density, angle }
# 
#   SEQ: protein sequence as read in signal plot
#   CHI: list of chi angles for each residue
#   DEGREES: array of the angle (in degrees) sampled in the signal plot
#   DENSITY: density at each DEGREE around a chi angle
#     $RINGER{DENSITY}{$chain}[$resi]{"X$chi"} = [ density, ... ]
#   DENSITY_MEAN and DENSITY_STDDEV: summarized density for the angle
#     Note that these are given an alt of {_} for consistency.
#     $RINGER{DENSITY_X}{$chain}[$resi]{"X$chi"}{_} = $density
#==============================================================================

sub Parse_Ringer_Data {

  my $CFG    = $_[0];
  my $title;
  my $dir    = $_[1];
  my $prefix = $_[2] || 'verbose';

  if ( Is_Structure_Config($CFG) ) { 
    $title  = $CFG->{TITLE};
    $dir    = $CFG->{DIRS }{RINGER};
    $prefix = $CFG->{FILES}{RINGER_BASE} || "verbose";
    if ( not defined $dir ) { 
      $dir = "$CFG->{DIR}/ringer_$title";
      $CFG->{FILES}{DIR_RINGER} = $dir;
    }
  } else { 
    if ( not defined $dir ) { 
      if ( not defined $title ) { 
        if ( -e "ringer" and -d "ringer" ) { 
          $dir = "ringer";
        } else { 
          confess "Phenix::Parse_Ringer_Data: No title or directory passed.\n";
        }
      } else { 
        $dir = "ringer_$title";
      }
    }
  }

  if ( not (-e $dir and -d $dir) ) { 
    confess "Phenix::Parse_Ringer_Data: Could not find ringer directory for title '$title', dir '$dir'.\n";
  }


  #------------------------------------ Primary Ringer Datastructure
  my %RINGER = (
    DIRS => { 
      RINGER => $dir,
    },
    FILES => { 
      out     => "ringer.out",
      ala     => "$prefix.ala_peaklist.txt",
      chi21   => "$prefix.chi2chi1_peaklist.txt",
      chi1    => "$prefix.ubchi1_peaklist.txt",
      chi2    => "$prefix.ubchi2_peaklist.txt",
      chi3    => "$prefix.chi3_peaklist.txt",
      chi4    => "$prefix.chi4_peaklist.txt",
      sig     => "$prefix.signal_plot.txt",
    },

    MULTI       => [],
    UNMATCHED   => [],
    LOW_DENSITY => [],
  );

  foreach my $file ( keys %{$RINGER{FILES}} ) { 
    $RINGER{PATHS}{$file} = "$dir/$RINGER{FILES}{$file}";
  }

  #------------------------------------ ringer.out file
  my $lines_out = Utility::Read_File($RINGER{PATHS}{out});

  #---------------- Multi Conformer Residues
  if ( $lines_out =~ /Removed multi-conformer residues:\s*\n(.*?\n)\n/s ) { 
    #or confess "Couldn't parse ringer.out for multi conf.\n";
    my @multi = map { [split "_", $_ ] }
                map { split /\,\s*/, $_ } 
                map { Utility::Trim($_) } 
                split "\n", 
                      $1;
    #print "Multiple conformer residues:\n";
    foreach my $multi ( @multi ) { 
      #print "$multi->[0], $multi->[1], $multi->[2]\n";
      push @{$RINGER{MULTI}}, { 
          resn  => $multi->[0], 
          chain => $multi->[1], 
          resi  => $multi->[2]
      };
    }
  }


  #---------------- Chi1 Conformation Doesn't Match Peak
  if ( $lines_out =~ /WARNING:  Following chi1 atoms do not match.*?position!\n(.*?\n)\n/s ) { 
    #or confess "Couldn't parse ringer.out for missed peaks.\n";
    my @peaks = map { [split /\s+/, $_] }
                map { Utility::Trim($_) } 
                split "\n", $1;

    shift @peaks; # throw away header
    foreach my $peak ( @peaks ) { 
      my ( $res, $built, $found ) = @$peak;
      my ( $resn, $chain, $resi ) = split "_", $res;
      push @{$RINGER{UNMATCHED}}, {
          resn => $resn,
          resi => $resi,
          chain => $chain,
          built => $built,
          found => $found,
      };
    }
  }

  #---------------- Density Less than 1.0 Sigma
  if ( $lines_out =~ /WARNING:  Following chi1 atoms have electron density.*?position!\n(.*?\n)~~~~~/s ) { 
    #or confess "Couldn't parse ringer.out for missed peaks.\n";
    my @lowden = map { [split /\s+/, $_] }
                 map { Utility::Trim($_) } 
                 split "\n", $1;

    shift @lowden; # throw away header
    foreach my $lowden ( @lowden ) { 
      my ( $res, $density) = @$lowden;
      my ( $resn, $chain, $resi ) = split "_", $res;
      push @{$RINGER{LOW_DENSITY}}, {
          resn    => $resn,
          resi    => $resi,
          chain   => $chain,
          density => $density,
      };
    }
  }



  #---------------- Miscellaneous Statistics
  $lines_out =~ /Total Sampled Side chains = (\d+)/i;           $RINGER{COUNT_TOTAL } = $1;
  $lines_out =~ /Side Chains Sampled for Alt Confs = (\d+)/i;   $RINGER{COUNT_READ  } = $1;
  $lines_out =~ /Ala side chains = (\d+)/i;                     $RINGER{COUNT_ALA   } = $1;
  $lines_out =~ /Methyl side chains = (\d+)/i;                  $RINGER{COUNT_METHYL} = $1;
  $lines_out =~ /UB Chi 1 side chains = (\d+)/i;                $RINGER{COUNT_CHI1  } = $1;
  $lines_out =~ /UB Chi 2 side chains = (\d+)/i;                $RINGER{COUNT_CHI2  } = $1;
  $lines_out =~ /Ring side chains = (\d+)/i;                    $RINGER{COUNT_RING  } = $1;
  $lines_out =~ /Chi 3 side chains= (\d+)/i;                    $RINGER{COUNT_CHI3  } = $1;
  $lines_out =~ /Chi 4 side chains= (\d+)/i;                    $RINGER{COUNT_CHI4  } = $1;




  #============================================== Chi Angles
  my %lines_chi  = ();
  foreach my $chi ( 1..4 ) { 
  
    #---- Read Data File 
    my $file = $RINGER{PATHS}{"chi$chi"};
    my @lines = Utility::Read_File($file);
 
    #===== Trash header
    shift @lines;

    #===== Parse Lines
    foreach my $line ( @lines ) { 
     
      #--- Get and parse residue specifier 
      my ( $residue, $peaks ) = split "\t", $line;
      $residue =~ /_(\w+)_(\w)_(\d+)$/ 
        or confess "Couldn't parse residue in file '$file' residue '$residue'\n";
      my ( $resn, $chain, $resi ) = ( $1, $2, $3 );

      #--- Parse peak data
      my @peaks_text = split /;\s*/, $peaks;
      $RINGER{PEAKS}{$chain}[$resi]{"X$chi"} = [];
      foreach my $peak ( @peaks_text ) { 
        $peak = Utility::Trim($peak);
        my ( $density, $angle ) = split /\s+/, $peak;
        push @{$RINGER{PEAKS}{$chain}[$resi]{"X$chi"}}, 
          { density => $density, angle => $angle, discrete => 1 };
      }

      #--- Retrieve, count, and sort the peaks
      my @peaks = @{$RINGER{PEAKS}{$chain}[$resi]{"X$chi"}};
      $RINGER{PEAK_COUNTS}{$chain}[$resi]{"X$chi"}{_}
        = scalar @peaks;
      # Peaks are sorted from highest to lowest
      @peaks = reverse sort { $a->{density} <=> $b->{density} } @peaks;


      #--- Test if each peak is surpassed by a close and higher peak
      foreach my $i ( 0..$#peaks ) { 
        my $peak_i   = $peaks[$i];
        next if not $peak_i->{discrete};
        my $angle_i  = $peak_i->{angle} % 360;

        foreach my $j ( $i+1 .. $#peaks ) { 
          my $peak_j   = $peaks[$j];
          my $angle_j  = $peak_j->{angle} % 360;

          my $d = Utility::Angular_Distance($angle_i, $angle_j);
          if ( $d < $DISCRETE_PEAK_CUTOFF ) { 
            $peak_j->{discrete} = 0;
            next;
          }
      }
      my $peak_count_discrete = grep { $_->{discrete} } @peaks;


      $RINGER{PEAK_COUNTS_DISCRETE}{$chain}[$resi]{"X$chi"}{_} = $peak_count_discrete;
    }
  }
     
  #============================================== Signal Plot
  ($RINGER{SEQ}, $RINGER{CHI}, $RINGER{DENSITY}, $RINGER{DEGREES})
    = Read_Ringer_Signal_Plot($RINGER{PATHS}{sig});

  my $Signal_Plot_Angular_Resolution = $RINGER{DEGREES}[1] - $RINGER{DEGREES}[0];

  my $MODES_CONFIG = { 
    circular  => 1,
    width     => $DISCRETE_PEAK_CUTOFF,
    binwidth  => $Signal_Plot_Angular_Resolution,
    length    => 360/$Signal_Plot_Angular_Resolution,
    start     => 0,
    stop      => 359,
    minmode   => $PEAK_DENSITY_MINIMUM,
  };

  };

  foreach my $chain ( keys %{$RINGER{DENSITY}} ) { 
    foreach my $resi ( 0..$#{$RINGER{DENSITY}{$chain}} ) { 
      my $resn = $RINGER{SEQ}{$chain}[$resi];
      my ( $aa, $ang_count );
      if ( defined $resn ) { 
        $aa = $Utility::AA{$resn};
        if ( defined $aa ) { 
          $ang_count = scalar @{$aa->{sc}};
        }
      } else { 
        $ang_count = 0;
      }

      #print "$chain, $resi, $resn, $ang_count\n";
      if ( not $ang_count ) { 
        $RINGER{PEAK_COUNTS         }{$chain}[$resi]{XX}{_} = 0;
        $RINGER{PEAK_COUNTS         }{$chain}[$resi]{AA}{_} = 0;
        $RINGER{PEAK_COUNTS_DISCRETE}{$chain}[$resi]{XX}{_} = 0;
        $RINGER{PEAK_COUNTS_DISCRETE}{$chain}[$resi]{AA}{_} = 0;
        $RINGER{DENSITY_MAX         }{$chain}[$resi]{XX}{_} = 0;
        $RINGER{DENSITY_MAX         }{$chain}[$resi]{AA}{_} = 0;
        $RINGER{DENSITY_MIN         }{$chain}[$resi]{XX}{_} = 0;
        $RINGER{DENSITY_MIN         }{$chain}[$resi]{AA}{_} = 0;
        $RINGER{DENSITY_MEAN        }{$chain}[$resi]{XX}{_} = 0;
        $RINGER{DENSITY_MEAN        }{$chain}[$resi]{AA}{_} = 0;
      }
      foreach my $chi ( 1..$ang_count ) { 

        # NOTE: Fixing PEAK COUNTS while we're here
        if (not defined $RINGER{PEAK_COUNTS}{$chain}[$resi]{"X$chi"}) { 
          $RINGER{PEAK_COUNTS}{$chain}[$resi]{"X$chi"}{_} = 1;
          $RINGER{PEAK_COUNTS_DISCRETE}{$chain}[$resi]{"X$chi"}{_} = 1;
        }

        my $density = $RINGER{DENSITY}{$chain}[$resi]{"X$chi"};
        if ( not defined $density ) { 
          #print "No density for $chain $resi $chi\n";
          $RINGER{DENSITY_MIN}{$chain}[$resi]{"X$chi"}{_} = 0;
          $RINGER{DENSITY_MAX}{$chain}[$resi]{"X$chi"}{_} = 0;
          ($RINGER{DENSITY_MEAN  }{$chain}[$resi]{"X$chi"}{_}, 
           $RINGER{DENSITY_STDDEV}{$chain}[$resi]{"X$chi"}{_}) = (0, 0);
        } else { 
          #print "Yes density for $chain $resi $chi\n";
          $RINGER{DENSITY_MIN}{$chain}[$resi]{"X$chi"}{_} = Utility::Min($density);
          $RINGER{DENSITY_MAX}{$chain}[$resi]{"X$chi"}{_} = Utility::Max($density);
          ($RINGER{DENSITY_MEAN  }{$chain}[$resi]{"X$chi"}{_}, 
           $RINGER{DENSITY_STDDEV}{$chain}[$resi]{"X$chi"}{_}) 
            = Utility::Mean_and_Stddev($density);
        }
      } #/CHI

      # Convenience property to map chi1 to aa for output
      if ( $ang_count > 0 ) { 
          $RINGER{PEAK_COUNTS_CHI1         }{$chain}[$resi]{AA}{_} 
        = $RINGER{PEAK_COUNTS              }{$chain}[$resi]{X1}{_};
          $RINGER{PEAK_COUNTS_DISCRETE_CHI1}{$chain}[$resi]{AA}{_} 
        = $RINGER{PEAK_COUNTS_DISCRETE     }{$chain}[$resi]{X1}{_};
      } else { 
          $RINGER{PEAK_COUNTS_CHI1         }{$chain}[$resi]{AA}{_} = 0;
          $RINGER{PEAK_COUNTS_DISCRETE_CHI1}{$chain}[$resi]{AA}{_} = 0;
      }
      $RINGER{PEAK_COUNTS_CHI1         }{$chain}[$resi]{AA}{'*'} = 
      $RINGER{PEAK_COUNTS_CHI1         }{$chain}[$resi]{AA}{'_'};
      $RINGER{PEAK_COUNTS_DISCRETE_CHI1}{$chain}[$resi]{AA}{'*'} =
      $RINGER{PEAK_COUNTS_DISCRETE_CHI1}{$chain}[$resi]{AA}{'_'};

    } #/RESI
  } #/CHAIN


  Collapse_Atom_Properties(\%RINGER, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'PEAK_COUNTS',         'sum');
  Collapse_Atom_Properties(\%RINGER, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'PEAK_COUNTS_DISCRETE','sum');
  Collapse_Atom_Properties(\%RINGER, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DENSITY_MIN',         'min');
  Collapse_Atom_Properties(\%RINGER, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DENSITY_MAX',         'max');
  Collapse_Atom_Properties(\%RINGER, $CFG->{PDB}{RESN}, $CFG->{PDB}{ALTS}, 'DENSITY_MEAN',        'avg');
 
  if ( Is_Structure_Config($CFG) ) { 
    Register_Analyses(
      $CFG, 
      [qw(ANGLE RESIDUE)], 
      'RINGER', 
      [ qw(PEAK_COUNTS PEAK_COUNTS_DISCRETE DENSITY_MEAN DENSITY_MIN DENSITY_MAX) ],
      ['PROTEIN'],
    );
    Register_Analyses (
      $CFG,
      [qw(RESIDUE)],
      'RINGER',
      [qw(PEAK_COUNTS_CHI1 PEAK_COUNTS_DISCRETE_CHI1)],
      ['PROTEIN'],
    );
    $CFG->{RINGER} = \%RINGER;
  }





  if ( wantarray ) { return  %RINGER }
  else             { return \%RINGER }
  

}

##############################################################################
#                                                 READ RINGER SIGNAL PLOT FILE
##############################################################################


sub Read_Ringer_Signal_Plot { 

  my $file  = $_[0];
  my @lines = Utility::Read_File($file);

  # Get and parse header
  my $header = shift @lines;
  $header =~ s/^\s*//;
  $header =~ s/\s*$//;
  my @header = map { [ split "_", $_ ] } split "\t", $header;
  #shift @header;

  # Parse data
  my @matrix = Utility::Transpose([map { [ split "\t", $_ ] } @lines]);


  my %SEQ = ();
  my %CHI = ();
  my %DEN = ();
  my @DEG = map { /^(\d+) .*$/; $1 } @{$matrix[0]};

  foreach my $i ( 0..$#header ) { 

    my ( $resn, $chain, $resi, $chi ) = ( @{$header[$i]} );
    $chi =~ s/chi/X/;

    #print "Found density for $chain $resi $resn $chi\n";

    $SEQ{$chain}[$resi] = $resn;
    $CHI{$chain}[$resi] = scalar @{$Utility::AA{lc $resn}{'sc'}};
    $DEN{$chain}[$resi]{$chi} = [ map { /^\d+ (.*)$/; $1 } @{$matrix[$i]} ];
  }
  return ( \%SEQ, \%CHI, \%DEN, \@DEG );
}






###############################################################################
###############################################################################
#                               DISTANCE ANALYSES
###############################################################################
###############################################################################

###############################################################################
#                                                            Get Coordinate Map
#------------------------------------------------------------------------------
# Convert a PDB structure (Chain, Residue, Atom, Alt) to a fast lookup list
# of atoms or a 2 level hierarchy of residue->atoms.
# atm_map is an ordered list of atoms
# res_map is an 2D matrix residues and its atoms
# res_set is an ordered list of hashes, to look up atoms by name in a residue
# coords is the list of xyz coordinates for every atom altloc
###############################################################################


sub Extract_Atom_Maps_Properties {
  my $CFG      = $_[0] or confess "Need CFG";
  my $PROPERTY = $_[1] or confess "Need Property";
  my $FILTER   = $_[2] or confess "Need filter subroutine";

  my @ATOM_MAP = ();
  my @RES_MAP  = ();
  my @RES_SET  = ();
  my @PROPS    = ();

  my $ri = 0;
  my $ai = 0;
  my $si = 0; 

  foreach my $residue ( @{$CFG->{RESIDUE_MAP}} ) { 
    my $resi = $residue->{resi};
    my $accepted = 0;
    foreach my $atom ( @{$residue->{atoms}} ) { 
      next if not &$FILTER($atom);
      $accepted = 1;
      $ATOM_MAP [$ai]       = $atom;
      $RES_MAP  [$ri]       = $atom;
      $RES_SET  [$ri][$si]  = $atom;
      $PROPS    [$ai]       = $atom->{$PROPERTY};
      $atom->{ri} = $ri;
      $atom->{ai} = $ai++;
      $atom->{si} = $si++;
    }
    $ri++ if $accepted;
  }
  return \@ATOM_MAP, \@RES_MAP, \@RES_SET, \@PROPS;
}







sub Get_Coordinate_Map { 

  my $CFG     = $_[0] or confess "Need CFG";
  my $FILTERS = $_[1] || {};
  $FILTERS->{type} //= ['P'];
  $FILTERS->{alt } //= ['_'];

  my $XYZ     = $CFG->{PDB}{COORDS};
  
  my @coords  = ();
  my @atm_map = ();
  my @res_map = ();
  my @res_set = ();

  my $atm_idx = 0;
  my $res_idx = 0;

  foreach my $chain ( sort keys %$XYZ ) { 
    foreach my $resi ( 0..$#{$XYZ->{$chain}} ) { 
      my $resn = $CFG->{PDB}{RESN}{$chain}[$resi];
      my $type = $CFG->{PDB}{TYPE}{$chain}[$resi];
      next if not defined $resn;
      next if not defined $type;
      #print "$chain $resi $resn $type\n";
      if ( defined $FILTERS->{resi} and scalar @{$FILTERS->{resi}} ) { 
        next if not Utility::Is_In_Set($FILTERS->{resi}, $resi);
      }
      if ( defined $FILTERS->{type} and scalar @{$FILTERS->{type}} ) { 
        next if not Utility::Is_In_Set($FILTERS->{type}, $type);
      }


      my @res = ();
      my $atm_lab = "";
      my $res_lab = "";
      my %res;
      foreach my $atom ( sort keys %{$XYZ->{$chain}[$resi]} ) { 
        foreach my $alt ( sort keys %{$XYZ->{$chain}[$resi]{$atom}} ) { 

          if ( defined $FILTERS->{alt} ) { 
            if ( not grep { $_ eq $alt } @{$FILTERS->{alt}} ) { 
              next;
            }
          }
          
          # Save Coordinates
          push @coords, $XYZ->{$chain}[$resi]{$atom}{$alt};

          # Create labels and save data structure
          $atm_lab = "$chain\_$resi\_$atom\_$alt";
          $res_lab = "$chain\_$resi" . ($resn ? "_$resn" : "");
          my $rolc = $Utility::AA{$resn}{ab} || $resn;
          my %atom = ( chain => $chain, resi => $resi, resn => $resn,
                       atom  => $atom,  alt  => $alt,   

                       atom_index    => $atm_idx,
                       atom_label    => $atm_lab,
                       #residue_name  => "$resn$resi",
                       #residue_name  => (ucfirst lc $resn) . $resi,
                       residue_name  => "$rolc$resi",
                       residue_index => $res_idx,
                       residue_label => $res_lab,
                     );
          %res  = ( chain => $chain, resi => $resi, resn => $resn,
                    #residue_name  => (ucfirst lc $resn) . $resi,
                    residue_name  => "$rolc$resi",
                    residue_index => $res_idx,
                    residue_label => $res_lab,
                  );

          push @atm_map, \%atom;
          push @res,  \%atom;
          $atm_idx++;
        }
      }
      next if (scalar @res) == 0;
      push @res_set, \@res;
      push @res_map, \%res;
      $res_idx++;
    }
  }

  return \@atm_map, \@res_map, \@res_set, \@coords;
}







###############################################################################
#                                                                Make Full DXM
#------------------------------------------------------------------------------
# I use this for both Difference Distance Matrix and Property Networks, so 
# what the hell.
###############################################################################

sub Make_Full_DXM {

  my $CFG = $_[0] or confess "Need CFG";
  my $FILTERS = $_[1] || {};
  $FILTERS->{type}   //= ['P'];
  $FILTERS->{alt }   //= ['_'];
  $FILTERS->{atoms } //= [   ];

  my $FILE_DXM = $CFG->{PATHS}{DXM};

  my %DXM = ();

  #=============== If DXM exists, just retrive it
  if ( -e $FILE_DXM ) { 
    print "Retrieving existing DXM file '$FILE_DXM'.\n";
    %DXM = %{retrieve($FILE_DXM)};
    
  #=============== Calculate DXM From Scratch
  } else { 
    print "Calculating new DXM set.\n";
    
    #------------- Get Coordinates and CA Atom Subset
    ($DXM{MAP}{ATOM}, $DXM{MAP}{RES}, $DXM{SET}{RES}, $DXM{COORDS})
      = Phenix::Get_Coordinate_Map($CFG, $FILTERS);  
      #= Phenix::Get_Coordinate_Map($CFG, '_', 'P');  
    $DXM{MAP}{CA} = $DXM{MAP}{RES};
    $DXM{SET}{CA} = [];
    foreach my $i ( 0..$#{$DXM{MAP}{RES}} ) { 
      $DXM{SET}{CA}[$i] = [ grep { $_->{atom} eq 'CA' or $_->{atom} eq 'FE' } @{$DXM{SET}{RES}[$i]} ];
    }
    print "    Atom Count: " . (scalar @{$DXM{MAP}{ATOM}}) . "\n";


    #------------- Calculate Distance Matrices
    print "  Creating distance matrix.\n";
    $DXM{D1M}{ATOM} = Phenix::Distance_Matrix($DXM{COORDS}, $DXM{COORDS});
    $DXM{D1M}{RES } = Phenix::Distance_Matrix_Groups($DXM{D1M}{ATOM}, $DXM{SET}{RES});
    $DXM{D1M}{CA  } = Phenix::Distance_Matrix_Groups($DXM{D1M}{ATOM}, $DXM{SET}{CA});

    #------------- Calculate Contact Matrices and Statistics
    print "  Calculating contact matrix.\n";
    foreach my $set ( qw(ATOM RES CA) ) { 
      $DXM{CONTACTS}{$set} = Phenix::Contact_Matrix($DXM{D1M}{$set});
      $DXM{LIST    }{$set} = [Utility::Flatten($DXM{D1M}{$set})];
      $DXM{STATS   }{$set} = Utility::Calculate_Distribution_Statistics($DXM{LIST}{$set}, {hist=>0});

      $DXM{MAX}{$set} = Utility::Max($DXM{MAX}{$set}, $DXM{STATS}{$set}{max});
      $DXM{MIN}{$set} = Utility::Min($DXM{MIN}{$set}, $DXM{STATS}{$set}{min});
    }

    #------------- Generate Labels
    print "  Making labels.\n";
    $DXM{LAB}{ATOM} = [ map { $_->{atom_label   } } @{$DXM{MAP}{ATOM}} ] ;
    $DXM{LAB}{RES } = [ map { $_->{resi         } } @{$DXM{MAP}{RES }} ] ;
    $DXM{LAB}{CA  } = [ map { $_->{resi         } } @{$DXM{MAP}{RES }} ] ;

    #$DXM{LABEL}{ATOM} = [ map { $_->{atom_label   } . 'Z' } @{$DXM{MAP}{ATOM}} ] ;
    #$DXM{LABEL}{RES } = [ map { $_->{residue_name } . 'X' } @{$DXM{MAP}{RES }} ] ;
    #$DXM{LABEL}{CA  } = [ map { $_->{residue_name } . 'Y' } @{$DXM{MAP}{RES }} ] ;
    
    $DXM{LABEL}{ATOM} = [ map { $_->{atom_label   } } @{$DXM{MAP}{ATOM}} ] ;
    $DXM{LABEL}{RES } = [ map { $_->{resi         } } @{$DXM{MAP}{RES }} ] ;
    $DXM{LABEL}{CA  } = [ map { $_->{resi         } } @{$DXM{MAP}{RES }} ] ;

    #------------- Set Output Matrix Properties
#    my $wm_atom = { col_labels  => $DXM{NAME}{ATOM},
#                    row_labels  => $DXM{NAME}{ATOM},
#                    overwrite   => 1,
#                    data_format => "%.2f",
#                    title       => "D1M_ATOM",
#                  };
#    my $wm_res  = { col_labels  => $DXM{NAME}{RES},
#                    row_labels  => $DXM{NAME}{RES},
#                    overwrite   => 1,
#                    data_format => "%.2f",
#                    title       => "D1M_RES",
#                  };
#    my $wm_ca   = { %$wm_res, 
#                    title => 'D1M_CA' };

    my $OUTPUT_LABEL = 'LABEL'; # Write Thr234
    #my $OUTPUT_LABEL = 'LAB';   # Residue Number

    my $wm_atom = { col_labels  => $DXM{$OUTPUT_LABEL}{ATOM},
                    row_labels  => $DXM{$OUTPUT_LABEL}{ATOM},
                    overwrite   => 1,
                    data_format => "%.2f",
                    title       => "D1M_ATOM",
                  };
    my $wm_res  = { col_labels  => $DXM{$OUTPUT_LABEL}{RES},
                    row_labels  => $DXM{$OUTPUT_LABEL}{RES},
                    overwrite   => 1,
                    data_format => "%.2f",
                    title       => "D1M_RES",
                  };
    my $wm_ca   = { %$wm_res, 
                    title => 'D1M_CA' };

    #------------- Write Distance Matrices
    print "  Writing D1M matrices.\n";
    Utility::Write_Matrix($CFG->{PATHS}{D1M_ATOM}, $DXM{D1M}{ATOM}, $wm_atom);
    Utility::Write_Matrix($CFG->{PATHS}{D1M_RES }, $DXM{D1M}{RES }, $wm_res );
    Utility::Write_Matrix($CFG->{PATHS}{D1M_CA  }, $DXM{D1M}{CA  }, $wm_ca  );

    #------------- Write Contact Matrices
    $wm_atom->{title} = 'CON_ATOM';
    $wm_res ->{title} = 'CON_RES';
    $wm_ca  ->{title} = 'CON_CA';
    print "  Writing Contact Matrices.\n";
    Utility::Write_Matrix($CFG->{PATHS}{CON_ATOM}, $DXM{CONTACTS}{ATOM}, $wm_atom);
    Utility::Write_Matrix($CFG->{PATHS}{CON_RES }, $DXM{CONTACTS}{RES }, $wm_res );
    Utility::Write_Matrix($CFG->{PATHS}{CON_CA  }, $DXM{CONTACTS}{CA  }, $wm_ca  );

    store(\%DXM, $FILE_DXM);

  }

  # Save %DXM into $CFG;
  $CFG->{DXM} = \%DXM;
}




###############################################################################
#                                                               Distance Matrix
#------------------------------------------------------------------------------
# Given a coordinate matrix (typically from Get_Coordinate_Map), calculate the
# distances between every atom pair.
###############################################################################

sub Distance_Matrix {

  my $coords1 = $_[0];
  my $coords2 = $_[1];
  
  my @distances = ();
  foreach my $i ( 0..$#$coords1 ) { 
    my $ci = $coords1->[$i];
    next if not defined $ci;
    foreach my $j ( 0..$#$coords2 ) {  
      my $cj = $coords2->[$j];
      next if not defined $cj;
      $distances[$i][$j] = Utility::Distance($ci, $cj);
    }
  }
  if ( wantarray ) { return  @distances }
  else             { return \@distances }
}

###############################################################################
#                                                        Distance Matrix Groups
#------------------------------------------------------------------------------
# Given a distance matrix (typically for atoms) and a set of group definitions
# (typically for residues from Get_Coordinate_Map, but could be others)
# Create a new group distance matrix with the minimum distance between any 
# pair of atoms from the two groups.
###############################################################################

sub Distance_Matrix_Groups { 
  my $distances_atoms = $_[0];
  my $groups = $_[1];

  my @distances_groups = ();

  foreach my $gi ( 0..$#$groups ) { 
    my $groupi = $groups->[$gi];
    foreach my $gj ( 0..$#$groups ) { 
      my $groupj = $groups->[$gj];
      my $distance;
      foreach my $ai ( 0..$#$groupi ) { 
        my $atomi = $groupi->[$ai];
        foreach my $aj ( 0..$#$groupj ) { 
          my $atomj = $groupj->[$aj];
          my $d = $distances_atoms->[$atomi->{atom_index}][$atomj->{atom_index}];
          if ( not defined $distance ) { $distance = $d }
          if ( $d < $distance        ) { $distance = $d }
          #if ( $atomi->{resi} == 65 and $atomj->{resi} == 119 ) { 
          #  print Dumper($atomi) . "\n";
          #  print Dumper($atomj) . "\n";
          #  printf "Distance: %.2f of %.2f\n", $d, $distance;
          #}
        }
      }
      $distances_groups[$gi][$gj] = $distance;
    }
  }
  if ( wantarray ) { return  @distances_groups }
  else             { return \@distances_groups } 
}

  

###############################################################################
#                                                    Distance Difference Matrix
#------------------------------------------------------------------------------
# Given two distance matrices, take the pairwise difference. Optionally 
# truncate to $min and $max distances.
###############################################################################


sub Distance_Difference_Matrix { 
  my $distances1 = $_[0];
  my $distances2 = $_[1];
  my $min        = $_[2];
  my $max        = $_[3];


  if ( scalar @$distances1 != scalar @$distances2 ) { 
    confess "Distance matrices are not the same length.\n";
  }

  my @D2M = ();
  foreach my $i ( 0..$#$distances1 ) { 
    foreach my $j ( 0..$#{$distances1->[$i]} ) { 
      $D2M[$i][$j] = 0;
      next if not defined $distances1->[$i][$j];
      next if not defined $distances2->[$i][$j];
      my $d2m = $distances2->[$i][$j] - $distances1->[$i][$j];
      if ( defined $min and ($d2m < $min) ) { $d2m = $min }
      if ( defined $max and ($d2m > $max) ) { $d2m = $max }
      $D2M[$i][$j] = $d2m;
    }
  }

  if ( wantarray ) { return  @D2M }
  else             { return \@D2M }
}

###############################################################################
#                                                                Contact Matrix
#------------------------------------------------------------------------------
# Given a distance matrix and a cutoff, return a matrix of atoms (or whatever, 
# as it only sees the distances) in contact.
###############################################################################

sub Contact_Matrix { 
  my $distances = $_[0];
  my $cutoff    = $_[1] || 5;
 
  my @contacts = (); 
  
  foreach my $i ( 0..$#$distances ) { 
    foreach my $j ( 0..$#{$distances->[$i]} ) { 
      my $contact;
      my $distance = $distances->[$i][$j];
      if ( not defined $distance   ) { $contact = 0 }
      elsif ( $distance <= $cutoff ) { $contact = 1 }
      else                           { $contact = 0 }
      $contacts[$i][$j] = $contact;
    }
  }
  if ( wantarray ) { return  @contacts }
  else             { return \@contacts }
}


###############################################################################
#                                                          Get Contacting Pairs
#------------------------------------------------------------------------------
# Given a contact matrix and optional subset filter, return a list of pairs
# of atoms that are in contact.
###############################################################################

sub Get_Contacting_Pairs { 
  my $contacts = $_[0];
  my $subset   = $_[1];
 
  my @pairs = (); 
  foreach my $i ( 0..$#$subset ) { 
    foreach my $j ( ($i+1)..$#$subset ) { 
      if ( $contacts->[$i][$j] ) { 
        push @pairs, [$i,$j];
      }
    }
  }
  if ( wantarray ) { return  @pairs } 
  else             { return \@pairs }

}







##############################################################################
##############################################################################
#                     DIFFRACTION IMAGE PROCESSING                            
##############################################################################
##############################################################################

##############################################################################
#                             DIFFDUMP                                       #
##############################################################################

our $DIFFDUMP_HOST = 'turing2.berkeley.edu';
our $DIFFDUMP_CMD = '/labs/sbgrid/programs/x86_64-linux/ccp4/6.2.0/ccp4-6.2.0/bin/diffdump';

sub diffdump { 

  my $file = $_[0];

  if ( not defined $PROGRAMS{DIFFDUMP} ) { 
    confess "Couldn't locate executable for diffdump.\n";
  }

  #======================================================= Handle Zipped Files
  # TODO
  if ( $file =~ /\.bz2$/ ) { 
  } elsif ( $file =~ /\.gzip$/ ) { 
  }

  my $dump = `$PROGRAMS{DIFFDUMP} $file`;

  my %D = ();

  $dump =~ /Format : (?<format>.+)
Manufacturer : (?<manufacturer>.+?)
Collection date : (?<date>.+?)
Exposure time : (?<exposure>[\d\.]+?) s
Detector S\/N : (?<sn>\d+?)
Wavelength : (?<wavelength>[\d\.]+?) Ang
Beam center : \((?<beamx>[\d\.]+?) mm,(?<beamy>[\d\.]+?) mm\)
Distance to detector : (?<distance>[\d\.]+?) mm
Image Size : \((?<imagex>\d+) px, (?<imagey>\d+) px\)
Pixel Size : \((?<pixelx>[\d\.]+) mm, (?<pixely>[\d\.]+) mm\)
Oscillation \(.*\) : (?<oscstart>[\d\.]+) -> (?<oscend>[\d\.]+) deg
Two Theta value: (?<twotheta>[-\d\.]+) deg
/s
  or carp "DFXN::DiffDump Couldn't parse diffdump output, '$dump'\n";


  foreach my $key ( qw(format manufacturer date exposure sn wavelength beamx beamy distance imagex imagey pixelx pixely oscstart oscend twotheta) ) { 
    $D{$key} = $+{$key};
  }

  $D{osc} = $D{oscend} - $D{oscstart};
  #$D{timestamp} = Date::Manip::ParseDate($D{date});

  if ( wantarray ) { return  %D }
  else             { return \%D }
}



##############################################################################
#                              SPOTFINDER                                    #
##############################################################################


our $SPOTFINDER_VERBOSE =  3;
our $SPOTFINDER_WIDTH   = 50;
our $SPOTFINDER_SERVER_LAUNCH = 1;

my $SPOTFINDER_SERVER_PID;
my $SPOTFINDER_SERVER_OUT;
my $SPOTFINDER_SERVER_ERR;
my $SPOTFINDER_SERVER_IN;

my $SPOTFINDER_FILE_OUT = "spotfinder.server.out";
my $SPOTFINDER_FILE_ERR = "spotfinder.server.err";

#my $SPOTFINDER_HOST   = "localhost";
my $SPOTFINDER_HOST   = "crystal.qb3.berkeley.edu";
my $SPOTFINDER_PORT   = 8126;
my $SPOTFINDER_CPUS   = 6;
my $SPOTFINDER_QUEUE  = 12;
my $SPOTFINDER_WAIT_SERVER = 5.0; # time to wait for the server to start up.
my $SPOTFINDER_WAIT_CLIENT = 0.1; # time to wait before launching another job.
my $SPOTFINDER_WAIT_POLL   = 1.0; # polling interval
my $SPOTFINDER_WAIT_STEPS  = 30;  # maximum tries

my $SPOTFINDER_SERVER_CMD = "distl.mp_spotfinder_server_read_file "
                        . "distl.port=$SPOTFINDER_PORT "
                        . "distl.processors=$SPOTFINDER_CPUS "
                        . "distl.bins.verbose=True "
                        #. "distl.verbose=True "
                        . ">$SPOTFINDER_FILE_OUT "
                        . "2>$SPOTFINDER_FILE_ERR "
                        ;

my $SPOTFINDER_CLIENT_CMD = "distl.thin_client FILE HOST PORT";

sub Spotfinder { 

  my @files = @_;

  print "Starting Spotfinder.\n"
    if $SPOTFINDER_VERBOSE >= 1;



  # START A SERVER IF NEEDED
  # NOTE: For now we will always start our own, because we may want
  # custom parameters. So we'll use a non-standard port.
  if ( $SPOTFINDER_SERVER_LAUNCH and not $SPOTFINDER_SERVER_PID ) { 
    if ( -e $SPOTFINDER_FILE_OUT ) { 
      unlink $SPOTFINDER_FILE_OUT or confess "Couldn't delete '$SPOTFINDER_FILE_OUT'. $!";
    }
    if ( -e $SPOTFINDER_FILE_ERR) { 
      unlink $SPOTFINDER_FILE_ERR or confess "Couldn't delete '$SPOTFINDER_FILE_ERR'. $!";
    }
    
    print "  Launching Spotfinder Server\n    $SPOTFINDER_SERVER_CMD\n" 
      if $SPOTFINDER_VERBOSE >= 1;
    eval { 
      $SPOTFINDER_SERVER_PID = open3($SPOTFINDER_SERVER_IN, 
                                     $SPOTFINDER_SERVER_OUT, 
                                     $SPOTFINDER_SERVER_ERR, 
                                     $SPOTFINDER_SERVER_CMD);
    };
    if ( $@ ) { 
      carp ("open3: $@");
    }

    print "  Launched with PID $SPOTFINDER_SERVER_PID\n"
        . "    Waiting for server to initialize\n"
        if $SPOTFINDER_VERBOSE >= 1;

    #==== First, poll log file creation

    my $log_found = 0;
    for ( 1..$SPOTFINDER_WAIT_STEPS ) { 
      if (-e $SPOTFINDER_FILE_OUT and (-s $SPOTFINDER_FILE_OUT > 100)) { 
        $log_found = 1;
        last ;
      }
      Time::HiRes::usleep($SPOTFINDER_WAIT_POLL * 1000000);
      print "." if $SPOTFINDER_VERBOSE >= 3;
    }
    if ( not $log_found ) { 
      confess ( "Server never made log file.\n");
    }

    #==== Consider launching a test thin_client to see if we get connection refused.
  }

  #==== Detect server launch failure.
  if ( not $SPOTFINDER_SERVER_PID ) { 
    my $error = readline($SPOTFINDER_SERVER_ERR);
    carp "Failed to launch spotfinder server, '$error'.";
  }

  
  # START A CIRCULAR QUEUE SUBMITTING JOBS TO THE SERVER
  my $q = 0;
  my %SPOTFINDER = ();
  
  my %Q = ( 
    ins   => [],
    outs  => [],
    errs  => [],
    pids  => [],
    files => [],
  );

  #my (@q_ins, @q_outs, @q_errs, @q_pids, @q_files) = ();

  print "Processing " . (scalar @files) . " files.\n"
    if $SPOTFINDER_VERBOSE >= 1;
  foreach my $i ( 0..$#files ) { 
    
    my $file = $files[$i];

    # We've (possibly) come around the queue, reap the job
    Reap_Spotfinder_Job(\%SPOTFINDER, \%Q, $q);

    # Execute a new job
    if    ( $SPOTFINDER_VERBOSE >= 2 ) { 
      printf "Loop %6d Q %3d File %s\n", $i, $q, $file;
    } elsif ( $SPOTFINDER_VERBOSE == 1 ) { 
      print ".";
      if ( not ( ($i+1) % $SPOTFINDER_WIDTH ) ) { 
        printf "  %s\n", ($i+1);
      }
    }

    my $SPOTFINDER_CLIENT_CMD = "distl.thin_client $file $SPOTFINDER_HOST $SPOTFINDER_PORT";
    #print "  $SPOTFINDER_CLIENT_CMD\n";
    $Q{files}[$q] = $file;
    eval { 
      $Q{pids}[$q] = open3($Q{ins}[$q], $Q{outs}[$q], $Q{errs}[$q], 
                    $SPOTFINDER_CLIENT_CMD);
    };
    if ( $@ ) { carp ("distl.thin_client: $@"); }
    Time::HiRes::usleep($SPOTFINDER_WAIT_CLIENT * 1000000);
    
    # Increment our place in the queue
    $q = ($q+1) % $SPOTFINDER_QUEUE;
  }

  print "\nDone Queueing\n" if $SPOTFINDER_VERBOSE >= 1;


  # REAP ANY JOBS NOT FINISHED YET
  foreach my $q ( 0..($SPOTFINDER_QUEUE-1) ) { 
    Reap_Spotfinder_Job(\%SPOTFINDER, \%Q, $q);
  }



  # Convert the hash to an ordered array
  my @SPOTFINDER = map { $SPOTFINDER{$_} } @files;

  if ( wantarray ) { return  @SPOTFINDER }
  else             { return \@SPOTFINDER }

}


sub Reap_Spotfinder_Job { 

  my $SF = $_[0]; # spotfinder data hash
  my $Q  = $_[1]; # spotfinder queue hash array
  my $q  = $_[2]; # index in the array

  return if not defined $Q->{pids}[$q];

  print "Waiting on pid $Q->{pids}[$q] (nohang)\n" if $SPOTFINDER_VERBOSE >= 3;
  my $kid;
  $kid = waitpid $Q->{pids}[$q], 0;
 
  #$| = 1;
  #for(1..100) { 
  #  $kid = waitpid $Q->{pids}[$q], 0;
  #  if ( not $kid ) { Time::HiRes::usleep($SPOTFINDER_WAIT_CLIENT) }
  #  else { last }
  #  #print ".";
  #}
  ##print "\n";

  if ( not $kid ) { 
    print "Client hung, terminating pid $Q->{pids}[$q]\n";
    print "kill it manually.\n";
  }

  my ( $result, $errors, $sf );

  #print "Reading results\n";
  if ( $Q->{errs}[$q] ) { 
    $result = join "", readline($Q->{outs}[$q]);
    $errors = join "", readline($Q->{errs}[$q]);
    printf "Loop %6s Q %3d File %s had an error.\n%s\n%s\n", '', $q, $Q->{files}[$q], $result, $errors;
    exit;
  } else { 
    $result = join "", readline($Q->{outs}[$q]);
    if ( $SPOTFINDER_VERBOSE >= 3 ) { 
      printf "Loop %6s Q %3d File %s reaped\n", '', $q, $Q->{files}[$q];
    } elsif ( $SPOTFINDER_VERBOSE == 2 ) { 
      print "r";
    }
    $sf = Parse_Spotfinder($result);
    $SF->{$Q->{files}[$q]} = $sf;
    #print Dumper($sf);
  }



  undef $Q->{ins  }[$q];
  undef $Q->{outs }[$q];
  undef $Q->{errs }[$q];
  undef $Q->{pids }[$q];
  undef $Q->{files}[$q];
  
  if ( $errors ) { return $errors }
  else           { return $result }
}







sub Parse_Spotfinder { 

  my $text = $_[0];

  # Detect common errrs
  if ( $text =~ /Connection refused/s ) { 
    confess "DFXN::Parse_Spotfinder client connection was refused.\n$text\n'";
  }

  my %S = ();
  my $binempty = "0 0 0 0 0 0 0 0 0 0 0 0";

  $text =~ /(Analysis of spots after ice removal.*?\n)\n/s;
  my $bins1 = ($1 || $binempty);

  $text =~ /(Analysis of spots after resolution filtering.*?\n)\n/s;
  my $bins2 = ($1 || $binempty);

  $text =~ /(Analysis of good Bragg spots.*?\n)\n/s;
  my $bins3 = ($1 || $binempty);


  #print "$bins1\n$bins2\n$bins3\n";

  # Split into lines
  my @binlines1 = split "\n", $bins1;
  my @binlines2 = split "\n", $bins2;
  my @binlines3 = split "\n", $bins3;
 
  # Take the summary line from the end 
  my @binsum1 = grep { defined } split /\s+/, $binlines1[-1];
  my @binsum2 = grep { defined } split /\s+/, $binlines2[-1];
  my @binsum3 = grep { defined } split /\s+/, $binlines3[-1];

  my @anal = qw(trash spots signal i bkg size isigi ecc skew);
  @S{ map { "ALL_$_" } (@anal, "text") } = (@binsum1, $bins1);
  @S{ map { "RES_$_" } (@anal, "text") } = (@binsum2, $bins2);
  @S{ map { "BRG_$_" } (@anal, "text") } = (@binsum3, $bins3);

  $text =~ /Image:\s+(?<image>.*?)
.*Spot Total :\s+(?<spots>\d+)
.*Resolution Total :\s+(?<spots2>\d+)
.*Good Bragg Candidates :\s+(?<bragg>\d+)
.*integrated signal.*?:\s+(?<signal>\d+)
.*Ice Rings : \s+(?<icerings>\d+)
.*Method 1 Resolution :\s+(?<rsln1>[\d\.]+)
.*Method 2 Resolution :\s+(?<rsln2>[\w\d\.]+)
.*Maximum unit cell :\s+(?<maxuc>[\d\.]+)
.*Saturation, Top \d+? Peaks :\s+(?<saturation>[\d\.]+) %
.*overloaded spots :\s+(?<overloads>\d+)
/s or confess "DFXN::Parse_Spotfinder Couldn't parse spotfinder output, '$text'.\n";

  foreach my $key ( qw(image spots spots2 bragg signal icerings rsln1 rsln2 maxuc saturation overloads) ) { 
    my $v = $+{$key};
    if ( $v eq 'None' ) { $v = 0; }
    $S{$key} = $v;
  }
  #$S{text} = $text;

  if ( wantarray ) { return  %S }
  else             { return \%S }
}



=pod

========== VERBOSE OUTPUT OF DISTL.THIN_CLIENT

scouras@crystal:~/crystals/Lipo_S3A_XT2_C5_x2_100K/DFXN> distl.thin_client Lipo_S3A_XT2_C5_x2_100K_HR_2_001.img localhost 8125

Limit: outer edge of (equal-volume) resolution shell (Angstroms)
Population: number of spots in the resolution shell
Missing: fraction of reciprocal space shell recorded in complete rotation, accounting for missing cone
Fract: fraction of shell recorded, accounting for truncated detector corner
PxlBkgrd: mean pixel background at center of scanbox windows in the shell
MnWndwSz: mean scanbox window size in pixels after discarding spot pixels
Integrated: total integrated signal within the shell in analog-digital units above background
MeanI: mean integrated signal for spots within the shell
MeanBkg: mean integrated background for spots within the shell
MeanSz: mean pixel count for spots within the shell
MeanIsigI: mean I over sig(I) for spots within the shell
MnEccen: mean eccentricity for spots within the shell; 0=circle, 1=parabola
MnSkew: mean skewness for spots within the shell; skew:= (maximum-center_of_gravity)/semimajor_axis

Analysis of spots after ice removal, but prior to resolution cutoff, for image 
/nfs/turing2/alber/scouras/crystals/Lipo_S3A_XT2_C5_x2_100K/DFXN/Lipo_S3A_XT2_C5_x2_100K_HR_2_001.img:
Limit  Population  Missing  Fract  PxlBkgrd  MnWndwSz  Integrated  MeanI  MeanBkg  MeanSz  MeanIsigI  MnEccen  MnSkew 
 4.80         150     1.00   1.00     341.3       598      407460   2716   8159.7      24     22.742    0.871   0.109 
 3.81         137     0.99   1.00     338.3       598      366115   2672   9642.9      29     20.326    0.832   0.106 
 3.33          68     0.99   1.00     316.5       598      188101   2766   8986.8      29     19.914    0.826   0.174 
 3.03          27     0.98   1.00     264.6       593       24663    913   5402.1      21     11.013    0.828   0.215 
 2.81          12     0.98   1.00     238.2       591        8742    728   4504.1      20      9.665    0.841   0.118 
 2.64           8     0.98   1.00     217.0       592        6496    812   4384.8      20     10.385    0.840   0.132 
 2.51           3     0.98   1.00     197.7       592        1128    376   2505.4      13      6.908    0.919   0.271 
 2.40           3     0.97   1.00     179.6       591        1632    544   2895.3      17      8.953    0.720   0.213 
 2.31           0     0.97   0.82     167.7       589           0      0      0.0       0      0.000    0.000   0.000 
 2.23           4     0.97   0.59     156.1       592       13240   3310   5801.1      38     23.725    0.879   0.178 
 2.16           1     0.97   0.46     147.6       587        2217   2217   1549.9      11     35.979    0.382   0.000 
 2.10           1     0.97   0.37     140.3       565         181    181   1615.8      12      4.233    0.882   0.000 
 2.04           2     0.96   0.29     135.3       650       13958   6979  12076.4      90     35.713    0.873   0.000 
        ----------                                      ----------                                                    
              416                                         1033932   2485   8312.1      26     19.953    0.845   0.128 

Analysis of spots after resolution filtering, but prior to spot quality heuristics, for image 
/nfs/turing2/alber/scouras/crystals/Lipo_S3A_XT2_C5_x2_100K/DFXN/Lipo_S3A_XT2_C5_x2_100K_HR_2_001.img:
Limit  Population  Missing  Fract  PxlBkgrd  MnWndwSz  Integrated  MeanI  MeanBkg  MeanSz  MeanIsigI  MnEccen  MnSkew 
 4.80         150     1.00   1.00     341.3       598      407460   2716   8159.7      24     22.742    0.871   0.109 
 3.81         137     0.99   1.00     338.3       598      366115   2672   9642.9      29     20.326    0.832   0.106 
 3.33          68     0.99   1.00     316.5       598      188101   2766   8986.8      29     19.914    0.826   0.174 
 3.03           3     0.98   1.00     264.6       593        2357    786   5227.3      18      9.939    0.820   0.368 
        ----------                                      ----------                                                    
              358                                          964032   2693   8859.8      27     21.173    0.847   0.122 

Analysis of good Bragg spots after quality heuristics, for image 
/nfs/turing2/alber/scouras/crystals/Lipo_S3A_XT2_C5_x2_100K/DFXN/Lipo_S3A_XT2_C5_x2_100K_HR_2_001.img:
Limit  Population  Missing  Fract  PxlBkgrd  MnWndwSz  Integrated  MeanI  MeanBkg  MeanSz  MeanIsigI  MnEccen  MnSkew 
 4.80         142     1.00   1.00     341.3       598      397701   2801   8223.2      24     23.297    0.873   0.094 
 3.81         124     0.99   1.00     338.3       598      317439   2560   9537.7      29     20.062    0.836   0.093 
 3.33          55     0.99   1.00     316.5       598      118195   2149   8338.4      27     17.892    0.826   0.142 
 3.03           1     0.98   1.00     264.6       593         919    919   5075.6      17     11.749    0.835   0.339 
        ----------                                      ----------                                                    
              322                                          834254   2591   8739.3      27     21.093    0.850   0.102 

Image: /nfs/turing2/alber/scouras/crystals/Lipo_S3A_XT2_C5_x2_100K/DFXN/Lipo_S3A_XT2_C5_x2_100K_HR_2_001.img
                                                     Spot Total :     433
                                      Method-2 Resolution Total :     358
                                          Good Bragg Candidates :     322
Total integrated signal, pixel-ADC units above local background :  834254
                                                      Ice Rings :       0
                                            Method 1 Resolution :    2.92
                                            Method 2 Resolution :    3.31
                                              Maximum unit cell :   191.1
                                       Saturation, Top 50 Peaks :   1.5 %
                                 In-resolution overloaded spots :       0

=cut















##############################################################################
##############################################################################

sub Cleanup {

  if ( $SPOTFINDER_SERVER_PID ) { 

    print "Exiting the spotfinder server gently.\n";
    my $SPOTFINDER_EXIT_CMD = "distl.thin_client EXIT $SPOTFINDER_HOST $SPOTFINDER_PORT";
    print "$SPOTFINDER_EXIT_CMD\n";
    print `$SPOTFINDER_EXIT_CMD`;
    sleep($SPOTFINDER_WAIT_SERVER);
    if ( kill 0, $SPOTFINDER_SERVER_PID ) { 
      #use Proc::Killfam;
      print "Killing the spotfinder server $SPOTFINDER_SERVER_PID\n";
      #killfam "SIGTERM", $SPOTFINDER_SERVER_PID;
    }

    #my $LINE = "--------------------";
    #print "$LINE\n"
    #    . "Spotfinder Server Output:\n"
    #    . (join "", readline($SPOTFINDER_SERVER_OUT))
    #    . "\n$LINE\n"
    #    . "Spotfinder Server Error:\n"
    #    . (join "", readline($SPOTFINDER_SERVER_ERR))
    #    . "\n$LINE\n"
    #    ;

    

    # Kill the server
    #print "Exiting the server gently.\n";
    #my $SPOTFINDER_EXIT_CMD = "distl.thin_client EXIT $SPOTFINDER_HOST $SPOTFINDER_PORT";
    #print "$SPOTFINDER_EXIT_CMD\n";
    #print `$SPOTFINDER_EXIT_CMD`;
    #sleep($SPOTFINDER_WAIT_SERVER);
    #if ( kill 0, $SPOTFINDER_SERVER_PID ) { 
    #  print "Terminating the server\n";
    #  kill "SIGTERM", $SPOTFINDER_SERVER_PID;
    #  sleep($SPOTFINDER_WAIT_SERVER);
    #}
    #if ( kill 0, $SPOTFINDER_SERVER_PID ) { 
    #  print "KILLING the server\n";
    #  kill "SIGKILL", $SPOTFINDER_SERVER_PID; 
    #  sleep($SPOTFINDER_WAIT_SERVER);
    #}
    #if ( kill 0, $SPOTFINDER_SERVER_PID ) { 
    #  print "Server still alive. :(\n";
    #}
    #print "Waiting on PID\n";
    #my $kid = waitpid $SPOTFINDER_SERVER_PID, 0;

  }
  exit;
}

END { &Cleanup() }
$SIG{INT} = \&Cleanup;




##############################################################################
##############################################################################
