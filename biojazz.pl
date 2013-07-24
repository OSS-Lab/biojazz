#!/usr/bin/perl -w
###############################################################################
## File:     biojazz.pl
#- Synopsys: Automated design or evolution of biochemical protein networks
#-           using a genetic algorithm.
#-#############################################################################
#- Detailed Description:
#- ---------------------
#- INVOCATION:
#-
#- biojazz.pl [option]...
#-
#- OPTIONS:
#-
#- --help            This help page.
#- --version         Print version/release info and quit.
#-
#- --verbosity=i     Verbosity level.  Defaults to 1.
#-
#- --config=s        Configuration file.
#-
#- --command=s       Execute given command.
#- --script=s        Execute given script.
#- --shell           Run an interactive shell after command/script execution.
#-
#- --seed=i          Random number generator seed value.
#- --tag=s           File tag and name of directory where results are stored.
#-
#- --cluster_type=s  Slave node creation.  Set to LOCAL to use localhost only,
#-                   SSH to use host_list in configuration_file, or PBS.
#- --cluster_size=i  Number of slave nodes to create.
#-
#- --score           Run the scoring function on the indicated genome.
#- --genome          Genome object file (*.obj) to score.
#-
#- --mrate=f         Mutation rate, overrides value for mutation_rate given in
#-                   the configuration file.
#- --inumg=i         Initial number of genomes. Overrides the value for
#-                   inum_genomes given in the configuration file.
#- --generation=i    specify which generation need to score
#-
#- --rescore         Run the rescore function to rescore (very) genome from
#-                   recorded obj files.
#
#-#############################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

#######################################################################################
# STANDARD PACKAGES
#######################################################################################
use Getopt::Long;
use English;			# Use english names for global system variables

use 5.008_000;            # require perl version 5.8.0 or higher
use Class::Std 0.0.8;     # require Class::Std version 0.0.8 or higher

#######################################################################################
# APPLICATION PACKAGES
#######################################################################################
use FindBin qw($Bin);

BEGIN {
    if (!defined $ENV{ANC_HOME} || !$ENV{ANC_HOME}) {
        print "ERROR: ANC_HOME environment variable is not defined\n";
        exit(1);
    }
    if (!defined $ENV{FACILE_HOME} || !$ENV{FACILE_HOME}) {
        print "ERROR: FACILE_HOME environment variable is not defined\n";
        exit(1);
    }
}
use lib "$ENV{ANC_HOME}/base";
use lib "$Bin/modules";
use lib "./custom";

use Utils;
use Globals qw (
$VERSION
$RELEASE_DATE
$verbosity
$TAG
$config_ref
$WORKSPACE
);

use BioJazz;

#######################################################################################
# MAIN PROGRAM
#######################################################################################

#======================================================================================
# VERSION AND COPYRIGHT
#======================================================================================

# ANC version
$VERSION = "ALPHA-02";
$RELEASE_DATE = "2012/06/15";

printn "##############################################################################";
printn "# BioJazz -- A biochemical network evolution/design tool";
printn "# Copyright (c) 2005-2012";
printn "# All rights reserved";
printn "# Author: Julien F. Ollivier";
printn "# Version: $VERSION";
printn "# Release Date: $RELEASE_DATE";
printn "##############################################################################";

#======================================================================================
# Don't buffer STDOUT
#======================================================================================

$OUTPUT_AUTOFLUSH = 1;

#======================================================================================
# GLOBALS
#======================================================================================

#======================================================================================
# CMD-LINE ARGUMENT PROCESSING AND DEFAULTS
#======================================================================================
use vars qw($HELP $SHELL $SCRIPT $COMMAND $SEED $GENOME $SCORE $GENERATION $RESCORE);
use vars qw($version_flag);

GetOptions(
    "help"            => \$HELP,
    "version"         => \$version_flag,
    "verbosity=i"     => \$verbosity,
    "config=s"        => \$config_ref->{config_file},
    "command=s"       => \$COMMAND,
    "script=s"        => \$SCRIPT,
    "shell"           => \$SHELL,
    "tag=s"           => \$TAG,
    "seed=i"          => \$SEED,
    "cluster_type=s"  => \$config_ref->{cluster_type},
    "cluster_size=i"  => \$config_ref->{cluster_size},
    "genome=s"        => \$GENOME,
    "generation=i"    => \$GENERATION,
    "score"           => \$SCORE,
    "rescore"         => \$RESCORE,
    "inumg=i"         => \$config_ref->{inum_genomes},
    "mrate=f"         => \$config_ref->{mutation_rate},
);

exit if ($version_flag);

# Print out the header; this also serves as help
if (defined $HELP) {
    my $help_tag = "#-";
    my $OUT = `grep -E '^$help_tag' $PROGRAM_NAME`;
    $OUT =~ s/#-/#/g;
    printn $OUT;
    exit;
}

if (!defined $TAG) {
    $TAG = `date +%F-%T`;
    chomp($TAG);
}
print "BioJazz running with tag \"$TAG\"\n";

#======================================================================================
# ECHO ENVIRONEMENT
#======================================================================================
printn "ANC_HOME=$ENV{ANC_HOME}";
printn "FACILE_HOME=$ENV{FACILE_HOME}";
printn "WORKSPACE=$WORKSPACE";

#======================================================================================
# READ CONFIG FILE
#======================================================================================
if (!defined $config_ref->{config_file}) {
    printn "WARNING: no configuration file was specified";
} else {
    print "BioJazz using configuration file $config_ref->{config_file}\n";
    read_config($config_ref, $config_ref->{config_file}, "NOCLOBBER"); # NOCLOBBER gives cmd-line args priority
}

#======================================================================================
# RUN COMMANDS AND SCRIPTS
#======================================================================================
evolve(seed => defined $SEED ? $SEED : -1) if !$SHELL && !$COMMAND && !$GENOME && !$SCORE && !$GENERATION && !$RESCORE && defined $config_ref->{config_file};

load_genome($GENOME) if defined $GENOME || $SCORE && defined $config_ref->{config_file};

score_genome() if $SCORE && defined $config_ref->{config_file};

score_generation(generation_num => $GENERATION) if defined $GENERATION && defined $config_ref->{config_file};

rescore_genomes() if $RESCORE && defined $config_ref->{config_file};


if (defined $COMMAND) {
    eval("$COMMAND");  warn if $@;
}

if (defined $SCRIPT) {
    interpreter("BIOJAZZ", $SCRIPT);
}

if (defined $SHELL) {
    interpreter("BIOJAZZ");
}

exit;
