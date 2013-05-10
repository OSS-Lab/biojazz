#-#####################################################################################
#- File:     BioJazz.pm
#- Synopsys: Top-level routines that can be called by the BioJazz interpreter.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

#######################################################################################
# TO-DO LIST
#######################################################################################

#######################################################################################
# PACKAGE INTERFACE
#######################################################################################
package BioJazz;

use strict;

use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK);
@ISA = qw(Exporter);
@EXPORT = qw(
$ref
$scoring_ref
create_workspace
evolve
load_genome
score_genome
load_history
save_history
collect_history_from_genomes
collect_history_from_logfile
export_history
);
#@EXPORT_OK = qw(
#		$xxx
#	    );

#######################################################################################
# MODULES USED
#######################################################################################
use FindBin qw($Bin);
use Storable qw(store retrieve);

use Utils;
use Globals qw ($verbosity $TAG $config_ref);

use GenAlg;
use History;

#######################################################################################
# PACKAGE GLOBALS
#######################################################################################
use vars qw(
$ref
$scoring_ref
$history_ref
);

#######################################################################################
# FUNCTIONS
#######################################################################################
#--------------------------------------------------------------------------------------
# Function: create_workspace
# Synopsys: Initialize and create a workspace
#--------------------------------------------------------------------------------------
sub create_workspace {
    # creat the workspace name
    my $workspace_name = shift || "bjazz";

    # creat workspace and cp configure files, custom files and module files into the workspace
    printn "BioJazz is creating your workspace...";
    system("mkdir -p $workspace_name");
    system("mkdir -p $workspace_name/config");
    system("mkdir -p $workspace_name/custom");
    system("mkdir -p $workspace_name/test/modules");
    system("mkdir -p $workspace_name/test/custom");
    system("cp -p $Bin/config/*.cfg $workspace_name/config");
    system("cp -p $Bin/custom/*.pm $workspace_name/custom");
    system("ln -s $Bin/Makefile $workspace_name/Makefile");
    printn "\nDone.  Please cd to $workspace_name before continuing\n";
    exit;
}

#--------------------------------------------------------------------------------------
# Function: evolve
# Synopsys: Initialization, configuration, run evolution.
#--------------------------------------------------------------------------------------
sub evolve {

    my %args = (
        seed => -1,
        @_,
    );

    # pass seed from arguments to seed ????
    check_args(\%args, 1);
    my $seed = $args{seed};

    # START TIME
    my $start_time = `date`;
    chomp($start_time);
    print "Start time: $start_time\n";

    # SIMULATION AND MODEL PARAMETERS
    if (!defined $config_ref->{inum_genomes}) {
        $config_ref->{inum_genomes} = 1;
    }
    print "BioJazz running with initial number of genomes equal to $config_ref->{inum_genomes}\n";

    if (!defined $config_ref->{mutation_rate}) {
        $config_ref->{mutation_rate} = 0.01;
    }
    print "BioJazz running with mutation rate of $config_ref->{mutation_rate}\n";

    if (!defined $config_ref->{cluster_type}) {
        $config_ref->{cluster_type} = "LOCAL";
    }
    print "BioJazz using cluster_type $config_ref->{cluster_type}\n";

    if (!defined $config_ref->{cluster_size}) {
        $config_ref->{cluster_size} = 1;
    }
    print "BioJazz using cluster_size $config_ref->{cluster_size}\n";

    if (!defined $config_ref->{nice}) {
        $config_ref->{nice} = 15;
    }
    print "BioJazz using nice $config_ref->{nice}\n";

    # CLEAN UP PREVIOUS RUNS
    if (defined $config_ref->{remove_old_files} && $config_ref->{remove_old_files} == 1) {
        printn "Removing old files... are you sure (y/n)? ";
        my $answer = <>; chomp($answer);
        if ($answer ne "y") {
            printn "no, then fix the configuration file";
            exit(1);
        }
        `echo $config_ref->{work_dir}/$TAG/*.log     | xargs rm -f`;
        `echo $config_ref->{work_dir}/$TAG/matlab/*  | xargs rm -f`;
        `echo $config_ref->{work_dir}/$TAG/obj/*     | xargs rm -f`;
        `echo $config_ref->{work_dir}/$TAG/source*   | xargs rm -rf`;
        `echo $config_ref->{work_dir}/$TAG/stats/*   | xargs rm -rf`;
        printn "Done removing files.";
    } else {
        printn "Keeping old files in $config_ref->{work_dir}/$TAG";
    }

    # CREATE OUTPUT DIRECTORIES
    printn "Creating output directories";
    `mkdir -p $config_ref->{work_dir}/$TAG`;
    `mkdir -p $config_ref->{work_dir}/$TAG/matlab`;
    `mkdir -p $config_ref->{work_dir}/$TAG/obj`;
    `mkdir -p $config_ref->{work_dir}/$TAG/report`;
    `mkdir -p $config_ref->{work_dir}/$TAG/stats`
    my $timestamp = `date +%F-%T`; chomp($timestamp);
    my $source_dir = "source_$timestamp";
    `mkdir -p $config_ref->{work_dir}/$TAG/$source_dir`;

    # COPY THE PACKAGE, CONFIG AND SOURCE CODE INTO WORKING FOR REFERENCE
    #system("cp -p $config_ref->{package}.pm $config_ref->{work_dir}/$TAG/$source_dir");
    system("cp -p $config_ref->{config_file} $config_ref->{work_dir}/$TAG/$source_dir");
    system("cp -p $Bin/*.pl $Bin/modules/*.pm $config_ref->{work_dir}/$TAG/$source_dir");
    system("cp -p ./custom/*.pm $config_ref->{work_dir}/$TAG/$source_dir");
    #system("bzr status $Bin > $config_ref->{work_dir}/$TAG/$source_dir/bzr_status");

    # START (creates nodes, forks node manager process)
    my $ga_ref = GenAlg->new({
            config_ref => $config_ref,
        });

    # SEED FOR REPEATABLE RANDOM NUMBER SEQUENCE
    # (n.b. seems important to set it after fork)
    #  srand(5334245);
    if (defined $seed && ($seed >= 0)) {
        srand($seed);
        printn("Seed was passed as parameter");
        `echo $start_time -- parameter -- tag=$TAG -- $seed >> seed`;
    } else {
        $seed = int rand (0xFFFFFFFF);
        srand($seed);
        printn("Seed was randomly generated");
        `echo $start_time -- random -- tag=$TAG -- $seed >> seed`;
    }
    printn("Using seed $seed");

    printn "First random number is ".rand;

    # RUN THE GA ROUTINES
    $ga_ref->evolve();

    # TERMINATION
    $ga_ref = undef;

    # CLOSE FILES

    # END TIME
    my $end_time = `date`;
    chomp($end_time);
    print "Start time: $start_time\n";
    print "End time:   $end_time\n";
}

#--------------------------------------------------------------------------------------
# Function: load_genome
# Synopsys: 
#--------------------------------------------------------------------------------------
sub load_genome {
    my $genome = shift;

    # retrieve the genome object
    return $ref = retrieve($genome);
}

#--------------------------------------------------------------------------------------
# Function: score_genome
# Synopsis: scoring genome with config_ref and scoring_ref which defined in configure files
#--------------------------------------------------------------------------------------
sub score_genome {
    eval("use $config_ref->{scoring_class};");
    if ($@) {print $@; return;}

    $scoring_ref = $config_ref->{scoring_class}->new({
            node_ID => 999,
            config_file => $config_ref->{config_file},
            work_dir => "$config_ref->{work_dir}/scratch/$TAG",
            matlab_startup_options => "-nodesktop -nosplash",  # need jvm
        });
    $config_ref = $scoring_ref->get_config_ref();

    $config_ref->{plot_input} = 1;
    $config_ref->{plot_output} = 1;
    $config_ref->{plot_phase} = 1;
    $config_ref->{plot_species} = 0 || $config_ref->{plot_species};
    $config_ref->{export_graphviz} = "network,collapse_states,collapse_complexes";
    $config_ref->{sprint_history} = 1;
    $config_ref->{sprint_transcript} = 1;
    $config_ref->{save_transcript} = 1;

    # set_elite_flag :: need to revisit to modify according to different selection method
    $ref->set_elite_flag(0);
    return $scoring_ref->score_genome($ref);
}

#--------------------------------------------------------------------------------------
# Function: xxx_history
# Synopsys: 
#--------------------------------------------------------------------------------------
sub load_history {
    my $history_obj = "$config_ref->{work_dir}/$TAG/obj/HISTORY.obj";
    printn "loading history from $history_obj";
    $history_ref = retrieve("$history_obj");
}
sub save_history {
    my $history_obj = "$config_ref->{work_dir}/$TAG/obj/HISTORY.obj";
    printn "saving history to $history_obj";
    store($history_ref, "$history_obj");
}
sub collect_history_from_genomes {
    my $obj_dir = "$config_ref->{work_dir}/$TAG/obj";

    $history_ref = History->new({});
    $history_ref->collect_history_from_genomes(
        genomes_dir => $obj_dir,
        max_generations => -1,
    );
}
sub collect_history_from_logfile {
    my $logfile = shift;

    $history_ref = History->new({});
    $history_ref->collect_history_from_logfile(
        logfile => $logfile,
    );
}
# subroutine to export history which could be used to analyse the simulation results?
sub export_history {
    $history_ref->export(
        filename => "$config_ref->{work_dir}/history.$TAG.xls",
        genome_attribute_names => $config_ref->{genome_attribute_names},
        # For kimura selection method there is no need to do population output
        # should put the analysis at the stage after evolution. From csv files 
        # generated "on-the-fly"
        #population_attribute_names => $config_ref->{population_attribute_names},
    );
}

#--------------------------------------------------------------------------------------
# Function: xxx
# Synopsys: 
#--------------------------------------------------------------------------------------
sub xxx {
}

#######################################################################################
# INITIALIZATION
#######################################################################################

# Package BEGIN must return true value
return 1;
