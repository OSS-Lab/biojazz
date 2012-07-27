#-#####################################################################################
#- File:     Template.pm
#- Synopsys: Use this package as a starting point for an application-specific
#-           scoring function.  The example score_genome() below shows how to call
#-           a few of the important routines to translate the genome into ODE form,
#-           simulate it in Matlab, and read out simulation results.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Template;
use Class::Std::Storable;
use base qw(Scoring);
{
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG);

    use Stimulus;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
#    my %A_of :ATTR(get => 'A', set => 'A', init_arg => 'A');  # constructor must supply initialization
#    my %B_of :ATTR(get => 'B', set => 'B', default => 'yyy'); # default value is yyy

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: score_genome
    # Synopsys: This is a template for the application-specific score_genome routine
    #           which is user implemented.  This example shows how to perform various
    #           operations and returns a random score.
    #--------------------------------------------------------------------------------------
    sub score_genome {
	my $self = shift; my $obj_ID = ident $self;
	my $genome_model_ref = shift;  # the genome to be scored

	confess "ERROR: internal error, $genome_model_ref not a GenomeModel" if !$genome_model_ref->isa('GenomeModel');

	my $config_ref = $self->get_config_ref();
	my $genome_name = $genome_model_ref->get_name();
	my $work_dir = $self->get_work_dir();
	my $local_dir = $self->get_local_dir();
	my $matlab_work = $self->get_matlab_work();

	printn "Template::score_genome: scoring genome $genome_name";

	my $score = undef;

	#---------------------------------------------------------
	# PARSE GENOME
	#---------------------------------------------------------
	my $genome_iref = $genome_model_ref->parse();
	$genome_model_ref->check();
	printn $genome_iref->sprint(colour_flag => 0);

	#---------------------------------------------------------
	# CREATE/PARSE I/O GENES
	#---------------------------------------------------------
	my $io_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
	    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
	    regulated_concentration => $config_ref->{regulated_concentration_min},
	    UNUSED => "1010",
	    domains => [
		{
		    allosteric_flag => 0,
		    RT_transition_rate => 1e2,
		    TR_transition_rate => 1e-2,
		    RT_phi => 0.55,
		    protodomains => [
			{
			    type => "bsite",
			    substrate_polarity => 1,
			    binding_profile => "101", # chosen to match first gene
			    kf_profile => "0011",
			    kb_profile => "111011",
			    kp_profile => "100011",
			    Keq_ratio => 1.0e-1,
			    kf_polarity_mask => "0000",
			    kb_polarity_mask => "110011",
			    kf_conformation_mask => "1111",
			    kb_conformation_mask => "001100",
			    kp_conformation_mask => "010101",
			    UNUSED => "1100000001",
			},
		       ],
		    UNUSED => "101",
		},
	       ],
	});
	printn "io_sequence=".$io_sequence_ref->get_sequence();
	$genome_model_ref->get_genome_parser_ref->parse(
	    dont_clear_flag => 1,
	    sequence_ref => $io_sequence_ref,
	    prefix => "L",
	   );

	my $stimulus_ref = staircase_equation(
	    #	 my $stimulus_ref = ramp_equation(
	    NODE => "LG_0000_x",
	    PERIOD => 100.0,
	    STRENGTH => 1.0,
	    CONCENTRATION => 1e-3,
	    DUTY => 75,
	    RFTIME =>25,
	    STEPS => 5,
	    DELAY => 10,
	   );
	my ($lg_source_staircase, $lg_sink_staircase) = @{$stimulus_ref->{equations}};
	my $event_list = join " ", @{$stimulus_ref->{events}};
	my $values_list = join " ", @{$stimulus_ref->{values}};
	printn "Stimulus:\n". join "\n",($lg_source_staircase, $lg_sink_staircase, $event_list, $values_list);

	#---------------------------------------------------------
	# TRANSLATE GENOME
	#---------------------------------------------------------
	$genome_model_ref->translate();
	# hack first protodomain to make it an msite
	my @pd_refs = $genome_model_ref->get_protodomain_parser_ref()->get_object_instances();
	$pd_refs[0]->get_translation_ref()->{type} = "msite";
	# hack second protodomain to make it phosphorylate the first msite
	$pd_refs[1]->get_translation_ref()->{substrate_polarity} = 0;

	#---------------------------------------------------------
	# BUILD NETWORK
	#---------------------------------------------------------
	my $genome_ref = $genome_model_ref->get_parser_ref();
	$genome_ref->build_network();
	# PROTODOMAIN CONNECTIVITY
	printn "Protodomains: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{protodomains}});
	printn $genome_ref->get_adjacency_matrix_ref()->{protodomains}->[0]->sprint_matrix();
	printn $genome_ref->get_connectivity_matrix_ref()->{protodomains}->sprint_matrix();
	# GENE CONNECTIVITY
	printn "Genes: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{genes}}) if $verbosity >= 1;
	printn $genome_ref->get_adjacency_matrix_ref()->{genes}->[0]->sprint_matrix();
	printn $genome_ref->get_connectivity_matrix_ref()->{genes}->sprint_matrix();
	# PRUNE NETWORK
	$genome_ref->prune_isolated_genes();

	#---------------------------------------------------------
	# GENERATE ANC/FACILE MODEL
	#---------------------------------------------------------
	my $anc_model = $genome_model_ref->get_genome_parser_ref()->export_anc(
	    max_external_iterations => $config_ref->{max_external_iterations},
	    max_internal_iterations => $config_ref->{max_internal_iterations},
	    max_complex_size => $config_ref->{max_complex_size},
	    max_species => $config_ref->{max_species},
	    max_csite_bound_to_msite_number => $config_ref->{max_csite_bound_to_msite_number},
	    default_steric_factor => $config_ref->{default_steric_factor},
	    equations => [$lg_source_staircase, $lg_sink_staircase],
	    export_graphviz => "network,collapse_states,collapse_complexes",
	    matlab_ode_solver => $config_ref->{solver},
	    matlab_odeset_options => ("odeset('InitialStep', $config_ref->{InitialStep}, ".
				      "'AbsTol', $config_ref->{AbsTol}, ".
				      "'RelTol', $config_ref->{RelTol}, ".
				      "'MaxStep', $config_ref->{MaxStep})"),
	    t_final => $config_ref->{t_final},
	    t_vector =>"[t0:$config_ref->{sampling_interval}:tf]",
	   );
	burp_file("$matlab_work/$genome_name.mod", $anc_model);
	system("$ENV{ANC_HOME}/anc.pl --report=species $matlab_work/$genome_name.mod");

	$self->anc_process_species_report("$matlab_work/$genome_name.species.rpt");
	my @anc_species = $self->anc_get_species();
	printn "ANC NUM SPECIES: ".scalar(@anc_species) if $verbosity >= 1;
	printn "ANC SPECIES: @anc_species" if $verbosity >= 2;

	# check that there was at least one species in the ANC model
	if (@anc_species == 0) {
	    printn "WARNING: no species in ANC file -- setting score to 0";
	    $genome_model_ref->set_score(0);
	    return;
	}

	#---------------------------------------------------------
	# RUN FACILE
	#---------------------------------------------------------
	$self->facile_run(
	    EQN_FILE => "$matlab_work/$genome_name.eqn",
	    SIM_TYPE => "matlab",
	   );

	#---------------------------------------------------------
	# RUN MATLAB SIM
	#---------------------------------------------------------
	printn "Template::score_genome: running matlab driver";
	my $matlab_ref = $self->get_matlab_ref();
	$matlab_ref->cmd("clear all; ${genome_name}Driver");
	$matlab_ref->wait_on("Facile.*done");

	#---------------------------------------------------------
	# RUN MATLAB CUSTOM SCORING FUNCTION
	#---------------------------------------------------------
#	# haven't tested this, but something like:
#	$matlab_ref->cmd("custom/MatlabScoringScript");
#	$score = $matlab_ref->matlab_get_variable("score");
	
	#---------------------------------------------------------
	# READ RAW RESULTS FROM MATLAB AND COMPUTE SCORE
	#---------------------------------------------------------
	# display name and state of first protein at 5s
	my $G_name = $anc_species[0];
	printn "G_name = $G_name";
	my $G_value = $self->matlab_get_state(complex => $G_name, t => 5.0);
	printn "G_value = $G_value";

	# get and display full system's state vector
	my @y = $self->matlab_get_state_vector(t => 8.0);
	printn "y = @y";
	# get and display state vector differential
	my @delta_y = @{$self->matlab_get_state_delta(t1 => 1.0, t2 => 9.0)->{delta}};
	printn "delta_y = @delta_y";

	# find maximum concentration of protein
	my $max_G = $self->matlab_get_max_value($G_name);
	printn "max_G = $max_G";
	$self->matlab_report_max_values();

	# find final concentration of protein
	my $final_G = $self->matlab_get_final_value($G_name);
	printn "final_G = $final_G";
	$self->matlab_report_final_values();

	# plot stimulus and protein waveforms
	$self->matlab_plot_complex(figure => 100,
				   complex => "LG_0000_x",
				   title_prefix => "BOGUS",
				  );
	$self->matlab_plot_complex(figure => 101,
				   complex => $G_name,
				   title_prefix => "BOGUS",
				  );
	system("sleep 5");

	if (defined $config_ref->{plot_species} && $config_ref->{plot_species}) {
	    $self->matlab_plot_all_complexes();
	}

	# generate some random numbers for stats and the score
	printn "random number = ".rand;
	printn "random number = ".rand;
	printn "random number = ".rand;
	printn "random number = ".rand;

	sleep(10);

	$score = int 100*rand;

	$genome_model_ref->set_score($score);

	$genome_model_ref->set_stats_ref({
	    stat1 => int 100*rand,
	    stat2 => int 100*rand,
	});

	#---------------------------------------------------------
	# MOVE FILES from LOCAL_DIR to WORK_DIR
	#---------------------------------------------------------
	if (defined $local_dir) {
	    my $file_glob = "$matlab_work/${genome_name}*";
	    my @files = glob($file_glob);
	    if (@files) {
		printn "Moving @files to $work_dir/matlab";
		system("mv @files $work_dir/matlab");
	    }
	}
    }
}


#--------------------------------------------------------------------------------------
# Function: run_testcases
# Synopsys: This routine should be used to test your scoring function on a random
#           or hand-crafted genome.  Here, we create and save a configuration file,
#           then we create a random genome and score it using the above function.
#--------------------------------------------------------------------------------------
sub run_testcases {
    $verbosity = 1;

    $TAG = "test";
    srand(33433);

    my $config_file = <<END;
#----------------------------------------
# CPU AND CLUSTER SETTINGS
#----------------------------------------
nice = 15
vmem = 750000

#----------------------------------------
# WORKSPACE AND CUSTOM SCORING MODULE
#----------------------------------------
scoring_class = Template
work_dir = template


#----------------------------------------
# GENOME PARAMS
#----------------------------------------
# Genome class
radius = 0
kf_max = 1e5
kf_min = 0.1
kb_max = 10
kb_min = 0.001
kp_max = 1000
kp_min = 1

# Gene class
regulated_concentration_width = 4
gene_unused_width = 32
regulated_concentration_max = 1e-3
regulated_concentration_min = 1e-12

# Domain class
RT_transition_rate_width = 5
TR_transition_rate_width = 6
RT_phi_width = 7
domain_unused_width = 8
RT_transition_rate_max = 1e2
RT_transition_rate_min = 1e-2
TR_transition_rate_max = 1e4
TR_transition_rate_min = 1e-4
RT_phi_max = 0.99
RT_phi_min = 0.01

# ProtoDomain class
binding_profile_width = 3
kf_profile_width = 4
kb_profile_width = 5
kp_profile_width = 6
Keq_profile_width = 7
protodomain_unused_width = 10
Keq_ratio_max = 1e2
Keq_ratio_min = 1e-2

#----------------------------------------
# ANC PARAMS
#----------------------------------------
max_external_iterations = 1
max_internal_iterations = 0
max_complex_size = 4
max_species = 40
max_csite_bound_to_msite_number = 2
default_steric_factor = 1e-3
export_graphviz = nothing
#export_graphviz = network,collapse_states,collapse_complexes

#----------------------------------------
# FACILE/MATLAB SETTINGS
#----------------------------------------
solver = ode23s

sampling_interval = 0.1
t_final = 100

# MATLAB odeset params
InitialStep = 1e-8
AbsTol = 1e-9
RelTol = 1e-3
MaxStep = 10.0

#----------------------------------------
# SIMULATION/SCORING PARAMS
#----------------------------------------
plot_input = 0
plot_output = 0
plot_species = 1
plot_phase = 0
plot_min = -1

END

    burp_file("test/custom/Template.cfg", $config_file);

    my $scoring_ref = Template->new({
	node_ID => 99,
	config_file => "test/custom/Template.cfg",
	work_dir => "test/custom",
	matlab_startup_options => "-nodesktop -nosplash",
    });

    printn $scoring_ref->_DUMP();

    my $config_ref = {};
    read_config($config_ref, "test/custom/Template.cfg");

    use GenomeModel;
    my $genome_model_ref = GenomeModel->new({
	name => "Template",
	Genome => {
	    radius => $config_ref->{radius},
	    kf_max => $config_ref->{kf_max},
	    kf_min => $config_ref->{kf_min},
	    kb_max => $config_ref->{kb_max},
	    kb_min => $config_ref->{kb_min},
	    kp_max => $config_ref->{kp_max},
	    kp_min => $config_ref->{kp_min},
	    Gene => {
		gene_start_code => "10000001",
		soft_linker_code => "0101",
		regulated_concentration_width => $config_ref->{regulated_concentration_width},
		unused_width => $config_ref->{gene_unused_width},
		regulated_concentration_max => $config_ref->{regulated_concentration_max},
		regulated_concentration_min => $config_ref->{regulated_concentration_min},
		Domain => {
		    hard_linker_code => "1001",
		    RT_transition_rate_width => $config_ref->{RT_transition_rate_width},
		    TR_transition_rate_width => $config_ref->{TR_transition_rate_width},
		    RT_phi_width => $config_ref->{RT_phi_width},
		    unused_width => $config_ref->{domain_unused_width},
		    RT_transition_rate_max => $config_ref->{RT_transition_rate_max},
		    RT_transition_rate_min => $config_ref->{RT_transition_rate_min},
		    TR_transition_rate_max => $config_ref->{TR_transition_rate_max},
		    TR_transition_rate_min => $config_ref->{TR_transition_rate_min},
		    RT_phi_max => $config_ref->{RT_phi_max},
		    RT_phi_min => $config_ref->{RT_phi_min},
		    ProtoDomain => {
			binding_profile_width => $config_ref->{binding_profile_width},
			kf_profile_width => $config_ref->{kf_profile_width},
			kb_profile_width => $config_ref->{kb_profile_width},
			kp_profile_width => $config_ref->{kp_profile_width},
			Keq_profile_width => $config_ref->{Keq_profile_width},
			unused_width => $config_ref->{protodomain_unused_width},
			Keq_ratio_max => $config_ref->{Keq_ratio_max},
			Keq_ratio_min => $config_ref->{Keq_ratio_min},
		    },
		},
	    },
	},
    });

    $genome_model_ref->generate_random_genome(2000);
    $scoring_ref->score_genome($genome_model_ref);
    printn $genome_model_ref->_DUMP();

    # save the genome object
    use Storable qw(store retrieve);
    store($genome_model_ref, "test/custom/Template.obj");
}


# Package BEGIN must return true value
return 1;
