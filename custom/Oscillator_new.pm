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

package Oscillator_new;
use Class::Std::Storable;
use base qw(Scoring);
{
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG $WORKSPACE);

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
        my $stats_ref = $genome_model_ref->get_stats_ref();
        if (!defined $stats_ref) {
            printn "WARNING: stats_ref is not defined for $genome_name";
            $stats_ref = {};
            $genome_model_ref->set_stats_ref($stats_ref);
        }
        my $history_ref = $genome_model_ref->get_history_ref();

        #---------------------------------------------------------
        # CREATE I/O GENES
        #---------------------------------------------------------
#MAKES "GENES" FOR TRANSCRIPTION FACTOR, LIGANDS, OTHER INPUT/OUTPUT SEPARATE FROM EVOLUTION

#LIGAND "GENE"
        my $lg_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
                START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                regulated_concentration => $config_ref->{LG_init}, # all-zeroes
                UNUSED => "0000",
                domains => [
                    {
                        allosteric_flag => 0,
                        RT_transition_rate => $config_ref->{RT_transition_rate_min}, # all-zeroes
                        TR_transition_rate => $config_ref->{TR_transition_rate_min}, # all-zeroes
                        RT_phi => $config_ref->{RT_phi_min}, # all-zeroes
                        protodomains => [
                            {
                                type => "bsite",
                                substrate_polarity => 0,
                                binding_profile => $config_ref->{lg_binding_profile},
                                kf_profile => "0",
                                kb_profile => "0",
                                kp_profile => "0",
                                Keq_ratio => $config_ref->{Keq_ratio_min},
                                kf_polarity_mask => "0",
                                kb_polarity_mask => "0",
                                kf_conformation_mask => "0",
                                kb_conformation_mask => "0",
                                kp_conformation_mask => "0",
                                UNUSED => "0",
                            },
                        ],
                        UNUSED => "0",
                    },
                ],
            });
        printn "lg_sequence=".$lg_sequence_ref->get_sequence() if $verbosity >= 2;

        #TRANSCRIPTION FACTOR "GENE"
        my $tg_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
                START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                regulated_concentration => $config_ref->{TG_init}, # all-zeroes
                UNUSED => "0000",
                domains => [
                    {
                        allosteric_flag => 0,
                        RT_transition_rate => $config_ref->{RT_transition_rate_min}, # all-zeroes
                        TR_transition_rate => $config_ref->{TR_transition_rate_min}, # all-zeroes
                        RT_phi => $config_ref->{RT_phi_min}, # all-zeroes
                        protodomains => [
                            {
                                type => "msite",
                                substrate_polarity => 0,
                                binding_profile => $config_ref->{tg_binding_profile},
                                kf_profile => "0",
                                kb_profile => "0",
                                kp_profile => "0",
                                Keq_ratio => $config_ref->{Keq_ratio_min},
                                kf_polarity_mask => "0",
                                kb_polarity_mask => "0",
                                kf_conformation_mask => "0",
                                kb_conformation_mask => "0",
                                kp_conformation_mask => "0",
                                UNUSED => "0",
                            },
                        ],
                        UNUSED => "0",
                    },
                ],
            });
        printn "tg_sequence=".$tg_sequence_ref->get_sequence() if $verbosity >= 2;

        #---------------------------------------------------------
        # STIMULUS EQUATIONS
        #---------------------------------------------------------
#	my $stimulus_ref = staircase_equation(
#	    #	 my $stimulus_ref = ramp_equation(
#	    NODE => "LG0000_x",
#	    PERIOD => 1000.0,
#	    STRENGTH => 1.0,
#	    CONCENTRATION => 10,
#	    DUTY => 50, #Percent time that stimulus is in stimulated state
#	    RFTIME =>0,#Total time for complete rise/fall
#	    STEPS => 1, #Number of steps in pulse rise/fall
#	    DELAY => 250,
#	   );

        my $stimulus_ref = impulse_equation(
            #	 my $stimulus_ref = ramp_equation(
            NODE => "LG0000_x",
            PERIOD => 1000.0,
            CONCENTRATION => 10,
            IMPULSE_LENGTH => 50,
            DELAY => 250,
        );


        my @lg_equations = @{$stimulus_ref->{equations}};
        my $event_list = join " ", @{$stimulus_ref->{events}};
        my $values_list = join " ", @{$stimulus_ref->{values}};
        printn "Stimulus:\n". join "\n",(@lg_equations, $event_list, $values_list);

        #---------------------------------------------------------
        # PARSE/TRANSLATE GENOME AND I/O GENES
        #---------------------------------------------------------
        my $genome_iref = $genome_model_ref->parse(
            [
                sequence_ref => $lg_sequence_ref,
                prefix => "L",
            ],
            [
                sequence_ref => $tg_sequence_ref,
                prefix => "T",
            ],
        );
        my $parse_successful = $stats_ref->{parse_successful} = $genome_model_ref->check();

        my $history = $genome_model_ref->sprint_history(10);
        printn $history if $verbosity >= 2 || $config_ref->{sprint_history};

        if ($parse_successful) {
            my $transcript = $genome_iref->sprint(colour_flag => 0);
            printn $transcript if $verbosity >= 2 || $config_ref->{sprint_transcript};
            burp_file("$matlab_work/$genome_name.tsc", "$history\n$transcript") if $config_ref->{save_transcript};
            $genome_model_ref->translate();

            #---------------------------------------------------------
            # BUILD/PRUNE NETWORK
            #---------------------------------------------------------
            my $genome_ref = $genome_model_ref->get_parser_ref();
            $genome_ref->build_network();

            # REPORT PROTODOMAIN CONNECTIVITY
            printn "Protodomains: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{protodomains}});
            printn $genome_ref->get_adjacency_matrix_ref()->{protodomains}->[0]->sprint_matrix();
            printn $genome_ref->get_connectivity_matrix_ref()->{protodomains}->sprint_matrix();
            # REPORT GENE CONNECTIVITY
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
                equations => [@lg_equations],
                export_graphviz => "network,collapse_states,collapse_complexes",
                matlab_ode_solver => $config_ref->{solver},
                matlab_odeset_options => ("odeset('InitialStep', $config_ref->{InitialStep}, ".
                    "'AbsTol', $config_ref->{AbsTol}, ".
                    "'RelTol', $config_ref->{RelTol}, ".
                    "'MaxStep', $config_ref->{MaxStep})"),
                t_final => $config_ref->{t_final},
                t_vector =>"[t0:$config_ref->{sampling_interval}:tf]",
                event_times => "$event_list",
            );
            burp_file("$matlab_work/$genome_name.mod", $anc_model);
            system("$ENV{ANC_HOME}/anc.pl --report=species $matlab_work/$genome_name.mod");

            $self->anc_process_species_report("$matlab_work/$genome_name.species.rpt");
            my @anc_species = $self->anc_get_species();
            printn "ANC NUM SPECIES: ".scalar(@anc_species) if $verbosity >= 1;
            printn "ANC SPECIES: @anc_species" if $verbosity >= 2;

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
            printn "Template::score_genome: running custom Matlab script" if $verbosity >= 1;
            $matlab_ref->cmd("score = MatlabScoringScript(t,LG0000_x, TG0000_1);");
            $score = $self->matlab_get_variable(name => "score");
            $genome_model_ref->set_score($score);

            #---------------------------------------------------------
            # PLOT RESULTS
            #---------------------------------------------------------
            if (defined $config_ref->{plot_input} && $config_ref->{plot_input}) {
                $self->matlab_plot_complex(figure => 900,
                    complex => "LG0000_x",
                    title_prefix => "$genome_name",
                );
            }
            if (defined $config_ref->{plot_output} && $config_ref->{plot_output}) {
                $self->matlab_plot_complex(figure => 901,
                    complex => "TG0000_0",
                    title_prefix => "$genome_name",
                );
                $self->matlab_plot_complex(figure => 902,
                    complex => "TG0000_1",
                    title_prefix => "$genome_name",
                );
            }
            if (defined $config_ref->{plot_species} && $config_ref->{plot_species}) {
                $self->matlab_plot_all_complexes();
            }
            $self->matlab_cmd("disp('Done plotting')\n");
            $self->matlab_wait_on("Done plotting");
            system("sleep 1");

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

            #	$genome_model_ref->set_stats_ref({
            #	    stat1 => int 100*rand,
            #	    stat2 => int 100*rand,
            #	});
        } else {
            # failed to parse, so set score to 0
            $genome_model_ref->set_score(0);
        }

        printn "final score=".$genome_model_ref->get_score();

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
#           then we create a hand-crafted genome and score it using the above function.
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
vmem = 450000

#----------------------------------------
# WORKSPACE AND CUSTOM SCORING MODULE
#----------------------------------------
scoring_class = Oscillator_new
work_dir = oscillator_new

#----------------------------------------
# GENOME PARAMS
#----------------------------------------
# Scaling: All concentrations in uM and rates in s^-1.
#          Hence all 2nd-order rates are in uM^-1 s^-1.

# Genome class
radius = 1
kf_max = 1e3    # uM^-1 s^-1
kf_min = 1e-3
kb_max = 1e3
kb_min = 1e-3
kp_max = 1e3
kp_min = 1e-3

# Gene class
regulated_concentration_width = 10
gene_unused_width = 4
regulated_concentration_max = 1e3    # 1mM
regulated_concentration_min = 1e-3   # 1nM  ~ 1 molecule in prokaryote

# Domain class
RT_transition_rate_width = 10
TR_transition_rate_width = 10
RT_phi_width = 10
domain_unused_width = 4
RT_transition_rate_max = 1e2
RT_transition_rate_min = 1e-2
TR_transition_rate_max = 1e2
TR_transition_rate_min = 1e-2
RT_phi_max = 1.0
RT_phi_min = 0.0

# ProtoDomain class
binding_profile_width = 8
kf_profile_width = 20
kb_profile_width = 20
kp_profile_width = 10
Keq_profile_width = 10
protodomain_unused_width = 4
Keq_ratio_max = 1e2
Keq_ratio_min = 1e-2

#----------------------------------------
# ANC PARAMS
#----------------------------------------
max_external_iterations = -1
max_internal_iterations = -1
max_complex_size = 4
max_species = 160
max_csite_bound_to_msite_number = 1
default_max_count = 2          # this prevents polymerization (see ANC manual)
default_steric_factor = 1e3
#export_graphviz = nothing
export_graphviz = network,collapse_states,collapse_complexes

#----------------------------------------
# GENERIC SIMULATION/SCORING PARAMETERS
#----------------------------------------
plot_input = 1
plot_output = 1
plot_species = 0
plot_phase = 1
plot_min = -1

sprint_transcript = 1
sprint_history = 1

# MATLAB solver
solver = ode23s

# MATLAB odeset params
InitialStep = 1e-8
AbsTol = 1e-9
RelTol = 1e-3
MaxStep = 1.0

# MATLAB sampling interval
sampling_interval = 0.1

#----------------------------------------
# APPLICATION-SPECIFIC PARAMETERS
#----------------------------------------
t_final = 1000

lg_binding_profile = 01000011
tg_binding_profile = 10100011

LG_init = 1e-3
TG_init = 200

END

    burp_file("test/modules/Oscillator_new.cfg", $config_file);

    my $scoring_ref = Oscillator_new->new({
            node_ID => 97,
            config_file => "test/modules/Oscillator_new.cfg",
            work_dir => "test/custom",
            matlab_startup_options => "-nodesktop -nosplash",
        });

    printn $scoring_ref->_DUMP();

    my $config_ref = {};
    read_config($config_ref, "test/modules/Oscillator_new.cfg");

    use GenomeModel;
    my $genome_model_ref = GenomeModel->new({
            name => "Oscillator_new",
            Genome => {
                radius => $config_ref->{radius},
                kf_max => $config_ref->{kf_max},
                kf_min => $config_ref->{kf_min},
                kb_max => $config_ref->{kb_max},
                kb_min => $config_ref->{kb_min},
                kp_max => $config_ref->{kp_max},
                kp_min => $config_ref->{kp_min},
                Gene => {
                    regulated_concentration_width => $config_ref->{regulated_concentration_width},
                    unused_width => $config_ref->{gene_unused_width},
                    regulated_concentration_max => $config_ref->{regulated_concentration_max},
                    regulated_concentration_min => $config_ref->{regulated_concentration_min},
                    Domain => {
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

    # CONFIGURE/CREATE GENOME
    my $lg_binding_profile = $config_ref->{lg_binding_profile};
    my $tg_binding_profile = $config_ref->{tg_binding_profile};
    my $sequence_ref = $genome_model_ref->get_genome_parser_ref()->create_sequence({
            PRE_JUNK => undef, POST_JUNK => "0000",
            genes => [
                {
#PHOSPHATASE
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 2.0, # uM
                    UNUSED => "0000",
                    domains => [
                        {
                            allosteric_flag => 1,
                            RT_transition_rate => 0.01,
                            TR_transition_rate => 1.0,
                            RT_phi => 1.0,
                            protodomains => [
                                {
                                    type => "bsite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement($lg_binding_profile)->sprint(),
                                    kf_profile => "00000000000000000000",
                                    kb_profile => "11111000000000000000",
                                    kp_profile => "11111",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111111111111111111",
                                    kb_conformation_mask => "00000111110000000000",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                                {
                                    type => "csite",
                                    substrate_polarity => 1,
                                    binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
                                    kf_profile => "11111111111111111111",
#				kf_profile => "11111111111111111110",
                                    kb_profile => "11111000000000000000",
                                    kp_profile => "1110000011",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0", #XOR'd with original kf profile when phosphorylated
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111111111111100000", #XOR'd with original kf profile when T
                                    kb_conformation_mask => "11111000000000000000",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                            ],
                            UNUSED => "0",
                        },
                    ],
                },
                {
#KINASE
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 1.0, # uM
                    UNUSED => "0000",
                    domains => [
                        {
                            allosteric_flag => 0,
                            RT_transition_rate => 1.0,
                            TR_transition_rate => 1.0,
                            RT_phi => 0.0,
                            protodomains => [
                                {
                                    type => "csite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
                                    kf_profile => "11111111111111111111",  # gives kf=1000, kb=31.6, kp=186 hence Km=~0.2
                                    kb_profile => "11111000000000000000",
                                    kp_profile => "1110000011",
                                    Keq_ratio => 2.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "0",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                            ],
                            UNUSED => "0",
                        },
                    ],
                },
            ],
        });
    printn "sequence=".$sequence_ref->get_sequence();
    $genome_model_ref->set_sequence_ref($sequence_ref);

    # save the genome object
    use Storable qw(store retrieve);
    store($genome_model_ref, "test/custom/Oscillator_new.obj");

    $scoring_ref->score_genome($genome_model_ref);
    printn $genome_model_ref->_DUMP();

    sleep 60;
}


# Package BEGIN must return true value
return 1;
