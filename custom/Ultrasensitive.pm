#-#####################################################################################
#- File:     Ultrasensitive.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Ultrasensitive;
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
    # Synopsys: Score for ultrasensitivity.
    #--------------------------------------------------------------------------------------
    sub score_genome {
        my $self = shift;
        my $genome_model_ref = shift;

        confess "ERROR: internal error, $genome_model_ref not a GenomeModel" if !$genome_model_ref->isa('GenomeModel');

        my $config_ref = $self->get_config_ref();
        my $genome_name = $genome_model_ref->get_name();
        my $work_dir = $self->get_work_dir();
        my $local_dir = $self->get_local_dir();
        my $matlab_work = $self->get_matlab_work();

        my $stats_ref = $genome_model_ref->get_stats_ref();
        if (!defined $stats_ref) {
            printn "WARNING: stats_ref is not defined for $genome_name";
            $stats_ref = {};
            $genome_model_ref->set_stats_ref($stats_ref);
        }
        my $history_ref = $genome_model_ref->get_history_ref();

        printn "Ultrasensitive::score_genome scoring genome $genome_name";

        #---------------------------------------------------------
        # INIT SCORING
        #---------------------------------------------------------
        my $elite_flag = $genome_model_ref->get_elite_flag();
        if ($elite_flag) {
            printn "Ultrasensitive::score_genome elite individual already scored, previous score=$stats_ref->{score}";
            return if ($config_ref->{rescore_elite} == 0);
            printn "Ultrasensitive::score_genome re-scoring elite individual";

            # Should clear stats in the evolution to prevent stats loss when rescore genomes.
            # $genome_model_ref->clear_stats(preserve => ["created_kinase", "created_phosphatase"]);
            $stats_ref->{score_count} = 1;
            $stats_ref->{score} = 0;
        } else {
            printn "Ultrasensitive::score_genome scoring non-elite individual...";
 
            # Should clear stats in the evolution to prevent stats loss when rescore genomes.
            # $genome_model_ref->clear_stats(preserve => ["created_kinase", "created_phosphatase"]);
            $stats_ref->{score_count} = 1;
            $stats_ref->{score} = 0;
        }

        #---------------------------------------------------------
        # CREATE I/O GENES
        #---------------------------------------------------------
        my $lg_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
                START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                regulated_concentration => $config_ref->{regulated_concentration_min}, # all-zeroes
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
        # STIMULUS/SAMPLING EQUATIONS
        #---------------------------------------------------------
        # if random delay, delay from 1/4 period to 1/2 period inclusive
        my $stimulus_sub_ref = \&{$config_ref->{stimulus}};
        my $stimulus_ref = undef;
        if ($config_ref->{stimulus} eq "ss_ramp_equation") {
            $stimulus_ref = &$stimulus_sub_ref(
                NODE => "LG0000",
                DELAY => $config_ref->{LG_delay},	
                RANGE => $config_ref->{LG_range},
                STRENGTH => $config_ref->{LG_strength},
                RAMP_TIME => $config_ref->{LG_ramp_time},
                STEPS => $config_ref->{LG_steps},
            );
        } elsif ($config_ref->{stimulus} eq "staircase_equation") {
            $stimulus_ref = &$stimulus_sub_ref(
                NODE => "LG0000",
                PERIOD => $config_ref->{LG_period},
                DELAY => $config_ref->{LG_delay},
                STRENGTH => $config_ref->{LG_strength},
                CONCENTRATION => $config_ref->{LG_concentration},
                DUTY => $config_ref->{LG_duty},
                RFTIME => $config_ref->{LG_rftime},
                STEPS => $config_ref->{LG_steps},
            );
        } else {
            confess "ERROR: unknown stimulus subroutine";
        }
        my ($lg_source_eqn, $lg_sink_eqn) = @{$stimulus_ref->{equations}};
        my @stimulus_event_times = @{$stimulus_ref->{events}};
        my @stimulus_values_list = @{$stimulus_ref->{values}};
        if ($verbosity >= 2) {
            printn "Stimulus:";
            printn $lg_source_eqn;
            printn $lg_sink_eqn;
            printn join " ", @stimulus_event_times;
            printn join " ", @stimulus_values_list;
        }

        my @input_vector = @stimulus_values_list;

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

        #############################################################################
        my $network_connectivity = $stats_ref->{network_connectivity} = 0;
        if ($parse_successful) {
            my $transcript = $genome_iref->sprint(colour_flag => 0);
            printn $transcript if $verbosity >= 2 || $config_ref->{sprint_transcript};
            burp_file("$matlab_work/$genome_name.tsc", "$history\n$transcript") if $config_ref->{save_transcript};
            $genome_model_ref->translate();

            #---------------------------------------------------------
            # BUILD/PRUNE NETWORK
            #---------------------------------------------------------
            # BUILD NETWORK
            my $genome_ref = $genome_model_ref->get_parser_ref();
            $genome_ref->build_network();

            # REPORT PROTODOMAIN CONNECTIVITY
            printn "Protodomains: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{protodomains}});
            printn $genome_ref->get_adjacency_matrix_ref()->{protodomains}->[0]->sprint_matrix() if $verbosity >= 2;
            printn $genome_ref->get_connectivity_matrix_ref()->{protodomains}->sprint_matrix() if $verbosity >= 2;
            # REPORT GENE CONNECTIVITY
            printn "Genes: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{genes}}) if $verbosity >= 1;
            printn $genome_ref->get_adjacency_matrix_ref()->{genes}->[0]->sprint_matrix() if $verbosity >= 2;
            printn $genome_ref->get_connectivity_matrix_ref()->{genes}->sprint_matrix() if $verbosity >= 2;

            # PRUNE GENES
            $stats_ref->{num_pruned_genes} = scalar $genome_ref->prune_isolated_genes();

            my $gene_ref = $genome_model_ref->get_gene_parser_ref();
            my $lg_gene_ref = $gene_ref->lookup_object_instance_by_name("LG0000");
            my $tg_gene_ref = $gene_ref->lookup_object_instance_by_name("TG0000");
            my $protodomain_ref = $genome_model_ref->get_protodomain_parser_ref();
            my ($lg_protodomain_ref)= $protodomain_ref->grep_instances_by_name("LPD");
            my ($tg_protodomain_ref)= $protodomain_ref->grep_instances_by_name("TPD");

            #---------------------------------------------------------
            # SCORING: 2 pts -- LG/TG connected to anything
            #---------------------------------------------------------
            if ($lg_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "LG is connected" if $verbosity >= 1;
            }
            if ($tg_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "TG is connected" if $verbosity >= 1;
            }

            #---------------------------------------------------------
            # SCORING: 90 + 400 pts -- LG/TG subnets
            #---------------------------------------------------------
            if ($network_connectivity == 2) { # LG/TF connected
                my (@lg_subnet, @tg0_subnet, @tg1_subnet);   # protodomain subnets
                @lg_subnet = $genome_ref->get_connected(key => "protodomains", ref => $lg_protodomain_ref);
                @tg0_subnet = $genome_ref->get_connected(key => "protodomains", ref => $tg_protodomain_ref, state => 0);
                @tg1_subnet = $genome_ref->get_connected(key => "protodomains", ref => $tg_protodomain_ref, state => 1);

                printn "LG protodomain connects to ".join ",",   (map {$_->[2]} @lg_subnet) if $verbosity >= 1;
                printn "TG/0 protodomain connects to ".join ",", (map {$_->[2]} @tg0_subnet) if $verbosity >= 1;
                printn "TG/1 protodomain connects to ".join ",", (map {$_->[2]} @tg1_subnet) if $verbosity >= 1;

                # max 90 points for subnet size
                my $lg_subnet_size = (@lg_subnet > 30) ? 30 : @lg_subnet;
                my $tg0_subnet_size = (@tg0_subnet > 30) ? 30 : @tg0_subnet;
                my $tg1_subnet_size = (@tg1_subnet > 30) ? 30 : @tg1_subnet;
                $network_connectivity += ($lg_subnet_size + $tg0_subnet_size + $tg1_subnet_size);

                #########################################################################################
                # score -- LG/TG connected to each other
                if (grep /LPD/, (map {$_->[2]} @tg0_subnet)) {
                    printn "TG/0 fans out to LG" if $verbosity >= 1;
                    $network_connectivity += 100;
                }
                if (grep /LPD/, (map {$_->[2]} @tg1_subnet)) {
                    printn "TG/1 fans out to LG" if $verbosity >= 1;
                    $network_connectivity += 100;
                }
                if (grep /TPD\d*\/0/, (map {$_->[2]} @lg_subnet)) {
                    printn "LG fans out to TG/0" if $verbosity >= 1;
                    $network_connectivity += 100;
                }
                if (grep /TPD\d*\/1/, (map {$_->[2]} @lg_subnet)) {
                    printn "LG fans out to TG/1" if $verbosity >= 1;
                    $network_connectivity += 100;
                }

                $network_connectivity = 500 if $network_connectivity >= 400;  # max LG/TG connectivity score
            }

            my @adjacent_kinases = $genome_ref->find_adjacent_csites($tg_protodomain_ref, 0);
            my @adjacent_phosphatases = $genome_ref->find_adjacent_csites($tg_protodomain_ref, 1);
            $stats_ref->{num_adjacent_kinases} = scalar(@adjacent_kinases);
            $stats_ref->{num_adjacent_phosphatases} = scalar(@adjacent_phosphatases);
            printn "Found ".@adjacent_kinases." adjacent kinases";
            printn "Found ".@adjacent_phosphatases." adjacent phosphatases";

            # to speed up evolution, we can cheat: if network_connectivity >= 200, and  < 500,
            # this means LG is connected to either TG/0, TG/1, but not both.
            # provided there are at least 2 connected proteins,
            # let's flip some bits to create the necessary enzymes out
            # of the adjacent proteins.  however, we will do this only
            # once per lineage by setting the persistent created_kinase/phosphatase flags
            my @tg_adjacent_genes = $genome_ref->get_adjacent(key => "genes", ref => $tg_gene_ref);
            my $tg_adjacent_genes = $stats_ref->{tg_adjacent_genes} = scalar @tg_adjacent_genes;
            if ($tg_adjacent_genes >= 2) {
                printn "TG has > 2 connected proteins" if $verbosity >= 1;
            }
            if ($network_connectivity >= 200 && $network_connectivity < 500 && $tg_adjacent_genes >= 2) {
                printn "Editing genome to create kinase/phosphatase";
                confess "ERROR: unexpected condition" if (@adjacent_kinases && @adjacent_phosphatases);
                #---------------------------------------------------------
                # PICK AN ADJACENT PROTEIN TO MAKE INTO KINASE/PHOSPHATASE
                #---------------------------------------------------------

                my @adjacent_protodomains = union(
                    [map {$_->[0]} $genome_ref->get_adjacent(key => "protodomains", ref => $tg_protodomain_ref, state => 0)],
                    [map {$_->[0]} $genome_ref->get_adjacent(key => "protodomains", ref => $tg_protodomain_ref, state => 1)],
                );
                @adjacent_protodomains = simple_difference(
                    \@adjacent_protodomains,
                    [$lg_protodomain_ref]
                );
                @adjacent_protodomains = simple_difference(
                    \@adjacent_protodomains,
                    [$tg_protodomain_ref]
                );

                my $csite_field_value = $tg_protodomain_ref->untranslate_field("type", "csite");
                my $modified_flag = 0;
                if (@adjacent_kinases < 1 && !defined $genome_model_ref->get_stat("created_kinase")) {
                    my @protodomains = (@adjacent_phosphatases ? # preserve one phosphatase
                        simple_difference(\@adjacent_protodomains, [$adjacent_phosphatases[int rand @adjacent_phosphatases]]) :
                        @adjacent_protodomains);
                    if (@protodomains >= 1) {
                        my $index = int rand @protodomains;
                        my $protodomain_ref = $protodomains[$index];
                        printn "Editing genome to create kinase (".$protodomain_ref->get_name().")";
                        $protodomain_ref->set_field("type", $csite_field_value);
                        $protodomain_ref->set_field("substrate_polarity", 0);
                        @adjacent_kinases = (@adjacent_kinases, $protodomain_ref); # update list
                        $stats_ref->{created_kinase} += 1;  # !!! store name of containing gene instead ???
                        $genome_model_ref->add_history("Edited genome protodomain ".$protodomain_ref->get_name()." to create kinase csite");
                        $modified_flag = 1;
                    }
                }
                if (@adjacent_phosphatases < 1 && !defined $genome_model_ref->get_stat("created_phosphatase")) {
                    my @protodomains = (@adjacent_kinases ? # preserve one kinase
                        simple_difference(\@adjacent_protodomains, [$adjacent_kinases[int rand @adjacent_kinases]]) :
                        @adjacent_protodomains);
                    if (@protodomains >= 1) {
                        my $index = int rand @protodomains;
                        my $protodomain_ref = $protodomains[$index];
                        printn "Editing genome to create phosphatase (".$protodomain_ref->get_name().")";
                        $protodomain_ref->set_field("type", $csite_field_value);
                        $protodomain_ref->set_field("substrate_polarity", 1);
                        @adjacent_phosphatases = (@adjacent_phosphatases, $protodomain_ref); # update list
                        $stats_ref->{created_phosphatase} += 1;
                        $genome_model_ref->add_history("Edited genome protodomain ".$protodomain_ref->get_name()." to create phosphatase csite");
                        $modified_flag = 1;
                    }
                }

                if ($modified_flag) {
                    #---------------------------------------------------------
                    # RE-PARSE/TRANSLATE GENOME AND I/O GENES
                    #---------------------------------------------------------
                    $genome_iref = $genome_model_ref->parse(
                        [
                            sequence_ref => $lg_sequence_ref,
                            prefix => "L",
                        ],
                        [
                            sequence_ref => $tg_sequence_ref,
                            prefix => "T",
                        ],
                    );
                    $genome_model_ref->check();
                    $genome_model_ref->translate();
                }

                #---------------------------------------------------------
                # SCORING: if we have adjacent kinases/phosphatases, we
                #          should now have a score that reflects this
                #---------------------------------------------------------
                if (@adjacent_kinases && @adjacent_phosphatases) {
                    $network_connectivity = 500;
                }
            }

            if ($network_connectivity >= 500) {
                $stats_ref->{network_connected_flag} = 1;

                # exclude proteins not in LG/TG subnet from export
                my (@lg_subnet, @tg_subnet);   # gene subnets
                @lg_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $lg_gene_ref);
                @tg_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $tg_gene_ref);
                printn "LG protein connects to ".join ",",   (map {$_->get_name} @lg_subnet) if $verbosity >= 1;
                printn "TG protein connects to ".join ",", (map {$_->get_name} @tg_subnet) if $verbosity >= 1;

                my @proteins_not_in_subnet = simple_difference([$genome_model_ref->get_genes()], [union(\@lg_subnet, \@tg_subnet)]);
                map {$_->set_export_flag(0)} @proteins_not_in_subnet;

                #---------------------------------------------------------
                # GENERATE ANC/FACILE MODEL
                #---------------------------------------------------------
                my $anc_model = $genome_model_ref->get_genome_parser_ref()->export_anc(
                    max_external_iterations => $config_ref->{max_external_iterations},
                    max_internal_iterations => $config_ref->{max_internal_iterations},
                    max_complex_size => $config_ref->{max_complex_size},
                    max_species => $config_ref->{max_species},
                    max_csite_bound_to_msite_number => $config_ref->{max_csite_bound_to_msite_number},
                    default_max_count => $config_ref->{default_max_count},
                    default_steric_factor => $config_ref->{default_steric_factor},
                    export_graphviz => ref $config_ref->{export_graphviz} ? (join ",",@{$config_ref->{export_graphviz}}) : $config_ref->{export_graphviz},
                    equations => [$lg_source_eqn, $lg_sink_eqn],
                    matlab_ode_solver => $config_ref->{solver},
                    matlab_odeset_options => ("odeset('InitialStep', $config_ref->{InitialStep}, ".
                        "'AbsTol', $config_ref->{AbsTol}, ".
                        "'RelTol', $config_ref->{RelTol}, ".
                        "'MaxStep', $config_ref->{MaxStep})"),
                    t_final => $config_ref->{LG_timeout},
                    t_vector =>"[t0:$config_ref->{sampling_interval}:tf]",
                    ode_event_times => (join " ", @stimulus_event_times),
                    SS_timescale => $config_ref->{SS_timescale},
                );
                burp_file("$matlab_work/$genome_name.mod", $anc_model);
                system("$ENV{ANC_HOME}/anc.pl --report=species $matlab_work/$genome_name.mod");
                my @facile_model = slurp_file("$matlab_work/$genome_name.eqn");

                $self->anc_process_species_report("$matlab_work/$genome_name.species.rpt");
                my @anc_species = $self->anc_get_species();
                $stats_ref->{num_anc_species} = @anc_species;
                printn "ANC NUM SPECIES: ".scalar(@anc_species) if $verbosity >= 1;
                printn "ANC SPECIES: @anc_species" if $verbosity >= 2;

                #---------------------------------------------------------
                # RUN FACILE
                #---------------------------------------------------------
                my $sampling_interval = $config_ref->{sampling_interval};
                $self->facile_run(
                    EQN_FILE => "$matlab_work/$genome_name.eqn",
                    SIM_TYPE => "matlab",
                );

                ###############################################################################
                #---------------------------------------------------------
                # SCORE COMPLEXITY
                #---------------------------------------------------------
#		my $line_count = `wc -l $matlab_work/$genome_name.eqn`;
#		$line_count =~ /^\s*(\S+)\s*/; # extract the line count
#		$stats_ref->{complexity} = $1;
                my $num_protodomains = @{[$anc_model =~ /ProtoDomain :/g]};
                my $num_domains = @{[$anc_model =~ /\sDomain :/g]};
                my $num_proteins = @{[$anc_model =~ /Protein :/g]};
                my $num_rules = @{[$anc_model =~ /CanBindRule :/g]};
                printn "ANC model complexity: $num_protodomains + $num_domains + $num_proteins + $num_rules";
                $stats_ref->{complexity} = $num_protodomains + $num_domains + $num_proteins + $num_rules;
                $stats_ref->{complexity_score} = n_hill($stats_ref->{complexity}, 100, 1);

                #---------------------------------------------------------
                # CHECK ANC/FACILE MODEL
                #---------------------------------------------------------
                # check that TF_0 and TF_1 are products of at least 1 reaction each
                my $num_reactions_tg_1 = grep (/(->|<-).* TG00001/, @facile_model);
                my $num_reactions_tg_0 = grep (/(->|<-).* TG00000/, @facile_model);
                $stats_ref->{num_reactions_tg_0} = $num_reactions_tg_0;
                $stats_ref->{num_reactions_tg_1} = $num_reactions_tg_1;
                $network_connectivity += 100 * ($num_reactions_tg_1 > 1 ? 1 : $num_reactions_tg_1);
                $network_connectivity += 100 * ($num_reactions_tg_0 > 1 ? 1 : $num_reactions_tg_0);

                # check that number of species is less than maximum
                if ($stats_ref->{num_anc_species} < $config_ref->{max_species}) {
                    $network_connectivity += 100;   # didn't go over max no. of species
                }
            }
            my $ANC_ok_flag = $stats_ref->{ANC_ok_flag} = ($network_connectivity >= 800) ? 1 : 0;

            #########################################################################################
            #---------------------------------------------------------
            # NETWORK SIMULATION
            #---------------------------------------------------------
            # sim_flag indicates that network was successfully simulated
            # and that calculated results are valid
            $stats_ref->{sim_flag} = 0;
            if ($ANC_ok_flag) {
                $stats_ref->{sim_flag} = 1;

                #---------------------------------------------------------
                # RUN MATLAB SIM
                #---------------------------------------------------------
                printn "Ultrasensitive::score_genome: running matlab driver...";
                my $matlab_ref = $self->get_matlab_ref();
                $matlab_ref->cmd("clear all; ${genome_name}Driver");
                $matlab_ref->wait_on("Facile.*done");

                my @event_times = $self->matlab_get_variable(name => "event_times");
                my @event_flags = $self->matlab_get_variable(name => "event_flags");
                my $timeout_flag = $stats_ref->{timeout_flag} = (
                    ($event_times[$#event_times] == $config_ref->{LG_timeout}) ||
                    (grep {$_ != 1.0} @event_flags) != 0
                ) ? 1 : 0;

                #---------------------------------------------------------
                # PLOT RESULTS
                #---------------------------------------------------------
                if (defined $config_ref->{plot_input} && $config_ref->{plot_input}) {
                    $self->matlab_plot_complex(figure => 900,
                        complex => "LG0000",
                        title_prefix => "$genome_name",
                    );
                }
                if (defined $config_ref->{plot_output} && $config_ref->{plot_output}) {
                    $self->matlab_plot_complex(figure => 901,
                        complex => "TG00000",
                        title_prefix => "$genome_name",
                    );
                    $self->matlab_plot_complex(figure => 902,
                        complex => "TG00001",
                        title_prefix => "$genome_name",
                    );
                }
                if (defined $config_ref->{plot_species} && $config_ref->{plot_species}) {
                    $self->matlab_plot_all_complexes();
                }
                $self->matlab_cmd("disp('Done plotting')\n");
                $self->matlab_wait_on("Done plotting");

                #########################################################################
                #---------------------------------------------------------
                # SCORE STEADY STATE
                #---------------------------------------------------------
                printn "computing steady-state slopes...";

                my $steady_state_threshold = $config_ref->{steady_state_threshold};
                $stats_ref->{steady_state_score} = n_hill(
                    $event_times[1], $steady_state_threshold,1);

                #---------------------------------------------------------
                # REPORT RESULT VECTOR
                #---------------------------------------------------------
                printn "RESULT VECTOR: INPUT = LG OUTPUT = TG00001 DELAY = $config_ref->{LG_delay}";
                my (@pos_output_vector, @neg_output_vector, @expected_output_vector);

                my @sampling_times = @event_times;
                for (my $i=0; $i < @sampling_times; $i++) {
                    my $t = $sampling_times[$i]; 
                    $pos_output_vector[$i] = $self->matlab_get_state(complex => "TG00001", t => $t);
                    $neg_output_vector[$i] = $self->matlab_get_state(complex => "TG00000", t => $t);
                    $expected_output_vector[$i] = (
                        $config_ref->{TG_init} *
                        p_hill($input_vector[$i], $config_ref->{hill_k}, $config_ref->{hill_n})); # positive hill function

                    if ($config_ref->{round_values_flag}) {
                        $pos_output_vector[$i] = $self->matlab_round_value(
                            value => $pos_output_vector[$i],
                            AbsTol => $config_ref->{AbsTol},
                            RelTol => $config_ref->{RelTol},
                        );
                        $neg_output_vector[$i] = $self->matlab_round_value(
                            value => $neg_output_vector[$i],
                            AbsTol => $config_ref->{AbsTol},
                            RelTol => $config_ref->{RelTol},
                        );
                        $expected_output_vector[$i] = $self->matlab_round_value(
                            value => $expected_output_vector[$i],
                            RelTol => $config_ref->{RelTol},
                        );
                    }

                    printf("RESULT VECTOR:  t=%-6.2f input vector:  %8.4g pos_output_vector: %8.6g neg_output_vector: %8.6g (expected = %8.4g)\n",
                        $t, $input_vector[$i], $pos_output_vector[$i], $neg_output_vector[$i], $expected_output_vector[$i]);
                }

                if (!$timeout_flag) {
                    #---------------------------------------------------------
                    # SELECT BEST OUTPUT VECTOR
                    #---------------------------------------------------------
                    my $i_dy2_bottom = ($config_ref->{LG_steps}-1)/2; # index at bottom of middle step
                    my $i_dy2_top    = ($config_ref->{LG_steps}+1)/2; # index at top of middle step
                    my $i_dy2n_bottom = $#sampling_times - ($config_ref->{LG_steps}-1)/2;   # index at bottom of middle step
                    my $i_dy2n_top    = $#sampling_times - ($config_ref->{LG_steps}+1)/2;   # index at top of middle step
                    #printn "XXX $i_dy2_bottom $i_dy2_top $i_dy2n_bottom $i_dy2n_top";
                    my $t_bottom = $sampling_times[$i_dy2_bottom];
                    my $t_top = $sampling_times[$i_dy2_top];
                    my $t_bottom_n = $sampling_times[$i_dy2n_bottom];
                    my $t_top_n = $sampling_times[$i_dy2n_top];

                    my $pos_t1 = $pos_output_vector[$i_dy2_bottom];
                    my $pos_t2 = $pos_output_vector[$i_dy2_top];
                    my $neg_t1 = $neg_output_vector[$i_dy2n_bottom];
                    my $neg_t2 = $neg_output_vector[$i_dy2n_top];

                    if ($config_ref->{round_values_flag}) {
                        ($pos_t1, $pos_t2, $neg_t1, $neg_t2) = map {
                        $self->matlab_round_value(value => $_, AbsTol => $config_ref->{AbsTol}, RelTol => $config_ref->{RelTol})
                        } ($pos_t1, $pos_t2, $neg_t1, $neg_t2);
                    }
                    my $pos_delta = ($pos_t2 - $pos_t1)/($pos_t2 + $pos_t1)*2; # relative measure
                    my $neg_delta = ($neg_t2 - $neg_t1)/($neg_t2 + $neg_t1)*2;

                    my $pos_output_select = ($pos_delta >= $neg_delta) ? 1 : 0;
                    my $output_complex = $pos_output_select ? "TG00001" : "TG00000";
                    my $delta = $stats_ref->{delta} = $pos_output_select ? $pos_delta : $neg_delta;
                    my $delta_score = $stats_ref->{delta_score} = p_hill($stats_ref->{delta}, $config_ref->{delta_threshold}, 1);

                    if ($delta <= 0) {
                        printn "\nWARNING: pos_delta and neg_delta are both zero or negative\n";
                        $stats_ref->{sim_flag} = 0;
                    }
                    printn "pos_output_select = $pos_output_select";
                    #		printn "t_bottom = $t_bottom t_top = $t_top t_bottom_n = $t_bottom_n t_top_n = $t_top_n";
                    printn "pos_delta = $pos_delta, neg_delta = $neg_delta";
                    my @output_vector = $pos_output_select ? @pos_output_vector : @neg_output_vector;

                    #---------------------------------------------------------
                    # PHASE PLOT
                    #---------------------------------------------------------
                    if (defined $config_ref->{plot_phase} && $config_ref->{plot_phase}) {
                        my $output_node = $pos_output_select ? "TG00001" : "TG00000";
                        $self->matlab_plot_phase(
                            figure => 904,
                            X_complex => "LG0000",
                            Y_complex => $output_node,
                            title_prefix => "$TAG $genome_name",
                            axis_ref => [0, $config_ref->{LG_range},
                                0, $config_ref->{TG_init}],
                            filename => "png/$genome_name.phase.png",
                        );
                    }

                    #---------------------------------------------------------
                    # COMPUTE ERROR, STATS
                    #---------------------------------------------------------
                    my $mean_squared_err = $stats_ref->{mean_squared_err} = mean_squared_error(\@output_vector, \@expected_output_vector);
                    my $max_mean_squared_err = $stats_ref->{max_mean_squared_err} = ($config_ref->{TG_init}) ** 2; # since it's the maximum MEAN error^2
                    my $min_mean_squared_err = $stats_ref->{min_mean_squared_err} = 1e-6 * $stats_ref->{max_mean_squared_err};
                    $stats_ref->{mean_squared_err_score} = (log($max_mean_squared_err)-log($min_mean_squared_err + $mean_squared_err))/log($max_mean_squared_err/$min_mean_squared_err);

                    ######################################################################################################
                    # ultrasensitivity measure
                    my @output_vector_slopes = map {$output_vector[$_] - $output_vector[$_-1]} (1..$#output_vector);
                    my $LG_steps = $config_ref->{LG_steps};
                    confess "ERROR: LG_steps must be odd" if (int $LG_steps/2) == ($LG_steps/2);
                    #		my $dy1  = max_numeric(0, $output_vector_slopes[$i_dy2_bottom-1]);
                    my $dy2  = max_numeric(0, $output_vector_slopes[$i_dy2_top-1]);
                    #		my $dy3  = max_numeric(0, $output_vector_slopes[$i_dy2_top]);
                    #		my $dy1n = max_numeric(0, -$output_vector_slopes[$i_dy2n_bottom]);
                    my $dy2n = max_numeric(0, -$output_vector_slopes[$i_dy2n_top]);
                    #		my $dy3n = max_numeric(0, -$output_vector_slopes[$i_dy2n_top-1]);

                    my ($dy1_min, $dy1_max)   = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2_bottom-1],
                        t2 => $sampling_times[$i_dy2_bottom],
                    );
                    my ($dy1n_min, $dy1n_max) = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2n_bottom],
                        t2 => $sampling_times[$i_dy2n_bottom+1],
                    );
                    my ($dy3_min, $dy3_max)   = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2_top],
                        t2 => $sampling_times[$i_dy2_top+1],
                    );
                    my ($dy3n_min, $dy3n_max) = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2n_top-1],
                        t2 => $sampling_times[$i_dy2n_top],
                    );
                    printn "dy1_min/max = ($dy1_min, $dy1_max) dy1n_min/max = ($dy1n_min, $dy1n_max)";
                    printn "dy3_min/max = ($dy3_min, $dy3_max) dy3n_min/max = ($dy3n_min, $dy3n_max)";
                    my $dy1 =  $dy1_max  - $dy1_min;
                    my $dy1n = $dy1n_max - $dy1n_min;
                    my $dy3 =  $dy3_max  - $dy3_min;
                    my $dy3n = $dy3n_max - $dy3n_min;

                    printn "dy1 = $dy1 dy2 = $dy2 dy3 = $dy3";
                    printn "dy1n = $dy1n dy2n = $dy2n dy3n=$dy3n";

                    my $max_dy = $config_ref->{TG_init};

                    my $amplitude = 0;
                    my $ultrasensitivity_dy1 = 0;
                    my $ultrasensitivity_dy3 = 0;
                    if ($delta > 0 && $max_dy != 0 && $dy2 > 0 && $dy2n > 0) {
                        $amplitude = ($dy2 + $dy2n)/2/$max_dy;
                        my $mean_dy1 = ($dy1 + $dy1n)/2;
                        my $mean_dy2 = ($dy2 + $dy2n)/2;
                        my $mean_dy3 = ($dy3 + $dy3n)/2;
                        # adding 0.01*mean_dy2 means you score at most 100
                        #		    $ultrasensitivity_dy1 = $mean_dy2/($mean_dy1 + 0.01*$mean_dy2);
                        #		    $ultrasensitivity_dy3 = $mean_dy2/($mean_dy3 + 0.01*$mean_dy2);
                        $ultrasensitivity_dy1 = $mean_dy2/($mean_dy1);
                        $ultrasensitivity_dy3 = $mean_dy2/($mean_dy3);
                    }
                    $stats_ref->{amplitude} = $amplitude;
                    #		$stats_ref->{amplitude_score} = p_hill($amplitude, $config_ref->{amplitude_threshold}, 1);
                    $stats_ref->{amplitude_score} = $amplitude;
                    $stats_ref->{ultrasensitivity_dy1} = $ultrasensitivity_dy1;
                    $stats_ref->{ultrasensitivity_dy3} = $ultrasensitivity_dy3;
                    my $w_u1 = $config_ref->{w_u1};
                    my $w_u3 = $config_ref->{w_u3};
                    $stats_ref->{ultrasensitivity_score} = (
                        (p_hill($ultrasensitivity_dy1, $config_ref->{ultrasensitivity_threshold}, 1)**$w_u1) *
                        (p_hill($ultrasensitivity_dy3, $config_ref->{ultrasensitivity_threshold}, 1)**$w_u3)
                    )**(1/($w_u1+$w_u3));
                    #		$stats_ref->{ultrasensitivity_score} = (
                    #		    p_hill($ultrasensitivity_dy1, $config_ref->{ultrasensitivity_threshold}, 1));

                    # check if any concentrations are negative
                    foreach my $sample (@output_vector) {
                        if ($sample < 0) {
                            printn "\nWARNING: detected negative concentrations in output vector\n";
                            $stats_ref->{sim_flag} = 0;
                        }
                    }

                    if (($stats_ref->{mean_squared_err_score} < 0) || ($stats_ref->{mean_squared_err_score} > $stats_ref->{max_mean_squared_err})) {
                        #		if (($stats_ref->{mean_squared_err_score} < 0) || ($stats_ref->{mean_squared_err_score} > 1)) {
                        # numerics got messed up, set score to zero
                        printn "\nWARNING: computed mean_squared_error_score out of bounds\n";
                        $stats_ref->{sim_flag} = 0;
                    }
                }
            }
        }  # if $parse_successful
        $stats_ref->{network_connectivity} = $network_connectivity;
        $stats_ref->{network_score} = p_hill($stats_ref->{network_connectivity}, 500, 1);

        #---------------------------------------------------------
        # FINAL SCORE
        #---------------------------------------------------------
        my $final_score = 0;

        if ($parse_successful) {
            my $w_n = $config_ref->{w_n};
            my $w_c = $config_ref->{w_c};
            my $w_s = $config_ref->{w_s};
            my $w_a = $config_ref->{w_a};
            my $w_u = $config_ref->{w_u};

            # is the input connected to the output?
            my $g0  = $stats_ref->{network_connected_flag} ? 1 : 0;
            my $g0n = !$g0 ? 1 : 0;
            # ANC network is OK (i.e. max species not reached) and was therefore simulated?
            # Also, no timeout occurred while simulating to steady-state?
            # No strange numerical problems?
            my $g1  = $stats_ref->{sim_flag} && !$stats_ref->{timeout_flag} && $stats_ref->{sim_flag} ? 1 : 0;
            # is the output amplitude large enough that the output can be trusted artifact-free?
            my $g3  = $g1 && ($stats_ref->{delta_score} > 0.5) ? 1 : 0;
            my $g3n = !$g3 ? 1 : 0;

            my $network_score = $stats_ref->{network_score} || 0;
            my $complexity_score = $stats_ref->{complexity_score} || 0;
            my $steady_state_score = $stats_ref->{steady_state_score} || 0;
            my $amplitude_score = $stats_ref->{amplitude_score} || 0;
            my $ultrasensitivity_score = $stats_ref->{ultrasensitivity_score} || 0;

            # don't optimize network_score once the network is connected
            $final_score =  ($network_score * $g0n + $g0)**$w_n;
            # optimize complexity only if the network is connected
            $final_score *= (1e-3    + $complexity_score * $g0)**$w_c;
            # optimize amplitude if ANC output ok and no timeout during simulation
            $final_score *= (1e-12  + $amplitude_score * $g1)**$w_a;

            $final_score = $final_score**(1/($w_n + $w_c + $w_a));  # re-normalization
        }


        # prevent neg've scores
        if ($final_score < 0) {
            printn "\nWARNING: negative score detected, setting to 0\n";
            $final_score = 0;
        }

        $stats_ref->{score} = ($stats_ref->{score} * ($stats_ref->{score_count} - 1) + $final_score) / $stats_ref->{score_count};
        #    $stats_ref->{score} = $final_score;

        $genome_model_ref->set_score($stats_ref->{score});

        printn $genome_model_ref->sprint_stats();
        printn "final_score = $final_score cumulative_score = $stats_ref->{score} count = $stats_ref->{score_count}";

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


sub run_testcases {
    $verbosity = 1;

    $TAG = "test";
    srand(33433);

    my $config_file = <<END;
#----------------------------------------
# CPU AND CLUSTER SETTINGS
#----------------------------------------
nice = 10
vmem = 2000000

#----------------------------------------
# WORKSPACE AND CUSTOM SCORING MODULE
#----------------------------------------
scoring_class = Ultrasensitive
work_dir = ultrasensitive

#----------------------------------------
# GENOME PARAMS
#----------------------------------------

# Scaling: all concentrations in uM, all 2nd-order rates in uM^-1 s^-1

# Genome class
radius = 3
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
binding_profile_width = 10
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
max_species = 80
max_csite_bound_to_msite_number = 1
default_max_count = 2          # this prevents polymerization (see ANC manual)
default_steric_factor = 1e3    # in micro-mol/L
#export_graphviz = nothing
export_graphviz = network,collapse_states,collapse_complexes

#----------------------------------------
# FACILE/MATLAB SETTINGS
#----------------------------------------
solver = ode23s
#solver = stoch

sampling_interval = 1
SS_timescale = 500.0

# MATLAB odeset params
InitialStep = 1e-8
AbsTol = 1e-9
RelTol = 1e-3
MaxStep = 500.0

#----------------------------------------
# SIMULATION/SCORING PARAMS
#----------------------------------------
plot_input = 1
plot_output = 1
plot_species = 0
plot_phase = 1
plot_min = -1

round_values_flag = 0

steady_state_threshold = 1000   # IC settling time
steady_state_score_threshold = 0.5

delta_threshold = 0.01          # relative measure of amplitude used to filter out integration noise
amplitude_threshold = 0.01      # absolute measure of amplitude
ultrasensitivity_threshold = 5  # ratio of 2nd step over 1st step

w_n = 1.0
w_c = 0.1
#w_s = 0.1
w_a = 1.0
w_u = 1.0
w_u1 = 1.0
w_u3 = 0.1

LG_range = 10          # uM (about 6 molecules in 1e-18L vol ???)
LG_delay = ~
LG_strength = 4.0      # in Hz
LG_ramp_time = 3000
LG_steps = 3

LG_timeout = 20000

stimulus = ss_ramp_equation

hill_n = 40
hill_k = 5

TG_init = 1  # uM
cell_volume = 1e-18             # 1e-18L --> sub-cellular volume

lg_binding_profile = 0101001110
tg_binding_profile = 1010001110

END

    burp_file("test/custom/Ultrasensitive.cfg", $config_file);

    my $scoring_ref = Ultrasensitive->new({
            node_ID => 98,
            config_file => "test/custom/Ultrasensitive.cfg",
            work_dir    => "test/custom",
            matlab_startup_options => "-nodesktop -nosplash",
        });

    printn $scoring_ref->_DUMP();

    my $config_ref = {};
    read_config($config_ref, "test/custom/Ultrasensitive.cfg");

    use GenomeModel;
    my $genome_model_ref = GenomeModel->new({
            name => "Ultrasensitive",
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
            PRE_JUNK => undef,   # undef == ""
            POST_JUNK => "0000",
            genes => [
                {
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 1.0, # uM
                    UNUSED => "0000",
                    domains => [
                        {
                            allosteric_flag => 0,
                            RT_transition_rate => 1.00,
                            TR_transition_rate => 1.00,
                            RT_phi => 1.0,
                            protodomains => [
                                {
                                    type => "bsite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement($lg_binding_profile)->sprint(),
                                    kf_profile => "00000000000000000000",
                                    kb_profile => "11000000110000001000",
                                    kp_profile => "00011111000111110011",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                                {
                                    type => "bsite",
                                    substrate_polarity => 0,
                                    binding_profile => "0101001011",
                                    kf_profile => "11111100111111001110",
                                    kb_profile => "11000000110000001000",
                                    kp_profile => "11111100111111001110",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                            ],
                            UNUSED => "0",
                        },
                    ],
                },
                {
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 1.0, # uM
                    UNUSED => "0000",
                    domains => [
                        {
                            allosteric_flag => 0,
                            RT_transition_rate => 1.00,
                            TR_transition_rate => 1.00,
                            RT_phi => 1.0,
                            protodomains => [
                                {
                                    type => "msite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement("0101001011")->sprint(),
                                    kf_profile => "00000000000000000000",
                                    kb_profile => "11000000110000001000",
                                    kp_profile => "00011111000111110011",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                                {
                                    type => "csite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
                                    kf_profile => "11111100111111001110",
                                    kb_profile => "11000000110000001000",
                                    kp_profile => "11111100111111001110",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                            ],
                            UNUSED => "0",
                        },
                    ],
                },
                {
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 0.1, # uM
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
                                    substrate_polarity => 1,
                                    binding_profile => "0111001011",
                                    kf_profile => "11111100111111001110",
                                    kb_profile => "11000000110000001000",
                                    kp_profile => "11111100111111001110",
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
    store($genome_model_ref, "test/custom/Ultrasensitive.obj");

    $scoring_ref->score_genome($genome_model_ref);
    printn $genome_model_ref->_DUMP();
    sleep 20;
}


# Package BEGIN must return true value
return 1;

