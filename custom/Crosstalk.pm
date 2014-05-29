#
#===============================================================================
#
#         FILE: Crosstalk.pm
#
#  DESCRIPTION: To scoring two or three signaling pathways (typically with
#               two/thress input LG0, LG1, LG2 and two/three output TG0, TG1, TG2
#               those input and output can not mutate/evolve. and try to evolve
#               the network to compute all signals
#
#        FILES: Scoring.pm, Stimulus.pm, MatlabDriver.pm, GenomeModel.pm
#         BUGS: Newly created
#        NOTES: Need to implement different types of signaling dynamics: adaptation
#               , ultrasensitivity, (damped) oscillation?
#               Ideally, if we could associate the dynamics with real biological
#               funcitons, it will be more meaningful.
#       AUTHOR: Song Feng
# ORGANIZATION: LifeWorks, OSSLAB
#      VERSION: 1.0
#      CREATED: 01/15/2014 10:26:20
#     REVISION: ---
#===============================================================================

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

    
#===  FUNCTION  ================================================================
#         NAME: score_genome
#      PURPOSE: scoring the genome to select crosstalk networks
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
    
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
            printn "Ultrasensitive::score_genome re-scoring elite individual" if $verbosity > 1;

            # Should clear stats in the evolution to prevent stats loss when rescore genomes.
            $genome_model_ref->clear_stats();
            $stats_ref->{score} = 0;
        } else {
            printn "Ultrasensitive::score_genome scoring non-elite individual..." if $verbosity > 1;
 
            # Should clear stats in the evolution to prevent stats loss when rescore genomes.
            $genome_model_ref->clear_stats();
            $stats_ref->{score} = 0;
        }

        #---------------------------------------------------------
        # CREATE I/O GENES
        #---------------------------------------------------------
        my $l0g_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
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
                                binding_profile => $config_ref->{l0g_binding_profile},
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
        printn "l0g_sequence=".$l0g_sequence_ref->get_sequence() if $verbosity > 1;
        my $t0g_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
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
                                binding_profile => $config_ref->{t0g_binding_profile},
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
        printn "t0g_sequence=".$t0g_sequence_ref->get_sequence() if $verbosity > 1;

        my $l1g_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
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
                                binding_profile => $config_ref->{l1g_binding_profile},
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
        printn "l1g_sequence=".$l1g_sequence_ref->get_sequence() if $verbosity > 1;
        my $t1g_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
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
                                binding_profile => $config_ref->{t1g_binding_profile},
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
        printn "t1g_sequence=".$t1g_sequence_ref->get_sequence() if $verbosity > 1;

        #---------------------------------------------------------
        # STIMULUS/SAMPLING EQUATIONS
        #---------------------------------------------------------
        # basically generate two square input signals for each two LG and overlap one of them but seperate the other ones.
        # eventually we could generate such signal to calculate the output comparing input.
        my $stimulus_sub_ref = \&{$config_ref->{stimulus}};
        my $stimulus0_ref = undef;
        my $stimulus1_ref = undef;
        if ($config_ref->{stimulus} eq "clamping_equation") {
            $stimulus0_ref = &$stimulus_sub_ref(
                NODE => "L0G0000",
                DELAY => $config_ref->{L0G_delay},
                PERIOD => $config_ref->{L0G_period},
                STRENGTH => $config_ref->{L0G_strength},
                CONCENTRATION => $config_ref->{L0G_range},
                DUTY => $config_ref->{L0G_duty},
            );
            $stimulus1_ref = &$stimulus_sub_ref(
                NODE => "L1G0000",
                DELAY => $config_ref->{L1G_delay},
                PERIOD => $config_ref->{L1G_period},
                STRENGTH => $config_ref->{L1G_strength},
                CONCENTRATION => $config_ref->{L1G_range},
                DUTY => $config_ref->{L1G_duty},
            );
        } else {
            confess "ERROR: unknown stimulus subroutine";
        }
 
        my ($l0g_source_eqn, $l0g_sink_eqn) = @{$stimulus0_ref};
        my ($l1g_source_eqn, $l1g_sink_eqn) = @{$stimulus1_ref};
        if ($verbosity > 1) {
            printn "Stimulus:";
            printn $l0g_source_eqn;
            printn $l0g_sink_eqn;
            printn $l1g_source_eqn;
            printn $l1g_sink_eqn;
        }

        #---------------------------------------------------------
        # PARSE/TRANSLATE GENOME AND I/O GENES
        #---------------------------------------------------------
        my $genome_iref = $genome_model_ref->parse(
            [
                sequence_ref => $l0g_sequence_ref,
                prefix => "L0",
            ],
            [
                sequence_ref => $t0g_sequence_ref,
                prefix => "T0",
            ],
            [
                sequence_ref => $l1g_sequence_ref,
                prefix => "L1",
            ],
            [
                sequence_ref => $t1g_sequence_ref,
                prefix => "T1",
            ],
        );
        my $parse_successful = $stats_ref->{parse_successful} = $genome_model_ref->check();

        my $history = $genome_model_ref->sprint_history(10);
        printn $history if $verbosity > 1 || $config_ref->{sprint_history};


        ##########################################################
        # Evaluate the connectivity, make sure all output T0 and T1 
        # have at least one kinase and one phophatase
        ##########################################################

        my $network_connectivity = $stats_ref->{network_connectivity} = 0;
        if ($parse_successful) {
            my $transcript = $genome_iref->sprint(colour_flag => 0);
            printn $transcript if $verbosity > 2 || $config_ref->{sprint_transcript};
            burp_file("$matlab_work/$genome_name.tsc", "$history\n$transcript") if $config_ref->{save_transcript};
            $genome_model_ref->translate();

            #---------------------------------------------------------------
            # BUILD/PRUNE Network
            #---------------------------------------------------------------
            my $genome_ref = $genome_model_ref->get_parser_ref();
            $genome_ref->build_network();
            
            # REPORT PROTODOMAIN CONNECTIVITY
            printn "Protodomains: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{protodomains}}) if $verbosity > 1;
            printn $genome_ref->get_adjacency_matrix_ref()->{protodomains}->[0]->sprint_matrix() if $verbosity > 2;
            printn $genome_ref->get_connectivity_matrix_ref()->{protodomains}->sprint_matrix() if $verbosity > 2;
            # REPORT GENE CONNECTIVITY
            printn "Genes: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{genes}}) if $verbosity > 1;
            printn $genome_ref->get_adjacency_matrix_ref()->{genes}->[0]->sprint_matrix() if $verbosity > 2;
            printn $genome_ref->get_connectivity_matrix_ref()->{genes}->sprint_matrix() if $verbosity > 2;

            # PRUNE GENES
            $stats_ref->{num_pruned_genes} = scalar $genome_ref->prune_isolated_genes();

            my $gene_ref = $genome_model_ref->get_gene_parser_ref();
            my $l0g_gene_ref = $gene_ref->lookup_object_instance_by_name("L0G0000");
            my $t0g_gene_ref = $gene_ref->lookup_object_instance_by_name("T0G0000");
            my $l1g_gene_ref = $gene_ref->lookup_object_instance_by_name("L1G0000");
            my $t1g_gene_ref = $gene_ref->lookup_object_instance_by_name("T1G0000");
            my $protodomain_ref = $genome_model_ref->get_protodomain_parser_ref();
            my ($l0g_protodomain_ref)= $protodomain_ref->grep_instances_by_name("L0PD");
            my ($t0g_protodomain_ref)= $protodomain_ref->grep_instances_by_name("T0PD");
            my ($l1g_protodomain_ref)= $protodomain_ref->grep_instances_by_name("L1PD");
            my ($t1g_protodomain_ref)= $protodomain_ref->grep_instances_by_name("T1PD");

            #---------------------------------------------------------
            # SCORING: 2 pts -- LG/TG connected to anything
            #---------------------------------------------------------
            if ($l0g_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "LG is connected" if $verbosity > 1;
            }
            if ($t0g_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "TG is connected" if $verbosity > 1;
            }
            if ($l1g_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "LG is connected" if $verbosity > 1;
            }
            if ($t1g_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "TG is connected" if $verbosity > 1;
            }

            #---------------------------------------------------------
            # SCORING: 90 + 400 pts -- LG/TG subnets
            #---------------------------------------------------------
            my (@l0g_subnet, @t0g_subnet, @l1g_subnet, @t1g_subnet);   # gene subnets
            my (@t0g_adjacent_kinases, @t0g_adjacent_phosphatases, @l0g_adjacent_protodomains);
            my (@t1g_adjacent_kinases, @t1g_adjacent_phosphatases, @l1g_adjacent_protodomains);
            if ($network_connectivity == 4) { # LG/TF connected
                my (@l0g_pd_subnet, @t0g0_pd_subnet, @t0g1_pd_subnet);   # protodomain subnets
                my (@l1g_pd_subnet, @t1g0_pd_subnet, @t1g1_pd_subnet);   # protodomain subnets
                @l0g_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $l0g_protodomain_ref);
                @t0g0_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $t0g_protodomain_ref, state => 0);
                @t0g1_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $t0g_protodomain_ref, state => 1);

                @l1g_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $l1g_protodomain_ref);
                @t1g0_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $t1g_protodomain_ref, state => 0);
                @t1g1_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $t1g_protodomain_ref, state => 1);

                printn "L0G protodomain connects to ".join ",", (map {$_->[2]} @l0g_pd_subnet) if $verbosity > 1;
                printn "T0G/0 protodomain connects to ".join ",", (map {$_->[2]} @t0g0_pd_subnet) if $verbosity > 1;
                printn "T0G/1 protodomain connects to ".join ",", (map {$_->[2]} @t0g1_pd_subnet) if $verbosity > 1;

                printn "L1G protodomain connects to ".join ",", (map {$_->[2]} @l1g_pd_subnet) if $verbosity > 1;
                printn "T1G/0 protodomain connects to ".join ",", (map {$_->[2]} @t1g0_pd_subnet) if $verbosity > 1;
                printn "T1G/1 protodomain connects to ".join ",", (map {$_->[2]} @t1g1_pd_subnet) if $verbosity > 1;

                # max 2 * 90 points for subnet size
                my $l0g_pd_subnet_size = (@l0g_pd_subnet > 30) ? 30 : @l0g_pd_subnet;
                my $t0g0_pd_subnet_size = (@t0g0_pd_subnet > 30) ? 30 : @t0g0_pd_subnet;
                my $t0g1_pd_subnet_size = (@t0g1_pd_subnet > 30) ? 30 : @t0g1_pd_subnet;
                $network_connectivity += ($l0g_pd_subnet_size + $t0g0_pd_subnet_size + $t0g1_pd_subnet_size);

                my $l1g_pd_subnet_size = (@l1g_pd_subnet > 30) ? 30 : @l1g_pd_subnet;
                my $t1g0_pd_subnet_size = (@t1g0_pd_subnet > 30) ? 30 : @t1g0_pd_subnet;
                my $t1g1_pd_subnet_size = (@t1g1_pd_subnet > 30) ? 30 : @t1g1_pd_subnet;
                $network_connectivity += ($l1g_pd_subnet_size + $t1g0_pd_subnet_size + $t1g1_pd_subnet_size);

                @t0g_adjacent_kinases = $genome_ref->find_adjacent_csites($t0g_protodomain_ref, 0);
                @t0g_adjacent_phosphatases = $genome_ref->find_adjacent_csites($t0g_protodomain_ref, 1);
                $stats_ref->{num_adjacent_kinases_t0g} = scalar(@t0g_adjacent_kinases);
                $stats_ref->{num_adjacent_phosphatases_t0g} = scalar(@t0g_adjacent_phosphatases);
                printn "Found ".@t0g_adjacent_kinases." adjacent kinases for T0G";
                printn "Found ".@t0g_adjacent_phosphatases." adjacent phosphatases for T0G";

                @t1g_adjacent_kinases = $genome_ref->find_adjacent_csites($t1g_protodomain_ref, 0);
                @t1g_adjacent_phosphatases = $genome_ref->find_adjacent_csites($t1g_protodomain_ref, 1);
                $stats_ref->{num_adjacent_kinases_t1g} = scalar(@t1g_adjacent_kinases);
                $stats_ref->{num_adjacent_phosphatases_t1g} = scalar(@t1g_adjacent_phosphatases);
                printn "Found ".@t1g_adjacent_kinases." adjacent kinases for T1G";
                printn "Found ".@t1g_adjacent_phosphatases." adjacent phosphatases for T1G";

                @l0g_adjacent_protodomains = union(
                    [map {$_->[0]} $genome_ref->get_adjacent(key => "protodomains", ref => $l0g_protodomain_ref)],
                );
                @l0g_adjacent_protodomains = simple_difference(
                    \@l0g_adjacent_protodomains,
                    [$l0g_protodomain_ref]
                );
                $stats_ref->{num_receptive_protodomains_l0g} = scalar (@l0g_adjacent_protodomains);

                @l1g_adjacent_protodomains = union(
                    [map {$_->[0]} $genome_ref->get_adjacent(key => "protodomains", ref => $l1g_protodomain_ref)],
                );
                @l1g_adjacent_protodomains = simple_difference(
                    \@l1g_adjacent_protodomains,
                    [$l1g_protodomain_ref]
                );
                $stats_ref->{num_receptive_protodomains_l1g} = scalar (@l1g_adjacent_protodomains);

                # now use the gene subnet to determine the connectivity
                @l0g_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $l0g_gene_ref);
                @t0g_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $t0g_gene_ref);
                printn "L0G protein connects to ".join ",", (map {$_->get_name} @l0g_subnet) if $verbosity > 1;
                printn "T0G protein connects to ".join ",", (map {$_->get_name} @t0g_subnet) if $verbosity > 1;

                @l1g_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $l1g_gene_ref);
                @t1g_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $t1g_gene_ref);
                printn "L1G protein connects to ".join ",", (map {$_->get_name} @l1g_subnet) if $verbosity > 1;
                printn "T1G protein connects to ".join ",", (map {$_->get_name} @t1g_subnet) if $verbosity > 1;

                #########################################################################################
                # score -- LG/TG connected to each other
                if (grep /L0G/, (map {$_->get_name()} @t0g_subnet)) {
                    printn "T0G0000 fans out to L0G0000" if $verbosity > 1;
                    $network_connectivity += 100;
                }
                if (grep /T0G/, (map {$_->get_name()} @l0g_subnet)) {
                    printn "L0G0000 fans out to T0G0000" if $verbosity > 1;
                    $network_connectivity += 100;
                }
                if (scalar(@t0g_adjacent_kinases) > 0) {
                    $network_connectivity += 100;
                }
                if (scalar(@t0g_adjacent_phosphatases) > 0) {
                    $network_connectivity += 100;
                }

                if (grep /L1G/, (map {$_->get_name()} @t1g_subnet)) {
                    printn "T1G0000 fans out to L1G0000" if $verbosity > 1;
                    $network_connectivity += 100;
                }
                if (grep /T1G/, (map {$_->get_name()} @l1g_subnet)) {
                    printn "L1G0000 fans out to T1G0000" if $verbosity > 1;
                    $network_connectivity += 100;
                }
                if (scalar(@t1g_adjacent_kinases) > 0) {
                    $network_connectivity += 100;
                }
                if (scalar(@t1g_adjacent_phosphatases) > 0) {
                    $network_connectivity += 100;
                }


                $network_connectivity += 100 if $network_connectivity >= 900;  # max LG/TG connectivity score
            }
 
            if ($network_connectivity >= 1000) {
                $stats_ref->{network_connected_flag} = 1;
                # exclude proteins not in LG/TG subnet from export
                my @proteins_not_in_subnet = simple_difference([$genome_model_ref->get_genes()], [union(\@l0g_subnet, \@t0g_subnet, \@l1g_subnet, \@t1g_subnet)]);
                map {$_->set_export_flag(0)} @proteins_not_in_subnet;
                $stats_ref->{num_protein_out_subnet} = scalar @proteins_not_in_subnet;

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
                    matlab_solver_options => ('matlab_solver_options{InitialStep} = ' . "$config_ref->{InitialStep};\n" .
                        'matlab_solver_options{AbsTol} = ' . "$config_ref->{AbsTol};\n" . 
                        'matlab_solver_options{RelTol} = ' . "$config_ref->{RelTol};\n" . 
                        'matlab_solver_options{MaxStep} = ' . "$config_ref->{MaxStep}"),
                    t_final => $config_ref->{LG_timeout},
                    t_vector =>"[t0:$config_ref->{sampling_interval}:tf]",
                    SS_timescale => $config_ref->{SS_timescale},
                );
                burp_file("$matlab_work/$genome_name.mod", $anc_model);
                system("$ENV{ANC_HOME}/anc.pl --report=species $matlab_work/$genome_name.mod");
                my @facile_model = slurp_file("$matlab_work/$genome_name.eqn");

                $stats_ref->{species_report_flag} = $self->anc_process_species_report("$matlab_work/$genome_name.species.rpt");
                if ($stats_ref->{species_report_flag} == 0) {
                    my @anc_species = $self->anc_get_species();
                    $stats_ref->{num_anc_species} = @anc_species;
                    printn "ANC NUM SPECIES: ".scalar(@anc_species) if $verbosity > 1;
                    printn "ANC SPECIES: @anc_species" if $verbosity > 2;

                    #---------------------------------------------------------
                    # OUTPUT KINASE AND PHOSPHATASE
                    #---------------------------------------------------------
                    my @adjacent_kinase_names = map {$_->get_name()} @tg_adjacent_kinases;
                    my @kinase_gene_names = map {$_->get_upper_ref()->get_upper_ref()->get_name()} @tg_adjacent_kinases;
                    my @adjacent_phosphatase_names = map {$_->get_name()} @tg_adjacent_phosphatases;
                    my @phosphatase_gene_names = map {$_->get_upper_ref()->get_upper_ref()->get_name()} @tg_adjacent_phosphatases;

                    my $K1 = 0.0;
                    if (scalar @adjacent_kinase_names > 0) {
                        for (my $i = 0; $i < @adjacent_kinase_names; $i++) {
                            my $pd_name = $adjacent_kinase_names[$i];
                            my $gene_name = $kinase_gene_names[$i];
                            my $protein_concentration = 0;
                            if ($anc_model =~ /Init : \{\s+structure\s?=>\s?$gene_name,\s+IC\s?=>\s?(\S+),/g) {
                                $protein_concentration = $1 + 0;
                                $stats_ref->{$gene_name} = $protein_concentration;
                            }
                            my @K1s = ();
                            while ($anc_model =~ /CanBindRule : \{\s+name\s?=>\s?\S$pd_name\s?(TPD\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)\S,\n.*\n.*\n.*\s+kf\s?=>\s?(\S+),\s+kb\s?=>\s?(\S+),\s+kp\s?=>\s?(\S+),/g) {
                                my $rule_name = 'K1_'.$pd_name.'_'."$1".'_'."$2"."$3"."$4"."$5";
                                my $rule_rate = ($7 + $8) / $6;
                                $stats_ref->{$rule_name} = $rule_rate;
                                push(@K1s, $rule_rate);
                            }
                            if (scalar @K1s > 0) {
                                $K1 = $K1s[0] / $config_ref->{TG_init};
                                for (my $i = 1; $i < @K1s; $i++) {
                                    $K1 *= ($K1s[$i] / $config_ref->{TG_init});
                                }
                            } else {
                                die "didn't find the rate of phosphorylation rule";
                            }
                            $K1 = $K1**(1/(scalar @K1s));
                        }
                    }
                    $stats_ref->{tg_K1} = $K1;

                    my $K2 = 0.0;
                    if (scalar @adjacent_phosphatase_names > 0) {
                        for (my $i = 0; $i < @adjacent_phosphatase_names; $i++) {
                            my $pd_name = $adjacent_phosphatase_names[$i];
                            my $gene_name = $phosphatase_gene_names[$i];
                            my $protein_concentration = 0;
                            if ($anc_model =~ /Init : \{\s+structure\s?=>\s?$gene_name,\s+IC\s?=>\s?(\S+),/g) {
                                $protein_concentration = $1 + 0;
                                $stats_ref->{$gene_name} = $protein_concentration;
                            }
                            my @K2s = ();
                            while ($anc_model =~ /CanBindRule : \{\s+name\s?=>\s?\S$pd_name\s?(TPD\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)\S,\n.*\n.*\n.*\s+kf\s?=>\s?(\S+),\s+kb\s?=>\s?(\S+),\s+kp\s?=>\s?(\S+),/g) {
                                my $rule_name = 'K2_'.$pd_name.'_'."$1".'_'."$2"."$3"."$4"."$5";
                                my $rule_rate = ($7 + $8) / $6;
                                $stats_ref->{$rule_name} = $rule_rate;
                                push(@K2s, $rule_rate);
                            }
                            if (scalar @K2s > 0) {
                                $K2 = $K2s[0] / $config_ref->{TG_init};
                                for (my $i = 1; $i < @K2s; $i++) {
                                    $K2 *= ($K2s[$i] / $config_ref->{TG_init});
                                }
                            } else {
                                die "didn't find the rate of phosphorylation rule";
                            }
                            $K2 = $K2**(1/(scalar @K2s));
                        }
                    }
                    $stats_ref->{tg_K2} = $K2;


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
                    # Basically compute the number of genes, domains, protodomains, rules
                    # and put those values in account as how complex is the network
                    #---------------------------------------------------------
                    my $num_protodomains = @{[$anc_model =~ /ReactionSite :/g]};
                    my $num_domains = @{[$anc_model =~ /AllostericStructure :/g]};
                    my $num_proteins = @{[$anc_model =~ /\sStructure :/g]};
                    my $num_rules = @{[$anc_model =~ /CanBindRule :/g]};
                    $stats_ref->{num_rules} = $num_rules;
                    printn "ANC model complexity: $num_protodomains + $num_domains + $num_proteins + $num_rules" if $verbosity >= 1;
                    $stats_ref->{complexity} = $num_protodomains + $num_domains + $num_proteins + $num_rules;
                    my $complexity_threshold = defined $config_ref->{complexity_threshold} ? $config_ref->{complexity_threshold} : 250;
                    $stats_ref->{complexity_score} = n_hill($stats_ref->{complexity}, $complexity_threshold, 1);
                    #---------------------------------------------------------
                    # CHECK ANC/FACILE MODEL
                    #---------------------------------------------------------
                    # check that TF_0 and TF_1 are products of at least 1 reaction each
                    my $num_reactions_tg_1 = grep (/(->|<-).* TG00001/, @facile_model);
                    my $num_reactions_tg_0 = grep (/(->|<-).* TG00000/, @facile_model);
                    $stats_ref->{num_reactions_tg_0} = $num_reactions_tg_0;
                    $stats_ref->{num_reactions_tg_1} = $num_reactions_tg_1;
                    $network_connectivity += 200 * ($num_reactions_tg_0 > 1 ? 1 : $num_reactions_tg_0);
                    $network_connectivity += 200 * ($num_reactions_tg_1 > 1 ? 1 : $num_reactions_tg_1);


                    # check that number of species is less than maximum
                    if (!defined $config_ref->{max_species} || $config_ref->{max_species} < 0 || $stats_ref->{num_anc_species} < $config_ref->{max_species}) {
                        $network_connectivity += 100;
                    }
                    #############################################################################
                    #----------------------------------------------------------
                    # Score expression cost
                    # compute and add up number of protodomains times concentration
                    # of each gene.
                    # ---------------------------------------------------------
                    my @genes = $genome_model_ref->get_genes();
                    my $expression_cost = 0;
                    foreach my $gene_instance_ref (@genes) {
                        my $pd_num = 0;
                        my @domains = $gene_instance_ref->get_domains();
                        foreach my $domain_ref (@domains) {
                            $pd_num += scalar $domain_ref->get_protodomains();
                        }
                        $expression_cost += $pd_num * ($gene_instance_ref->get_translation_ref()->{regulated_concentration});
                    }
                    $expression_cost -= $config_ref->{TG_init};
                    my $expression_threshold = defined $config_ref->{expression_threshold} ? $config_ref->{expression_threshold} : 500;
                    $stats_ref->{expression_score} = n_hill($expression_cost, $expression_threshold, 1);
                }
            }





           
        }


        return 1;
    } ## --- end sub score_genome

}



