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




        return 1;
    } ## --- end sub score_genome

}



