#===============================================================================
#
#         FILE: Multistability.pm
#
#  DESCRIPTION: This is the module to scoring the multistable response accroding
#               to the ramping up with 5 steps from input signal
#
#        FILES: Dependent on Scoring.pm and Stimuli.pm
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Song Feng 
# ORGANIZATION: LifeWorks
#      VERSION: 1.0
#      CREATED: 12/04/2013 15:44:11
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use diagnostics;

package Multistability;
use Class::Std::Storable;
use base qw(Scoring);
{
    use Carp;
    use Utils;
    use Global qw($verbosity $TAG);

    use Stimulus;

    #==========================================================================
    # Class Attributes
    #==========================================================================

    #==========================================================================
    # Instance Methods
    #==========================================================================

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


        return 1;
    } ## --- end sub score_genome
}
 

