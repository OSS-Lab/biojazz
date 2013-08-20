#-#####################################################################################
#- File:     GeneInstance.pm
#- Synopsys: Instance of Gene
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package GeneInstance;
use Class::Std::Storable;
use base qw(ParserInstance);
{
    use Carp;

    use Utils;
    use BitString;

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # FUNCTIONS
    #######################################################################################
    my %export_flag_of :ATTR(get => 'export_flag', set => 'export_flag', default => 1);

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: get_domains
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_domains {
        my $self = shift;

        return @{$self->get_field_ref(["domains"])};
    }

    #--------------------------------------------------------------------------------------
    # Function: translate_field
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub translate_field {
        my $self = shift;
        my $field_name = shift;

        confess "ERROR: unexpected repeat sequence" if ref $field_name;

        my $transcription_ref = $self->get_transcription_ref();

        my $field_ref = $transcription_ref->{$field_name};
        my $field_sequence = $field_ref->[0];
        my $field_int_value = bin2dec($field_sequence);

        my $new_value;
        SWITCH: {
            if (($field_name eq "START_CODE") ||
                ($field_name eq "STOP_CODE") ||
                ($field_name eq "UNUSED")
            ) {
                # no translation necessary
                $new_value = $field_sequence;
                last SWITCH;
            }
            if (($field_name eq "regulated_concentration")
            ) {
                my ($max, $min);
                eval "\$max = \$self->get_parent_ref->get_${field_name}_max";
                eval "\$min = \$self->get_parent_ref->get_${field_name}_min";
                confess "ERROR: internal error -- max not defined" if !defined $max;
                confess "ERROR: internal error -- min not defined" if !defined $min;
                $new_value = loglinear($min, $max, (2**length($field_sequence)) - 1, $field_int_value);

                last SWITCH;
            }
            confess "ERROR: translate_field() -- unknown field $field_name";
        }			# SWITCH

        return $new_value;
    }

    #--------------------------------------------------------------------------------------
    # Function: export_anc
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub export_anc {
        my $self = shift;
        my $max_count = shift;

        $max_count = -1 if !defined $max_count;

        confess "ERROR: sequence not translated" if (!$self->get_translation_valid_flag());

        my $name = $self->get_name();
        my $regulated_concentration = $self->get_translation_ref->{regulated_concentration};

        my $str .= (join "\n", map {$_->[0]->export_anc()} @{$self->get_transcription_ref->{domains}});

        my $elements = join(",",map {$_->[0]->get_name()} @{$self->get_transcription_ref->{domains}});

        $str .= <<END;
Structure : {  # IC of $name = $regulated_concentration
    name => "$name",
    elements => [$elements],
    max_count => $max_count,
}

END
        return $str;
    }

}


sub run_testcases {

}


# Package BEGIN must return true value
return 1;

