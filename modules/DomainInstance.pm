#-#####################################################################################
#- File:     DomainInstance.pm
#- Synopsys: Instance of Domain
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package DomainInstance;
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

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: get_protodomains
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_protodomains {
	my $self = shift;
	
	return @{$self->get_field_ref(["protodomains"])};
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
	    if (($field_name eq "allosteric_flag") ||
		($field_name eq "UNUSED")
	       ) {
		# no translation necessary
		$new_value = $field_sequence;
		last SWITCH;
	    }
	    if (($field_name eq "RT_transition_rate") ||
		($field_name eq "TR_transition_rate")
	       ) {
		my ($max, $min);
		eval "\$max = \$self->get_parent_ref->get_${field_name}_max";
		eval "\$min = \$self->get_parent_ref->get_${field_name}_min";
		confess "ERROR: internal error -- max not defined" if !defined $max;
		confess "ERROR: internal error -- min not defined" if !defined $min;
		$new_value = loglinear($min, $max, (2**length($field_sequence)) - 1, $field_int_value);

		last SWITCH;
	    }
	    if ($field_name eq "RT_phi") {
		my ($max, $min);
		eval "\$max = \$self->get_parent_ref->get_${field_name}_max";
		eval "\$min = \$self->get_parent_ref->get_${field_name}_min";
		confess "ERROR: internal error -- max not defined" if !defined $max;
		confess "ERROR: internal error -- min not defined" if !defined $min;
		$new_value = linear($min, $max, (2**length($field_sequence)) - 1, $field_int_value);

		last SWITCH;
	    }
	    confess "ERROR: translate_field() -- unknown field $field_name";
	}			# SWITCH

	# special treatment for allosteric flag -- propagate value down
	# to ProtoDomainInstance attributes in the DomainInstance
	if ($field_name eq "allosteric_flag") {
	    map {$_->[0]->set_allosteric_flag($new_value)} @{$self->get_transcription_ref()->{protodomains}};
	}

	return $new_value;
    }


    #--------------------------------------------------------------------------------------
    # Function: export_anc
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub export_anc {
	my $self = shift;

	confess "ERROR: sequence not translated" if (!$self->get_translation_valid_flag());


	my $name = $self->get_name();
	my $allosteric_flag = $self->get_translation_ref->{allosteric_flag};
	my $RT_transition_rate = $self->get_translation_ref->{RT_transition_rate};
	my $TR_transition_rate = $self->get_translation_ref->{TR_transition_rate};
	my $RT_phi = $self->get_translation_ref->{RT_phi};
	
	my $str = (join "\n", map {$_->[0]->export_anc()} @{$self->get_transcription_ref->{protodomains}});

	my $elements = join(",",map {$_->[0]->get_name()} @{$self->get_transcription_ref->{protodomains}});

	$str .= <<END;
AllostericStructure : {
    name => "$name",
    allosteric_flag => $allosteric_flag,
    allosteric_transition_rates => [$RT_transition_rate,$TR_transition_rate],
    Phi => $RT_phi,
    elements => [$elements],
}
END
	
	return $str;
    }

#    #--------------------------------------------------------------------------------------
#    # Function: xxx
#    # Synopsys: 
#    #--------------------------------------------------------------------------------------
#    sub xxx {
#	my $self = shift;
#    }
}


sub run_testcases {

}


# Package BEGIN must return true value
return 1;

