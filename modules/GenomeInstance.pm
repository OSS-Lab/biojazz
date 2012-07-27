#-#####################################################################################
#- File:     GenomeInstance.pm
#- Synopsys: Instance of Genome
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package GenomeInstance;
use Class::Std::Storable;
use base qw(ParserInstance);
{
    use Carp;

    use Utils;

    use Globals;
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
    # Function: get_genes
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_genes {
	my $self = shift;
	
	return @{$self->get_field_ref(["genes"])};
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
	my $field_int_value = bin2dec($field_sequence) if $field_name !~ /JUNK/;

	my $new_value;
      SWITCH: {
	    if (($field_name eq "PRE_JUNK") ||
		($field_name eq "POST_JUNK")
	       ) {
		# no translation necessary
		$new_value = $field_sequence;
		last SWITCH;
	    }
	    if (0) {  # keep code for future use....
		$new_value = $field_int_value;

		my $scaling_factor;
		eval "\$scaling_factor = \$self->get_parent_ref->get_${field_name}_scaling_factor";
		confess "ERROR: internal error" if !defined $scaling_factor;
		$new_value *= $scaling_factor;

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
	my $self = shift; my $obj_ID = ident $self;
	my $default_max_count = shift;

	$default_max_count = -1 if !defined $default_max_count;

	confess "ERROR: sequence not translated" if (!$self->get_translation_valid_flag());

	my $name = $self->get_name();

	my $str = <<END;
#-------------------------
MODEL:  # $name
#-------------------------
END

	my @genes = map {$_->[0]} @{$self->get_transcription_ref->{genes}};
	@genes = grep {$_->get_export_flag()} @genes;

	# APPEND PROTEINS, DOMAINS, PROTODOMAINS
	$str .= (join "\n", map {$_->export_anc($default_max_count)} @genes);

	# APPEND INIT SECTION
	my @element_names = map {$_->get_name()} @genes;
	my @initial_concentrations = map {$_->get_translation_ref()->{regulated_concentration}} @genes;

	my $INIT = join "\n", map {"Init : {\n".
				   "	structure => $element_names[$_],\n".
				   "	IC => $initial_concentrations[$_],\n".
				   "}"} (0..$#element_names);

	$str .= <<END;

# Initial concentrations
$INIT

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

