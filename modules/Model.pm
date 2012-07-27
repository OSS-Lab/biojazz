#-#####################################################################################
#- File:     Model.pm
#- Synopsys: A model generator comprises a Sequence and the Parser required to
#-           transcribe/translate this Sequence into model attributes.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Model;
use Class::Std::Storable;
use base qw(Named ClassData);
{
    use Carp;
    use Utils;

    use Sequence;
    use Parser;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    Model->set_class_data("PARSER_CLASS", "Parser");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %sequence_ref_of :ATTR(get => 'sequence_ref', set => 'sequence_ref');
    my %parser_ref_of   :ATTR(get => 'parser_ref',   set => 'parser_ref');

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
#    #--------------------------------------------------------------------------------------
#    # Function: XXX
#    # Synopsys: 
#    #--------------------------------------------------------------------------------------
#    sub XXX {
#	my $class = shift;
#    }

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: START
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;
	my $class = ref $self;

	my $sequence_ref = $sequence_ref_of{$obj_ID} = Sequence->new();

	my $parser_class = $class->get_class_data("PARSER_CLASS");
	my $parser_ref = $parser_class->new({
	    %{$arg_ref->{$parser_class} || {}},
	    name => $parser_class,
	});
	$parser_ref_of{$obj_ID} = $parser_ref;
    }

    #--------------------------------------------------------------------------------------
    # Function: parse
    # Synopsys: Parse sequence_ref and any additional sequences provided
    #--------------------------------------------------------------------------------------
    sub parse {
	my $self = shift; my $obj_ID = ident $self;

	my $parser_ref = $self->get_parser_ref();
	my $iref = $parser_ref->parse(sequence_ref => $sequence_ref_of{$obj_ID});
	foreach my $arg_ref (@_) {
	    $parser_ref->parse(dont_clear_flag => 1, @$arg_ref);
	}
	return $iref;
    }

    #--------------------------------------------------------------------------------------
    # Function: translate
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub translate {
	my $self = shift;

	$self->get_parser_ref()->translate();
    }

    #--------------------------------------------------------------------------------------
    # Function: clear_parser_instances
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub clear_parser_instances {
	my $self = shift;

	$self->get_parser_ref()->clear_parser_instances();
    }
}


sub run_testcases {
    my $model_ref = Model->new({
	name => "model1",
    });
    printn $model_ref->_DUMP();
}


# Package BEGIN must return true value
return 1;

