#-#####################################################################################
#- File:     Parser.pm   (!!!  better name --> Scribe/RNAP/Expressor ???)
#- Synopsys: A Parser comprises the information required to  transcribe a Sequence (of bits)
#            into individual fields, and the mutation rates of these fields.  It also
#            comprises a Set of sub-sequence Parsers.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-
#- The structure_ref contains the parsing information and its entries comprise a
#- field name with an associated regular expression.
#-
#- The mutator_ref give the mutation rates of the fields defined in structure_ref.
#-
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Parser;
use Class::Std::Storable;
use base qw(Instantiable Set);
{
    use Carp;
    use WeakRef;

    use Utils;
    use Globals;

    use Sequence;
    use ParserInstance;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    Parser->set_class_data("ELEMENT_CLASS", "Parser");
    Parser->set_class_data("INSTANCE_CLASS", "ParserInstance");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    #
    # sequence_structure_ref = [
    #                           [attrib, REGEX, repeat]
    #                           [attrib, REGEX, repeat]
    #                           ....
    #                          ]
    # attrib = scalar | functional_sequence_ref
    # REGEX = a regular expression
    # repeat = if attrib is a functional_sequence_ref
    #
    my %structure_ref_of        :ATTR(get => 'structure_ref', set => 'structure_ref');
    my %linker_ref_of        :ATTR(get => 'linker_ref', set => 'linker_ref');
    my %mutation_rate_ref_of    :ATTR(get => 'mutation_rate_ref', set => 'mutation_rate_ref');

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
    # Function: BUILD
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

	confess "ERROR: obsolete init_arg sequence_ref" if exists $arg_ref->{sequence_ref};

	$structure_ref_of{$obj_ID} = (exists $arg_ref->{structure_ref}) ? $arg_ref->{structure_ref} : [];
	$mutation_rate_ref_of{$obj_ID} = (exists $arg_ref->{mutation_rate_ref}) ? $arg_ref->{mutation_rate_ref} : {};
    }

    sub START {
        my ($self, $obj_ID, $arg_ref) = @_;

	# now that structure_ref is initialized by sub-class,
	# can fill in appropriate defaults for mutation_rate_ref
	for (my $i=0; $i < @{$structure_ref_of{$obj_ID}}; $i++) {
	    my $field = $structure_ref_of{$obj_ID}->[$i][0];
	    if (!defined $mutation_rate_ref_of{$obj_ID}{$field}) {
		$mutation_rate_ref_of{$obj_ID}{$field} =  1.0;
	    }
	    if (defined $structure_ref_of{$obj_ID}->[$i][2]) {  # is there a repeat sequence?
		my $linker_field = "${field}_linker";
		if (!defined $mutation_rate_ref_of{$obj_ID}{$linker_field}) {
		    $mutation_rate_ref_of{$obj_ID}{$linker_field} =  1.0;
		}
	    }
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: set_field_mutation_rate
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub set_field_mutation_rate {
	my $self = shift; my $obj_ID = ident $self;
	my $field = shift;
	my $rate = shift;

	$mutation_rate_ref_of{$obj_ID}{$field} = $rate;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_field_mutation_rate
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_field_mutation_rate {
	my $self = shift; my $obj_ID = ident $self;
	my $field = shift;

	return $mutation_rate_ref_of{$obj_ID}{$field};
    }

    #--------------------------------------------------------------------------------------
    # Function: clear_parser_instances
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub clear_parser_instances {
	my $self = shift;

	printn "clear_parser_instances: clearing all instances of ".$self->get_name() if $verbosity >= 2;

	$self->clear_object_instances();
	foreach my $subseq_parser_ref ($self->get_elements()) {
	    $subseq_parser_ref->clear_parser_instances();
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: STORABLE_freeze_pre
    # Synopsys: Hook provided by Class::Std::Storable used to clear parser instances.
    #--------------------------------------------------------------------------------------
     sub STORABLE_freeze_pre: CUMULATIVE {
 	my ($self, $clone_flag) = @_;

	# we clear any parser instances because they contain weakrefs
	# to their parents which are a pain to freeze and thaw
	$self->clear_parser_instances();
    }

    #--------------------------------------------------------------------------------------
    # Function: parse (aka transcribe)
    # Synopsys: This routine will parse the given sequence starting at the given position.
    #           It does so by creating a new instance, and calling the instance's parse()
    #           routine, will recursively parses any substructures.
    #           If there is a successful match, all (sub)instances created will be
    #           registered.
    #--------------------------------------------------------------------------------------
    sub parse {
	my $self = shift; my $obj_ID = ident $self;
	my %args = (
	    sequence_ref => undef,
	    start_pos => 0,
	    dont_clear_flag => 0,
	    prefix => "",
	    @_,
	   );
	check_args(\%args, 4);

	my $sequence_ref = $args{sequence_ref};
	my $start_pos = $args{start_pos};
	my $dont_clear_flag = $args{dont_clear_flag};
	my $prefix = $args{prefix};

	confess "ERROR: sequence_ref is not a Sequence" if (!ref $sequence_ref) || (!$sequence_ref->isa('Sequence'));

	$self->clear_parser_instances() if (!$dont_clear_flag);

	my $instance_ref = $self->new_object_instance({
	    UNREGISTERED => 1,
	    sequence_ref => $sequence_ref,
	    start_pos => $start_pos,
	    prefix => $prefix,
	});

	my @parse_result = $instance_ref->parse($prefix);

	if (@parse_result) {
	    # since parsing was successful, we register the instances (including top, i.e. instance_ref)
	    foreach my $subseq_instance_ref (@{$parse_result[2]}) {
		# create a name based on locus prior to registering
		my $name = $subseq_instance_ref->get_name();
		$name = join("", ($name =~ /[A-Z]/g));  # concat all uppercase chars
		my $locus = sprintf("%04d",$subseq_instance_ref->get_locus());
		$subseq_instance_ref->set_name("__$name$locus");  # double underscore will be removed upon registration
		$subseq_instance_ref->get_parent_ref()->register_instance($subseq_instance_ref);
	    }
	    return $instance_ref;
	} else {
	    # all created instances are destroyed on exit since they were not registered
	    return undef;
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: translate
    # Synopsys: Translate all instances of object.
    #--------------------------------------------------------------------------------------
    sub translate {
	my $self = shift;

	return map {$_->translate()} $self->get_object_instances();
    }

    #--------------------------------------------------------------------------------------
    # Function: untranslate_field
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub untranslate_field {
	my $self = shift;
	my $class = ref $self;

	confess "ERROR: sub-class $class must implement an untranslate_field() routine";
    }

    #--------------------------------------------------------------------------------------
    # Function: untranslate
    # Synopsys: Given field translated values, produce the corresponding transcription.
    #--------------------------------------------------------------------------------------
    sub untranslate {
	my $self = shift; my $obj_ID = ident $self;
	my $translation_ref = shift;

	my $class = ref $self;
	my $structure_ref = $structure_ref_of{$obj_ID};
	my $transcription_ref = {};

	# undo the translation
	for (my $i=0; $i < @$structure_ref; $i++) {
	    my $field = $structure_ref->[$i][0];
	    confess "ERROR: missing field $field when untranslating $class" if !exists $translation_ref->{$field};
	    if (defined $structure_ref->[$i][2]) {  # is there a repeat sequence?
		my @field_values = @{$translation_ref->{$field}};
		for (my $j=0; $j < @field_values; $j++) {
		    if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
			my $subseq_ref = $structure_ref->[$i][1];
			push @{$transcription_ref->{$field}}, $subseq_ref->untranslate($field_values[$j]);
		    } else {
			push @{$transcription_ref->{$field}}, $self->untranslate_field($field, $field_values[$j]);
		    }
		}
	    } else {
		my $field_value = $translation_ref->{$field};
		if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
		    my $subseq_ref = $structure_ref->[$i][1];
		    $transcription_ref->{$field} = $subseq_ref->untranslate($field_value);
		} else {
		    $transcription_ref->{$field} = $self->untranslate_field($field, $field_value);
		}
	    }
	}
	
	return $transcription_ref;
    }

    #--------------------------------------------------------------------------------------
    # Function: untranscribe
    # Synopsys: Reverse the transcription to get sequence, using linkers when necessary
    #--------------------------------------------------------------------------------------
    sub untranscribe {
	my $self = shift; my $obj_ID = ident $self;
	my $transcription_ref = shift;

	my $class = ref $self;
	my $structure_ref = $structure_ref_of{$obj_ID};
	my $sequence = "";

	for (my $i=0; $i < @$structure_ref; $i++) {
	    my $field = $structure_ref->[$i][0];
	    if (defined $structure_ref->[$i][2]) {  # is there a repeat sequence?
		my @field_values = @{$transcription_ref->{$field}};
		my @sequences;
		for (my $j=0; $j < @field_values; $j++) {
		    if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
			my $subseq_ref = $structure_ref->[$i][1];
			push @sequences, $subseq_ref->untranscribe($field_values[$j])->get_sequence();
		    } else {
			push @sequences, $field_values[$j];
		    }
		}
		my $linker_ref = $self->get_linker_ref();
		confess "ERROR: undefined linker_ref" if !defined $linker_ref;
		my $linker = $linker_ref->{$field};
		confess "ERROR: undefined linker for $field field" if !defined $linker;
		$sequence .= join $linker, @sequences;
	    } else {
		my $field_value = $transcription_ref->{$field};
		if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
		    my $subseq_ref = $structure_ref->[$i][1];
		    $sequence .= $subseq_ref->untranscribe($field_value)->get_sequence();
		} else {
		    $sequence .= $field_value;
		}
	    }
	}

	return Sequence->new({sequence => $sequence});
    }

    #--------------------------------------------------------------------------------------
    # Function: create_sequence
    # Synopsys: Given field translated values, produce the corresponding transcription.
    #--------------------------------------------------------------------------------------
    sub create_sequence {
	my $self = shift; my $obj_ID = ident $self;
	my $translation_ref = shift;

	my $untranslated_ref = $self->untranslate($translation_ref);
	my $untranscribed_ref = $self->untranscribe($untranslated_ref);
	return $untranscribed_ref;
    }
}

sub run_testcases {
    use Globals;
    $verbosity = 3;

    my $sequence = "000000110101000111111000011110000011001011010001110100";
    my $seq_ref = Sequence->new({
	#                 [  ][ ][]  []               [  ][ ][]  []  []  []
	sequence => $sequence,
    });

    my $test_parser_ref = Parser->new({
	name => 'Test',
	structure_ref => [
	    ["A", '0110'],
	    ["B", '\G[0-1]{3}'],
	    ["elements_ref", '\G[0-1]{2}', '\G01'],
	   ],
    });

    my $test_parser_I0_ref = $test_parser_ref->parse(sequence_ref => $seq_ref);     # will find first match at offset 5
    printn $test_parser_I0_ref->_DUMP();
    printn $test_parser_I0_ref->sprint(context_flag => 1);
    my $test_parser_I1_ref = $test_parser_ref->parse(sequence_ref => $seq_ref, start_pos => 6, dont_clear_flag => 1);  # will find second match at offset 33
    printn $test_parser_I1_ref->_DUMP();
    printn $test_parser_I1_ref->sprint(context_flag => 1);
    my $test_parser_I2_ref = $test_parser_ref->parse(sequence_ref => $seq_ref, start_pos => 39, dont_clear_flag => 1); # last matching position at offset 39
    printn $test_parser_I2_ref->_DUMP();
    printn $test_parser_I2_ref->sprint(context_flag => 1);
    my $test_parser_I3_ref = $test_parser_ref->parse(sequence_ref => $seq_ref, start_pos => 40, dont_clear_flag => 1); # won't match

    # print the test parser
    printn $test_parser_ref->_DUMP();

    my $protodomain_parser_ref = Parser->new({
	name => 'ProtoDomain',
	structure_ref => [
	    ["A", '..'],
	    ["B", '..', '\G(?!11)..'],  # scoop up bits 2 by 2 until terminating '11'
	   ],
    });

    my $domain_parser_ref = Parser->new({
	name => 'Domain',
	structure_ref => [
	    ["S", '0110'],
	    ["PD", $protodomain_parser_ref, '\G(?!110)...'], # until terminating 110
	    ["T", '\G...'],
	   ],
	mutation_rate_ref => {
	    S => 0.0,
	},
	elements_ref => [$protodomain_parser_ref],
    });

    printn "PARSING DOMAIN....";
    my ($domain_parser_I0_ref) = $domain_parser_ref->parse(sequence_ref => $seq_ref, prefix => "X");

    printn "DUMPING PARSERS";
    printn $domain_parser_ref->_DUMP();
    printn $protodomain_parser_ref->_DUMP();

    printn "DUMPING DOMAIN INSTANCE";
    printn $domain_parser_I0_ref->_DUMP();
    printn $domain_parser_I0_ref->sprint(context_flag => 1);
    printn $domain_parser_I0_ref->sprint(context_flag => 1, hierarchical_flag => 0, colour_flag => 0);

    printn "DUMPING PROTODOMAIN INSTANCES";
    foreach my $fsi_ref ($protodomain_parser_ref->get_object_instances()) {
	printn $fsi_ref->_DUMP();
    }

    printn "FIELD MANIPULATION";
    my $before = $domain_parser_I0_ref->get_sequence_ref->sprint();
    printn "BEFORE: $before";
    $domain_parser_I0_ref->set_field("S", "1001");
    $domain_parser_I0_ref->set_field(["PD", 0], "*" x $domain_parser_I0_ref->get_field_length(["PD", 0]));
    printn "AFTER:  ".$domain_parser_I0_ref->get_sequence_ref->sprint();
    printn "FIELD MUTATION";
    $domain_parser_I0_ref->mutate_field(["PD", 1], 1.0);
    printn "BEFORE: $before";
    printn "AFTER:  ".$domain_parser_I0_ref->get_sequence_ref->sprint();

    printn "GLOBAL MUTATION";
    # reset to orig sequence and clear parsing instances
    $seq_ref->set_sequence($sequence);
    $domain_parser_I0_ref = $domain_parser_ref->parse(sequence_ref => $seq_ref);
    $domain_parser_I0_ref->mutate();
    printn "BEFORE: $sequence";
    printn "AFTER:  ".$seq_ref->sprint();

    printn "STORABLE TEST";
    use Storable;
    my $ice_ref = Storable::freeze($test_parser_ref);
    my $water_ref = Storable::thaw($ice_ref);
    printn $test_parser_ref->_DUMP();
    printn $water_ref->_DUMP();
    map {printn $_->_DUMP()} ($water_ref->get_object_instances());
}


# Package BEGIN must return true value
return 1;

