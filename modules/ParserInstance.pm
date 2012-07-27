#-#####################################################################################
#- File:     ParserInstance.pm   (!!! better name mRNA ???)
#- Synopsys: Instance of Parser, created when Parser::parse() is called with a given
#-           Sequence object as argument, and containing all the extracted information.
#-           Contains routines to get, set and mutate specific fields using the
#-           rates specified in parent object.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package ParserInstance;
use Class::Std::Storable;
use base qw(Instance);
{
    use Carp;
    use Term::ANSIColor; # constants for generating color in text output.
    use WeakRef;

    use Utils;
    use Globals;

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %sequence_ref_of             :ATTR(get => 'sequence_ref', set => 'sequence_ref', init_arg => 'sequence_ref');

    my %start_pos_of                :ATTR(get => 'start_pos', set => 'start_pos', init_arg => 'start_pos');
    my %locus_of                    :ATTR(get => 'locus', set => 'locus');
    my %length_of                   :ATTR(get => 'length', set => 'length');
    my %transcription_ref_of        :ATTR(get => 'transcription_ref', set => 'transcription_ref');
    my %transcription_valid_flag_of :ATTR(get => 'transcription_valid_flag', set => 'transcription_valid_flag');

    my %translation_ref_of          :ATTR(get => 'translation_ref', set => 'translation_ref');
    my %translation_valid_flag_of   :ATTR(get => 'translation_valid_flag', set => 'translation_valid_flag');

    # keep a reference to upper level parser instance
    my %upper_ref_of                :ATTR(get => 'upper_ref', set => 'upper_ref', init_arg => 'upper_ref', default => "");

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    sub START {
        my ($self, $obj_ID, $arg_ref) = @_;

	# check sequence_ref
	my $sequence_ref = $sequence_ref_of{$obj_ID};
	confess "ERROR: sequence_ref must be a Sequence object" if !$sequence_ref->isa('Sequence');

	# check that upper parser points to the same Sequence
	if (ref $upper_ref_of{obj_ID}) {
	    my $upper_sequence_ID = ident($upper_ref_of{obj_ID}->get_sequence_ref());
	    confess "ERROR: internal error -- upper parser doesn't point to same Sequence" if $upper_sequence_ID != (ident $sequence_ref);
	}

	# weaken sequence_ref
	weaken($sequence_ref_of{$obj_ID});

	# weaken upper_ref to break reference loop
	if ($upper_ref_of{$obj_ID}) {
	    weaken($upper_ref_of{$obj_ID});
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: parse
    # Synopsys: Recursively parse the sequence, extracting fields and creating instances
    #           for sub-sequences as necessary.
    #--------------------------------------------------------------------------------------
    sub parse {
	my $self = shift;
	my $obj_ID = ident $self;
	my $prefix = shift || "";

	my $structure_ref = $self->get_structure_ref();

	my $start_pos = $start_pos_of{$obj_ID};
	my $sequence_ref = $self->get_sequence_ref();
	my $sequence = $sequence_ref->get_sequence();

	# init/clear parse results
	my $transcription_ref = $transcription_ref_of{$obj_ID} = {};
	delete $locus_of{$obj_ID};
	delete $length_of{$obj_ID};
	delete $transcription_valid_flag_of{$obj_ID};

	my $pos = pos($sequence) = $start_pos;

	my @subseq_instances = ($self);

	for (my $i=0; $i < @$structure_ref; $i++) {
	    my $field = $structure_ref->[$i][0];
	    confess "ERROR: you cannot repeat fields in structure_ref" if exists $transcription_ref->{$field};
	    my $field_RE = $structure_ref->[$i][1];
	    my $repeat_RE = defined $structure_ref->[$i][2] ? $structure_ref->[$i][2] : undef;
	    my $subseq_flag = ref $structure_ref->[$i][1] ? 1 : 0;
	    my $subseq_ref = $subseq_flag ? $structure_ref->[$i][1] : undef;
	    confess "ERROR: reference must be a Parser object" if $subseq_flag && !$subseq_ref->isa('Parser');

	    my $RE = $field_RE;
	    my $success_flag;
	    while (1) {
		if ($subseq_flag && ($RE eq $field_RE)) {
		    printn "Matching subfield $field" if $verbosity >= 2;
		    my $subseq_instance_ref = $subseq_ref->new_object_instance({
			prefix => $prefix,
			UNREGISTERED => 1,
			sequence_ref => $sequence_ref,
			start_pos => $pos,
			upper_ref => $self,
		    });
		    $success_flag = my @subseq_parse_result = $subseq_instance_ref->parse($prefix);
		    last if !$success_flag;
		    printn "Matched subfield $field" if $verbosity >= 2;

		    push @subseq_instances, @{$subseq_parse_result[2]};

		    $pos = pos($sequence) = $subseq_parse_result[0];
		    $locus_of{$obj_ID} = $subseq_instance_ref->get_locus() if !defined $locus_of{$obj_ID};

		    if (!defined $repeat_RE) {
			$transcription_ref->{$field} = [$subseq_instance_ref,
							$subseq_instance_ref->get_locus(),
							$subseq_instance_ref->get_length()];  # subseq ref, locus, length
			last;
		    } else {
			if ($RE eq $field_RE) {
			    push @{$transcription_ref->{$field}}, [$subseq_instance_ref,
								   $subseq_instance_ref->get_locus(),
								   $subseq_instance_ref->get_length()];  # subseq ref, locus, length
			    $RE = $repeat_RE;
			} else {
			    my $linker_field_name = "${field}_linker";
			    push @{$transcription_ref->{$linker_field_name}}, [$1, $-[0], $pos - $-[0]];
			    $RE = $field_RE;
			}
		    }
		} else {
		    pos($sequence) = $pos;
		    $success_flag = $sequence =~ /($RE)/g;
		    last if !$success_flag;

		    $pos = pos $sequence;
		    $locus_of{$obj_ID} = $-[0] if !defined $locus_of{$obj_ID};

		    if (!defined $repeat_RE) {
			printn "Matched field $field=$1 at $-[0] using RE $structure_ref->[$i][1]" if $verbosity >= 2;
			$transcription_ref->{$field} = [$1, $-[0], $pos - $-[0]];  # sequence, locus, length
			last;
		    } else {
			if ($RE eq $field_RE) {
			    printn "Matched field $field=$1 at $-[0] using RE $structure_ref->[$i][1]" if $verbosity >= 2;
			    push @{$transcription_ref->{$field}}, [$1, $-[0], $pos - $-[0]];
			    $RE = $repeat_RE;
			} else {
			    my $linker_field_name = "${field}_linker";
			    push @{$transcription_ref->{$linker_field_name}}, [$1, $-[0], $pos - $-[0]];
			    $RE = $field_RE;
			}
		    }
		}
	    }
	    if (!$success_flag && ((!defined $repeat_RE) || ($RE eq $field_RE))) {
		printn "NOTE: couldn't match regexp for field $field using RE $field_RE" if $verbosity >= 2;
		# init/clear parse results
		delete $transcription_ref_of{$obj_ID};
		delete $locus_of{$obj_ID};
		delete $length_of{$obj_ID};
		$transcription_valid_flag_of{$obj_ID} = 0;
		return ();
	    }
	}

	my $length = $length_of{$obj_ID} = $pos - $locus_of{$obj_ID};

	$transcription_valid_flag_of{$obj_ID} = 1;
	return ($pos, $length, \@subseq_instances);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_sequence
    # Synopsys: Get the sequence of this parser instance.
    #--------------------------------------------------------------------------------------
    sub get_sequence {
	my $self = shift;
	my $obj_ID = ident $self;

	return $self->get_sequence_ref()->get_subseq($locus_of{$obj_ID}, $length_of{$obj_ID});
    }

    #--------------------------------------------------------------------------------------
    # Function: get_stop_locus
    # Synopsys: Get location of where of ParserInstance transcription stops.
    #--------------------------------------------------------------------------------------
    sub get_stop_locus {
	my $self = shift;
	my $obj_ID = ident $self;

	return $locus_of{$obj_ID} + $length_of{$obj_ID} - 1;
    }

    #--------------------------------------------------------------------------------------
    # Function: set_field
    # Synopsys: Using locus and length info from transcription, set corresponding sub-sequence
    #           to given value. Note that the actual sequence is changed, not the contents of
    #           transcription_ref, so sprint() will not show the change until the sequence
    #           is re-parsed. For a repeated field, supply index also e.g.
    #                $self->set_field(fA, 10011);      # non-repeated field
    #                $self->set_field([fR,3], 10011);  # repeated field
    #--------------------------------------------------------------------------------------
    sub set_field {
	my $self = shift;
	my $obj_ID = ident $self;
	my $field_ref = shift;
	my $value = shift;

	$transcription_valid_flag_of{$obj_ID} = 0;   # mark as dirty

	my $repeated_flag = ref $field_ref eq "ARRAY" ? 1 : 0;
	my $field = $repeated_flag ? $field_ref->[0] : $field_ref;
	my $sequence_ref = $self->get_sequence_ref();
	my $transcription_ref = $transcription_ref_of{$obj_ID};

	confess "ERROR: no such field" if !defined $transcription_ref->{$field};

	if ($repeated_flag) {  # repeated field?
	    confess "ERROR: this is not a repeated field" if ref $transcription_ref->{$field}[0] ne "ARRAY";
	    my $index = $field_ref->[1];
	    confess "ERROR: index out of range" if !defined $transcription_ref->{$field}[$index];
	    my $locus = $transcription_ref->{$field}[$index][1];
	    my $length = $transcription_ref->{$field}[$index][2];
	    return $sequence_ref->set_subseq($value, $locus, $length);
	} else {
	    confess "ERROR: this is a repeated field" if ref $transcription_ref->{$field}[0] eq "ARRAY";
	    my $field = $field_ref;
	    my $locus = $transcription_ref->{$field}[1];
	    my $length = $transcription_ref->{$field}[2];
	    return $sequence_ref->set_subseq($value, $locus, $length);
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: get_field_info
    # Synopsys: Using locus and length info from transcription, extract a field reference
    #           (if applicable), sequence (from actual sequence not transcription), locus,
    #           and length.  For a repeated field, supply index also e.g.
    #                $self->get_field_info(fA);      # non-repeated field
    #                $self->get_field_info([fR,3]);  # repeated field, specific element
    #                $self->get_field_info([fR]);    # repeated field, all elements
    #
    # Returns (subseq_ref, subseq, locus, length) or corresponding lists if all elements.
    #--------------------------------------------------------------------------------------
    sub get_field_info {
	my $self = shift;
	my $obj_ID = ident $self;
	my $arg_ref = shift;

	my $repeated_flag = ref $arg_ref eq "ARRAY" ? 1 : 0;
	my $field = $repeated_flag ? $arg_ref->[0] : $arg_ref;
	my $sequence_ref = $self->get_sequence_ref();
	my $transcription_ref = $transcription_ref_of{$obj_ID};

	confess "ERROR: no such field" if !defined $transcription_ref->{$field};

	if ($repeated_flag) {  # repeated field?
	    confess "ERROR: field $field is not a repeated field" if ref $transcription_ref->{$field}[0] ne "ARRAY";
	    my $index = $arg_ref->[1];
	    if (defined $index) {
		confess "ERROR: index $index out of range in field $field" if !defined $transcription_ref->{$field}[$index];
		my $field_ref = $transcription_ref->{$field}[$index];
		my $subseq_ref = ref $field_ref->[0] ? $field_ref->[0] : undef;
		my $locus = $field_ref->[1];
		my $length = $field_ref->[2];
		my $sequence = $sequence_ref->get_subseq($locus, $length);
		return ($subseq_ref, $sequence, $locus, $length);
	    } else {
		my @field_refs = @{$transcription_ref->{$field}};  # array of field entries
		my @subseq_refs = map {ref $_->[0] ? $_->[0] : undef} @field_refs;
		my @loci = map {$_->[1]} @field_refs;
		my @lengths = map {$_->[2]} @field_refs;
		my @sequences = map {$sequence_ref->get_subseq($loci[$_], $lengths[$_])} (0..$#field_refs);
		return (\@subseq_refs, \@sequences, \@loci, \@lengths);
	    }		
	} else {
	    confess "ERROR: this is a repeated field" if ref $transcription_ref->{$field}[0] eq "ARRAY";
	    my $field_ref = $transcription_ref->{$field};
	    my $subseq_ref = ref $field_ref->[0] ? $field_ref->[0] : undef;
	    my $locus = $field_ref->[1];
	    my $length = $field_ref->[2];
	    my $sequence = $sequence_ref->get_subseq($locus, $length);
	    return ($subseq_ref, $sequence, $locus, $length);
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: get_field_ref
    # Synopsys: Using locus and length info from transcription, get a field sub-sequence
    #           reference (if exists).
    #           For a repeated field, supply index also e.g.
    #                $self->get_field_ref(fA);      # non-repeated field
    #                $self->get_field_ref([fR,3]);  # repeated field
    #--------------------------------------------------------------------------------------
    sub get_field_ref {
	my $self = shift;
	my @info = $self->get_field_info(@_);
	return $info[0];
    }

    #--------------------------------------------------------------------------------------
    # Function: get_field_sequence
    # Synopsys: Using locus and length info from transcription, get a field sequence.
    #           For a repeated field, supply index also e.g.
    #                $self->get_field_sequence(fA);      # non-repeated field
    #                $self->get_field_sequence([fR,3]);  # repeated field
    #--------------------------------------------------------------------------------------
    sub get_field_sequence {
	my $self = shift;
	my @info = $self->get_field_info(@_);
	return $info[1];
    }

    #--------------------------------------------------------------------------------------
    # Function: get_field_locus
    # Synopsys: Using locus and length info from transcription, get a field value.
    #           For a repeated field, supply index also e.g.
    #                $self->get_field_size(fA);      # non-repeated field
    #                $self->get_field_size([fR,3]);  # repeated field
    #--------------------------------------------------------------------------------------
    sub get_field_locus {
	my $self = shift;
	my @info = $self->get_field_info(@_);
	return $info[2];
    }
    #--------------------------------------------------------------------------------------
    # Function: get_field_length
    # Synopsys: Using locus and length info from transcription, get a field value.
    #           For a repeated field, supply index also e.g.
    #                $self->get_field_length(fA);      # non-repeated field
    #                $self->get_field_length([fR,3]);  # repeated field
    #--------------------------------------------------------------------------------------
    sub get_field_length {
	my $self = shift;
	my @info = $self->get_field_info(@_);
	return $info[3];
    }

    #--------------------------------------------------------------------------------------
    # Function: get_field_stop
    # Synopsys: Get location of where of ParserInstance field stops.
    #--------------------------------------------------------------------------------------
    sub get_field_stop {
	my $self = shift;

	my @info = $self->get_field_info(@_);

	return $info[2] + $info[3] - 1;
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate_field
    # Synopsys: Using locus and length info from transcription, mutate corresponding sub-sequence.
    #           For a repeated field, supply index also e.g.
    #                $self->mutate_field(fA);      # non-repeated field
    #                $self->mutate_field([fR,3]);  # repeated field
    #--------------------------------------------------------------------------------------
    sub mutate_field {
	my $self = shift;
	my $obj_ID = ident $self;
	my $arg_ref = shift;
	my $probability = shift;

	my $repeated_flag = ref $arg_ref eq "ARRAY" ? 1 : 0;
	my $field = $repeated_flag ? $arg_ref->[0] : $arg_ref;
	my $sequence_ref = $self->get_sequence_ref();
	my $transcription_ref = $transcription_ref_of{$obj_ID};

	confess "ERROR: no such field" if !defined $transcription_ref->{$field};

	my $num_bits;
	if ($repeated_flag) {  # repeated field?
	    confess "ERROR: this is not a repeated field" if ref $transcription_ref->{$field}[0] ne "ARRAY";
	    my $index = $arg_ref->[1];
	    confess "ERROR: index out of range" if !defined $transcription_ref->{$field}[$index];
	    my $locus = $transcription_ref->{$field}[$index][1];
	    my $length = $transcription_ref->{$field}[$index][2];
	    $num_bits = $sequence_ref->mutate_subseq($probability, $locus, $length);
	} else {
	    confess "ERROR: this is a repeated field" if ref $transcription_ref->{$field}[0] eq "ARRAY";
	    my $locus = $transcription_ref->{$field}[1];
	    my $length = $transcription_ref->{$field}[2];
	    $num_bits = $sequence_ref->mutate_subseq($probability, $locus, $length);
	}
	printn "mutated_field: mutated $num_bits bits in field $field" if $verbosity >= 2;
	return $num_bits;
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate
    # Synopsys: Recursively mutate a sequence field-by-field, using the given mutation
    #           rates for each field.  WILL MUTATE REPEAT SEQUENCES ALSO!!!
    #           Also takes a global scaling factor.
    #--------------------------------------------------------------------------------------
    sub mutate {
	my $self = shift;
	my $obj_ID = ident $self;
	my $mutation_rate_scaling_factor = shift;
	$mutation_rate_scaling_factor = !defined $mutation_rate_scaling_factor ? 1.0 : $mutation_rate_scaling_factor;

	return 0 if ($mutation_rate_scaling_factor == 0);  # short circuit

	confess "ERROR: sequence not parsed" if (!$self->get_transcription_valid_flag());

	my $structure_ref = $self->get_structure_ref();
	my $mutation_rate_ref = $self->get_mutation_rate_ref();
	my $transcription_ref = $transcription_ref_of{$obj_ID};

	my $num_bits = 0;
	for (my $i=0; $i < @$structure_ref; $i++) {
	    my $field = $structure_ref->[$i][0];
	    my $mutation_rate = $mutation_rate_ref->{$field};
	    if (defined $structure_ref->[$i][2]) {  # is there a repeat sequence?
		for (my $j=0; $j < @{$transcription_ref->{$field}}; $j++) {
		    # mutate the field itself
		    if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
			$num_bits += $transcription_ref->{$field}[$j][0]->mutate($mutation_rate_scaling_factor * $mutation_rate);
		    } else {
			$num_bits += $self->mutate_field([$field, $j], $mutation_rate_scaling_factor * $mutation_rate);
		    }
		    # mutate the linker
		    my $linker_field = "${field}_linker";
		    my $linker_mutation_rate = $mutation_rate_ref->{$linker_field};
		    if (@{$transcription_ref->{$field}} > 1 && $j < @{$transcription_ref->{$linker_field}}) {
			$num_bits += $self->mutate_field([$linker_field, $j], $mutation_rate_scaling_factor * $linker_mutation_rate);
		    }
		}
	    } else {
		if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
		    $num_bits += $transcription_ref->{$field}[0]->mutate($mutation_rate_scaling_factor * $mutation_rate);
		} else {
		    $num_bits += $self->mutate_field($field, $mutation_rate_scaling_factor * $mutation_rate);
		}
	    }
	}

	printn "mutated: mutated $num_bits bits in ".$self->get_name() if $verbosity >= 2;
	return $num_bits;
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint
    # Synopsys: Recursively sprint the transcribed sequence.
    #--------------------------------------------------------------------------------------
    sub sprint {
	my $self = shift;
	my $obj_ID = ident $self;
	my %arg_ref = (
	    hierarchical_flag => 1,
	    colour_flag  => 1,
	    level => 0,
	    pos => $start_pos_of{$obj_ID},
	    context_flag => 0,
	    @_,
	   );
	check_args(\%arg_ref, 5);

	my $hierarchical_flag = $arg_ref{hierarchical_flag};
	my $colour_flag = $arg_ref{colour_flag};
	my $level = $arg_ref{level};
	my $indent = $hierarchical_flag ? (" " x (2*$level)) : "";
	my $newline = $hierarchical_flag ? "\n" : "";

	my $structure_ref = $self->get_structure_ref();

	my $start_pos = $start_pos_of{$obj_ID};
	my $transcription_ref = $transcription_ref_of{$obj_ID};

	my @field_colours = ("reset", "reverse");
	my $field_colour_index = 0;

	my @level_colours = ("green", "red", "yellow", "black");

	my $pos = $arg_ref{pos};
	my ($locus, $length);

	my $result = $hierarchical_flag ? "${indent}name => ".$self->get_name()."\n" : "";

	for (my $i=0; $i < @$structure_ref; $i++) {
	    my $field = $structure_ref->[$i][0];
	    if (defined $structure_ref->[$i][2]) {  # is there a repeat sequence?
		for (my $j=0; $j < @{$transcription_ref->{$field}}; $j++) {
		    # change to level colour
		    $result .= color($level_colours[$level]) if $colour_flag;
		    # print leader for untranscribed stuff
		    $locus = $transcription_ref->{$field}[$j][1];
		    $length = $transcription_ref->{$field}[$j][2];
		    $result .= ("." x ($locus - $pos)).($locus - $pos > 0 ? $newline : "");
		    $pos = $locus;
		    # print the field name with indentation and reset the colour
		    $result .= "$indent$field => " if $hierarchical_flag;
		    $result .= color("reset") if $colour_flag;

		    if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
			# print the field...
			$result .= "$newline".$transcription_ref->{$field}[$j][0]->sprint(%arg_ref, level => $level + 1, pos => $pos);
			$pos = $locus + $length;

			# ...and also print linker
			my $linker_field = "${field}_linker";
			$result .= color($level_colours[$level]) if $colour_flag;
			if (@{$transcription_ref->{$field}} > 1 && $j < @{$transcription_ref->{$linker_field}}) {
			    $result .= "$indent$linker_field => " if $hierarchical_flag;
			    $result .= "$transcription_ref->{$linker_field}[$j][0]$newline";
			    $pos = $transcription_ref->{$linker_field}[$j][1] + $transcription_ref->{$linker_field}[$j][2];
			}
			$result .= color("reset") if $colour_flag;
		    } else {
			# set level and field colour, print the field, reset colour
			$result .= color($field_colours[($field_colour_index++) % scalar(@field_colours)]) if $colour_flag;
			$result .= color($level_colours[$level]) if $colour_flag;
		 	$result .= "$transcription_ref->{$field}[$j][0]$newline";
			$pos = $locus + $length;
			$result .= color("reset") if $colour_flag;

			# ...and also print linker with level colour and indent
			my $linker_field = "${field}_linker";
			$result .= color($level_colours[$level]) if $colour_flag;
			if (@{$transcription_ref->{$field}} > 1 && $j < @{$transcription_ref->{$linker_field}}) {
			    $result .= "$indent$linker_field => "if $hierarchical_flag;
			    $result .= "$transcription_ref->{$linker_field}[$j][0]$newline";
			    $pos = $transcription_ref->{$linker_field}[$j][1] + $transcription_ref->{$linker_field}[$j][2];
			}
			$result .= color("reset") if $colour_flag;
		    }
		}
	    } else {
		# change to level colour
		$result .= color($level_colours[$level]) if $colour_flag;
		# print leader for untranscribed stuff
		$locus = $transcription_ref->{$field}[1];
		$length = $transcription_ref->{$field}[2];
		$result .= ("." x ($locus - $pos)).($locus - $pos > 0 ? $newline : "");
		$pos = $locus;
		# print the field name with indentation and reset the colour
		$result .= "$indent$field => " if $hierarchical_flag;
		$result .= color("reset") if $colour_flag;

		if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
		    $result .= "$newline".$transcription_ref->{$field}[0]->sprint(%arg_ref, level => $level + 1, pos => $pos);
		    $pos = $locus + $length;
		} else {
		    # set level and field colour, print the field, reset colour
		    $result .= color($field_colours[$field_colour_index++ % @field_colours]) if $colour_flag;
		    $result .= color($level_colours[$level]) if $colour_flag;
		    $result .= "$transcription_ref->{$field}[0]$newline";
		    $pos = $locus + $length;
		    $result .= color("reset") if $colour_flag;
		}
	    }
	}
	
	# if level 0, print out context
	if ($arg_ref{context_flag}) {
	    my $sequence_ref = $self->get_sequence_ref();
	    my $sequence = $sequence_ref->get_sequence();
	    $result .= color($level_colours[$level]) if $colour_flag;
	    $result .= ("." x (length($sequence) - $pos)).(length($sequence) - $pos > 0 ? $newline : "") if $level == 0;
	    $result .= color("reset") if $colour_flag;
	}
	
	return $result;
    }


    #--------------------------------------------------------------------------------------
    # Function: translate_field
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub translate_field {
	my $self = shift;

	confess "ERROR: sub-class must implement translate_field() routine";
    }

    #--------------------------------------------------------------------------------------
    # Function: translate
    # Synopsys: Recursively translate a sequence field-by-field.
    #--------------------------------------------------------------------------------------
    sub translate {
	my $self = shift;
	my $obj_ID = ident $self;

	confess "ERROR: sequence not parsed" if (!$self->get_transcription_valid_flag());

	my $structure_ref = $self->get_structure_ref();
	my $transcription_ref = $transcription_ref_of{$obj_ID};

	# init/clear translation results
	my $translation_ref = $translation_ref_of{$obj_ID} = {};

	for (my $i=0; $i < @$structure_ref; $i++) {
	    my $field = $structure_ref->[$i][0];
	    if (defined $structure_ref->[$i][2]) {  # is there a repeat sequence?
		for (my $j=0; $j < @{$transcription_ref->{$field}}; $j++) {
		    if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
			push @{$translation_ref->{$field}}, $transcription_ref->{$field}[$j][0]->translate();
		    } else {
			push @{$translation_ref->{$field}}, $self->translate_field([$field, $j]);
		    }
		}
	    } else {
		if (ref $structure_ref->[$i][1]) { # is this a sub-sequence?
		    $translation_ref->{$field} = $transcription_ref->{$field}[0]->translate();
		} else {
		    $translation_ref->{$field} = $self->translate_field($field);
		}
	    }
	}

	$translation_valid_flag_of{$obj_ID} = 1;
	return $translation_ref;
    }
}

sub run_testcases {

}


# Package BEGIN must return true value
return 1;

