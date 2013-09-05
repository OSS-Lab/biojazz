#-#####################################################################################
#- File:     Sequence.pm
#- Synopsys: A sequence of bits.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Sequence;
use Class::Std::Storable;
use base qw();
{
    use Carp;

    use Utils;
    use Globals;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %sequence_of :ATTR(get => 'sequence', set => 'sequence', init_arg => 'sequence', default => "");

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
    # Function: as_str
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    # Convert object to a string (automatically in string contexts)...
    # (n.b. this function gets called when printing out the object reference)
    sub as_str : STRINGIFY {
        my $self = shift;
        my $class = ref $self;

        return "$class=SCALAR(".$self->get_sequence().")";
    }

    #--------------------------------------------------------------------------------------
    # Function: get_length
    # Synopsys: Get sequence length.
    #--------------------------------------------------------------------------------------
    sub get_length {
        my $self = shift; my $obj_ID = ident $self;

        return length($sequence_of{$obj_ID});
    }

    #--------------------------------------------------------------------------------------
    # Function: get_subseq
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_subseq {
        my $self = shift; my $obj_ID = ident $self;
        my $locus = shift || 0;
        my $length = shift; $length = !defined $length ? length($sequence_of{$obj_ID}) : $length;

        return substr($sequence_of{$obj_ID}, $locus, $length);
    }

    #--------------------------------------------------------------------------------------
    # Function: set_subseq
    # Synopsys: Sets part of the sequence to given sub-sequence.
    #           (By default, replaces entire sequence with given).
    #--------------------------------------------------------------------------------------
    sub set_subseq {
        my $self = shift; my $obj_ID = ident $self;
        my $subseq = shift; $subseq = !defined $subseq ? "" : $subseq;
        my $locus = shift || 0;
        my $length = shift; $length = !defined $length ? length($subseq) : $length;

        return substr($sequence_of{$obj_ID}, $locus, $length, $subseq);
    }

    #--------------------------------------------------------------------------------------
    # Function: splice_subseq
    # Synopsys: Will splice in given sequence of at locus, replacing given length of
    #           original sequence (inserting by default).
    #--------------------------------------------------------------------------------------
    sub splice_subseq {
        my $self = shift; my $obj_ID = ident $self;
        my $subseq = shift; $subseq = !defined $subseq ? "" : $subseq;
        my $locus = shift || length($sequence_of{$obj_ID});  # append by default
        my $length = shift || 0;                             # insert by default

        return substr($sequence_of{$obj_ID}, $locus, $length, $subseq);
    }

    #--------------------------------------------------------------------------------------
    # Function: generate_random_sequence
    # Synopsys: Generate a sequence of a given length in bits.
    #--------------------------------------------------------------------------------------
    sub generate_random_sequence {
        my $self = shift;
        my $sequence_length = shift;

        printn "generate_random_sequence: generating sequence of length $sequence_length" if ($verbosity >= 2);
        confess "ERROR: generate_random_sequence -- need positive argument" if ($sequence_length <= 0);

        my @sequence = ();

        while ($sequence_length--) {
            my $bit = int(rand(2));
            push @sequence, $bit;
        }

        # return one big string
        return $sequence_of{ident $self} = join "", @sequence;
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate_subseq
    # Synopsys: Mutate given length of bits from starting at locus, with given probability.
    #           Returns number of bits mutated.
    #--------------------------------------------------------------------------------------
    sub mutate_subseq {
        my $self = shift; my $obj_ID = ident $self;
        my $probability = shift || 0.0;
        my $locus = shift || 0;
        my $length = shift; $length = (!defined $length) ? length($sequence_of{$obj_ID}) - $locus : $length;

        my $stop = $locus + $length;
        printn "mutate_subseq: mutating $length sequence bits starting at $locus with probability $probability" if ($verbosity >= 3);

        if (length($sequence_of{$obj_ID}) < $locus + $length) {
            confess "ERROR: mutate_subseq -- can't mutate past end of sequence";
        }

        my $subseq = $self->get_subseq($locus, $length);

        my $num_bits = 0;
        for (my $i=0; $i < $length; $i++) {
            if ((rand 1) < $probability) {
                my $value = substr($subseq, $i, 1);
                confess "ERROR: mutate_subseq -- sequence appears corrupted" if ($value ne "1") && ($value ne "0");
                my $new_value = $value eq "1" ? "0" : "1";
                substr($subseq, $i, 1) = $new_value;
                printn "mutate_subseq: mutating bit at position ".($locus+$i).", $value -> $new_value" if ($verbosity >= 3);
                $num_bits++;
            }
        }
        printn "mutate_subseq: mutated $num_bits sequence bits" if ($verbosity >= 3);

        $self->set_subseq($subseq, $locus, $length);

        return $num_bits;
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint
    # Synopsys: Alias for get_subseq.
    #--------------------------------------------------------------------------------------
    sub sprint {
        my $self = shift;

        return $self->get_subseq(@_);
    }

    #--------------------------------------------------------------------------------------
    # Function: load_sequence
    # Synopsys: Load a sequence from a file
    #--------------------------------------------------------------------------------------
    sub load_sequence {
        my $self = shift;
        my $filename = shift;

        open (SEQUENCE, "<$filename") or die "ERROR: can't read sequence file\n";

        my $sequence = "";

        my $bit;
        while (read SEQUENCE, $bit, 1) {
            confess "ERROR: file $filename is not a sequence" if (($bit ne '0') && ($bit ne '1'));
            $sequence .= $bit;
        }

        $sequence_of{ident $self} = ($sequence);
    }

    #--------------------------------------------------------------------------------------
    # Function: save_sequence
    # Synopsys: Save a sequence to a file
    #--------------------------------------------------------------------------------------
    sub save_sequence {
        my $self = shift;
        my $filename = shift;

        die "ERROR: save_sequence is not implemented";
    }
}


sub run_testcases {
    use Globals;
    $verbosity = 3;

    my $s_ref = Sequence->new({sequence => "11110000"});
    printn $s_ref->sprint();
    printn $s_ref->get_subseq(2,4);
    printn $s_ref->set_subseq("0011",2,4);
    printn $s_ref->sprint();
    printn $s_ref->mutate_subseq(1.0,2,4);
    printn $s_ref->sprint();
    printn $s_ref;  # as_str test
    printn $s_ref->splice_subseq("1010");
    printn $s_ref->sprint;

    printn "STORABLE TEST";
    use Storable;
    my $ice_ref = Storable::freeze($s_ref);
    my $water_ref = Storable::thaw($ice_ref);
    printn $s_ref->_DUMP();
    printn $water_ref->_DUMP();
}


# Package BEGIN must return true value
return 1;

