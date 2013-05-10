#-#####################################################################################
#- File:     BitVector.pm
#- Synopsys: 
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package BitVector;
use Class::Std::Storable;
use base qw();
{
    use Carp;
    use Utils;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %vector_of :ATTR(get => 'vector', set => 'vector', init_arg => 'vector');
    my %width_of :ATTR(get => 'width', set => 'width', init_arg => 'width', default => -1);

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # OPERATOR METHODS
    #######################################################################################

    #######################################################################################
    # Function: flip_lr (CONSTRUCTOR)
    # Synopsys: Bit flip (R/L) the vector and return result.
    #######################################################################################
    sub lr_flip {
        my $self = shift;
        my $class = ref $self || $self;
        # if class method call assume a string was supplied and create object on the fly
        $self = $class->new({vector => shift}) if !ref $self;
        my $obj_ID = ident $self;

        my @vector = @{$vector_of{$obj_ID}};
        my $width = $width_of{$obj_ID};

        confess "ERROR: internal error -- incorrect width" if $width != @vector;

        # flipping ?
        # what is the map operator and the # symbol
        my @flipped = map {$vector[$#vector - $_] } (0..$#vector);

        return $class->new({vector => \@flipped, width => $width});
    }

    #######################################################################################
    # Function: NOT (CONSTRUCTOR)
    # Synopsys: Returns the 1s-complement of the vector.
    #######################################################################################
    sub NOT {
        my $self = shift;
        my $class = ref $self || $self;
        # if class method call assume a string was supplied and create object on the fly
        $self = $class->new({vector => shift}) if !ref $self;
        my $obj_ID = ident $self;

        my @vector = @{$vector_of{$obj_ID}};
        my $width = $width_of{$obj_ID};

        confess "ERROR: internal error -- incorrect width" if $width != @vector;

        my @NOT = map {$_ eq "1" ? "0" : "1"} @vector;

        return $class->new({vector => \@NOT, width => $width});
    }

    #######################################################################################
    # Function: XOR (CONSTRUCTOR)
    # Synopsys: Returns the bitwise XOR of the (at least 2) arguments.
    #######################################################################################
    sub XOR {
        my $class = shift;
        my ($x_ref, $y_ref) = (shift, shift);

        # if strings were supplied as arguments, create objects on the fly
        $x_ref = $class->new({vector => $x_ref}) if !ref $x_ref;
        $y_ref = $class->new({vector => $y_ref}) if !ref $y_ref;

        my $x_ID = ident $x_ref;
        my $y_ID = ident $y_ref;

        my @x_vector = @{$vector_of{$x_ID}};
        my @y_vector = @{$vector_of{$y_ID}};
        my $x_width = $width_of{$x_ID};
        my $y_width = $width_of{$y_ID};

        confess "ERROR: internal error -- incorrect width (x)" if $x_width != @x_vector;
        confess "ERROR: internal error -- incorrect width (y)" if $y_width != @y_vector;
        confess "ERROR: different sizes" if $y_width != $x_width;

        my $width = $x_width;
        my @xor = map {int($x_vector[$_]) ^ int($y_vector[$_])} (0..$width-1);

        my $result = $class->new({vector => \@xor, width => $width});

        if (@_) {  # XOR again if more arguments
            $result = $class->XOR($result, @_);
        }
        return $result;
    }

    #######################################################################################
    # Function: hamming_distance
    # Synopsys: Return the hamming distance between two vectors.
    #######################################################################################
    sub hamming_distance {
        my $class = shift;
        my ($x_ref, $y_ref) = @_;

        # if strings were supplied as arguments, create objects on the fly
        $x_ref = $class->new({vector => $x_ref}) if !ref $x_ref;
        $y_ref = $class->new({vector => $y_ref}) if !ref $y_ref;

        my $xor_ref = BitVector->XOR($x_ref, $y_ref);

        return $xor_ref->ones();
    }

    #######################################################################################
    # Function: zeroes/ones
    # Synopsys: Returns the number of zeroes/ones in the vector.
    #######################################################################################
    sub zeroes {
        my $self = shift;
        my $class = ref $self || $self;
        # if class method call assume a string was supplied and create object on the fly
        $self = $class->new({vector => shift}) if !ref $self;

        my $vector = $self->sprint();

        my $count = $vector =~ tr/0//;
        return $count;
    }
    sub ones {
        my $self = shift;
        my $class = ref $self || $self;
        # if class method call assume a string was supplied and create object on the fly
        $self = $class->new({vector => shift}) if !ref $self;

        my $vector = $self->sprint();

        my $count = $vector =~ tr/1//;
        return $count;
    }

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
    sub START {
        my ($self, $obj_ID, $arg_ref) = @_;

        my $vector = $vector_of{$obj_ID};
        $vector = join "", @{$vector} if ref $vector;  # convert to string if list given
        my $width = $width_of{$obj_ID};
        my $init_length = length($vector);

        # use vector length if width argument supplied
        $width = $width_of{$obj_ID} = $init_length if $width == -1;

        # check initializers
        confess "ERROR: vector is longer than specified length" if ($init_length > $width);

        # zero-pad
        $vector = ("0" x ($width - $init_length)) . $vector if ($init_length < $width);

        # convert to list
        $vector_of{$obj_ID} = [split "", $vector];
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub sprint {
        my $self = shift;
        return join "", @{$vector_of{ident $self}};
    }
}


sub run_testcases {
    my $ref1 = BitVector->new({vector => "010110", width => 16});
    printn $ref1->sprint()." (ref1)";
    printn $ref1->lr_flip->sprint();
    printn $ref1->NOT->sprint();

    my $ref2 = BitVector->XOR($ref1, $ref1->lr_flip);
    printn $ref2->sprint()." (ref2)";
    my $ref3 = BitVector->XOR($ref2, $ref1);
    printn $ref3->sprint()." (ref3)";

    printn $ref1->zeroes();
    printn $ref2->zeroes();
    printn $ref3->zeroes();

    printn $ref1->ones();
    printn $ref2->ones();
    printn $ref3->ones();


    my $ref4 = BitVector->new({vector => "1111111111101001"});
    printn $ref4->sprint()." (ref4)";

    my $hdist = BitVector->hamming_distance($ref3, $ref4);
    printn $hdist;

    my $ref5 = $ref1->lr_flip()->NOT();
    printn $ref5->sprint()." (ref5)";

    # test with string arguments
    printn (BitVector->lr_flip("11100")->sprint());
    printn (BitVector->NOT("10100")->sprint());
    printn (BitVector->ones("10100"));
    printn (BitVector->zeroes("10100"));
    printn (BitVector->hamming_distance("10100", "00100"));
    printn (BitVector->XOR("0001", "1111", "1000")->sprint());
    printn (BitVector->XOR("1100011111111111000000000000000000000000000000011", "1000011111111111000000000000000000000000000000001")->sprint());
}


# Package BEGIN must return true value
return 1;

