#-#####################################################################################
#- File:     BindingProfile.pm
#- Synopsys: 
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package BindingProfile;
use Class::Std::Storable;
use base qw(BitVector);
{
    use Carp;
    use Utils;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # OPERATOR METHODS
    #######################################################################################
    #######################################################################################
    # Function: binding_complement
    # Synopsys: Returns L/R flipped, 1s-complement of the vector.
    #######################################################################################
    sub binding_complement {
	my $self = shift;
	my $class = ref $self || $self;
	# if class method call assume a string was supplied and create object on the fly
	$self = $class->new({vector => shift}) if !ref $self;
	my $obj_ID = ident $self;
		
	return $self->lr_flip()->NOT();
    }

    #--------------------------------------------------------------------------------------
    # Function: mismatch
    # Synopsys: A hamming distance-based measure of the mismatch between two BindingProfiles
    #           when the first is compared to L-R flipped, inverted version of the second.
    #--------------------------------------------------------------------------------------
    sub mismatch {
	my $class = shift;
	my $x_ref = shift;
	my $y_ref = shift;

	# if strings were supplied, create objects on the fly
	$x_ref = $class->new({vector => $x_ref}) if !ref $x_ref;
	$y_ref = $class->new({vector => $y_ref}) if !ref $y_ref;

	# apply masks
	while (@_) {
	    my $x_mask_ref = shift;
	    my $y_mask_ref = shift;
	    $x_mask_ref = BindingProfile->new({vector => $x_mask_ref}) if !ref $x_mask_ref;
	    $y_mask_ref = BindingProfile->new({vector => $y_mask_ref}) if !ref $y_mask_ref;
	    $x_ref = BindingProfile->XOR($x_ref, $x_mask_ref);
	    $y_ref = BindingProfile->XOR($y_ref, $y_mask_ref);
	}

	return BindingProfile->hamming_distance($x_ref, $y_ref->binding_complement());
    }

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################

}


sub run_testcases {
    my $ref1 = BindingProfile->new({vector => "111000"});
    my $ref2 = BindingProfile->new({vector => "100000"});
    printn $ref1->sprint();
    printn $ref2->sprint();
    printn ($ref2->lr_flip()->NOT()->sprint());
    printn (BindingProfile->mismatch($ref1, $ref2));

    my $ref3 = BindingProfile->new({vector => "101010"});
    my $ref4 = BindingProfile->new({vector => "101010"});
    printn $ref3->sprint();
    printn $ref4->sprint();
    printn ($ref4->lr_flip()->NOT()->sprint());
    printn (BindingProfile->mismatch($ref3, $ref4));
}


# Package BEGIN must return true value
return 1;

