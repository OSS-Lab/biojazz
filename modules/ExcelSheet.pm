#-#####################################################################################
#- File:     ExcelSheet.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package ExcelSheet;
use Class::Std::Storable;
use base qw(Named);
{
    use Carp;
    use Utils;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %column_labels_of :ATTR(get => 'column_labels', set => 'column_labels');
    my %matrix_ref_of    :ATTR(get => 'matrix_ref', set => 'matrix_ref');
    my %title_of         :ATTR(get => 'title', set => 'title', init_arg => 'title');

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

        $column_labels_of{$obj_ID} = defined $arg_ref->{column_labels} ? $arg_ref->{column_labels} : [];
        $matrix_ref_of{$obj_ID} = [];
    }

    #--------------------------------------------------------------------------------------
    # Function: set_element
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub set_element {
        my $self = shift;  my $obj_ID = ident $self;
        my ($row, $col) = (shift, shift);
        my $value = shift;

        my $matrix_ref = $matrix_ref_of{$obj_ID};

        $matrix_ref->[$row] = [] if !defined $matrix_ref->[$row];
        $matrix_ref->[$row][$col] = $value;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_element
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_element {
        my $self = shift;  my $obj_ID = ident $self;
        my ($row, $col) = (shift, shift);

        my $matrix_ref = $matrix_ref_of{$obj_ID};

        return defined $matrix_ref->[$row] ? $matrix_ref->[$row][$col] : undef;
    }

    #--------------------------------------------------------------------------------------
    # Function: set_row
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub set_row {
        my $self = shift;  my $obj_ID = ident $self;
        my $row = shift;

        my @values = @_;

        my $matrix_ref = $matrix_ref_of{$obj_ID};
        $matrix_ref->[$row] = \@values;
    }

    #--------------------------------------------------------------------------------------
    # Function: set_column
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub set_column {
        my $self = shift;  my $obj_ID = ident $self;
        my $col = shift;

        my @values = @_;

        map {$self->set_element($_, $col, $values[$_])} (0..@values-1);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_row
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_row {
        my $self = shift;  my $obj_ID = ident $self;
        my $row = shift;

        my $matrix_ref = $matrix_ref_of{$obj_ID};

        return defined $matrix_ref->[$row] ? @{$matrix_ref->[$row]} : ();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_num_rows
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_num_rows {
        my $self = shift;  my $obj_ID = ident $self;

        my $matrix_ref = $matrix_ref_of{$obj_ID};

        return scalar(@{$matrix_ref});
    }

}


sub run_testcases {

    printn "NO TESTCASES!!!";
}


# Package BEGIN must return true value
return 1;

