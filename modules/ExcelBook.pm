#-#####################################################################################
#- File:     ExcelBook.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package ExcelBook;
use Class::Std::Storable;
use base qw(ClassData Named Set);
{
    use Carp;

    use Spreadsheet::WriteExcel;

    use Utils;
    use ExcelSheet;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    ExcelBook->set_class_data("ELEMENT_CLASS", "ExcelSheet");

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
    # Function: add_sheet
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub add_sheet {
        my $self = shift;

        my $sheet_ref = ExcelSheet->new({@_});
        $self->add_element($sheet_ref);
        return $sheet_ref;
    }


    #--------------------------------------------------------------------------------------
    # Function: write
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub write {
        my $self = shift;
        my $name = $self->get_name();

        printn "Writing spreadsheet....";
        my $spreadsheet_ref = Spreadsheet::WriteExcel->new("$name");

        my @sheets = $self->get_elements();
        for (my $i=0; $i < @sheets; $i++) {
            my $sheet_ref = $sheets[$i];
            my @column_labels = @{$sheet_ref->get_column_labels()};
            my $worksheet_ref = $spreadsheet_ref->add_worksheet($sheet_ref->get_name());
            # write first line containing title
            my $sheet_title = $sheet_ref->get_title();
            $worksheet_ref->write(0, 0, $sheet_title);
            # write second line containing column labels
            for (my $k = 0; $k < @column_labels; $k++) {
                $worksheet_ref->write(1, $k, $column_labels[$k]) if defined $column_labels[$k];
            }
            # remaining lines consist of the matrix itself
            my $num_rows = $sheet_ref->get_num_rows();
            for (my $j = 0; $j < $num_rows; $j++) {
                my @row = $sheet_ref->get_row($j);
                for (my $k=0; $k < @row; $k++) {
                    $worksheet_ref->write($j+2, $k, $row[$k]) if defined $row[$k] && !ref $row[$k];
                    $worksheet_ref->write($j+2, $k, join ",",@{$row[$k]}) if defined $row[$k] && (ref $row[$k]) eq "ARRAY";
                }
            }
        }
    }
}


sub run_testcases {

    my $book_ref = ExcelBook->new({
            name => "test/modules/ExcelBook.xls",
        });

    my $fitness_sheet_ref = $book_ref->add_sheet(
        name => "fitness",
        title => "SHEET TITLE",
        column_labels => [qw(GENERATION FITNESS)],
    );

    my $other_sheet_ref = $book_ref->add_sheet(
        name => "other",
        title => "OTHER SHEET TITLE",
        column_labels => [qw(GENERATION OTHER)],
    );

    $fitness_sheet_ref->set_column(0, (0..10));

    $fitness_sheet_ref->set_element(0,1,0.1);
    $fitness_sheet_ref->set_element(5,1,0.2);

    $other_sheet_ref->set_element(6,4,0.4);

    $book_ref->write();
}


# Package BEGIN must return true value
return 1;

