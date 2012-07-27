#-#####################################################################################
#- File:     BitString.pm
#- Synopsys: Package for bit string manipulation and conversion
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

#######################################################################################
# TO-DO LIST
#######################################################################################

#######################################################################################
# Package interface
#######################################################################################
package BitString;

use strict;
use Carp;

use Exporter;
use vars qw(@ISA @EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(
	     dec2bin
	     bin2dec
);


#######################################################################################
# Function: dec2bin/bin2dec
# Synopsys: Convert from bit string to decimal and back.  Max 32 bits!!
#######################################################################################
# From Perl Cookbook, recipe 2.4
#######################################################################################
sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

sub bin2dec {
    my $arg = shift;
    confess "ERROR: bin2dec argument is more than 32 bits" if length($arg) > 32;
    return unpack("N", pack("B32", substr("0" x 32 . $arg, -32)));
}

