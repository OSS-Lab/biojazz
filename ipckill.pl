#!/usr/bin/perl -w

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

my $username = `whoami`;
chomp $username;

open(FH, "ipcs -s |");

while (<FH>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /\S+\s*(\S+)\s*$username/) {
	my $id = $1;
	print "found user sema $line\n";
	system ("./ipckill_c 0 $id");
    }
}

close(FH);

open(FH, "ipcs -m |");

while (<FH>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /\S+\s*(\S+)\s*$username/) {
	my $id = $1;
	print "found user shmem $line\n";
	system ("./ipckill_c 1 $id");
    }
}

close(FH);



# use ipcs -s/-m to get list of semas and shms

