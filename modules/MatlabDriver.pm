#-#####################################################################################
#- File:     MatlabDriver.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package MatlabDriver;
use Class::Std::Storable;
use base qw(Named);
{
    use Carp;
    use Utils;

    use Globals;

    use IPC::Open2;
    use FileHandle;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %READ_FH_OF  :ATTR(get => 'READ_FH', set => 'READ_FH',);
    my %WRITE_FH_OF :ATTR(get => 'WRITE_FH', set => 'WRITE_FH',);
    my %LOG_FH_OF   :ATTR(get => 'LOG_FH', set => 'LOG_FH',);

    my %logfile_of  :ATTR(get => 'logfile', set => 'logfile', init_arg => 'logfile');
    my %host_of     :ATTR(get => 'host', set => 'host', init_arg => 'host', default => 'localhost');
    my %vmem_of     :ATTR(get => 'vmem', set => 'vmem', init_arg => 'vmem', default => 'unlimited');
    my %options_of  :ATTR(get => 'options', set => 'options');

    my %echo_of     :ATTR(get => 'echo', set => 'echo', init_arg => 'echo', default => 0);

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################

    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

	$options_of{$obj_ID} = $arg_ref->{options} || "-nodesktop -nosplash";
    }

    #--------------------------------------------------------------------------------------
    # Function: START
    # Synopsys: Fork matlab process on given host and initialize I/O filehandles.
    #--------------------------------------------------------------------------------------
    sub START {
        my ($self, $obj_ID, $arg_ref) = @_;

	my $logfile = $logfile_of{$obj_ID};
	my $options = $options_of{$obj_ID};
	my $host = $host_of{$obj_ID};
	my $vmem = $vmem_of{$obj_ID};

	my $matlab_command;
	if ($host eq "localhost") {
	    $matlab_command = "ulimit -S -v $vmem; ulimit -a; matlab $options";
	} else {
	    $matlab_command = "ssh -X $host \"cd $ENV{PWD}; ulimit -S -v $vmem; ulimit -a; matlab $options\"";
	}

	my $LOG_FH;
	open ($LOG_FH, ">$logfile") or die "Couldn't open $logfile for writing.\n";
	$LOG_FH->autoflush(1);

	# start the matlab process
	my ($READ_FH, $WRITE_FH);
	my $date = `date`; chomp($date);
	printn "MatlabDriver::START: starting ".$self->get_name()." on $host ($date)";
	printn "MatlabDriver::START: logfile is $logfile";
	printn "MatlabDriver::START: executing \"$matlab_command\"";
	my $shell_pid = open2($READ_FH, $WRITE_FH, $matlab_command) or die ("Couldn't start MATLAB process: $!");
	printn "MatlabDriver::START: MATLAB shell pid=$shell_pid";

	# store filehandles
	$READ_FH_OF{$obj_ID} = $READ_FH;
	$WRITE_FH_OF{$obj_ID} = $WRITE_FH;
	$LOG_FH_OF{$obj_ID} = $LOG_FH;

	$self->wait_on("www.mathworks.com");
    }

    sub DEMOLISH {
	my ($self, $obj_ID) = @_;

	printn "MatlabDriver::DEMOLISH: Closing ".$self->get_name();
	$self->cmd("exit");      # exit from matlab
    }

    #--------------------------------------------------------------------------------------
    # Function: cmd
    # Synopsys: Wrapper for sending a command to MATLAB process.
    #--------------------------------------------------------------------------------------
     sub cmd {
 	my $self = shift;
 	my $obj_ID = ident $self;

	my $LOG_FH = $LOG_FH_OF{$obj_ID};
	my $WRITE_FH = $WRITE_FH_OF{$obj_ID};

 	foreach my $string (@_) {
	    print "$string\n" if $echo_of{$obj_ID};
 	    print $LOG_FH "$string\n";     # echo commands in logfile
 	    print $WRITE_FH "$string\n";
 	}
     }

    #--------------------------------------------------------------------------------------
    # Function: wait_on
    # Synopsys: Read from filehandle until expected output, return line that matched.
    #--------------------------------------------------------------------------------------
    sub wait_on {
	my $self = shift; my $obj_ID = ident $self;
	my $cue = shift;

	my $READ_FH = $READ_FH_OF{$obj_ID};
	my $LOG_FH = $LOG_FH_OF{$obj_ID};

	my @problem = ();
	while (<$READ_FH>) {
	    print $_ if $echo_of{$obj_ID};
	    print $LOG_FH $_;
	    if ($_ =~ /$cue/) {
		return ($_, @problem);
	    }
	    if ($_ =~ /^sim time/i) {
		print $_ if $verbosity >= 3;
	    }
	    if ($_ =~ /Error/i) {
		printn "ERROR: matlab_wait_on -- matlab appears to have reported an error";
		printn "  --> $_";
		push @problem, $_;
	    }
	    if ($_ =~ /Warning/i) {
		printn "ERROR: matlab_wait_on -- matlab appears to have reported an warning";
		printn "  --> $_";
		push @problem, $_;
	    }
	    if ($_ =~ /out of memory/i) {
		printn "ERROR: matlab_wait_on -- matlab appears to have run out of memory";
		printn "  --> $_";
		push @problem, $_;
	    }
	    if ($_ =~ /Undefined function or variable/i) {
		printn "ERROR: matlab_wait_on -- matlab appears to be missing a function or variable";
		printn "  --> $_";
		push @problem, $_;
	    }
	}

	printn "ERROR: matlab_wait_on -- no more output";
	exit;
    }
}


sub run_testcases {

    $verbosity = 3;

    my $d1_ref = MatlabDriver->new({
	name => "driver1",
	logfile => "test/modules/MatlabDriver.localhost1.log",
	host => "localhost",
	vmem => 2000000,
	echo => 1,
    });
	
    $d1_ref->cmd("3 + 4");
    $d1_ref->wait_on("7");

    my $d2_ref = MatlabDriver->new({
	name => "driver2",
	logfile => "test/modules/MatlabDriver.localhost2.log",
	host => "localhost",
	vmem => 2000000,
	echo => 1,
    });
	
    $d2_ref->cmd("5 + 20");
    $d2_ref->wait_on("25");

}


# Package BEGIN must return true value
return 1;

