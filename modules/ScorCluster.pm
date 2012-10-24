#!/usr/bin/perl -w
#-#####################################################################################
#- File:     ScorCluster.pm
#- Synopsys: Create and manage a cluster of genome scoring nodes.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################


use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package ScorCluster;
use Class::Std::Storable;
use base qw();
{
    use Carp;
    use Linux::Pid;

    use Utils;
    use ScorNode;
    use Globals qw($verbosity $TAG);

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    use vars qw($ctrl_C_flag);  # do NOT make this a 'my' variable!!
    $ctrl_C_flag = 0;
    use vars qw($sigusr1_flag);
    $sigusr1_flag = 0;

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %config_file_of :ATTR(get => 'config_file', set => 'config_file', init_arg => 'config_file');
    my %cluster_size_of :ATTR(get => 'cluster_size', set => 'cluster_size', init_arg => 'cluster_size');
    my %cluster_type_of :ATTR(get => 'cluster_type', set => 'cluster_type', init_arg => 'cluster_type');
    my %host_list_of :ATTR(get => 'host_list', set => 'host_list', init_arg => 'host_list');
    my %nice_of :ATTR(get => 'nice', set => 'nice', init_arg => 'nice');
    my %work_dir_of :ATTR(get => 'work_dir', set => 'work_dir', init_arg => 'work_dir');
    my %local_dir_of        :ATTR(get => 'local_dir', set => 'local_dir'); # custom init_arg
    my %scoring_class_of :ATTR(get => 'scoring_class', set => 'scoring_class', init_arg => 'scoring_class');
    my %node_list_of :ATTR(get => 'node_list', set => 'node_list');

    my %ppid_of :ATTR(get => 'ppid', set => 'ppid');
    my %cpid_of :ATTR(get => 'cpid', set => 'cpid');

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

	$local_dir_of{$obj_ID} = $arg_ref->{local_dir} if exists $arg_ref->{local_dir};
    }

    #--------------------------------------------------------------------------------------
    # Function: START
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub START {
        my ($self, $obj_ID, $arg_ref) = @_;

	my $config_file = $config_file_of{$obj_ID};
	my $cluster_size = $cluster_size_of{$obj_ID};
	my $cluster_type = $cluster_type_of{$obj_ID};
	my $host_list = $host_list_of{$obj_ID};
	my $nice = $nice_of{$obj_ID};
	my $work_dir = $work_dir_of{$obj_ID};
	my $local_dir = $local_dir_of{$obj_ID};
	my $scoring_class = $scoring_class_of{$obj_ID};
	
	$node_list_of{$obj_ID} = [];

	printn "ScorCluster::new() -- New ScorCluster cluster_size=$cluster_size, work_dir=$work_dir, scoring_class=$scoring_class, cluster_type=$cluster_type";
	# create tagged work_dir
	system ("mkdir -p $work_dir/$TAG");
	system ("mkdir -p $local_dir/$TAG") if defined $local_dir;
	# create the nodes
	my $timestamp = `date +%F_%T`; chomp($timestamp);
	$timestamp =~ s/[-:]//g;
	for (my $i = 0; $i < $cluster_size; $i++) {
	    push @{$node_list_of{$obj_ID}}, ScorNode->new({
		node_ID => $i,
		cluster_type => $cluster_type,
		nice => $nice,
		config_file => $config_file,
		work_dir => "$work_dir/$TAG",
		local_dir => (defined $local_dir) ? "$local_dir/$TAG" : undef,
		logfile => "ScorNode.$i.$timestamp.log",
	    });
	}

	return $self;
    }

    #--------------------------------------------------------------------------------------
    # Function: spawn_rrobin
    # Synopsys: Forks a child process whose function is to tickle the ScorNodes in a
    #           round-robin fashion.
    #--------------------------------------------------------------------------------------
    sub spawn_rrobin {
	my $self = shift; my $obj_ID = ident $self;
	my $cluster_type = $cluster_type_of{$obj_ID};
	my $cluster_size = $cluster_size_of{$obj_ID};
	my $scoring_class = $scoring_class_of{$obj_ID};
	my $nice = $nice_of{$obj_ID};
	my $pid;

      FORK: {
	    if ($pid = fork) {
		# parent here
		# child process pid is available in $pid
		$ppid_of{$obj_ID} = $$;
		$cpid_of{$obj_ID} = $pid;
		ScorNode->set_cpid($cpid_of{$obj_ID});

		# NOTES ON SIGINT
		# i)   When the user hits CTRL-C, both parent and child receive it
		# ii)  On normal exit, no one gets CTRL-C
		# iii) A SIGINT handler is necessary to call exit() and for DEMOLISH
		#      to be executed
		# iv)  Multiple SIGINTs may kill process prematurely?
		$SIG{INT} = sub {    # trapping CTRL-C ensures graceful exit via DEMOLISH/END BLOCKS etc.
		    printn "ScorCluster: Parent process $$ exiting from CTRL-C\n";
		    $ctrl_C_flag = 1;
		    exit;
		};
		$SIG{USR1} = sub {
		    printn "ScorCluster: Parent process $$ received SIGUSR1\n";
		    $sigusr1_flag = 1;
		};

		# init nodes
		printn "ScorCluster::start -- parent process started (pid=$$, child=$pid)";
		for (my $i = 0; $i < $cluster_size; $i++) {
		    printn "node $i status at start (parent): " . $node_list_of{$obj_ID}->[$i]->get_node_status() if $verbosity >=2;
		    $node_list_of{$obj_ID}->[$i]->init($scoring_class);
		}
		sleep $cluster_size + 1;
		printn "ScorCluster::start -- parent process returning";
		return $cpid_of{$obj_ID};
	    } elsif (defined $pid) { # $pid is zero here if defined
		# child here
		# parent process pid is available with getppid
		$ppid_of{$obj_ID} = getppid();
		$cpid_of{$obj_ID} = $$;
		ScorNode->set_cpid($cpid_of{$obj_ID});

		$SIG{INT} = sub {    # trapping CTRL-C ensures graceful exit via DEMOLISH/END BLOCKS etc.
		    printn "ScorCluster: Child process $$ received CTRL-C (ignoring)\n";
		    $ctrl_C_flag = 1;
		    # take no other action, with for SIGUSR1
		};
		$SIG{USR1} = sub {   # SIGUSR1 is exit signal from parent
		    printn "ScorCluster: Child process $$ received SIGUSR1 (shutting down)\n";
		    $sigusr1_flag = 1;
		    exit;
		};

		sleep 1;
		printn "ScorCluster::start -- child process started (pid=$$, parent=$ppid_of{$obj_ID})";

		my @host_list;
		if ($cluster_type eq "SSH") {
		    @host_list = @{$host_list_of{$obj_ID}};
		} else {
		    @host_list = ();
		}

		for (my $i = 0; $i < $cluster_size; $i++) {
		    printn "node $i status at start (child): " . $node_list_of{$obj_ID}->[$i]->get_node_status() if $verbosity >=2;
		    $node_list_of{$obj_ID}->[$i]->spawn((@host_list && $i != 0) ? $host_list[($i-1) % @host_list] : undef);
		}
		printn "ScorCluster::start -- child process done spawning, started multitask";
		while ($ppid_of{$obj_ID} == Linux::Pid::getppid()) { # while parent still alive
		    for (my $i = 0; $i < $cluster_size; $i++) {
			$node_list_of{$obj_ID}->[$i]->tickle();
		    }
		    sleep 1;
		}
		printn "ScorCluster::start -- parent died unexpectedly, cleaning up and exiting child...";
		ScorNode->clean_up_all(); # child process must clean up shared memory in case parent didn't
		exit;
	    } elsif ($! =~ /No more process/) {
		# EAGAIN, supposedly recoverable fork error
		sleep 5;
		redo FORK;
	    } else {
		# weird fork error
		die "Can't fork: $!\n";
	    }
	}

	# init processes

	printn "ScorCluster::new() -- Done.";
	return $self;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_free_node
    # Synopsys: Poll each node's status flag until one is found which is ready.
    #--------------------------------------------------------------------------------------
    sub get_free_node {
	my $self = shift; my $obj_ID = ident $self;

	printn "ScorCluster::get_free_node -- looking for free node...";
	while (1) {
	    for (my $i = 0; $i < $cluster_size_of{$obj_ID}; $i++) {
		if (defined $node_list_of{$obj_ID}->[$i]) {
		    printn "node $i status: " . $node_list_of{$obj_ID}->[$i]->get_node_status() if $verbosity >=3;
		    if ($node_list_of{$obj_ID}->[$i]->get_node_status()) {
			printn "ScorCluster::get_free_node -- node $i is free";
			return $node_list_of{$obj_ID}->[$i];
		    }
		} else {
		    printn "node $i status: undefined";
		    printn "ERROR: get_free_node -- undefined node";
		    exit(1);
		}
	    }
	    sleep 1;
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: get_busy_node_id_list
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_busy_node_id_list {
	my $self = shift; my $obj_ID = ident $self;
	
	my @busy_node_list = grep {$_->get_node_status() != 1} @{$node_list_of{$obj_ID}};
	@busy_node_list = map {$_->get_node_ID()} @busy_node_list;
	
	return @busy_node_list;
    }

    sub DEMOLISH {
	my $self = shift; my $obj_ID = ident $self;
	printn "ScorCluster::DEMOLISH called (pid $$, ctrl_C_flag=$ctrl_C_flag, sigusr1_flag=$sigusr1_flag)";
	if ($$ == $ppid_of{$obj_ID}) {
	    printn "---------------------------------------";
	    printn "Shutting down the cluster";
	    printn "---------------------------------------";
	    # We provide a SIGUSR1 signal to child to signal
	    # it to shutdown, then wait for the processing to occur.
	    # If user pressed CTRL-C, both processes got the signal, but
	    # the child's handler ignores it.
	    my $cpid = $cpid_of{$obj_ID};
	    printn "ScorCluster::DEMOLISH am parent (pid $$), shutting down child and waiting...";
	    kill ("USR1", $cpid);  # shutdown child process
	    my $child_done_filename = $work_dir_of{$obj_ID}."/$TAG/$cpid.cleanup";
	    my $timeout = 0;
	    while ((! -e $child_done_filename) && $timeout < 120) {
		sleep 1;
		$timeout++;
	    }
	    unlink $child_done_filename;
	    printn "ScorCluster::DEMOLISH am parent (pid $$), done cleanup...";
	    printn "---------------------------------------";
	    printn "Shutdown complete";
	    printn "---------------------------------------";
	} elsif ($$ == $cpid_of{$obj_ID}) {
	    printn "ScorCluster::DEMOLISH am child (pid $$)";
	    my $child_done_filename = $work_dir_of{$obj_ID}."/$TAG/$$.cleanup";
	    $node_list_of{$obj_ID} = undef;       # shutdown the ScorNodes
	    burp_file($child_done_filename, "")   # signal that ScorNodes have been DEMOLISHED
	} else {
	    printn "ERROR: neither parent nor child -- did you call spawn_rrobin?";
	}
    }
}


sub run_testcases {
    $TAG = "ScorCluster";
    system("rm -f test/modules/ScorCluster/*.log");  # remove old files

    my $config_file = <<END;
#----------------------------------------
# CPU AND CLUSTER SETTINGS
#----------------------------------------
nice = 15
vmem = 2000000
local_dir = scoring/localdir

#----------------------------------------
# WORKSPACE AND CUSTOM SCORING MODULE
#----------------------------------------
scoring_class = Scoring
work_dir = scoring

END

    burp_file("test/modules/ScorCluster.cfg", $config_file);

    my $ref = ScorCluster->new({
	config_file => "test/modules/ScorCluster.cfg",
	cluster_size => 2,
	cluster_type => "LOCAL",
	host_list => [],
	nice => 10,
	work_dir => "test/modules",
	local_dir => "test/modules/localdir",
	scoring_class => "Scoring",
    });
    $ref->spawn_rrobin();

    while ($ref->get_busy_node_id_list()) {
	sleep 1;
    }

    $ref = undef;
    printn;
    printn "LOGFILE test/modules/ScorCluster/ScorNode.0.*.log:";

    system("cat test/modules/ScorCluster/ScorNode.0.*.log");
}


# Package BEGIN must return true value
return 1;

