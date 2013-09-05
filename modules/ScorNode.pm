#-#####################################################################################
#- File:     ScorNode.pm
#- Synopsys: Create (i.e. spawn) and manage a single genome scoring node using Expect.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#- A shared memory FIFO is implemented to queue print/expect commands to the
#- node.  Commands are queued to the FIFO using the node_print() and
#- node_expect() routines.  One or more commands are pulled off this queue
#- and processed when tickle() is called.
#
#- The FIFO (and a status flag) are in shared memory because it is expected
#- that a parent process will queue node print/expect commands, while
#- a child process tickles the nodes periodically to unqueue and process them.
#-#####################################################################################
#
#-#####################################################################################
use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package ScorNode;
use Class::Std::Storable;

use base qw();
{
    use Carp;
    use English;

    use Utils;
    use Globals qw($verbosity $TAG);

    use Expect;                       # for interaction w/ processes expecting a terminal
    $Expect::Multiline_Matching = 1;  # makes ^ and $ match as expected

    use IPC::Shareable;  # share variables with other process

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    my $ppid = $$;
    my $cpid;
    my @node_IO;
    # key is IPC_PRIVATE by default, so only child processes can share
    # also, setting destroy option should cause parent process to destroy the shared memory
    my $node_IO_handle = tie @node_IO, 'IPC::Shareable', { create => 1, destroy => 1 };
    confess "ERROR: could not initialize shared memory" if !defined $node_IO_handle;

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %node_ID_of         :ATTR(get => 'node_ID', set => 'node_ID', init_arg => 'node_ID');
    my %cluster_type_of    :ATTR(get => 'cluster_type', set => 'cluster_type', init_arg => 'cluster_type');
    my %nice_of            :ATTR(get => 'nice', set => 'nice', init_arg => 'nice');
    my %work_dir_of        :ATTR(get => 'work_dir', set => 'work_dir', init_arg => 'work_dir');
    my %local_dir_of        :ATTR(get => 'local_dir', set => 'local_dir'); # custom init_arg
    my %logfile_of         :ATTR(get => 'logfile', set => 'logfile', init_arg => 'logfile');
    my %config_file_of     :ATTR(get => 'config_file', set => 'config_file', init_arg => 'config_file');

    my %shell_ref_of       :ATTR(get => 'shell_ref', set => 'shell_ref');
    my %last_problem_of    :ATTR(get => 'last_problem', set => 'last_problem', default => 0);

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: clean_up_all
    # Synopsys: Class method to explicitely clean up shared memory
    #--------------------------------------------------------------------------------------
    sub clean_up_all {
        my $class = shift;
        if (tied @node_IO) {
            printn "ScorNode::clean_up_all: Cleaning up tied variable ($$)";
            IPC::Shareable->clean_up_all();
            return 1;
        } else {
            printn "ScorNode::clean_up_all: Not tied ($$)";
            return 0;
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: set_cpid
    # Synopsys: Set the cpid after a fork.
    #--------------------------------------------------------------------------------------
    sub set_cpid {
        my $class = shift;

        $cpid = shift;
    }

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

        my $node_ID = $node_ID_of{$obj_ID};
        my $work_dir = $work_dir_of{$obj_ID};
        my $local_dir = $local_dir_of{$obj_ID};
        my $logfile = $logfile_of{$obj_ID};
        my $cluster_type = $cluster_type_of{$obj_ID};

        printn "Creating ScorNode ID=$node_ID logfile=".(defined $local_dir ? $local_dir : $work_dir)."/$logfile cluster_type=$cluster_type";

        $self->set_node_status(0);

        system("mkdir -p $work_dir");
        system("mkdir -p $local_dir") if defined $local_dir;
        system("rm -f $work_dir/$logfile");
        system("rm -f $local_dir/$logfile") if defined $local_dir;
        @{$node_IO[$node_ID]->{command_queue}} = ();
    }

    #--------------------------------------------------------------------------------------
    # Function: DEMOLISH
    # Synopsys: Report pid, move logfile from local_dir to work_dir.
    #--------------------------------------------------------------------------------------
    sub DEMOLISH {
        my $self = shift; my $obj_ID = shift;
        printn "ScorNode::DEMOLISH called (pid $$)";

        if (!defined $cpid || $$ == $cpid) {  # child process or didn't fork?
            $shell_ref_of{$obj_ID}->print(chr(0x03)."\n");   # send CTRL-C to terminal
            $shell_ref_of{$obj_ID}->print("exit\n");         # exit the shell (ssh, qsub, etc.)
            $shell_ref_of{$obj_ID}->expect(10, "Scoring::DEMOLISH: done");  # wait for DEMOLISH

        # commented out this line: works but grabs special characters (the ctrl-c?) 
        # that cause problems when cat'ing the file
        #$shell_ref_of{$obj_ID}->expect(1, "DUMMY");                     # flush to logfile

            $shell_ref_of{$obj_ID} = undef; # destroy Expect object

            my $work_dir = $work_dir_of{$obj_ID};
            my $local_dir = $local_dir_of{$obj_ID};
            my $logfile = $logfile_of{$obj_ID};
            if (defined $local_dir) {
                printn "Moving $local_dir/$logfile to $work_dir/$logfile";
                system("mv $local_dir/$logfile $work_dir/$logfile");
            }
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: spawn
    # Synopsys: Spawn shell process and initialize associated tty.
    #--------------------------------------------------------------------------------------
    sub spawn {
        my $self = shift; my $obj_ID = ident $self;
        my $host = shift;

        my $cluster_type = $self->get_cluster_type();
        my $node_ID = $node_ID_of{$obj_ID};

        # create the pseudo-term
        my $shell_ref = $shell_ref_of{$obj_ID} = new Expect;

        # the spawn command does a fork/exec and returns FHs to object
        if ($cluster_type eq "LOCAL" || $node_ID == 0) {
            printn "ScorNode::spawn spawning local node";
            $shell_ref->spawn("bash ") or die "Couldn't start program $!\n";
        } elsif ($cluster_type eq "SSH") {
            confess "ERROR: host argument is undefined" if !defined $host;
            sleep 1;	   # prevent race condition on Xauthority file
            printn "ScorNode::spawn spawning node on $host";
            $shell_ref->spawn("ssh -X $host") or die "Couldn't start program $!\n";
        } elsif ($cluster_type eq "PBS") {
            printn "ScorNode::spawn spawning PBS node";
            $shell_ref->raw_pty(1);
            my $delay = $node_ID * 4;
            #$shell_ref->spawn("sleep $delay; /usr/pbs/bin/qsub -I genalg.pbs") or die "Couldn't start program $!\n";
            $shell_ref->spawn("/usr/pbs/bin/qsub -I genalg.pbs") or die "Couldn't start program $!\n";
            $shell_ref->stty(qw(raw -echo));
        } else {
            printn "ERROR: ScorNode::spawn() -- invalid node type";
            exit(1);
        }

        my $work_dir = $work_dir_of{$obj_ID};
        my $local_dir = $local_dir_of{$obj_ID};
        my $logfile = defined $local_dir ? "$local_dir/$logfile_of{$obj_ID}" : "$work_dir/$logfile_of{$obj_ID}";

        $shell_ref->log_stdout(0);
        $shell_ref->log_file($logfile);
        $shell_ref->autoflush(1);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_node_status
    # Synopsys: Get the node's ready flag.
    #--------------------------------------------------------------------------------------
    sub get_node_status {
        my $self = shift;

        $node_IO_handle->shlock();
        my $ready = $node_IO[$self->get_node_ID()]->{ready};
        $node_IO_handle->shunlock();

        return $ready;
    }

    #--------------------------------------------------------------------------------------
    # Function: set_node_status
    # Synopsys: Set the node's status flag.  0=BUSY, 1=READY.
    #--------------------------------------------------------------------------------------
    sub set_node_status {
        my $self = shift;
        my $status = shift;

        $node_IO_handle->shlock();
        $node_IO[$self->get_node_ID()]->{ready} = $status;
        $node_IO_handle->shunlock();
    }

    #--------------------------------------------------------------------------------------
    # Function: node_expect
    # Synopsys: Queue an expect command to given node and set status flag to BUSY(0).
    #--------------------------------------------------------------------------------------
    sub node_expect {
        my $self = shift;
        my $timeout = shift;
        my $regexp = shift;

        $self->set_node_status(0);

        $node_IO_handle->shlock();
        push @{$node_IO[$self->get_node_ID()]->{command_queue}}, ["expect", $timeout, $regexp];
        $node_IO_handle->shunlock();
    }

    #--------------------------------------------------------------------------------------
    # Function: node_print
    # Synopsys: Queue a print command to given node and set status flag to BUSY(0)
    #--------------------------------------------------------------------------------------
    sub node_print {
        my $self = shift;
        my $arg = shift;

        $self->set_node_status(0);

        $node_IO_handle->shlock();
        push @{$node_IO[$self->get_node_ID()]->{command_queue}}, ["print", $arg];
        $node_IO_handle->shunlock();
    }

    #--------------------------------------------------------------------------------------
    # Function: init
    # Synopsys: Wait for bash shell to be ready and start perl shell.
    #--------------------------------------------------------------------------------------
    sub init {
        my $self = shift; my $obj_ID = ident $self;
        my $scoring = shift;

        confess "ERROR: specify scoring class as argument" if !defined $scoring;
        confess "ERROR: TAG is not defined" if !defined $TAG;

        my $project_dir = `pwd`;
        my $work_dir = $work_dir_of{$obj_ID};
        my $local_dir = $local_dir_of{$obj_ID};

        my $node_ID = $self->get_node_ID();

        # change prompt to something known
        $self->node_print("PS1='[ScorNode $node_ID]'\n");

        $self->node_expect(undef, "\\[ScorNode $node_ID]");
        $self->node_print("echo my shell is ready\n");

        $self->node_expect(undef, "^my shell is ready");
        $self->node_expect(undef, "[ScorNode $node_ID]");

        $self->node_print("stty ocrnl -onlcr\n"); # prevents line wrap-around without carriage return
        $self->node_expect(undef, "[ScorNode $node_ID]");

        $self->node_print("date\n");
        $self->node_expect(undef, "[ScorNode $node_ID]");

        $self->node_print("hostname\n");
        $self->node_expect(undef, "[ScorNode $node_ID]");

        $self->node_print("cd $project_dir\n");
        $self->node_expect(undef, "[ScorNode $node_ID]");

        my $nice = $self->get_nice();

        use FindBin qw($Bin);  # need application path
        my $command = ("nice -$nice perl -I$Bin/modules -Icustom -I$ENV{ANC_HOME}/base -MUtils -MGenomeModel -M$scoring ".
            "-e \'\$|=1; print \"PERL_SHELL> \";while(<STDIN>){printn \$_;eval \"\$_\";".
            "print \"\$@\" if \$@;print \"PERL_SHELL> \";}'"
        );
        $self->node_print("$command\n");
        $self->node_expect(undef, '^PERL_SHELL>');

        # initialize Globals
        $self->node_print("use Storable qw(store retrieve);\n");
        $self->node_expect(undef, '^PERL_SHELL>');
        $self->node_print("use Globals qw(\$verbosity \$TAG);\n");
        $self->node_expect(undef, '^PERL_SHELL>');
        $self->node_print("\$verbosity=$verbosity; \$TAG=\"$TAG\";\n");
        $self->node_expect(undef, '^PERL_SHELL>');

        # get node to run custom initialization
        my $config_file = $self->get_config_file();
        $self->node_print("\$scoring_ref = $scoring->new({".
            "config_file => \"$config_file\", ".
            "node_ID => $node_ID, ".
            "work_dir => \"$work_dir\", ".
            (defined $local_dir ? "local_dir => \"$local_dir\", " : "").
            "});\n");
        $self->node_expect(undef, '^PERL_SHELL>');

        # this signals multitasker to set ready flag
        $self->node_print("NODE_READY");
    }

    #--------------------------------------------------------------------------------------
    # Function: tickle
    # Synopsys: Pull a print/expect command off the queue, and process it.
    #           Back-to-back "print" commands are processed immediately by echoing
    #           them to the terminal.
    #           For an "expect" command, the terminal is polled to see whether the
    #           expected pattern has been output.  If so, it will continue processing
    #           back-to-back "print" and "expect" commands until it encounters a
    #           timeout condition on one of the "expect" commands.
    #           If a "print" command with the argument NODE_READY is issued, the node's
    #           status is set to READY(1).  NODE_READY is queued by client as a final
    #           command to signal when processing is finished, since it causes the node
    #           to become READY again -- i.e. it can be used for another job.
    #           The node status is always returned.
    #--------------------------------------------------------------------------------------
    sub tickle {
        my $self = shift; my $obj_ID = ident $self;

        my $shell_ref = $self->get_shell_ref();
        my $node_ID = $self->get_node_ID();

        printn "ScorNode::tickle: running on node $node_ID" if $verbosity >=3;

        TICKLE : while (@{$node_IO[$node_ID]->{command_queue}}) {
            $node_IO_handle->shlock();
            my $command_ref = shift @{$node_IO[$node_ID]->{command_queue}};
            $node_IO_handle->shunlock();

            printn "ScorNode::tickle: node $node_ID command is: ".join(",",@$command_ref) if $verbosity >=3;

            if ($command_ref->[0] eq "print") {
                if ($command_ref->[1] ne "NODE_READY") {
                    $shell_ref->print($command_ref->[1]);
                } else {
                    if (@{$node_IO[$node_ID]->{command_queue}} != 0) {
                        confess "ERROR: node $node_ID -- received ready command on non-empty queue";
                    }
                    printn "ScorNode::tickle: setting node $node_ID status to 1" if $verbosity >=3;
                    $self->set_node_status(1);
                }
                # keep control, processing all pending back-to-back print commands
                next TICKLE;
            }

            if ($command_ref->[0] eq "expect") {
                my $timeout = $command_ref->[1];
                my $regexp = $command_ref->[2];
                my $next_timeout;
                if (defined $timeout) {
                    $next_timeout = $timeout - 1;
                } else {
                    $next_timeout = undef;
                }

                # always use timeout of 0 -- i.e. just polling
                $timeout = 0;

                my ($matched_pattern_position,
                    $error,
                    $expect_match,
                    $expect_before,
                    $expect_after) = $shell_ref->expect($timeout, -re => $regexp);
                #	    printn "EXPECT DEBUG: node_ID=$node_ID timeout=$timeout regexp=$regexp position=$matched_pattern_position error=$error";
                #	    printn "EXPECT DEBUG: before=$expect_before";
                #	    printn "EXPECT DEBUG: match=$expect_match";
                #	    printn "EXPECT DEBUG: after=$expect_after;
                if (defined $error && $error !~ /TIMEOUT/) {
                    printn "\nERROR: tickle (node_ID=$node_ID) -- Expect signalled an error condition (e.g. EOF or process died) errno=".$shell_ref->error();
                    printn "\nERROR: tickle (node_ID=$node_ID) -- expect.before: ".$shell_ref->before();
                    printn "\nERROR: tickle (node_ID=$node_ID) -- expect.regexp:  $regexp";
                };

                # last_problem stores the location of the last problem found
                # in the current (and possibly growing) expect_before string
                # ... this allows us to prevent multiple reporting of the same problem
                my $old_last_problem = $last_problem_of{$obj_ID};
                if ($expect_before =~ /(^.*?(ERROR).*?$)/mi) {
                    my $length_prematch = length $PREMATCH;
                    if ($length_prematch > $old_last_problem) {
                        $last_problem_of{$obj_ID} = $length_prematch if $length_prematch > $last_problem_of{$obj_ID};
                        print "\nWARNING: tickle (matched $2) -- node $node_ID appears to have a problem";
                        printn " -- expect.before:\n\t==> " . $1;
                        my $next_line = $POSTMATCH;  # get POSTMATCH from regexp and extract first non-empty line
                        my @next_lines = split /[\n\r\l]+/, $next_line;
                        @next_lines = grep {$_} @next_lines;  # remove empty lines
                        printn "\t==> $next_lines[0]";
                    }
                }
                if ($expect_before =~ /(^.*?(WARNING).*?$)/mi) {
                    my $length_prematch = length $PREMATCH;
                    if ($length_prematch > $old_last_problem) {
                        $last_problem_of{$obj_ID} = $length_prematch if $length_prematch > $last_problem_of{$obj_ID};
                        print "\nWARNING: tickle (matched $2) -- node $node_ID appears to have a problem";
                        printn " -- expect.before:\n\t==> " . $1;
                        my $next_line = $POSTMATCH;  # get POSTMATCH from regexp and extract first non-empty line
                        my @next_lines = split /[\n\r\l]+/, $next_line;
                        @next_lines = grep {$_} @next_lines;  # remove empty lines
                        printn "\t==> $next_lines[0]";
                    }
                }
                if ($expect_before =~ /(^.*?(Undefined).*?$)/m) {
                    my $length_prematch = length $PREMATCH;
                    if ($length_prematch > $old_last_problem) {
                        $last_problem_of{$obj_ID} = $length_prematch if $length_prematch > $last_problem_of{$obj_ID};
                        print "\nWARNING: tickle (matched $2) -- node $node_ID appears to have a problem";
                        printn " -- expect.before:\n\t==> " . $1;
                        my $next_line = $POSTMATCH;  # get POSTMATCH from regexp and extract first non-empty line
                        my @next_lines = split /[\n\r\l]+/, $next_line;
                        @next_lines = grep {$_} @next_lines;  # remove empty lines
                        printn "\t==> $next_lines[0]";
                    }
                }
                if ($expect_before =~ /(^.*? (line) \d.*?$)/m) {
                    my $length_prematch = length $PREMATCH;
                    if ($length_prematch > $old_last_problem) {
                        $last_problem_of{$obj_ID} = $length_prematch if $length_prematch > $last_problem_of{$obj_ID};
                        print "\nWARNING: tickle (matched $2) -- node $node_ID appears to have a problem";
                        printn " -- expect.before:\n\t==> " . $1;
                        my $next_line = $POSTMATCH;  # get POSTMATCH from regexp and extract first non-empty line
                        my @next_lines = split /[\n\r\l]+/, $next_line;
                        @next_lines = grep {$_} @next_lines;  # remove empty lines
                        printn "\t==> $next_lines[0]";
                    }
                }
                if ($expect_before =~ /(^.*?(Couldn't open).*?$)/m) {
                    my $length_prematch = length $PREMATCH;
                    if ($length_prematch > $old_last_problem) {
                        $last_problem_of{$obj_ID} = $length_prematch if $length_prematch > $last_problem_of{$obj_ID};
                        print "\nWARNING: tickle (matched $2) -- node $node_ID appears to have a problem";
                        printn " -- expect.before:\n\t==> " . $1;
                        my $next_line = $POSTMATCH;  # get POSTMATCH from regexp and extract first non-empty line
                        my @next_lines = split /[\n\r\l]+/, $next_line;
                        @next_lines = grep {$_} @next_lines;  # remove empty lines
                        printn "\t==> $next_lines[0]";
                    }
                }

                if (!defined $matched_pattern_position) {
                    # didn't match anything, so put command back to front of queue if timeout not reached
                    if (!defined $next_timeout || $next_timeout > 0) {
                        $node_IO_handle->shlock();
                        unshift @{$node_IO[$node_ID]->{command_queue}}, ["expect", $next_timeout, $regexp];
                        $node_IO_handle->shunlock();
                    }
                    # only get one chance, return control
                    last TICKLE;
                } else {
                    # match was successful, so the expect_before string will be cleared
                    $last_problem_of{$obj_ID} = 0;
                    # we got something, continue processing since we didn't wait for it
#		    next TICKLE;
                    last TICKLE;   # why last not next???
                }
            }
            printn "ERROR: ScorNode::tickle -- unknown command ".join ",", @$command_ref;
            exit(1);
        }
        return $self->get_node_status();
    }
}

sub run_testcases {
    # have to set the tag to something
    $TAG = "ScorNode";

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

    burp_file("test/modules/ScorNode.cfg", $config_file);

    my $node_ref = ScorNode->new({
            node_ID => 0,
            cluster_type => "LOCAL",
            nice => "10",
            work_dir => "test/modules",
            local_dir => "test/modules/localdir",
            logfile => "ScorNode.0.log",
            config_file => "test/modules/ScorNode.cfg",
        });
    $node_ref->spawn();
    $node_ref->init("Scoring");
    while(!$node_ref->tickle()) {
        sleep 1;
    }
    $node_ref = undef;
    printn;
    printn "LOGFILE test/modules/ScorNode.0.log:";
    system("cat test/modules/ScorNode.0.log");
}

# Package BEGIN must return true value
return 1;

