#-#####################################################################################
#- File:     Scoring.pm
#- Synopsys: Base package for user-defined Scoring classes.
#-           See class ScoringFunctionTemplate for an example of how to implement
#-           an application-specific Scoring class.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Scoring;
use Class::Std::Storable;
use base qw();
{
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG $WORKSPACE);

    use MatlabDriver;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %node_ID_of     :ATTR(get => 'node_ID', set => 'node_ID', init_arg => 'node_ID');
    my %config_file_of :ATTR(get => 'config_file', set => 'config_file', init_arg => 'config_file');
    my %config_ref_of  :ATTR(get => 'config_ref', set => 'config_ref');
    my %work_dir_of        :ATTR(get => 'work_dir', set => 'work_dir', init_arg => 'work_dir');
    my %local_dir_of        :ATTR(get => 'local_dir', set => 'local_dir'); # custom init_arg
    my %matlab_work_of :ATTR(get => 'matlab_work', set => 'matlab_work');
    my %logfile_of         :ATTR(get => 'logfile', set => 'logfile');
    my %matlab_ref_of  :ATTR(get => 'matlab_ref', set => 'matlab_ref');
    my %anc_ref_of     :ATTR(get => 'anc_ref', set => 'anc_ref');

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: matlab_round_value (class or instance method)
    # Synopsys: Applies optional AbsTol/RelTol to a variable's value.
    #           resulting value is rounded appropriately given RelTol and also
    #           a cutoff of AbsTol is applied below which the value is forced to 0.0.
    #--------------------------------------------------------------------------------------
    sub matlab_round_value {
        my $class = shift;
        my %args = (
            value => undef,
            AbsTol => -1,
            RelTol => -1,
            @_,
        );
        check_args(\%args, 3);

        my $value = $args{value};
        my $AbsTol = $args{AbsTol};
        my $RelTol = $args{RelTol};

        $value = ($value < $AbsTol) ? 0.0 : $value       if ($AbsTol != -1); # apply AbsTol
        $value = round2sig($value, 1-log_10($RelTol))    if ($RelTol != -1); # apply RelTol, e.g. RelTol=1e-3 yields 4 sig. digits

        return $value;
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

        $config_ref_of{$obj_ID} = {};
        $anc_ref_of{$obj_ID} = {};
    }

    #--------------------------------------------------------------------------------------
    # Function: START
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub START {
        my ($self, $obj_ID, $arg_ref) = @_;

        # get configuration
        read_config($config_ref_of{$obj_ID}, $config_file_of{$obj_ID}, "NOCLOBBER");

        my $node_ID = $node_ID_of{$obj_ID};
        my $work_dir = $work_dir_of{$obj_ID};
        my $local_dir = $local_dir_of{$obj_ID};
        my $matlab_work = $matlab_work_of{$obj_ID} = defined $local_dir ? "$local_dir/matlab" : "$work_dir/matlab";
        system("mkdir -p $work_dir/matlab");
        system("mkdir -p $local_dir/matlab") if defined $local_dir;
        my $logfile = $logfile_of{$obj_ID} = "matlab.$node_ID.log";
        my $matlab_ref = $matlab_ref_of{$obj_ID} = MatlabDriver->new({
                name => "matlab($node_ID)",
                logfile => "$matlab_work/$logfile",
                host => "localhost",
                vmem => $config_ref_of{$obj_ID}->{vmem},
                echo => 0,
                options => $arg_ref->{matlab_startup_options} || undef,
            });
        # ADD CUSTOM DIRECTORY TO PATH
        $matlab_ref->cmd("path(path,'$WORKSPACE/custom')");

        # CHANGE WORKING DIR
        $matlab_ref->cmd("cd $matlab_work; format long;");

        # check initializers
        # ...
    }

    #--------------------------------------------------------------------------------------
    # Function: CTRL-C handler
    # Synopsys: Trap CTRL-C and exit cleanly.
    #           This is necessary to ensure DEMOLISH is called when a CTRL-C is received.
    #--------------------------------------------------------------------------------------
    $SIG{INT} = sub {		# trapping CTRL-C ensures graceful exit via DEMOLISH/END BLOCKS etc.
        printn "Scoring: process $$ exiting from CTRL-C\n";
        exit;
    };

    #--------------------------------------------------------------------------------------
    # Function: DEMOLISH
    # Synopsys: Move logfiles from local_dir to work_dir.
    #--------------------------------------------------------------------------------------
    sub DEMOLISH {
        my $self = shift; my $obj_ID = shift;

        printn "Scoring::DEMOLISH: called";

        my $work_dir = $work_dir_of{$obj_ID};
        my $local_dir = $local_dir_of{$obj_ID};
        my $matlab_work = $matlab_work_of{$obj_ID};
        my $logfile = $logfile_of{$obj_ID};
        $matlab_ref_of{$obj_ID} = undef; # shut down matlab
        if (defined $local_dir) {
            printn "Moving $matlab_work/$logfile to $work_dir/matlab";
            system("mv $matlab_work/$logfile $work_dir/matlab");
            printn "Done moving.";
        }

        printn "Scoring::DEMOLISH: done";
    }

    #--------------------------------------------------------------------------------------
    # Function: anc_process_species_report
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub anc_process_species_report {
        my $self = shift; my $obj_ID = ident $self;
        my $filename = shift;

        $anc_ref_of{$obj_ID} = {};

        my @file = slurp_file($filename);
        my $stats_line = shift @file;
        $stats_line =~ /(\S+)\s+(\S+)/;
        my $num_complexes = $anc_ref_of{$obj_ID}->{num_complexes} = $1;
        my $num_species = $anc_ref_of{$obj_ID}->{num_species} = $2;

        @{$anc_ref_of{$obj_ID}->{species}} = ();
        $anc_ref_of{$obj_ID}->{complexes} = {};

        foreach my $line (@file) {
            my @split_line = split(/\s+/, $line);
            push @{$anc_ref_of{$obj_ID}->{species}}, @split_line[1..$#split_line];
            $anc_ref_of{$obj_ID}->{complexes}{$split_line[0]} = [@split_line[1..$#split_line]];
        }

        my $species_report_status = 0;
        if ($num_complexes != keys %{$anc_ref_of{$obj_ID}->{complexes}}) {
            printn "ERROR: ANC species report is messed up (1)";
            $species_report_status++;
        }
        if ($num_species != @{$anc_ref_of{$obj_ID}->{species}}) {
            $species_report_status++;
            printn "ERROR: ANC species report is messed up (2)";
        }
        return $species_report_status;
    }

    #--------------------------------------------------------------------------------------
    # Function: anc_get_species
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub anc_get_species {
        my $self = shift; my $obj_ID = ident $self;
        return @{$anc_ref_of{$obj_ID}->{species}};
    }

    #--------------------------------------------------------------------------------------
    # Function: facile_run
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub facile_run {
        my $self = shift; my $obj_ID = ident $self;

        my %args = (
            # default values
            EQN_FILE    => undef,
            SIM_TYPE    => "matlab",
            SOLVER      => "UNDEF",
            SOLVER_OPTIONS => "UNDEF",
            T_FINAL     => "UNDEF",
            T_SAMPLING    => "UNDEF",
            T_TICK      => "UNDEF",
            T_EVENTS    => "UNDEF",
            SPLIT       => 0,
            VERBOSE     => 0,
            @_,			# argument pair list overwrites defaults
        );
        check_args(\%args, 10);

        my $eqn_file = $args{EQN_FILE};
        my $sim_type = $args{SIM_TYPE};
        my $solver   = $args{SOLVER};
        my $solver_options   = $args{SOLVER_OPTIONS};
        my $tf       = $args{T_FINAL};
        my $t_sampling = $args{T_SAMPLING};
        my $tk       = $args{T_TICK};
        my $split_flag = $args{SPLIT};

        my $file_root = $eqn_file; $file_root =~ s/\..*//;

        # 	 open(EQN_INI, ">$work_dir/$rootname.eqn.ini") or die "ERROR: facile_run -- Couldn't open $work_dir/$rootname.eqn.ini for writing.\n";
        # 	 $| = 1;
        # 	 print EQN_INI "INIT:\n";
        # 	 my ($protein, $state, $concentration, $units);
        # 	 $units = "mol/L";
        # 	 my $concentration_Nmolecules;
        # 	 foreach $protein (sort keys %{$sdb->{protein_table}}) {
        # 	     if (defined $sdb->{protein_table}{$protein}{input_concentration}) {
        # 		 $concentration = $sdb->{protein_table}{$protein}{input_concentration};  # input conc. in mol/L (no scaling)
        # 		 $concentration_Nmolecules = $concentration * $cell_volume * $avogadro;
        # 		 $state = $sdb->{complex_table}{$protein}{state_list}->[0];
        # 		 printn "facile_run: initial input level for protein ${protein}_$state is $concentration $units";
        # 		 print EQN_INI "${protein}_$state = $concentration;  # ".sprintf("%.2f",$concentration_Nmolecules)." N\n";
        # 	     }
        # 	     elsif (defined $sdb->{protein_table}{$protein}{regulated_concentration}) {
        # 		 $concentration = $sdb->{protein_table}{$protein}{regulated_concentration} * $concentration_scaling_factor;  # regulated conc. is integer and needs scaling
        # 		 $concentration_Nmolecules = $concentration * $cell_volume * $avogadro;
        # 		 $state = $sdb->{complex_table}{$protein}{state_list}->[0];
        # 		 printn "facile_run: initial regulated level for protein ${protein}_$state is $concentration $units";
        # 		 print EQN_INI "${protein}_$state = $concentration;  # ".sprintf("%.2f",$concentration_Nmolecules)." N\n";
        # 	     } else {
        # 		 printn "ERROR: facile_run -- no initial concentration for $protein";
        # 		 exit;
        # 	     }
        # 	 }
        # 	 close(EQN_INI) or die "Couldn't close $work_dir/$rootname.eqn.ini\n";
        #     `cat $eqn_file $work_dir/$rootname.eqn.ini > $work_dir/$rootname.eqn.all`;
        #     `rm -f $eqn_file $work_dir/$rootname.eqn.ini`;   # remove the unused files


        #    my $facile_cmd = "$ENV{FACILE_HOME}/facile.pl -q ";
        my $facile_cmd = "$ENV{FACILE_HOME}/facile.pl ". ($args{VERBOSE} ? "--verbose " : "");
        if ($sim_type =~ "matlab") {
            $facile_cmd .= " --matlab";
            $facile_cmd .= " --solver $solver" if $solver ne "UNDEF";
            $facile_cmd .= " --events \"$args{T_EVENTS}\"" if ($args{T_EVENTS} ne "UNDEF");
            $facile_cmd .= " --solver_options \"$solver_options\"" if ($solver_options ne "UNDEF");
            $facile_cmd .= " --t_tick $tk " if $tk ne "UNDEF";
            $facile_cmd .= " --t_final $tf " if $tf ne "UNDEF";
            $facile_cmd .= " --t_sampling \"$t_sampling\" " if $t_sampling ne "UNDEF";
            $facile_cmd .= " --split " if $split_flag;
            $facile_cmd .= " $eqn_file";
        } elsif ($sim_type eq "easystoch") {
            # 	     $facile_cmd .= " --easystoch";
            # 	     if ($args{T_EVENTS} ne "") {
            # 		 $facile_cmd .= " --events \"$args{T_EVENTS}\"";
            # 	     }
            # 	     $facile_cmd .= " -C $cell_volume $eqn_file";
        } else {
            printn "ERROR: unsupported sim_type \"$sim_type\"";
            exit;
        }

        printn "facile_run: facile command is $facile_cmd";
        printn "facile_run: started facile on " . `date`;
        system("$facile_cmd");
        system("echo ${file_root}_r.m   |   xargs rm -f");
        system("echo ${file_root}_s.m   |   xargs rm -f");
        if ($?) {
            printn "ERROR: Facile reported an error ($?)";
            exit(1);
        }
        printn "facile_run: finished facile on " . `date`;
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_cmd
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub matlab_cmd {
        my $self = shift; my $obj_ID = ident $self;

        my $matlab_ref = $matlab_ref_of{$obj_ID};
        return $matlab_ref->cmd(@_);
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_wait_on
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub matlab_wait_on {
        my $self = shift; my $obj_ID = ident $self;

        my $matlab_ref = $matlab_ref_of{$obj_ID};
        return $matlab_ref->wait_on(@_);
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_variable
    # Synopsys: Returns the elements of an NxN variable as a list;
    #--------------------------------------------------------------------------------------
    sub matlab_get_variable {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            name => undef,
            @_,
        );
        check_args(\%args, 1);
        my $name = $args{name};

        my $matlab_ref = $matlab_ref_of{$obj_ID};

        printn "matlab_get_variable: getting value of $name" if $verbosity >= 3;
        $matlab_ref->cmd("fv_exists = exist('$name');");
        $matlab_ref->cmd("if (fv_exists)");
        $matlab_ref->cmd("str1=sprintf('GET VALUE $name :');\n");
        $matlab_ref->cmd("str2=sprintf(' %.15e', $name);\n");
        $matlab_ref->cmd("str=[str1, str2];");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("else");
        $matlab_ref->cmd("str=sprintf('GET VALUE $name : UNDEFINED');\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("end");
        # Wait for matlab to be done
        my ($line) = $matlab_ref->wait_on("GET VALUE $name");
        $line =~ s/.*GET VALUE $name :\s*//;
        my @values = split /\s+/, $line;

        return @values;
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_state
    # Synopsys: Returns the concentration at time t of a state variable.
    #           Finds closest prior value if exact time point is not in t vector.
    #--------------------------------------------------------------------------------------
    sub matlab_get_state {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            complex => undef,
            t => undef,
            @_,
        );
        check_args(\%args, 2);
        my $complex = $args{complex};
        my $t = $args{t};

        my $matlab_ref = $matlab_ref_of{$obj_ID};

        printn "matlab_get_state: $complex at t=$t" if $verbosity >= 3;
        $matlab_ref->cmd("fv_exists = exist('$complex');");
        $matlab_ref->cmd("if (fv_exists)");
        $matlab_ref->cmd("I = find(t<=$t,1,'last');\n");
        $matlab_ref->cmd("str=sprintf('GET VALUE $complex : %.15e', $complex(I));\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("else");
        $matlab_ref->cmd("str=sprintf('GET VALUE $complex : UNDEFINED');\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("end");
        # Wait for matlab to be done
        my ($line) = $matlab_ref->wait_on("GET VALUE $complex");
        $line =~ /.*\s(\S+)/;
        my $value = $1;

        return $value;
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_state_range
    # Synopsys: Returns the concentration range in time interval [t1,t2]
    #           of a state variable.
    #           Finds closest prior value if exact time point is not in t vector.
    #--------------------------------------------------------------------------------------
    sub matlab_get_state_range {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            complex => undef,
            t1 => undef,
            t2 => undef,
            @_,
        );
        check_args(\%args, 3);
        my $complex = $args{complex};
        my $t1 = $args{t1};
        my $t2 = $args{t2};

        my $matlab_ref = $matlab_ref_of{$obj_ID};

        printn "matlab_get_state_range: $complex from t1=$t1 to t2=$t2" if $verbosity >= 3;
        $matlab_ref->cmd("fv_exists = exist('$complex');");
        $matlab_ref->cmd("if (fv_exists)");
        $matlab_ref->cmd("I1 = find(t<=$t1,1,'last');\n");
        $matlab_ref->cmd("I2 = find(t<=$t2,1,'last');\n");
        $matlab_ref->cmd("range_max = max($complex(I1:I2));\n");
        $matlab_ref->cmd("range_min = min($complex(I1:I2));\n");
        $matlab_ref->cmd("str=sprintf('GET RANGE MIN $complex : %.15e', range_min);\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("str=sprintf('GET RANGE MAX $complex : %.15e', range_max);\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("else");
        $matlab_ref->cmd("str=sprintf('GET RANGE MIN $complex : UNDEFINED');\n");
        $matlab_ref->cmd("str=sprintf('GET RANGE MAX $complex : UNDEFINED');\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("end");
        # Wait for matlab to be done
        my ($min_line) = $matlab_ref->wait_on("GET RANGE MIN $complex");
        $min_line =~ /.*\s(\S+)/;
        my $range_min = $1;
        my ($max_line) = $matlab_ref->wait_on("GET RANGE MAX $complex");
        $max_line =~ /.*\s(\S+)/;
        my $range_max = $1;

        return ($range_min, $range_max);
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_state_vector
    # Synopsys: Returns concentration state vector at time t
    #           Finds closest prior value if exact time point is not in t vector.
    #--------------------------------------------------------------------------------------
    sub matlab_get_state_vector {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            t => undef,
            state_var => "y",
            @_,
        );
        check_args(\%args, 2);

        my $t = $args{t};
        my $state_var = $args{state_var};

        my $matlab_ref = $matlab_ref_of{$obj_ID};

        printn "matlab_get_state_vector: state vector $state_var at t=$t" if $verbosity >= 3;
        $matlab_ref->cmd("fv_exists = exist('$state_var');");
        $matlab_ref->cmd("if (fv_exists)");
        $matlab_ref->cmd("I = find(t<=$t,1,'last');");
        $matlab_ref->cmd("y_t=$state_var(I,:);");

        $matlab_ref->cmd("str1=sprintf('GET STATE at t=$t :');");
        $matlab_ref->cmd("str2=sprintf(' %.15e ', y_t);");
        $matlab_ref->cmd("str=[str1, str2];");
        $matlab_ref->cmd("disp(str)");
        $matlab_ref->cmd("else");
        $matlab_ref->cmd("str=sprintf('ERROR: cannot find state variable $state_var !!! ');");
        $matlab_ref->cmd("disp(str)");
        $matlab_ref->cmd("end");

        # Wait for matlab to be done
        my ($line) = $matlab_ref->wait_on("GET STATE at t=$t");

        $line =~ s/.*GET STATE at t=$t :\s*//;
        my @state_vector = split /\s+/, $line;

        if ($verbosity >= 3) {
            for (my $i = 0; $i < @state_vector; $i++) {
                printn "matlab_get_state_vector: $state_var(" . ($i+1) . ") = $state_vector[$i]";
            }
        }
        return @state_vector;
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_state_vector_range
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub matlab_get_state_vector_range {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            t1 => undef,
            t2 => undef,
            state_var => "y",
            @_,
        );
        check_args(\%args, 3);

        my $t1 = $args{t1};
        my $t2 = $args{t2};
        my $state_var = $args{state_var};

        my $matlab_ref = $matlab_ref_of{$obj_ID};

        printn "matlab_get_state_vector_range: state vector $state_var from t1=$t1 t2=$t2" if $verbosity >= 3;
        $matlab_ref->cmd("fv_exists = exist('$state_var');");
        $matlab_ref->cmd("if (fv_exists)");
        $matlab_ref->cmd("I1 = find(t<=$t1,1,'last')");
        $matlab_ref->cmd("I2 = find(t<=$t2,1,'last')");
        $matlab_ref->cmd("y_range=$state_var(I1:I2,:);");
        $matlab_ref->cmd("size(y)");
        $matlab_ref->cmd("size(y_range)");
        $matlab_ref->cmd("range_min = min(y_range);\n");
        $matlab_ref->cmd("range_max = max(y_range);\n");
        $matlab_ref->cmd("str1=sprintf('GET RANGE MIN : ');");
        $matlab_ref->cmd("str2=sprintf('%.15e ', range_min);");
        $matlab_ref->cmd("str=[str1, str2];");
        $matlab_ref->cmd("disp(str)");
        $matlab_ref->cmd("str1=sprintf('GET RANGE MAX : ');");
        $matlab_ref->cmd("str2=sprintf('%.15e ', range_max);");
        $matlab_ref->cmd("str=[str1, str2];");
        $matlab_ref->cmd("disp(str)");
        $matlab_ref->cmd("else");
        $matlab_ref->cmd("str=sprintf('ERROR: cannot find state variable $state_var !!! ');");
        $matlab_ref->cmd("disp(str)");
        $matlab_ref->cmd("end");

        # Wait for matlab to be done
        my ($min_line) = $matlab_ref->wait_on("GET RANGE MIN :");
        $min_line =~ s/.*GET RANGE MIN :\s*//;
        my @range_min = split " ", $min_line;
        my ($max_line) = $matlab_ref->wait_on("GET RANGE MAX :");
        $max_line =~ s/.*GET RANGE MAX :\s*//;
        my @range_max = split " ", $max_line;

        if ($verbosity >= 3) {
            for (my $i = 0; $i < @range_min; $i++) {
                printn "matlab_get_state_vector_range: min $state_var(" . ($i+1) . ",:) = $range_min[$i]";
                printn "matlab_get_state_vector_range: max $state_var(" . ($i+1) . ",:) = $range_max[$i]";
            }
        }

        return {
            range_min => \@range_min,
            range_max => \@range_max,
        };
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_state_delta
    # Synopsys: Returns state vector delta y(t2) - y(t1), and state change relative to
    #           dynamic range from t0 to t2.  Will return relative delta of 0
    #           if the average value of a state variable is 0.0.
    #           (Note: uses closest prior time value for t if exact time point
    #            is not in t vector)
    #           Will perform appropriate rounding if passed Abs/RelTol.
    #--------------------------------------------------------------------------------------
    sub matlab_get_state_delta {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            t0 => 0,
            t1 => undef,
            t2 => undef,
            state_var => "y",
            AbsTol => -1,
            RelTol => -1,
            @_,
        );
        check_args(\%args, 6);

        my $t0 = $args{t0};
        my $t1 = $args{t1};
        my $t2 = $args{t2};
        my $state_var = $args{state_var};

        printn "matlab_get_state_delta: state delta y @ dt = ($t2 - $t1)" if $verbosity >= 3;

        my @state_vector_t1 = $self->matlab_get_state_vector(t => $t1, state_var => $state_var);
        my @state_vector_t2 = $self->matlab_get_state_vector(t => $t2, state_var => $state_var);

        @state_vector_t1 = map {$self->matlab_round_value(value => $_,
        AbsTol => $args{AbsTol},
        RelTol => $args{RelTol})} @state_vector_t1;
        @state_vector_t2 = map {$self->matlab_round_value(value => $_,
        AbsTol => $args{AbsTol},
        RelTol => $args{RelTol})} @state_vector_t2;

        my $range_ref = $self->matlab_get_state_vector_range(
            state_var => $state_var,
            t1 => $t0,
            t2 => $t2,
        );

        my @range_min = @{$range_ref->{range_min}};
        my @range_max = @{$range_ref->{range_max}};
        my @state_vector_dynamic_range = map {$range_max[$_] - $range_min[$_]} (0..$#range_min);

        if (@state_vector_t1 != @state_vector_t2) {
            confess "ERROR: matlab_get_state_delta -- inconsistent state vector sizes";
        }

        my @state_vector_delta = ();
        my @state_vector_relative_delta = ();
        for (my $i = 0; $i < @state_vector_t1; $i++) {
            $state_vector_delta[$i] = $state_vector_t2[$i] - $state_vector_t1[$i];
            my $dynamic_range = $state_vector_dynamic_range[$i];
            $state_vector_relative_delta[$i] = ($dynamic_range == 0.0) ? 0.0 : ($state_vector_delta[$i] / $dynamic_range);
            if ($verbosity >= 3) {
                printn "matlab_get_state_delta: y_t1(" . ($i+1) . ") = $state_vector_t1[$i]";
                printn "matlab_get_state_delta: y_t2(" . ($i+1) . ") = $state_vector_t2[$i]";
                printn "matlab_get_state_delta: delta_y(" . ($i+1) . ") = $state_vector_delta[$i]";
                printn "matlab_get_state_delta: range_y(" . ($i+1) . ") = $state_vector_dynamic_range[$i]";
                printn "matlab_get_state_delta: relative_delta_y(" . ($i+1) . ") = $state_vector_relative_delta[$i]";
            }
        }
        return {
            delta => \@state_vector_delta,
            relative_delta => \@state_vector_relative_delta,
            dynamic_range => \@state_vector_dynamic_range,
            state_vector_t1 => \@state_vector_t1,
            state_vector_t2 => \@state_vector_t2,
        };
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_max_value
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub matlab_get_max_value {
        my $self = shift; my $obj_ID = ident $self;
        my $complex = shift;

        my $matlab_ref = $matlab_ref_of{$obj_ID};

        printn "matlab_get_max_value: $complex" if $verbosity >= 3;
        $matlab_ref->cmd("fv_exists = exist('$complex');");
        $matlab_ref->cmd("if (fv_exists)");
        $matlab_ref->cmd("str=sprintf('GET MAX VALUE $complex : %.15e', max($complex));\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("else");
        $matlab_ref->cmd("str=sprintf('GET MAX VALUE $complex : UNDEFINED');\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("end");
        # Wait for matlab to be done
        my ($line) = $matlab_ref->wait_on("GET MAX VALUE $complex");
        $line =~ /.*\s(\S+)/;
        my $value = $1;
        printn "matlab_get_max_value: $complex MAX VALUE = $value" if $verbosity >= 3;
        return $value;
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_report_max_values
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub matlab_report_max_values {
        my $self = shift; my $obj_ID = ident $self;

        my @anc_species = $self->anc_get_species();

        printn "matlab_report_max_values: getting max values from matlab...";
        foreach my $species_name (sort @anc_species) {
            printn "matlab_report_max_values: MAX value for $species_name is " . $self->matlab_get_max_value($species_name);
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_get_final_value
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub matlab_get_final_value {
        my $self = shift; my $obj_ID = ident $self;
        my $complex = shift;

        my $matlab_ref = $matlab_ref_of{$obj_ID};

        printn "matlab_get_final_value: for $complex" if $verbosity >= 3;
        $matlab_ref->cmd("fv_exists = exist('$complex');");
        $matlab_ref->cmd("if (fv_exists)");
        $matlab_ref->cmd("str=sprintf('FINAL VALUE $complex : %.15e', $complex(end));\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("else");
        $matlab_ref->cmd("str=sprintf('FINAL VALUE $complex : UNDEFINED');\n");
        $matlab_ref->cmd("disp(str)\n");
        $matlab_ref->cmd("end");
        # Wait for matlab to be done
        my ($line) = $matlab_ref->wait_on("FINAL VALUE $complex");
        $line =~ /.*\s(\S+)/;
        my $final_value = $1;
        return $final_value;
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_report_final_values
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub matlab_report_final_values {
        my $self = shift; my $obj_ID = ident $self;

        my @anc_species = $self->anc_get_species();

        printn "matlab_report_final_values: getting final values from matlab...";
        foreach my $species_name (sort @anc_species) {
            printn "matlab_report_final_values: final value for $species_name is " . $self->matlab_get_final_value($species_name);
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_plot_complex
    # Synopsys: Plot a specific complex
    #--------------------------------------------------------------------------------------
    sub matlab_plot_complex {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            figure => undef,
            complex => undef,
            title_prefix => "",
            plot_command => "plot",
            filename => "",
            @_);
        check_args(\%args, 5);

        my $figure = $args{figure};
        my $complex = $args{complex};
        my $title_prefix = $args{title_prefix};
        my $plot_command = $args{plot_command};
        my $filename = $args{filename};

        my $title = $complex;

        # matlab treats underscores as indication of subscript, so escape them
        $title =~ s/_/\\_/g;
        $title_prefix =~ s/_/\\_/g;

        my $matlab_ref = $matlab_ref_of{$obj_ID};
        my $command = "";
        my $rescoring = 0 || $config_ref_of{$obj_ID}->{rescoring};
        if ($rescoring==1) {
            $command = "io=figure(\'Visible\',\'off\');$plot_command(t, $complex); title(\'$title_prefix $title\')";
        } else {
            $command = "io=figure($figure);$plot_command(t, $complex); title(\'$title_prefix $title\')";
        }
        printn "matlab_plot_complex: Figure $figure -- $complex" if $verbosity >= 1;
        $matlab_ref->cmd("$command");
        $matlab_ref->cmd("saveas(io, \'$filename\', \'png\')") if $filename;
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_plot_all_complexes
    # Synopsys: Plot all complexes whose max conc gte min_value.
    #           A -ve min_value means don't plot anything.
    #--------------------------------------------------------------------------------------
    sub matlab_plot_all_complexes {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            min_value => 0,
            plot_command => "plot",
            figure => 1,
            @_,
        );
        check_args(\%args, 3);

        my $min_value = $args{min_value};
        my $plot_command = $args{plot_command};
        my $figure = $args{figure};

        return if ($min_value < 0);

        my @anc_species = $self->anc_get_species();

        foreach my $species_name (sort @anc_species) {
            my $max_value = $self->matlab_get_max_value($species_name);
            if ($max_value !~ /UNDEF/ && $max_value >= $min_value) {
                $self->matlab_plot_complex(
                    figure => $figure++,
                    complex => $species_name,
                    plot_command => $plot_command,
                );
            }
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: matlab_plot_phase
    # Synopsys: Plot phase plot.
    #--------------------------------------------------------------------------------------
    sub matlab_plot_phase {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            figure => undef,
            X_complex => undef,
            Y_complex => undef,
            title_prefix => "",
            plot_command => "plot",
            filename => "",
            axis_ref => "",
            @_);
        check_args(\%args, 7);

        my $figure = $args{figure};
        my $X_complex = $args{X_complex};
        my $Y_complex = $args{Y_complex};
        my $title_prefix = $args{title_prefix};
        my $plot_command = $args{plot_command};
        my $filename = $args{filename};
        my $axis_ref = $args{axis_ref};

        my $title = "PHASE PLOT ($Y_complex vs $X_complex)";

        # matlab treats underscores as indication of subscript, so escape them
        $title =~ s/_/\\_/g;
        $title_prefix =~ s/_/\\_/g;

        my $matlab_ref = $matlab_ref_of{$obj_ID};
        $matlab_ref->cmd("halfway = floor(size($X_complex,1)/2)");
        my $rescoring = 0 || $config_ref_of{$obj_ID}->{rescoring};
        if ($rescoring==1) {
            $matlab_ref->cmd("h=figure(\'Visible\',\'off\'); $plot_command($X_complex(1:halfway), $Y_complex(1:halfway));title(\'$title_prefix $title\')");
        } else {
            $matlab_ref->cmd("h=figure($figure); $plot_command($X_complex(1:halfway), $Y_complex(1:halfway));title(\'$title_prefix $title\')");
        }  
        $matlab_ref->cmd("hold on; $plot_command($X_complex(halfway+1:end), $Y_complex(halfway+1:end), \'r\'); hold off;");
        $matlab_ref->cmd("axis([".join(" ", @$axis_ref)."])") if $axis_ref;
        #$matlab_ref->cmd("hgsave(h, \'$filename\')") if $filename;
        $matlab_ref->cmd("saveas(h, \'$filename\', \'png\')") if $filename;
        printn "matlab_plot_phase: Figure $figure -- $Y_complex vs $X_complex" if $verbosity >= 1;
    }

    #--------------------------------------------------------------------------------------
    # Function: score_genome
    # Synopsys: This routine is application-specific and should be provided by a sub-class.
    #--------------------------------------------------------------------------------------
    sub score_genome {
        my $self = shift; my $obj_ID = ident $self;
        my $genome_model_ref = shift;

        printn "ERROR: you need to define a your own application specific Scoring class";
    }
}

#--------------------------------------------------------------------------------------
# Function: run_testcases
# Synopsys: 
#--------------------------------------------------------------------------------------
sub run_testcases {

    printn "NO TESTCASES!!!!";
}


# Package BEGIN must return true value
return 1;

