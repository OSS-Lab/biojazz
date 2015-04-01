#-#####################################################################################
#- File:     Stimulus.pm
#- Synopsys: Functions for stimulating network.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-
#-#####################################################################################

#######################################################################################
# TO-DO LIST
#######################################################################################

#######################################################################################
# Package interface
#######################################################################################
package Stimulus;

use strict;

use Exporter;
use vars qw(@ISA @EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(
clamping_equation
impulse_equation
wash_equation
staircase_equation
staircase_sample
ramp_equation
ss_ramp_equation
rand_ss_ramp_equation
);

use Utils;
use Globals;

#--------------------------------------------------------------------------------------
# Function: clamping_equation
# Synopsys: Generate a source/sink equation pair to clamp a node at a given concentration
#--------------------------------------------------------------------------------------
sub clamping_equation {
    my %args = (
        # default values
        NODE => undef,
        PERIOD => undef,
        DELAY => 0,
        STRENGTH => undef,
        CONCENTRATION    => undef,
        DUTY => 100,
        @_,        # argument pair list overwrites defaults
    );

    $args{PERIOD} = -1 if ($args{DUTY} == 100);

    check_args(\%args, 6);

    my $node = $args{NODE};
    my $period = $args{PERIOD};
    my $delay = $args{DELAY};
    my $strength = $args{STRENGTH};
    my $concentration = $args{CONCENTRATION};
    my $duty = $args{DUTY};

    if ($duty < 100) {
        return ("null -> $node; clamp_source_$node=\"0.5*$concentration*$strength*(square(2*pi*(t-$delay)/$period, $duty)+1)\"", 
            "$node -> null; clamp_sink_$node=$strength");
    } else {
        return ("null -> $node; clamp_source_$node=".$concentration*$strength,
            "$node -> null; clamp_sink_$node=$strength");
    }
}

#--------------------------------------------------------------------------------------
# Function: impulse_equation
# Synopsys: 
#--------------------------------------------------------------------------------------
sub impulse_equation {
    my %args = (
        # default values
        NODE => undef,
        PERIOD => undef,
        DELAY => undef,
        CONCENTRATION => undef,
        IMPULSE_LENGTH    => undef,
        @_,        # argument pair list overwrites defaults
    );

    check_args(\%args, 5);

    my $node = $args{NODE};
    my $period = $args{PERIOD};
    my $delay = $args{DELAY};
    my $concentration = $args{CONCENTRATION};
    my $impulse_length = $args{IMPULSE_LENGTH};

    if ($delay + $impulse_length > $period) {
        printn "impulse_equation: ERROR -- impulse period wrap-around";
        exit(1);
    }

    my $level = $concentration/$impulse_length;
    my $duty = $impulse_length/$period*100;

    my $leading_edge = $delay;
    my $trailing_edge = $delay + $impulse_length;

    return {
        equations => [
            "null -> $node; impulse_source_$node=\"0.5*$level*(square(2*pi*(t-$delay)/$period, $duty)+1)\"",
        ],
        events => [$leading_edge, $trailing_edge],
        values => [],
    };
}

#--------------------------------------------------------------------------------------
# Function: wash_equation
# Synopsys: 
#--------------------------------------------------------------------------------------
sub wash_equation {
    my %args = (
        # default values
        NODE => undef,
        PERIOD => -1,
        DELAY => 0,
        STRENGTH => undef,
        WASH_LENGTH    => undef,
        SINK_NODE => "null",
        @_,        # argument pair list overwrites defaults
    );

    $args{DELAY} = -1 if ($args{PERIOD} == -1);
    $args{WASH_LENGTH} = -1 if ($args{PERIOD} == -1);

    check_args(\%args, 6);

    my $node = $args{NODE};
    my $sink_node = $args{SINK_NODE};
    my $period = $args{PERIOD};
    my $delay = $args{DELAY};
    my $strength = $args{STRENGTH};
    my $wash_length = $args{WASH_LENGTH};
    my $duty;
    if ($period != -1) {
        if ($delay + $wash_length > $period) {
            printn "wash_equation: ERROR -- wash period wrap-around";
            exit;
        }
        $duty = $wash_length/$period*100;
        return ("$node -> $sink_node; wash_sink_$node=\"0.5*$strength*(square(2*pi*(t-$delay)/$period, $duty)+1)\"");
    } else {
        return ("$node -> $sink_node; wash_sink_$node=$strength");
    }
}

#--------------------------------------------------------------------------------------
# Function: staircase_equation
# Synopsys: Generate a source/sink equation pair staircase-style ramp up and down.
#           DUTY indicates time spent at high level (incl. rise/fall times), RFTIME
#           gives the rise/fall times, and STEPS gives how many steps used when changing.
#           Also returns a list of event times where step function changes
#           (first period only), and the expected steady-state value of the function
#           right BEFORE the event.
#           The ramp_equation will touch at the bottom/top corners of the corresponding
#           staircase function when rising/falling.
#--------------------------------------------------------------------------------------
sub staircase_equation {
    my %args = (
        # default values
        NODE => undef,
        PERIOD => undef,
        DELAY => 0,
        STRENGTH => undef,
        CONCENTRATION => undef,
        DUTY => 50,
        RFTIME => undef,
        STEPS => undef,   # of steps in each direction
        @_,        # argument pair list overwrites defaults
    );

    check_args(\%args, 8);

    my $node = $args{NODE};
    my $period = $args{PERIOD};
    my $delay = $args{DELAY};
    my $strength = $args{STRENGTH};
    my $concentration = $args{CONCENTRATION};
    my $duty = $args{DUTY};
    my $rftime = $args{RFTIME};
    my $steps = $args{STEPS};

    my $pulse_width = $duty * $period / 100.0;
    my $square_width = $pulse_width - $rftime;
    my $square_duty = ($square_width) / $period * 100.0;

    if ($steps < 1) {
        printn "ERROR: staircase_equation -- step parameter must be >= 1";
        exit(1);
    }

    if ($pulse_width < 2 * $rftime) {
        printn "ERROR: staircase_equation -- parameters are inconsistent (rise/fall time exceeds computed pulse width)";
        exit(1);
    }

    if ($square_duty <= 0) {
        printn "ERROR: staircase_equation -- parameters are inconsistent (rise/fall time too large?)";
        exit(1);
    }

    if (($steps < 2) && ($rftime > 0)) {
        printn "ERROR: staircase_equation -- parameters are inconsistent (cannot implement non-zero rise/fall time with indicated # of steps)";
        exit(1);
    }

    my $square_delay = ($steps < 2) ? 0 : $rftime / ($steps);

    my @stair_source_nodes;
    my @leading_edge_events;
    my @trailing_edge_events;
    my @leading_edge_values;
    my @trailing_edge_values;
    for (my $i=0; $i < $steps; $i++) {
        my $shift = $i * $square_delay + $delay;
        push @leading_edge_events, $shift;
        push @leading_edge_values, $i / $steps * $concentration;
        push @trailing_edge_events, $shift + $square_width;
        push @trailing_edge_values, $concentration * ($steps - $i) / $steps;
        # to prevent roundoff error, multiply by pi as very last step, else
        # we may get **old** value at event time instead of the new value
        push @stair_source_nodes, "0.5*$concentration*$strength*(1 + square((t-$shift)/$period*2*pi, $square_duty))";
    }
    my $stair_source_node = "(".join(" + ", @stair_source_nodes).")/$steps";

    my @events = (@leading_edge_events, @trailing_edge_events);
    # only events in the first period are included here. also, if the delay offset
    # is too large, some events will fall outside the period, so we must use
    # a modulus operator to fix this.  we also use Utils::fmod which can handle
    # fractional dividends
    @events = map {fmod($_,$period)} @events;
    my @order = sort {$events[$a] <=> $events[$b]} (0..$#events);
    @order = grep {$events[$_] != 0} @order;   # don't include 0.0 as an event
    @events = map {$events[$_]} @order;
    my @values = (@leading_edge_values, @trailing_edge_values);
    @values = map {$values[$_]} @order;

    printn "staircase_equation: event list is @events" if $verbosity >= 3;
    printn "staircase_equation: value list is @values" if $verbosity >= 3;

    if ($duty < 100) {
        # final value is value at last event +/- one step
        my $final_value = $values[@values-1];
        $final_value += (
            fmod(($period - $delay),$period) <= ($pulse_width - $rftime) &&
            fmod(($period - $delay), $period) > 0 ?
            ($concentration/$steps) :
            (-$concentration/$steps)
        );
        return {
            equations => [
                "null -> $node; clamp_source_$node=\"$stair_source_node\"",
                "$node -> null; clamp_sink_$node=$strength",
            ],
            events => \@events,
            values => \@values,
            final_value => $final_value,
        };
    } else {
        return {
            equations => [
                "null -> $node; clamp_source_$node=".$concentration*$strength,
                "$node -> null; clamp_sink_$node=$strength",
            ],
            events => [],
            values => [],
            final_value => $concentration,
        };
    }
}

#--------------------------------------------------------------------------------------
# Function: staircase_sample
# Synopsys: Call staircase_equation with dummy values
#--------------------------------------------------------------------------------------
sub staircase_sample {
    my %args = (
        # default values
        CONCENTRATION => undef,
        PERIOD => undef,
        DELAY => 0,
        DUTY => 50,
        RFTIME => undef,
        STEPS => undef,   # of steps in each direction
        @_,        # argument pair list overwrites defaults
    );

    check_args(\%args, 6);

    my $staircase_ref = staircase_equation(
        NODE => "DUMMY",
        STRENGTH => 1.0,
        %args,
    );
    return {
        events => $staircase_ref->{events},
        values => $staircase_ref->{values},
        final_value => $staircase_ref->{final_value},
    }
}

#--------------------------------------------------------------------------------------
# Function: ramp_equation
# Synopsys: Generate a source/sink equation pair to generate linear ramp up and down.
#           Arguments and return values are identical to staircase_equation
#           for compatibility, even though some arguments are not applicable.
#--------------------------------------------------------------------------------------
sub ramp_equation {
    my %args = (
        # default values
        NODE => undef,
        PERIOD => undef,
        DELAY => 0,
        STRENGTH => undef,
        CONCENTRATION    => undef,
        DUTY => 50,
        RFTIME => undef,
        STEPS => undef,   # N/A
        @_,               # argument pair list overwrites defaults
    );

    check_args(\%args, 8);

    my $node = $args{NODE};
    my $period = $args{PERIOD};
    my $delay = $args{DELAY};
    my $strength = $args{STRENGTH};
    my $concentration = $args{CONCENTRATION};
    my $duty = $args{DUTY};
    my $rftime = $args{RFTIME};
    my $steps = $args{STEPS};

    my $pulse_width = $duty * $period / 100.0;
    my $ramp_duty = ($rftime) / $period * 100.0;
    my $square_duty = ($pulse_width - 2*$rftime) / $period * 100.0;

    if ($pulse_width < 2 * $rftime) {
        printn "ERROR: ramp_equation -- parameters are inconsistent (rise/fall time exceeds computed pulse width)";
        exit(1);
    }

    # n.b. in the square() function, very important to *(pi) as very last step, else roundoff error introduced messes up the square()
    # function interval, such that square(falling edge time) = +1 instead of -1 --> this introduces an extra discontinuity!!
    # i.e. at event time, we will add old value of leading square with new value of lagging square, and this gives a spike with twice
    # the intended value at (and only at) the event time!!!
    my $ramp_source_node = "(";
    my $rise_ramp_delay = $delay;
    my $fall_ramp_delay = $delay + $pulse_width - $rftime;
    $ramp_source_node .= "0.5*(mod(t,$period)-$rise_ramp_delay)/$rftime*$concentration*$strength*(1 + square((t-$rise_ramp_delay)/$period*2*pi, $ramp_duty))";
    $ramp_source_node .= "+ 0.5*$concentration*$strength*(1 + square((t-$rise_ramp_delay-$rftime)/$period*2*pi, $square_duty))";
    $ramp_source_node .= "+ 0.5*(-mod(t,$period)+$fall_ramp_delay+$rftime)/$rftime*$concentration*$strength*(1 + square((t-$fall_ramp_delay)/$period*2*pi, $ramp_duty))";
    $ramp_source_node .= ")";

    my @events;
    if ($delay + $rftime != $fall_ramp_delay) {
        @events= ($delay, $delay + $rftime, $fall_ramp_delay, $fall_ramp_delay + $rftime);
    } else {
        @events= ($delay, $delay + $rftime, $fall_ramp_delay + $rftime);
    }
    # only events in the first period are included here. also, if the delay offset
    # is too large, some events will fall outside the period, so we must use
    # a modulus operator to fix this.  we also use Utils::fmod which can handle
    # fractional dividends
    @events = map {fmod($_, $period)} @events;
    my @order = sort {$events[$a] <=> $events[$b]} (0..$#events);
    @order = grep {$events[$_] != 0} @order;   # don't include 0.0 as an event
    @events = map {$events[$_]} @order;
    my @values = (0, $concentration, $concentration, 0);
    @values = map {$values[$_]} @order;

    printn "ramp_equation: event list is @events" if $verbosity >= 3;
    printn "ramp_equation: value list is @values" if $verbosity >= 3;

    if ($duty < 100) {
        return {
            equations => [
                "null -> $node; clamp_source_$node=\"$ramp_source_node\"",
                "$node -> null; clamp_sink_$node=$strength",
            ],
            events => \@events,
            values => \@values,
        }
    } else {
        return {
            equations => [
                "null -> $node; clamp_source_$node=".$concentration*$strength,
                "$node -> null; clamp_sink_$node=$strength"
            ],
            events => [],
            values => [],
        };
    }
}


#--------------------------------------------------------------------------------------
# Function: ss_ramp_equation
# Synopsys: Generate a non-periodic source/sink equation pair to generate linear ramp
#           up and down over several steps.  The system is brought to steady-state at
#           each step.
#--------------------------------------------------------------------------------------
sub ss_ramp_equation {
    my %args = (
        # default values
        NODE => undef,
        DELAY => "~",  # default is to wait for steady-state before applying stimulus
        STRENGTH => undef,
        RANGE    => undef,
        STEPS => undef,
        RAMP_TIME => undef,
        @_,               # argument pair list overwrites defaults
    );

    check_args(\%args, 6);

    my $node = $args{NODE};
    my $delay = $args{DELAY};
    my $strength = $args{STRENGTH};
    my $range = $args{RANGE};
    my $steps = $args{STEPS};
    my $ramp_time = $args{RAMP_TIME};

    my $step_size = $range / $steps;  # the concentration changing size
    my $step_time = $ramp_time / $steps;

    my @events = ($delay, map {"~"} (1..2*$steps));
    my @values = (0);
    push @values, map {$_*$step_size} (1..$steps);
    push @values, map {($steps-$_)*$step_size} (1..$steps);

    my $ramp_source_node = "(";
    for (my $i=0; $i < (@events-1)/2; $i++) {
        my $ii = $i + 1;
        my $jj = $i + 1 + (@events-1)/2;
        $ramp_source_node .= "+(event_flags($ii) && ~event_flags($jj))*min((t-event_times($ii))/$step_time, 1)*$step_size*$strength";
        $ramp_source_node .= "+event_flags($jj)*max(1-(t-event_times($jj))/$step_time, 0)*$step_size*$strength";
    }
    $ramp_source_node .= ")";


    printn "ramp_equation: event list is @events" if $verbosity >= 3;
    printn "ramp_equation: value list is @values" if $verbosity >= 3;

    return {
        equations => [
            "null -> $node; clamp_source_$node=\"$ramp_source_node\"",
            "$node -> null; clamp_sink_$node=$strength",
        ],
        events => \@events,
        values => \@values,
    }
}


#--------------------------------------------------------------------------------------
# Function: rand_ss_ramp_equation
# Synopsys: Generate a non-periodic source/sink equation pair to generate linear ramp
#           up and down over several steps but with random step size.
#           The system is brought to steady-state at each step.
#--------------------------------------------------------------------------------------
sub rand_ss_ramp_equation {
    my %args = (
        # default values
        NODE => undef,
        DELAY => "~",  # default is to wait for steady-state before applying stimulus
        STRENGTH => undef,
        RANGE    => undef,
        STEPS => undef,
        RAMP_TIME => undef,
        @_,               # argument pair list overwrites defaults
    );

    check_args(\%args, 6);

    my $node = $args{NODE};
    my $delay = $args{DELAY};
    my $strength = $args{STRENGTH};
    my $range = $args{RANGE};
    my $steps = $args{STEPS};
    my $ramp_time = $args{RAMP_TIME};

    my $step_size = $range;  # the concentration changing size
    my $step_size_exp = 3/($steps-1);
    my $step_time = $ramp_time / $steps || 1;

    my @events = ($delay, map {"~"} (1..2*$steps));
    my @values = (0.0001);
    for (my $i=0; $i < $steps; $i++) {
        $step_size = $values[$i] + ($range ** (rand() * 3));
        push @values, $step_size;
    }
    for (my $i=0; $i < $steps; $i++) {
        $step_size = max_numeric(($values[-1] - ($range ** (rand() * 3))), 0);
        push @values, $step_size;
    }

    my $ramp_source_node = "(";
    for (my $i=0; $i < (@events-1)/2; $i++) {
        my $ii = $i + 1;
        my $jj = $i + 1 + (@events-1)/2;
        my $jj_pre = $jj - 1;
        my $amp1 = $values[$ii]-$values[$i];
        my $amp2 = $values[$jj]-$values[$jj_pre];
        $ramp_source_node .= "+(event_flags($ii) && ~event_flags($jj))*min((t-event_times($ii))/$step_time, 1)*($amp1)*$strength";
        $ramp_source_node .= "+event_flags($jj)*max(1-(t-event_times($jj))/$step_time, 0)*($amp2)*$strength";
    }
    $ramp_source_node .= ")";


    printn "ramp_equation: event list is @events" if $verbosity >= 3;
    printn "ramp_equation: value list is @values" if $verbosity >= 3;

    return {
        equations => [
            "null -> $node; clamp_source_$node=\"$ramp_source_node\"",
            "$node -> null; clamp_sink_$node=$strength",
        ],
        events => \@events,
        values => \@values,
    }
}

