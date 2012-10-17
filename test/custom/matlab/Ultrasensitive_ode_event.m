%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:     ode_event.m
% Synopsys: ODE_EVENT(@odesolver, @odefun, tv, y0, opt, ode_events,
%                     SS_timescale, SS_RelTol, SS_AbsTol, varargin)
%           Splits an ODE system integration into multiple integration
%           intervals. Each interval is integrated until the specified
%           event time, or until a steady-state condition is attained.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DETAILED DESCRIPTION:
% ---------------------
% The user supplies a list of N event times via the parameter EVENTS.  The
% integration of the system will by split into (N + 1) intervals.  Each
% event may be either
%   i) a positive value indicating an absolute time,
%  ii) 0, to indicate that the system should be integrated to steady-state,
% iii) a negative value, indicating a time relative to the last event 
%      that occurred (e.g. reaching a steady-state).
% 
% Specifying non-zero event times is useful when dealing with time-varying
% parameters in a system of ODEs that change discontinuously at specific
% times.  Here, splitting the integration into multiple intervals instead
% of calling the solver once, avoids a situation whereby the solver could
% fail to detect the discontinuity.
%
% Specifying 0 for an event time indicates that the system should be
% integrated to steady-state during that interval.
%
% Input Arguments:
% ----------------
%
% odesolver - The solver to use, e.g. ode23s or ode15s.
%
% odefun - an ode file, just as would be input to the odesolver itself.
% 
% tv     - interval of integration and times at which the result is to be
%          sampled, usually just [t0 tf], but can be e.g. [t0:0.01:tf]
%
% y0     - a vector of initial conditions for the ode file at time t_start.
%
% opt    - options to pass to the odesolver. If no options are desired, the
%          input should be opt = [].
%
% events - A vector of event times at which the integration should
%          re-started. These event times should fall between the initial
%          and final times specified by the TV argument.  Given a set of
%          event times of the form events=[event1 event2 ... eventk ...
%          eventn] and also tv = [t0 tf], the integration will be split
%          into the following intervals:
%
%               [t0     event1]
%               [event1 event2]
%                ...
%               [eventn tf]
%
%          An event time of 0 indicates to the solver that the integration
%          over the current interval should continue until a steady-state
%          stopping condition is reached.  An negative event time means
%          that the integration should continue until a time of (-eventk)
%          has elapsed since the last event.  Thus eventk=100 indicates
%          that the kth event occurs at time t=100, while eventk=-50
%          indicates that the kth event occurs 50 time units after the last
%          event.
%
%          This encoding of event times requires the user to specify a
%          t0>=0. Consequently, a positive value must be supplied for tf,
%          with tf>t0.  
%
%          When eventk indicates that the integration should continue to
%          steady-state, a stopping condition parametrized by
%          SS_timescale, SS_RelTol and SS_AbsTol is evaluated at each
%          integration step to determine if a steady-state is reached.
%          The stopping condition is a component-wise check based on the
%          following formula:
%
%          dy/dt * SS_timescale < max(SS_RelTol * abs(y), SS_AbsTol)
%
%          This condition ensures that the significant digits of the
%          solution vector y, as specified by SS_RelTol, are unchanging
%          over the given timescale. Since the most significant digit
%          is always changing if y exponentially decays to zero, and also
%          because the digits below the Matlab's AbsTol threshold are not 
%          accurate, (c.f. odeset documentation for an explanation of
%          AbsTol) enforcement of the condition is relaxed as |y|
%          approachesSS_AbsTol. Once the stopping condition is true for all
%          solution components, integration is halted.
%
%          The integration will terminate if t=tf, or earlier if the
%          last event is to find the steady-state (i.e. eventn=0) and the
%          stopping condition occurs before tf. In the case where eventn=0,
%          only N intervals are integrated instead of (N+1). Therefore when
%          simulating to steady-state, tf acts as a timeout value for the
%          integration rather than specifying the final integration time.
%
%          The simplest situation involving finding the steady-state is to
%          specify a single steady-state event and a timeout.  For example,
%          specifying EVENT=[0] and TV=[0:0.1:1000] will simulate the system
%          to steady-state or to t=1000, whichever occurs first.
%
%          Two global variables are updated by ode_event and may be accessed
%          from within the user-specified ODEFUN to implement event-dependent
%          or time-varying equations.  The global variable event_flags is
%          a boolean vector with the same size as ODE_EVENTS and indicating
%          whether or not the corresponding event has occurred.  Similarly,
%          the global event_times indicates at which time the corresponding
%          event has occurred.
%
% SS_timescale - A vector of up to N+1 timescales used to detect the
%          reaching of a steady-state condition by slow processes.  Each
%          timescale should correspond to the time constant of the
%          slowest process in the corresponding interval.  If a partial
%          vector is supplied, the last element is extended as necessary.
%          If empty, all elements default to 100.
%
% SS_RelTol - A vector of up to N+1 steady-state tolerances which
%          specify an upper bound for the relative change in y over the
%          corresponding SS_timescale.
%
% SS_AbsTol - A vector of up to N+1 relative steady-state absolute 
%          tolerances for y below which enforcement of SS_RelTol
%          is relaxed.
%
% varargin - any other variables to pass to odefun, such as parameter
%             values.
%
% Returned Values:
% ----------------
%
% [t, y] - the solution to the ode file on the interval [t_start t_end].
%
% l      - a vector containing the length of each segment solved by odesolver,
%          ie. the number of points generated in [t, y] between events.
%          This is useful for working with only certain parts of the
%          solution. For example, you may give a model time to come to
%          equilibrium before perturbing the system and throw that time
%          away as a transient.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y, l] = ode_event(odesolver, odefun, tv, y0, opt, ...
    ode_events, SS_timescale, SS_RelTol, SS_AbsTol, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: ode_event
% Synopsys: Split integration into multiple intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the ouputs:
t = [];
y = [];
l = [];

% check arguments
if (length(tv) < 2)
    disp('ERROR: ode_event -- tv must have at least 2 elements');
    return
end
if (tv(1) < 0)   % check t0
    disp('ERROR: ode_event -- integration start time cannot be negative');
    return
end
if (tv(end) <= 0)  % check tf
    disp('ERROR: ode_event -- integration end time must be positive');
    return
end
if (ode_events(1) > 0 && (ode_events(1) <= tv(1)))
    disp('ERROR: ode_event -- events must be within specified integration vector, compare t0 and first event');
    return
end
if (ode_events(end) > 0 && (ode_events(end) >= tv(end)))
    disp('ERROR: ode_event -- events must be within specified integration vector, compare tf and last event');
    return
end

% extract and save the RelTol and AbsTol options
% n.b. AbsTol is only used as a default value for SS_AbsTol
AbsTol = opt.AbsTol;
if isempty(AbsTol)
    AbsTol = 1e-6;  % this is the Matlab default
end
% n.b. RelTol is only used as a default value for SS_RelTol
RelTol = opt.RelTol;
if isempty(RelTol)
    RelTol = 1e-3;  % this is the Matlab default
end

% initialize global event_flags and event_times variables
global event_flags;
global event_times;
event_flags = [];
event_flags(1:length(ode_events)) = 0;
event_times = ode_events;

% initialize SS_timescale
if (length(SS_timescale) == 0)
    % assign default timescale of 100
    SS_timescale(1:length(ode_events)+1) = 100;
elseif (length(SS_timescale) >= 1)
    % extend the last value to any missing intervals
    SS_timescale(end+1:length(ode_events)+1) = SS_timescale(end);
end
if (length(find(SS_timescale <= 0)) > 0)
    disp('ERROR: ode_event -- SS_timescale must be positive in each interval');
    return
end

% initialize SS_RelTol
if (length(SS_RelTol) == 0)
    % assign default value of RelTol
    SS_RelTol(1:length(ode_events)+1) = RelTol;
elseif (length(SS_RelTol) >= 1)
    % extend the last value to any missing intervals
    SS_RelTol(end+1:length(ode_events)+1) = SS_RelTol(end);
end
if (length(find(SS_RelTol <= 0)) > 0)
    disp('ERROR: ode_event -- SS_RelTol must be positive in each interval');
    return
end
if (length(find(SS_RelTol >= 1)) > 0)
    disp('ERROR: ode_event -- SS_RelTol must be < 1 in each interval');
    return
end

% initialize SS_AbsTol
if (length(SS_AbsTol) == 0)
    % assign default value of AbsTol
    SS_AbsTol(1:length(ode_events)+1) = AbsTol;
elseif (length(SS_AbsTol) >= 1)
    % extend the last value to any missing intervals
    SS_AbsTol(end+1:length(ode_events)+1) = SS_AbsTol(end);
end
if (length(find(SS_AbsTol <= 0)) > 0)
    disp('ERROR: ode_event -- SS_AbsTol must be positive in each interval');
    return
end

% initialize events variable
% n.b. indexing is +1 offset relative to event_flags and event_times
events=[tv(1) ode_events tv(end)];

% call the solver at each event time specified.
for i = 1:length(events)-1
    last_event_time = events(i);
    next_event_time = events(i+1);
    if next_event_time < 0
        next_event_time = last_event_time + abs(next_event_time);
        event_times(i) = next_event_time;
        events(i+1) = next_event_time;
    end
    timeout = events(find(events > last_event_time,1));
    if (isempty(timeout))
        disp('WARNING: reached last timeout value, cannot complete all integration intervals');
        break;
    end
    if (next_event_time > 0 && next_event_time > last_event_time)  % regular integration interval?
        TV = [last_event_time tv(find(tv > last_event_time & tv < next_event_time)) next_event_time];
        str=sprintf('ode_event: integrating from %f to %f', last_event_time, next_event_time);
        disp(str);
        % uncomment the following 2 lines to debug integration vector
        %    str=[sprintf('ode_event: current integration vector is ') sprintf('%.2f ', TV)];
        %    disp(str);
        [T, Y] = odesolver(odefun, TV, y0, opt, varargin{:});
        % initial condition for next interval is final value of this one
        y0 = Y(end,:);
        % append results of current interval to final result
        length_t = length(t);
        if length_t == 0
            y = Y;
            t = T;
        else
            y(end:end + length(T) - 1, :) = Y;
            t(end:end + length(T) - 1, 1) = T;
        end
        % save length of current interval
        l(i) = length(T);
        % update event_flags, marking current interval as passed
        event_flags(i) = 1;
        % no need to update event_times
    elseif (next_event_time == 0) % integrate to steady-state?
        TV = [last_event_time tv(find(tv > last_event_time & tv < timeout)) timeout];
        str=sprintf('ode_event: integrating from %f to steady-state (or to timeout at %f)', last_event_time, timeout);
        disp(str);
        
        check_steady_state_fh = @(t,y,flag,varargin) check_steady_state(...
            t,y,flag,odefun,SS_timescale(i),...
            SS_RelTol(i),SS_AbsTol(i),varargin{:});
        
        opt_steady_state = odeset(opt, 'OutputFcn', check_steady_state_fh);
                
        [T, Y] = odesolver(odefun, TV, y0, opt_steady_state, varargin{:});
        % initial condition for next interval is final value of this one
        y0 = Y(end,:);
        % append results of current interval to final result
        length_t = length(t);
        if length_t == 0
            y = Y;
            t = T;
        else
            y(end:end + length(T) - 1, :) = Y;
            t(end:end + length(T) - 1, 1) = T;
        end
        % save length of current interval
        l(i) = length(T);
        events(i+1) = T(end);
        event_times(i) = T(end);
        event_flags(i) = 1;
        str=sprintf('ode_event: integrated from %f to %f', last_event_time, events(i+1));
        disp(str);
        if (events(i+1) >= timeout)
            disp('WARNING: did not reach steady-state');
        end
        if (i == length(events)-2)  % steady-state event is last before tf ?
            break;
        end
    elseif (next_event_time <= last_event_time)
        str=sprintf('WARNING: skipping interval from t=%f to t=%f because it has already passed', last_event_time, next_event_time);
        disp(str);
        events(i+1) = last_event_time;
    else
        str=sprintf('ERROR: internal error while processing interval from t=%f to t=%f', last_event_time, next_event_time);
        disp(str);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: check_steady_state
% Synopsys: Function to pass as OutputFcn in odeset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = check_steady_state(t,y,flag, ...
    odefun,SS_timescale,SS_RelTol,SS_AbsTol,varargin)

% status variable expected by odesolver as return value
status = 0;

if strcmp(flag, 'init')
    %    disp('check_steady_state OutputFcn initialized');
elseif ~strcmp(flag, 'done')
    y_end = y(:,end); % for speed
    dy = odefun(t(:,end),y_end,varargin{:}) * SS_timescale;
    dy_threshold = max(SS_RelTol * abs(y_end), SS_AbsTol);
    ss_condition = abs(dy) < dy_threshold;
    if (ss_condition)
        str=sprintf('ode_event: steady-state stopping condition reached at t=%f', t(end));
        disp(str);
        fast_check_flag = 1;
        status = 1;
    end
end
end

