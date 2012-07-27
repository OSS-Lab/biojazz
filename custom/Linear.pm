#-#####################################################################################
#- File:     Linear.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Linear;
use Class::Std::Storable;
use base qw(Scoring);
{
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG);

    use Stimulus;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
#    my %A_of :ATTR(get => 'A', set => 'A', init_arg => 'A');  # constructor must supply initialization
#    my %B_of :ATTR(get => 'B', set => 'B', default => 'yyy'); # default value is yyy

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
    # Function: score_genome
    # Synopsys: Score for ultrasensitivity.
    #--------------------------------------------------------------------------------------
    sub score_genome {
	my $self = shift;
	my $genome_model_ref = shift;

	confess "ERROR: internal error, $genome_model_ref not a GenomeModel" if !$genome_model_ref->isa('GenomeModel');

	my $config_ref = $self->get_config_ref();
	my $genome_name = $genome_model_ref->get_name();
	my $matlab_work = $self->get_matlab_work();

	my $stats_ref = $genome_model_ref->get_stats_ref();
	if (!defined $stats_ref) {
	    printn "WARNING: stats_ref is not defined for $genome_name";
	    $stats_ref = {};
	    $genome_model_ref->set_stats_ref($stats_ref);
	}
	my $history_ref = $genome_model_ref->get_history_ref();

	printn "Linear::score_genome scoring $matlab_work/$genome_name.mod";

	#---------------------------------------------------------
	# INIT SCORING
	#---------------------------------------------------------
	my $elite_flag = $genome_model_ref->get_elite_flag();
	if ($elite_flag) {
	    printn "Linear::score_genome elite individual already scored, previous score=$stats_ref->{score}";
	    return if ($config_ref->{rescore_elite} == 0);
	    printn "Linear::score_genome re-scoring elite individual";
	    $genome_model_ref->clear_stats(preserve => ["score_count", "score"]);
	    $stats_ref->{score_count}++;
	} else {
	    printn "Linear::score_genome scoring non-elite individual...";
	    $genome_model_ref->clear_stats(preserve => []);
	    $stats_ref->{score_count} = 1;
	    $stats_ref->{score} = 0;
	}

	#---------------------------------------------------------
	# CREATE I/O GENES
	#---------------------------------------------------------
	my $lg_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
	    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
	    regulated_concentration => $config_ref->{regulated_concentration_min}, # all-zeroes
	    UNUSED => "0000",
	    domains => [
		{
		    allosteric_flag => 0,
		    RT_transition_rate => $config_ref->{RT_transition_rate_min}, # all-zeroes
		    TR_transition_rate => $config_ref->{TR_transition_rate_min}, # all-zeroes
		    RT_phi => $config_ref->{RT_phi_min}, # all-zeroes
		    protodomains => [
			{
			    type => "bsite",
			    substrate_polarity => 0,
			    binding_profile => $config_ref->{lg_binding_profile},
			    kf_profile => "0",
			    kb_profile => "0",
			    kp_profile => "0",
			    Keq_ratio => $config_ref->{Keq_ratio_min},
			    kf_polarity_mask => "0",
			    kb_polarity_mask => "0",
			    kf_conformation_mask => "0",
			    kb_conformation_mask => "0",
			    kp_conformation_mask => "0",
			    UNUSED => "0",
			},
		       ],
		    UNUSED => "0",
		},
	       ],
	});
	printn "lg_sequence=".$lg_sequence_ref->get_sequence() if $verbosity >= 2;
	my $tg_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
	    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
	    regulated_concentration => $config_ref->{TG_init}, # all-zeroes
	    UNUSED => "0000",
	    domains => [
		{
		    allosteric_flag => 0,
		    RT_transition_rate => $config_ref->{RT_transition_rate_min}, # all-zeroes
		    TR_transition_rate => $config_ref->{TR_transition_rate_min}, # all-zeroes
		    RT_phi => $config_ref->{RT_phi_min}, # all-zeroes
		    protodomains => [
			{
			    type => "msite",
			    substrate_polarity => 0,
			    binding_profile => $config_ref->{tg_binding_profile},
			    kf_profile => "0",
			    kb_profile => "0",
			    kp_profile => "0",
			    Keq_ratio => $config_ref->{Keq_ratio_min},
			    kf_polarity_mask => "0",
			    kb_polarity_mask => "0",
			    kf_conformation_mask => "0",
			    kb_conformation_mask => "0",
			    kp_conformation_mask => "0",
			    UNUSED => "0",
			},
		       ],
		    UNUSED => "0",
		},
	       ],
	});
	printn "tg_sequence=".$tg_sequence_ref->get_sequence() if $verbosity >= 2;

	#---------------------------------------------------------
	# STIMULUS EQUATIONS
	#---------------------------------------------------------
	my $step_time = $config_ref->{LG_rftime}/($config_ref->{LG_steps} - 1);
	my $step_times_per_period = $config_ref->{LG_period}/$step_time;
	my $LG_delay = ($config_ref->{LG_delay} eq "random") ? (int(rand($step_times_per_period)) * $step_time) : $config_ref->{LG_delay};

	my $stimulus_sub_ref = \&{$config_ref->{stimulus}};
	my $stimulus_ref = &$stimulus_sub_ref(
	    NODE => "LG0000_x",
	    PERIOD => $config_ref->{LG_period},
	    STRENGTH => $config_ref->{LG_strength},
	    CONCENTRATION => $config_ref->{LG_level},
	    DUTY => $config_ref->{LG_duty},
	    RFTIME => $config_ref->{LG_rftime},
	    STEPS => $config_ref->{LG_steps},
	    DELAY => $LG_delay,
	   );
	my ($lg_source_staircase, $lg_sink_staircase) = @{$stimulus_ref->{equations}};
	my $event_list = join " ", @{$stimulus_ref->{events}};
	my $values_list = join " ", @{$stimulus_ref->{values}};
	printn "Stimulus:\n". join "\n",($lg_source_staircase, $lg_sink_staircase, $event_list, $values_list) if $verbosity >= 2;

	#---------------------------------------------------------
	# PARSE/TRANSLATE GENOME AND I/O GENES
	#---------------------------------------------------------
	my $genome_iref = $genome_model_ref->parse(
	    [
		sequence_ref => $lg_sequence_ref,
		prefix => "L",
	       ],
	    [
		sequence_ref => $tg_sequence_ref,
		prefix => "T",
	       ],
	   );
	$genome_model_ref->check();
	printn $genome_iref->sprint(colour_flag => 0) if $verbosity >= 2 || $config_ref->{sprint_genome};
	$genome_model_ref->translate();

	#---------------------------------------------------------
	# PRUNE NETWORK
	#---------------------------------------------------------
	my $genome_ref = $genome_model_ref->get_parser_ref();
	$genome_ref->build_network();
	printn join(",", map {$_->get_name()} @{$genome_ref->get_gene_adjacency_matrix_nodes_ref()}) if $verbosity >= 1;
	printn $genome_ref->get_gene_adjacency_matrix_ref()->[0]->sprint_matrix() if $verbosity >= 2;
	$genome_ref->compute_high_order_adjacency_matrix(8);
	$genome_ref->compute_connectivity_matrix();
	printn $genome_ref->get_gene_connectivity_matrix_ref()->sprint_matrix() if $verbosity >= 2;
	$stats_ref->{num_pruned_genes} = scalar $genome_ref->prune_isolated_genes();
	$stats_ref->{num_genes} = $genome_model_ref->get_num_genes();
	
	#---------------------------------------------------------
	# SCORING: 2 pts -- LG/TG connected to anything
	#---------------------------------------------------------
	my $connectivity_score = 0;
	my $gene_ref = $genome_model_ref->get_gene_parser_ref();
	my $lg_gene_ref = $gene_ref->lookup_object_instance_by_name("LG0000");
	my $tg_gene_ref = $gene_ref->lookup_object_instance_by_name("TG0000");
	if ($lg_gene_ref->get_export_flag()) {
	    $connectivity_score++;
	    printn "LG is connected" if $verbosity >= 1;
	}
	if ($tg_gene_ref->get_export_flag()) {
	    $connectivity_score++;
	    printn "TG is connected" if $verbosity >= 1;
	}

	#---------------------------------------------------------
	# SCORING: 90 + 300 pts -- LG/TG subnets
	#---------------------------------------------------------
	my (@lg_subnet, @tg_subnet);
	if ($connectivity_score == 2) {	# LG/TF connected
	    @lg_subnet = $genome_ref->get_connected_genes($lg_gene_ref);
	    @tg_subnet = $genome_ref->get_connected_genes($tg_gene_ref);

	    printn "LG connects to ".join ",", (map {$_->get_name()} @lg_subnet) if $verbosity >= 1;
	    printn "TG connects to ".join ",", (map {$_->get_name()} @tg_subnet) if $verbosity >= 1;

	    # max 90 points for subnet size
	    my $lg_subnet_size = (@lg_subnet > 30) ? 45 : @lg_subnet;
	    my $tg_subnet_size = (@tg_subnet > 30) ? 45 : @tg_subnet;
	    $connectivity_score += ($lg_subnet_size + $tg_subnet_size);

	    # score -- LG/TG connected to each other
	    if (grep /LG/, (map {$_->get_name()} @tg_subnet)) {
		printn "LG and TG are connected to each other" if $verbosity >= 1;
		$connectivity_score += 100;
	    }
	    if (grep /TG/, (map {$_->get_name()} @lg_subnet)) {
		printn "TG and LG are connected to each other" if $verbosity >= 1;
		$connectivity_score += 100;
	    }

	    # TG must have at least two connected proteins
	    my @tg_adjacent_genes = $genome_ref->get_adjacent_genes($tg_gene_ref);
	    if (@tg_adjacent_genes >= 2) {
		printn "TG has > 2 connected proteins" if $verbosity >= 1;
		$connectivity_score += 100;
	    }
	}

	#---------------------------------------------------------
	# PICK AN ADJACENT PROTEIN TO MAKE INTO KINASE/PHOSPHATASE
	#---------------------------------------------------------
	if ($connectivity_score >= 300) {
	    my ($tg_protodomain_ref) = map {$_->get_protodomains()} $tg_gene_ref->get_domains();
	    my @adjacent_kinases = $genome_ref->find_adjacent_csites($tg_protodomain_ref,0);
	    my @adjacent_phosphatases = $genome_ref->find_adjacent_csites($tg_protodomain_ref,1);

	    my $csite_field_value = $tg_protodomain_ref->untranslate_field("type", "csite");

	    my ($lg_protodomain_ref) = map {$_->get_protodomains()} $lg_gene_ref->get_domains();
	    my @adjacent_protodomains = simple_difference([$genome_ref->get_adjacent_protodomains($tg_protodomain_ref)], [$lg_protodomain_ref]);
	    if (@adjacent_kinases < 1 && !defined $genome_model_ref->get_stat("created_kinase")) {
		my @protodomains = simple_difference(\@adjacent_protodomains, \@adjacent_phosphatases);
		if (@protodomains >= 1) {
		    printn "editing genome to create kinase";
		    my $index = int rand @protodomains;
		    my $protodomain_ref = $protodomains[$index];
		    confess "ERROR: unexpected condition" if ($protodomain_ref->get_translation_ref()->{type} eq "csite");
		    $protodomain_ref->set_field("type", $csite_field_value);
		    $protodomain_ref->set_field("substrate_polarity", 0);
		    @adjacent_kinases = (@adjacent_kinases, $protodomain_ref);
		    $stats_ref->{created_kinase} = 1;
		    $genome_model_ref->add_history("Edited genome protodomain ".$protodomain_ref->get_name()." to create kinase csite");
		}
	    }
	    if (@adjacent_phosphatases < 1 && !defined $genome_model_ref->get_stat("created_phosphatase")) {
		my @protodomains = simple_difference(\@adjacent_protodomains, \@adjacent_kinases);
		if (@protodomains >= 1) {
		    printn "editing genome to create phosphatase";
		    my $index = int rand @protodomains;
		    my $protodomain_ref = $protodomains[$index];
		    confess "ERROR: unexpected condition" if ($protodomain_ref->get_translation_ref()->{type} eq "csite");
		    $protodomain_ref->set_field("type", $csite_field_value);
		    $protodomain_ref->set_field("substrate_polarity", 1);
		    @adjacent_phosphatases = (@adjacent_phosphatases, $protodomain_ref);
		    $stats_ref->{created_phosphatase} = 1;
		    $genome_model_ref->add_history("Edited genome protodomain ".$protodomain_ref->get_name()." to create phosphatase csite");
		}
	    }

	    #---------------------------------------------------------
	    # SCORING: 100 + 100 pts -- connected kinase/phosphatase
	    #---------------------------------------------------------
	    $connectivity_score += 100 * (@adjacent_kinases > 1 ? 1 : @adjacent_kinases);
	    $connectivity_score += 100 * (@adjacent_phosphatases > 1 ? 1 : @adjacent_phosphatases);
	    $stats_ref->{num_adjacent_kinases} = scalar(@adjacent_kinases);
	    $stats_ref->{num_adjacent_phosphatases} = scalar(@adjacent_phosphatases);
	}
	
	if ($connectivity_score >= 500 &&
	    $stats_ref->{num_adjacent_kinases} != 0 &&
	    $stats_ref->{num_adjacent_phosphatases} != 0) {

	    #---------------------------------------------------------
	    # RE-PARSE/TRANSLATE GENOME AND I/O GENES
	    #---------------------------------------------------------
	    $genome_iref = $genome_model_ref->parse(
		[
		    sequence_ref => $lg_sequence_ref,
		    prefix => "L",
		   ],
		[
		    sequence_ref => $tg_sequence_ref,
		    prefix => "T",
		   ],
	       );
	    $genome_model_ref->check();
	    $genome_model_ref->translate();

	    # *****  NETWORK PERTURBATION HERE *****
	    #    printn "PERTURBATION (BEFORE):  ".$sdb->{protodomain_table}{xxx}{yyy};
	    #	delete_gene(...)
	    #    printn "PERTURBATION (AFTER):  ".$sdb->{protodomain_table}{xxx}{yyy};
	    # **************************************

	    # exclude proteins not in LG/TG subnet from export
	    my @proteins_not_in_subnet = simple_difference([$genome_model_ref->get_genes()], [union(\@lg_subnet, \@tg_subnet)]);
	    map {$_->set_export_flag(0)} @proteins_not_in_subnet;

	    #---------------------------------------------------------
	    # GENERATE ANC/FACILE MODEL
	    #---------------------------------------------------------
	    my $anc_model = $genome_model_ref->get_genome_parser_ref()->export_anc(
		max_external_iterations => $config_ref->{max_external_iterations},
		max_internal_iterations => $config_ref->{max_internal_iterations},
		max_complex_size => $config_ref->{max_complex_size},
		max_species => $config_ref->{max_species},
		max_csite_bound_to_msite_number => $config_ref->{max_csite_bound_to_msite_number},
		default_steric_factor => $config_ref->{default_steric_factor},
		export_graphviz => ref $config_ref->{export_graphviz} ? (join ",",@{$config_ref->{export_graphviz}}) : $config_ref->{export_graphviz},
		equations => [$lg_source_staircase, $lg_sink_staircase],
	       );
	    burp_file("$matlab_work/$genome_name.mod", $anc_model);
	    system("$ENV{ANC_HOME}/anc.pl --report=species $matlab_work/$genome_name.mod");
	    my @facile_model = slurp_file("$matlab_work/$genome_name.eqn");

	    # now check that TF_0 and TF_1 are involved in at least 1 reaction each
	    my $num_reactions_tg_1 = grep (/(->|<-).* TG0000_1/, @facile_model);
	    my $num_reactions_tg_0 = grep (/(->|<-).* TG0000_0/, @facile_model);

	    $connectivity_score += 100 * ($num_reactions_tg_1 > 1 ? 1 : $num_reactions_tg_1);
	    $connectivity_score += 100 * ($num_reactions_tg_0 > 1 ? 1 : $num_reactions_tg_0);

	    $stats_ref->{num_reactions_tg_0} = $num_reactions_tg_0;
	    $stats_ref->{num_reactions_tg_1} = $num_reactions_tg_1;
	}
#	printn "($node_ID, $generation, $individual) connectivity_score=$connectivity_score, new_lg_subnet_size=$new_lg_subnet_size new_lg_subnet_size=$new_lg_subnet_size new_tf_subnet_size=$new_tf_subnet_size, num_adjacent_kinases=$num_adjacent_kinases, num_adjacent_phosphatases=$num_adjacent_phosphatases, num_reactions_tf_1=$num_reactions_tf_1, num_reactions_tf_0=$num_reactions_tf_0";

	if ($connectivity_score > 700 &&
	    $stats_ref->{num_reactions_tg_0} != 0 &&
	    $stats_ref->{num_reactions_tg_1} != 0) {

	    #---------------------------------------------------------
	    # RUN FACILE
	    #---------------------------------------------------------
	    $self->facile_run(
		SOLVER => $config_ref->{solver},
		SOLVER_OPTIONS => ("odeset('InitialStep', $config_ref->{InitialStep}, ".
				   "'AbsTol', $config_ref->{AbsTol}, ".
				   "'RelTol', $config_ref->{RelTol}, ".
				   "'MaxStep', $config_ref->{MaxStep})"),
		T_FINAL => $config_ref->{LG_period},
		T_SAMPLING =>'[t0:0.1:tf]',
		T_EVENTS => "$event_list",
		EQN_FILE => "$matlab_work/$genome_name.eqn",
	       );

	    $self->anc_process_species_report("$matlab_work/$genome_name.species.rpt");
	    my @anc_species = $self->anc_get_species();
	    printn "ANC_SPECIES: @anc_species" if $verbosity >= 2;

	    #---------------------------------------------------------
	    # RUN MATLAB SIM
	    #---------------------------------------------------------
	    printn "Linear::score_genome: running matlab driver...";
	    my $matlab_ref = $self->get_matlab_ref();
	    $matlab_ref->cmd("clear all; ${genome_name}Driver");
	    $matlab_ref->wait_on("Facile.*done");

	    #---------------------------------------------------------
	    # PLOT RESULTS
	    #---------------------------------------------------------
	    if (defined $config_ref->{plot_input} && $config_ref->{plot_input}) {
		$self->matlab_plot_complex(figure => 900,
					   complex => "LG0000_x",
					   title_prefix => "$genome_name",
					  );
	    }
	    if (defined $config_ref->{plot_output} && $config_ref->{plot_output}) {
		$self->matlab_plot_complex(figure => 901,
					   complex => "TG0000_0",
					   title_prefix => "$genome_name",
					  );
		$self->matlab_plot_complex(figure => 902,
					   complex => "TG0000_1",
					   title_prefix => "$genome_name",
					  );
	    }
	    if (defined $config_ref->{plot_species} && $config_ref->{plot_species}) {
		$self->matlab_plot_all_complexes();
	    }
	    $self->matlab_cmd("disp('Done plotting')\n");
	    $self->matlab_wait_on("Done plotting");
	    system("sleep 1");

	    #---------------------------------------------------------
	    # SCORE COMPLEXITY
	    #---------------------------------------------------------
	    my $line_count = `wc -l $matlab_work/$genome_name.eqn`;
	    $line_count =~ /^\s*(\S+)\s*/;   # extract the line count
	    $stats_ref->{complexity} = $1;
	    $stats_ref->{complexity_score} = 1/($stats_ref->{complexity} + 1);

	    #---------------------------------------------------------
	    # SCORE STEADY STATE
	    #---------------------------------------------------------
	    printn "computing steady-state slope...";

	    my @event_list = split " ", $event_list;
	    my @steady_state_slope;
	    my @steady_state_event_list = (@event_list, $config_ref->{LG_period});
	    for (my $i = 0; $i < @steady_state_event_list; $i++) {
		my $steady_state_event = $steady_state_event_list[$i];
		my $state_var = ($config_ref->{solver} eq "stoch") ? "simdata" : "y";
		my @state_vector_delta = $self->matlab_get_state_delta(t1 => $steady_state_event - 1, t2 => $steady_state_event, state_var => $state_var);
		$steady_state_slope[$i] = max_numeric(abs(max_numeric(@state_vector_delta)), abs(min_numeric(@state_vector_delta)));
		printn "slope @ $steady_state_event = ".sprintf("%.5e",$steady_state_slope[$i]);
	    }
	    my $steady_state_slope = $stats_ref->{steady_state_slope} = max_numeric(@steady_state_slope);
	    # MATLAB PLOT: clear all; ss_th = 1; SS=[-2:0.01:2]; ss=10.^(SS); k = ss/ss_th; n=1; kn=k.^n; semilogx(ss,1./(1+kn))
	    my $steady_state_threshold = $config_ref->{steady_state_threshold};
	    $stats_ref->{steady_state_score} = n_hill($steady_state_slope, $steady_state_threshold, 1);

	    #---------------------------------------------------------
	    # REPORT RESULT VECTOR
	    #---------------------------------------------------------
	    printn "RESULT VECTOR: INPUT = LG OUTPUT = TG0000_1 DELAY=$LG_delay";
	    my (@input_vector, @output_vector, @expected_output_vector);

	    # compute sample list based on staircase equation events (even if actually using ramp)
	    my $staircase_sample_ref = staircase_sample(
		CONCENTRATION => $config_ref->{LG_level},
		PERIOD => $config_ref->{LG_period},
		DUTY => $config_ref->{LG_duty},
		RFTIME => $config_ref->{LG_rftime},
		STEPS => $config_ref->{LG_steps},
		DELAY => $LG_delay,
	       );
	    my @events = (@{$staircase_sample_ref->{events}}, $config_ref->{LG_period});
	    my @values = (@{$staircase_sample_ref->{values}}, $staircase_sample_ref->{final_value});

	    for (my $i=0; $i < @events; $i++) {
		my $t = $events[$i]; 
		$input_vector[$i] = $values[$i];
		$output_vector[$i] = $self->matlab_get_value(complex => "TG0000_1", t => $t);
		$expected_output_vector[$i] = $config_ref->{TG_init} * p_hill($input_vector[$i], $config_ref->{hill_k}, $config_ref->{hill_n}); # positive hill function
		printf("RESULT VECTOR:  t=%-6.2f input vector:  %4.0g output_vector: %12.8g (expected = $expected_output_vector[$i])\n",
		       $t, $input_vector[$i], $output_vector[$i]);
	    }

	    #---------------------------------------------------------
	    # PHASE PLOT
	    #---------------------------------------------------------
	    if (defined $config_ref->{plot_phase} && $config_ref->{plot_phase}) {
		# plot first half in blue, second half in red
		$self->matlab_cmd("halfway = floor(size(LG0000_x,1)/2)");
		$self->matlab_cmd("figure(904); plot(LG0000_x(1:halfway), TG0000_1(1:halfway));title(\'PHASE PLOT\')");
		$self->matlab_cmd("hold on; plot(LG0000_x(halfway+1:end), TG0000_1(halfway+1:end), 'r');");
	    }

	    #---------------------------------------------------------
	    # COMPUTE ERROR, STATS
	    #---------------------------------------------------------
	    if ("@output_vector" !~ /UNDEFINED/) {
		printn "computing correlation...";
		#	    $cos_angle = norm_projection(\@output_vector, \@expected_output_vector);
		#	    if (!defined $cos_angle) {$cos_angle = 0;};  # norm_projection can return undef
		$stats_ref->{correlation} = correlation(@output_vector, @expected_output_vector);
		$stats_ref->{magnitude} = magnitude(@output_vector);
		$stats_ref->{variance} = variance(@output_vector);
		$stats_ref->{stdev} = $stats_ref->{variance} ** (0.5);
	    }

	    $stats_ref->{mean_squared_err} = mean_squared_error(\@output_vector, \@expected_output_vector);
	    $stats_ref->{max_mean_squared_err} = ($config_ref->{TG_init}) ** 2; # since it's the maximum MEAN error^2
	    $stats_ref->{mean_squared_err_score}  = (1-$stats_ref->{mean_squared_err}/$stats_ref->{max_mean_squared_err});

	    # check if any concentrations are negative
	    foreach my $sample (@output_vector) {
		if ($sample < 0) {
		    printn "WARNING: detected negative concentrations in output vector";
		    $stats_ref->{mean_squared_err_score} = 0;
		    $stats_ref->{correlation} = 0;
		}
	    }

	    if (($stats_ref->{mean_squared_err_score} < 0) || ($stats_ref->{mean_squared_err_score} > 1)) {
		# numerics got messed up, set score to zero
		printn "WARNING: computed mean_squared_error_score out of bounds, setting to zero";
		$stats_ref->{correlation} = 0;
		$stats_ref->{mean_squared_err_score} = 0;
	    }


	    #	$steady_state_score = ($steady_state_score > 0.5) ? 1.0 : $steady_state_score;
	    #$stats_ref->{sim_score} =  0.1 * $stats_ref->{steady_state_score};
	    #$stats_ref->{sim_score} += 0.9 * $stats_ref->{mean_squared_err_score} if ($stats_ref->{steady_state_score} > 0.5);
	    $stats_ref->{sim_score} = $stats_ref->{mean_squared_err_score};

	    # round the sim score to prevent integration error from affecting result
	    $stats_ref->{sim_score} = round2sig($stats_ref->{sim_score}, 4);
	}

	$stats_ref->{connectivity_score} = $connectivity_score;
	my @stats_report = ();
	foreach my $stat qw(mean_squared_err_score
			    mean_squared_err
			    max_mean_squared_err
			    steady_state_slope
			    steady_state_score
			    sim_score
			    correlation
			    variance
			    stdev
			    magnitude
			    complexity
			    complexity_score) {
	    push @stats_report, "$stat = ".(defined $stats_ref->{$stat} ? $stats_ref->{$stat} : "UNDEF");
	}
	printn join ", ", @stats_report;

	#---------------------------------------------------------
	# FINAL SCORE
	#---------------------------------------------------------
	$stats_ref->{sim_score} = 0 if !defined $stats_ref->{sim_score};
	my $final_score = 0;

	$final_score += ($stats_ref->{sim_score} > 0) ? 0.1 : 0.1 * (1-exp(-$stats_ref->{connectivity_score}/300));  	# get full points for connectivity if sim_score > 0
	$final_score += 0.8 * $stats_ref->{sim_score};
	$final_score += ($stats_ref->{sim_score} > 0) ? 0.1 * $stats_ref->{complexity_score} : 0.0;   # complexity counts whenever sim works

	$final_score = ($final_score < 0) ? 0 : $final_score;  # prevent neg've scores

# 	if (@problems) {
# 	    printn "Warning: matlab had some problems, tossing out this solution....";
# 	    for (my $i=0; $i < @problems; $i++) {
# 		printn "Warning: (from matlab_wait_on) ==> $problems[$i]";
# 		last if $i > 10;
# 	    }
# 	    $final_score = -100;
# 	}

	$stats_ref->{score} = ($stats_ref->{score} * ($stats_ref->{score_count} - 1) + $final_score) / $stats_ref->{score_count};
	#    $stats_ref->{score} = $final_score;

	printn "final_score = $final_score cumulative_score = $stats_ref->{score} count = $stats_ref->{score_count}";

	$genome_model_ref->set_score($stats_ref->{score});
    }
}


sub run_testcases {
    $verbosity = 1;

    $TAG = "test";
    srand(33433);

    my $config_file = "ultrasensitive.cfg";

    my $scoring_ref = Linear->new({node_ID => 98,
					   config_file => $config_file,
					   matlab_work => "test/modules",
				       });

    printn $scoring_ref->_DUMP();

    my $config_ref = {};
    read_config($config_ref, $config_file);

    use GenomeModel;
    my $genome_model_ref = GenomeModel->new({
	name => "Linear",
	Genome => {
	    radius => $config_ref->{radius},
	    kf_max => $config_ref->{kf_max},
	    kf_min => $config_ref->{kf_min},
	    kb_max => $config_ref->{kb_max},
	    kb_min => $config_ref->{kb_min},
	    kp_max => $config_ref->{kp_max},
	    kp_min => $config_ref->{kp_min},
	    Gene => {
		regulated_concentration_width => $config_ref->{regulated_concentration_width},
		unused_width => $config_ref->{gene_unused_width},
		regulated_concentration_max => $config_ref->{regulated_concentration_max},
		regulated_concentration_min => $config_ref->{regulated_concentration_min},
		Domain => {
		    RT_transition_rate_width => $config_ref->{RT_transition_rate_width},
		    TR_transition_rate_width => $config_ref->{TR_transition_rate_width},
		    RT_phi_width => $config_ref->{RT_phi_width},
		    unused_width => $config_ref->{domain_unused_width},
		    RT_transition_rate_max => $config_ref->{RT_transition_rate_max},
		    RT_transition_rate_min => $config_ref->{RT_transition_rate_min},
		    TR_transition_rate_max => $config_ref->{TR_transition_rate_max},
		    TR_transition_rate_min => $config_ref->{TR_transition_rate_min},
		    RT_phi_max => $config_ref->{RT_phi_max},
		    RT_phi_min => $config_ref->{RT_phi_min},
		    ProtoDomain => {
			binding_profile_width => $config_ref->{binding_profile_width},
			kf_profile_width => $config_ref->{kf_profile_width},
			kb_profile_width => $config_ref->{kb_profile_width},
			kp_profile_width => $config_ref->{kp_profile_width},
			Keq_profile_width => $config_ref->{Keq_profile_width},
			unused_width => $config_ref->{protodomain_unused_width},
			Keq_ratio_max => $config_ref->{Keq_ratio_max},
			Keq_ratio_min => $config_ref->{Keq_ratio_min},
		    },
		},
	    },
	},
    });

    # CONFIGURE/CREATE GENOME
    my $lg_binding_profile = $config_ref->{lg_binding_profile};
    my $tg_binding_profile = $config_ref->{tg_binding_profile};
    my $sequence_ref = $genome_model_ref->get_genome_parser_ref()->create_sequence({
	PRE_JUNK => undef, POST_JUNK => "0000",
	genes => [
	    {
		START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
		regulated_concentration => 2.0, # uM
		UNUSED => "0000",
		domains => [
		    {
			allosteric_flag => 1,
			RT_transition_rate => 0.01,
			TR_transition_rate => 1.00,
			RT_phi => 1.0,
			protodomains => [
			    {
				type => "bsite",
				substrate_polarity => 0,
				binding_profile => BindingProfile->binding_complement($lg_binding_profile)->sprint(),
				kf_profile => "11000000",
				kb_profile => "11000000",
				kp_profile => "11111",
				Keq_ratio => 1.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "11111111",
				kb_conformation_mask => "00111100",
				kp_conformation_mask => "0",
				UNUSED => "0",
			    },
			    {
				type => "bsite",
				substrate_polarity => 0,
				binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
				kf_profile => "11111111",
				kb_profile => "11000000",
				kp_profile => "11100000",
				Keq_ratio => 1.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "11111100",
				kb_conformation_mask => "11000000",
				kp_conformation_mask => "0",
				UNUSED => "0",
			    },
			   ],
			UNUSED => "0",
		    },
		   ],
	    },
	    {
		START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
		regulated_concentration => 1.0, # uM
		UNUSED => "0000",
		domains => [
		    {
			allosteric_flag => 0,
			RT_transition_rate => 1.0,
			TR_transition_rate => 1.0,
			RT_phi => 0.0,
			protodomains => [
			    {
				type => "bsite",
				substrate_polarity => 0,
				binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
				kf_profile => "11111111",
				kb_profile => "11000000",
				kp_profile => "11100000",
				Keq_ratio => 2.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "0",
				kb_conformation_mask => "0",
 				kp_conformation_mask => "0",
				UNUSED => "0",
			    },
			   ],
			UNUSED => "0",
		    },
		   ],
	    },
	   ],
    });
    printn "sequence=".$sequence_ref->get_sequence();
    $genome_model_ref->set_sequence_ref($sequence_ref);

    # save the genome object
    use Storable qw(store retrieve);
    store($genome_model_ref, "test/modules/Linear.obj");

    $scoring_ref->get_config_ref()->{plot_input} = 1;
    $scoring_ref->get_config_ref()->{plot_output} = 1;
    $scoring_ref->get_config_ref()->{plot_phase} = 1;
    $scoring_ref->get_config_ref()->{export_graphviz} = 'network,collapse_states,collapse_complexes';
    $scoring_ref->get_config_ref()->{sprint_genome} = 1;

    $scoring_ref->get_config_ref()->{LG_level} = 30;
    $scoring_ref->get_config_ref()->{LG_period} = 1000;
    $scoring_ref->get_config_ref()->{LG_strength} = 0.5;
    $scoring_ref->get_config_ref()->{LG_duty} = 95;
    $scoring_ref->get_config_ref()->{LG_rftime} = 450;
    $scoring_ref->get_config_ref()->{LG_delay} = 0;
    $scoring_ref->get_config_ref()->{LG_steps} = 10;

    $scoring_ref->get_config_ref()->{TG_init} = 70;

    $scoring_ref->score_genome($genome_model_ref);
    printn $genome_model_ref->_DUMP();
    sleep 20;
}


# Package BEGIN must return true value
return 1;

