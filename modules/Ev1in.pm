#!/usr/bin/perl -w
#-#####################################################################################
#- File:     Ev1in.pm
#- Synopsys: GA routines for 1 input network
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-#####################################################################################

#######################################################################################
# Package interface
#######################################################################################
package Ev1in;

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

use Cwd;

use Exporter;
use vars qw(@ISA @EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(
	     create_genome_sequence
	     score_genome
);

#######################################################################################
# TO-DO LIST
#######################################################################################

#######################################################################################
# EXTERNAL PACKAGES
#######################################################################################
use Storable qw(dclone);        # For deep copying of data structures
use Data::Dumper;               # Supports many formats to dump data structures into text
$Data::Dumper::Indent = 1;

#######################################################################################
# INTERNAL PACKAGES
#######################################################################################
use Utils;
use Sigtev;
use Matlab;
use Database;
#use Complex;
use Genome;
use Fields;
use CompileNetwork;
use Export;
use Network;
use Matrix;
use Protein;
use Globals;
#use BitString;
use GenAlg;

#######################################################################################
# PACKAGE GLOBALS
#######################################################################################
use vars qw($initial_sequence);   # initial sequence
use vars qw($MATLAB_RD $MATLAB_WR);  # matlab i/o handles
use vars qw($figure_num);  # matlab i/o handles

######################################################################################
# Function: create_test_ligand
# Synopsys: Create simple proteins to use as input ligands for testcases.
#######################################################################################
sub create_test_ligand {
    my ($ligand, $profile) = @_;

    create_simple_protein($ligand, 1.0, {
	bsite_profile => $profile,
	bs_control_tag => 0x0,
	cs_control_tag => 0x0,
	ms_control_tag => 0x0,
	ms_R => 0.11,
	bs_R => -1,  # make sure nothing can modulate
	cs_R => -1,  # make sure nothing can modulate
	protodomain_type => "bsite",
	ligand_polarity => 0,
#	csite_type => "null",
#	msite_type => "null",
	csite_default_state => "off",
	bsite_default_state => "on",
	kf_profile => 0x00,
	kb_profile => 0x00,
	kp_rate    => 0x00,
	kd_rate    => 0x00,
	site_type_conflict_resolution => 0,
    });
}

######################################################################################
# Function: create_test_tfactor
# Synopsys: Create simple proteins to use as output transcription factor.
#######################################################################################
sub create_test_tfactor {
    my ($tf, $profile) = @_;

    create_simple_protein($tf, 1.0, {
	bsite_profile => $profile,

	bs_control_tag => 0x0,
	cs_control_tag => 0x0,
	ms_control_tag => 0xF,
	ms_R => 0.1,
	bs_R => -1,  # make sure nothing can modulate
	cs_R => -1,  # make sure nothing can modulate

#	bs_control_tag => 0,
#	cs_control_tag => 0,
#	ms_control_tag => 0,
#	ms_R => 0,
#	bs_R => 1,  # make sure nothing can modulate
#	cs_R => 1,  # make sure nothing can modulate


	site_type_conflict_resolution => 0,
	protodomain_type => "msite",
	ligand_polarity => 0,
#	csite_type => "null",
#	msite_type => "phosphorylation",
	csite_default_state => "off",
	bsite_default_state => "on",
	kf_profile => 0x00,
	kb_profile => 0x00,
	kp_rate    => 0x00,
	kd_rate    => 0x00,
    });
}

#######################################################################################
# Function: node_init
# Synopsys: This routine is called when scoring nodes are created to do initialization
#######################################################################################
sub node_init {
    $genalg_config{config_file} = shift;
    read_config(\%genalg_config, $genalg_config{config_file}, "NOCLOBBER");
}


#######################################################################################
# Function: create_genome_sequence
# Synopsys: create_genome_sequence
#######################################################################################
sub create_genome_sequence {

    # given concentration_scaling_factor=1e-6, then
    # concentration_width=4  => max no. of molecules in 1e-18L is 15e-6   * 1e-18 * 6.022e23 ~ 9
    # concentration_width=10 => max no. of molecules in 1e-18L is 1023e-6 * 1e-18 * 6.022e23 ~ 616

    # given kinetic_profile_width = 20
    #   forward_rate_scaling_factor = 1e3  => forward rates can go as fast as 2e9 M^-1 s^-1   TOO FAST!!
    #   forward_rate_scaling_factor = 10   => forward rates can go as fast as 20e6 M^-1 s^-1

    set_field_widths({
	bsite_profile_width => 8,
	control_tag_width => 4,
	modR_seq_width => 2,
	kinetic_profile_width => 20,            # 1e6 dynamic range
	diffusion_rate_width => 8,
	concentration_width => 10,              # regulated conc. range from 0 to ~1000uM
	forward_rate_scaling_factor => 1e3,     # forward rate range is 1e3 -> 2e9 M^-1 s^-1
	reverse_rate_scaling_factor => 1e-5,    # reverse rate range is 10uHz -> 20Hz
	product_rate_scaling_factor => 1e-3,    # product rate range is 1 mHz -> 2kHz
	concentration_scaling_factor => 1e-6,   # all concentrations in micro-mol/L
	cell_volume => $genalg_config{cell_volume},  # cell volume in L
	kb_release => 10,                       # release reactions occur at 10Hz
    });

    parser_init();

    if ($genalg_config{initial_genome} =~ /random/) {
	generate_random_genome(5000);
    } else {
	printn "ERROR: initial genome must be random";
	exit(1);
    }
}

#######################################################################################
# Function: score_genome
# Synopsys: score_genome
#######################################################################################
sub score_genome {
    my $node_id = shift;
    my $generation = shift;;
    my $individual = shift;

    my $work_dir = cwd();
    $work_dir .= "/ev1in/$TAG/matlab";

    if (!defined $MATLAB_RD) {  # need to init matlab?
	($MATLAB_RD, $MATLAB_WR) = matlab_init_process(
	    LOGFILE => sprintf("ev1in/$TAG/matlab/ev1in.matlab.$node_id.log",
			       $generation,
			       $individual),
	   );
	matlab_printn ($MATLAB_WR, "cd $work_dir"); # change working dir
    }


    printn "score_genome: scoring generation=$generation, individual=$individual";

    if ($sdb->{genome}{history}->[$#{$sdb->{genome}{history}}] =~ /elite dup/) {
	printn "score_genome: elite individual, previous score=$sdb->{genome}{score}";
	return if ($genalg_config{rescore_elite} == 0);
	printn "score_genome: re-scoring elite individual";
	$sdb->{genome}{score_count}++;
    } else {
	printn "score_genome: mutated individual, scoring...";
	$sdb->{genome}{score_count} = 1;
	$sdb->{genome}{score} = 0;
    }
    $sdb->{genalg}{score} = $sdb->{genome}{score};  # save score because reset_database() deletes it

    # SCORE THE NETWORK
    reset_database(); parser_init(); parse_genome(); report_proteins();

#	my $lg_profile = create_random_bsite_profile($bsite_profile_width, 4);
#	my $tf_profile = create_random_bsite_profile($bsite_profile_width, 4);
    my $lg_profile = 0x03;
    my $tf_profile = 0x80;

    create_test_ligand("LG", $lg_profile);
    create_test_tfactor("TF", $tf_profile);

    compute_adjacency_matrix();

    export_ppi_sif(sprintf("ev1in/$TAG/sif/ppi.%04d.%04d.sif",$generation, $individual));
    prune_isolated_proteins(); compute_adjacency_matrix();

    my $connectivity_score = 0;
    my @new_lg_subnet;
    my @new_tf_subnet;
    my $new_lg_subnet_size = 0;
    my $new_tf_subnet_size = 0;

    my $num_adjacent_kinases = 0;
    my $num_adjacent_phosphatases = 0;
    my @adjacent_protein_list = ();
    my @adjacent_protodomain_list = ();
    my $new_protodomain_sequence;

    # 1st order scoring -- LG/TF connected to anything
    $connectivity_score++ if (grep(/LG/,@{$sdb->{network}{adjacency_matrix_protein_list}}));
    $connectivity_score++ if (grep(/TF/,@{$sdb->{network}{adjacency_matrix_protein_list}}));

    if ($connectivity_score == 2) {  # LG/TF connected
	# score -- size of subnets connected to LGx
	compute_high_order_adjacency_matrix(8);
	compute_connectivity_matrix();

#	printn "max subnet size is:" . max_numeric(@{msumrow($sdb->{network}{protein_connectivity_matrix})});

	@new_lg_subnet = network_get_connected_protein("LG");
	@new_tf_subnet = network_get_connected_protein("TF");

	# max 90 points for subnet size
	$new_lg_subnet_size = (@new_lg_subnet > 30) ? 45 : @new_lg_subnet;
	$new_tf_subnet_size = (@new_tf_subnet > 30) ? 45 : @new_tf_subnet;
	$connectivity_score += ($new_lg_subnet_size + $new_tf_subnet_size);

	# score -- LG/TF connected to each other
	if (grep /LG/, @new_tf_subnet) {
	    printn "LG and TF are connected to each other";
	    $connectivity_score += 100;
	}
	if (grep /TF/, @new_lg_subnet) {
	    printn "TF and LG are connected to each other";
	    $connectivity_score += 100;
	}

	# TF must have at least two connected proteins
	@adjacent_protein_list = network_get_adjacent_protein("TF");
	if (@adjacent_protein_list >= 2) {
	    printn "TF has > 2 connected proteins";
	    $connectivity_score += 100;
	}
    }

    # ONCE CONNECTED, MAKE DEFAULT STATE OF ALL PROTODOMAINS "ON", AND 
    # PICK AN ADJACENT PROTEIN TO MAKE INTO KINASE/PHOSPHATASE
    if ($connectivity_score >= 300) {
	if (!defined $sdb->{genome}{forced_default_on}) {
	    printn "editing genome to turn default states on";
	    # loop thru all protodomains and turn them on
	    foreach my $protodomain (keys %{$sdb->{protodomain_table}}) {
		next if ($protodomain =~ /TF/);
		next if ($protodomain =~ /LG/);
		$sdb->{protodomain_table}{$protodomain}{csite_default_state} = "on";
		$sdb->{protodomain_table}{$protodomain}{bsite_default_state} = "on";
		$new_protodomain_sequence = unparse_protodomain($protodomain);
 		replace_protodomain_sequence($protodomain, $new_protodomain_sequence);
	    }
	    push @{$sdb->{genome}{history}}, "G$generation I$individual edited genome to force defaults on";
	    $sdb->{genome}{forced_default_on} = 1;
	}

	$num_adjacent_kinases = find_adjacent_enzymes("TF",0);
	$num_adjacent_phosphatases = find_adjacent_enzymes("TF",1);

	@adjacent_protodomain_list = network_get_adjacent_protodomain("TF");
	if ($num_adjacent_kinases < 1 && !defined $sdb->{genome}{created_kinase}) {
	    printn "editing genome to create kinase";
	    for (my $i = 0; $i < 10; $i++) {
		my $index = int rand @adjacent_protodomain_list;
		my $protodomain = $adjacent_protodomain_list[$index];
		if ($sdb->{protodomain_table}{$protodomain}{protodomain_type} ne "csite") {
		    $sdb->{protodomain_table}{$protodomain}{protodomain_type} = "csite";
		    $sdb->{protodomain_table}{$protodomain}{ligand_polarity} = 0;
		    # unparse back to genome
		    $new_protodomain_sequence = unparse_protodomain($protodomain);
		    replace_protodomain_sequence($protodomain, $new_protodomain_sequence);
		    $num_adjacent_kinases++;
		    $sdb->{genome}{created_kinase} = 1;
		    push @{$sdb->{genome}{history}}, "G$generation I$individual edited genome to create kinase";
		    last;
		}
	    }
	}
	if ($num_adjacent_phosphatases < 1 && !defined $sdb->{genome}{created_phosphatase}) {
	    printn "editing genome to create phosphatase";
	    for (my $i = 0; $i < 10; $i++) {
		my $index = int rand @adjacent_protodomain_list;
		my $protodomain = $adjacent_protodomain_list[$index];
		if ($sdb->{protodomain_table}{$protodomain}{protodomain_type} ne "csite") {
		    $sdb->{protodomain_table}{$protodomain}{protodomain_type} = "csite";
		    $sdb->{protodomain_table}{$protodomain}{ligand_polarity} = 1;
		    # unparse back to genome
		    $new_protodomain_sequence = unparse_protodomain($protodomain);
		    replace_protodomain_sequence($protodomain, $new_protodomain_sequence);
		    $num_adjacent_phosphatases++;
		    $sdb->{genome}{created_phosphatase} = 1;
		    push @{$sdb->{genome}{history}}, "G$generation I$individual edited genome to create phosphatase";
		    last;
		}
	    }
	}

	# max out at one
	$connectivity_score += 100 * ($num_adjacent_kinases > 1 ? 1 : $num_adjacent_kinases);
	$connectivity_score += 100 * ($num_adjacent_phosphatases > 1 ? 1 : $num_adjacent_phosphatases);

	# report current score
	printn "lg connected to: " . join ",", @new_lg_subnet;
	printn "tf connected to: " . join ",", @new_tf_subnet;
	printn join ",", sort @{$sdb->{network}{adjacency_matrix_protein_list}};
	print_2d_array($sdb->{network}{protein_connectivity_matrix});
    }

    my $num_reactions_tf_0 = 0;
    my $num_reactions_tf_1 = 0;
    # CHECK REACTION TABLE TO MAKE SURE THERE IS AN EQN THAT LEADS TO DESIRED PRODUCT
    # must be fully connected, and adjacent to both kinase and phosphatase
    if ($connectivity_score >= 500 && $num_adjacent_kinases !=0 && $num_adjacent_phosphatases != 0) {
	# Before first sim, remove all genes not connected to input subnet, init matlab
	my @proteins_not_in_subnet = difference([get_parsed_genes()], [union(\@new_lg_subnet, \@new_tf_subnet)]);
#	    delete_gene(@proteins_not_in_subnet);

	printn "compiling the reaction network...";
	reset_database(); parse_genome();
	# recompute since may have deleted genes
	@proteins_not_in_subnet = difference([get_parsed_genes()], [union(\@new_lg_subnet, \@new_tf_subnet)]);
	create_test_ligand("LG", $lg_profile); create_test_tfactor("TF", $tf_profile);

    # *****  NETWORK PERTURBATION HERE *****
#    printn "PERTURBATION (BEFORE):  ".$sdb->{protodomain_table}{xxx}{yyy};
#	delete $sdb->{protein_table}{"P3946"};
#    printn "PERTURBATION (AFTER):  ".$sdb->{protodomain_table}{xxx}{yyy};
    # **************************************

	compute_adjacency_matrix(); prune_proteins(@proteins_not_in_subnet);
	compile_ligand_lists();
	compile_complex_table(MAX_ITERATIONS => -1,
			      MAX_COMPLEXES => 400,
			      MAX_COMPLEX_SIZE => 3,
			      MAX_DUPS => -1,
			      MAX_INTERNAL_ITERATIONS => -1,
			     );
	compile_reaction_table(COLOCALIZED_CONC => 1);

	# now check that TF_0 and TF_1 are product of at least 2 reactions
	$num_reactions_tf_1 = grep (/ ->.* TF:1/, (
	    keys %{$sdb->{reaction_table}{catalytic}},
	    keys %{$sdb->{reaction_table}{binding}},
	    keys %{$sdb->{reaction_table}{release}},
	   ));
	$num_reactions_tf_0 = grep (/ ->.* TF:0/, (
	    keys %{$sdb->{reaction_table}{catalytic}},
	    keys %{$sdb->{reaction_table}{binding}},
	    keys %{$sdb->{reaction_table}{release}},
	   ));
	$connectivity_score += 100 * ($num_reactions_tf_1 > 1 ? 1 : $num_reactions_tf_1);
	$connectivity_score += 100 * ($num_reactions_tf_0 > 1 ? 1 : $num_reactions_tf_0);
    }
    printn "($node_id, $generation, $individual) connectivity_score=$connectivity_score, new_lg_subnet_size=$new_lg_subnet_size new_lg_subnet_size=$new_lg_subnet_size new_tf_subnet_size=$new_tf_subnet_size, num_adjacent_kinases=$num_adjacent_kinases, num_adjacent_phosphatases=$num_adjacent_phosphatases, num_reactions_tf_1=$num_reactions_tf_1, num_reactions_tf_0=$num_reactions_tf_0";
	
    # IF CONNECTIVITY CONDITIONS MET, SIMULATE NETWORK AND SCORE
    my $sim_score = 0;
#    my $cos_angle = 0;
    my $correlation = 0;
    my $magnitude = 0;
    my $variance = 0;
    my $stdev = 0;
    my $steady_state_slope = 0;
    my $steady_state_score = 0;
    my $mean_squared_err = 0;
    my $max_mean_squared_err = 0;
    my $mean_squared_err_score = 0;
    my $complexity = 0;
    my $complexity_score = 0;
    
    my ($cue, @problems);

    if ($connectivity_score > 700 && $num_reactions_tf_0 != 0 && $num_reactions_tf_1 != 0) {
#    if (1) {
	my $reaction_sif_file = sprintf("ev1in/$TAG/sif/reaction.%04d.%04d.sif",$generation,$individual);
	export_reaction_sif($reaction_sif_file, "MORE_DETAIL");
	compute_adjacency_matrix();
#	export_ppi_sif(sprintf("ev1in/$TAG/sif/ppi.%04d.%04d.sif",$generation, $individual));

	my $output_tf = "TF_1";  # phosphorylated version of output TF

	my $staircase_period = 800;
#	my $staircase_strength = 0.2;  # 1s tau
	my $staircase_strength = 1.0;  # 1s tau
	my $staircase_level = $genalg_config{LG_range};
	my $staircase_duty = 75;
	my $staircase_rftime = 250;
	my $staircase_delay = 100;
	my $staircase_steps = 5;

	$sdb->{protein_table}{LG}{input_concentration} = 0.1 * $genalg_config{LG_range};   # in mol/L
	$sdb->{protein_table}{TF}{input_concentration} = $genalg_config{TF_range};   # in mol/L

	my ($lg_source_staircase, $lg_sink_staircase, $event_list) = matlab_staircase_equation(
#	my ($lg_source_staircase, $lg_sink_staircase, $event_list) = matlab_ramp_equation(
	    NODE => "LG_0",
	    PERIOD => $staircase_period,
	    STRENGTH => $staircase_strength,
	    CONCENTRATION => $staircase_level,
	    DUTY => $staircase_duty,
	    RFTIME => $staircase_rftime,
	    STEPS => $staircase_steps,
	    DELAY => $staircase_delay,
	   );

	my @event_list = split " ", $event_list;

	my $eqn_file = sprintf("ev1in/$TAG/matlab/ev1in_%04d_%04d.eqn",$generation,$individual);
	export_eqn($eqn_file, $lg_source_staircase, $lg_sink_staircase);

	my $solver = $genalg_config{solver};

	# multiplying by conversion factor converts no. molecules to mol/L
	my $cfac = ($solver eq "stoch") ? 1.0/($cell_volume * 6.022e23) : 1.0;

	printn "simulating network...";
	matlab_eqn_convert(
	    ROOTNAME => sprintf("ev1in_%04d_%04d", $generation, $individual),
	    EQN_FILE => sprintf("ev1in/$TAG/matlab/ev1in_%04d_%04d.eqn",$generation,$individual),
	    WORK_DIR => $work_dir,
	    SOLVER => $solver,
	    SOLVER_OPTIONS => "odeset('InitialStep', 1e-8, 'AbsTol', 1e-15, 'RelTol', 1e-3, 'MaxStep', 1.0)",
	    T_FINAL => $staircase_period,
	    T_VECTOR => "[t0:0.2:tf]",
	    T_TICK => 10,
	    T_EVENTS => $event_list,
	   );

	# compute a complexity score based on eqn file line count
	my $line_count = `wc -l $eqn_file.all`;
	$line_count =~ /^\s*(\S+)\s*/;   # extract the line count
	$complexity = $1;
	$complexity_score = 1/($complexity + 1);

	##############################
	# run simulation
	##############################
	my $rootname = sprintf("ev1in_%04d_%04d",$generation,$individual);
	matlab_run_sim(MATLAB_WR => $MATLAB_WR,
		       MATLAB_RD => $MATLAB_RD,
		       ROOTNAME => $rootname,
		       WORK_DIR => $work_dir,
		       SOLVER => $solver,
		       T_FINAL => $staircase_period,
		       PLOT_MIN_VALUE => $genalg_config{plot_min},
		       PLOT_FIGURE => 100,
#		       STOCH_COUT => 1,
		       STOCH_OUT_INTERVAL => 0.2,
#		       STOCH_TMIN => 30.0,
#		       STOCH_SFILTER => 688,
		      );
	##############################

	$figure_num = 1 if (!defined $figure_num);
	if (defined $genalg_config{plot_input} && $genalg_config{plot_input}) {
	    # input is always the same, so plot first time only
	    matlab_plot_complex($MATLAB_WR, 1, "LG_0") if ($figure_num == 1);
	}
	if (defined $genalg_config{plot_output} && $genalg_config{plot_output}) {
	    matlab_plot_complex($MATLAB_WR, 2, "TF_0", "GENOME $generation.$individual");
	    matlab_plot_complex($MATLAB_WR, 3, "TF_1", "GENOME $generation.$individual");
#	    matlab_plot_all_complexes($MATLAB_WR, $MATLAB_RD, 0.0, "stairs", 100);
	}
	matlab_printn ($MATLAB_WR, "disp('Done plotting')\n");
	matlab_wait_on($MATLAB_RD, "Done plotting");
	sleep 1;  # this gives X some time to finish plotting

	# score for steady state
	printn "computing steady-state slope...";
	my @steady_state_slope;
	my @steady_state_event_list = (@event_list, $staircase_period);
	for (my $i = 0; $i < @steady_state_event_list; $i++) {
	    my $steady_state_event = $steady_state_event_list[$i];
	    my $state_var = ($solver eq "stoch") ? "simdata" : "y";
	    my @state_vector_delta = matlab_get_state_delta($MATLAB_WR, $MATLAB_RD, $steady_state_event - 1, $steady_state_event, $state_var);
	    $steady_state_slope[$i] = $cfac * max_numeric(abs(max_numeric(@state_vector_delta)), abs(min_numeric(@state_vector_delta)));
	    printn "slope @ $steady_state_event = ".sprintf("%.5e",$steady_state_slope[$i]);
	}
	$steady_state_slope = max_numeric(@steady_state_slope);
	# MATLAB PLOT hill func: ss_th = -log10(1e-5); ss=[-20:0.01:0]; s=10.^(ss); k = -log10(s)/ss_th; n=4; kn=k.^n; semilogx(s,kn./(1+kn))
	my $steady_state_threshold = $genalg_config{steady_state_threshold};
	my $steady_state_hill_transition = -log_10($steady_state_threshold);  # this is -log10(1e-3/600) i.e. 1 molecule in 1e-18 L
	$steady_state_score = p_hill(-log_10($steady_state_slope), $steady_state_hill_transition, 4);
#	$steady_state_score = ($steady_state_score > 0.5) ? 1.0 : $steady_state_score;

	printn "RESULT VECTOR: INPUT = LG OUTPUT = ".simplify($output_tf);
	my (@input_vector, @output_vector, @expected_output_vector);

	# compute sample list based on staircase equation events (even if actually using ramp)
	my $sample_list = matlab_staircase_sample_vector(
	    PERIOD => $staircase_period,
	    DUTY => $staircase_duty,
	    RFTIME => $staircase_rftime,
	    STEPS => $staircase_steps,
	    DELAY => $staircase_delay,
	   );
	my @sample_list = split " ", "$sample_list $staircase_period";

	for (my $i=0; $i < @sample_list; $i++) {
	    my $t = $sample_list[$i]; 
	    $input_vector[0][$i] = matlab_get_value($MATLAB_WR, $MATLAB_RD, "LG_0", $t);
	    $output_vector[0][$i] = matlab_get_value($MATLAB_WR, $MATLAB_RD, "TF_1", $t);
	    $expected_output_vector[0][$i] = ($cfac * $input_vector[0][$i] > 0.5*$genalg_config{LG_range}) ? $genalg_config{TF_range}/$cfac : 0;   # thresholded output at 50% of input
	    printf("RESULT VECTOR:  t=%-6.2f input vector:  %4.0g output_vector: %12.8g (expected = $expected_output_vector[0][$i])\n",
		    $t, $input_vector[0][$i], $output_vector[0][$i]);
	}

	if (defined $genalg_config{plot_phase} && $genalg_config{plot_phase}) {
	    my $matlab_command = "figure(4); plot(LG_0(1:end/2), TF_1(1:end/2));title(\'PHASE PLOT\')";
	    matlab_printn($MATLAB_WR, "$matlab_command");
	    $matlab_command = "hold on; plot(LG_0(end/2+1:end), TF_1(end/2+1:end), 'r');";
	    matlab_printn($MATLAB_WR, "$matlab_command");
	}

	if ("@{$output_vector[0]}" !~ /UNDEFINED/) {
	    printn "computing correlation...";
#	    $cos_angle = norm_projection(@{$output_vector[0]}, @{$expected_output_vector[0]});
#	    if (!defined $cos_angle) {$cos_angle = 0;};  # norm_projection can return undef
	    $correlation = correlation(@{$output_vector[0]}, @{$expected_output_vector[0]});
	    $magnitude = magnitude(@{$output_vector[0]});
	    $variance = variance(@{$output_vector[0]});
	    $stdev = $variance ** (0.5);
	}

	$mean_squared_err = mean_squared_error($output_vector[0], $expected_output_vector[0]);
	$max_mean_squared_err = ($staircase_level/$cfac) ** 2;  # since it's the maximum MEAN error^2
	$mean_squared_err_score  = (1-$mean_squared_err/$max_mean_squared_err);

	# check if any concentrations are negative
	foreach my $sample (@{$output_vector[0]}) {
	    if ($sample < 0) {
		printn "NOTE: detected negative concentrations in output vector";
		$mean_squared_err_score = 0;
		$correlation = 0;
	    }
	}

	if (($mean_squared_err_score < 0) || ($mean_squared_err_score > 1)) {
	    # numerics got messed up, set score to zero
	    printn "NOTE: computed mean_squared_error_score out of bounds, setting to zero";
	    $correlation = 0;
	    $mean_squared_err_score = 0;
	}


	my $scaled_stdev = $stdev/1e-3;

#	$sim_score =  0.1 * $steady_state_score;
#	$sim_score += 0.9 * $mean_squared_err_score if ($steady_state_score == 1);
	$sim_score = $mean_squared_err_score;

    }

    $sdb->{genome}{stats}{steady_state_slope} = $steady_state_slope;
    $sdb->{genome}{stats}{steady_state_score} = $steady_state_score;
    $sdb->{genome}{stats}{sim_score} = $sim_score;
    $sdb->{genome}{stats}{correlation} = $correlation;
    $sdb->{genome}{stats}{variance} = $variance;
    $sdb->{genome}{stats}{stdev} = $stdev;
    $sdb->{genome}{stats}{magnitude} = $magnitude;
    $sdb->{genome}{stats}{connectivity_score} = $connectivity_score;
    $sdb->{genome}{stats}{max_mean_squared_err} = $max_mean_squared_err;
    $sdb->{genome}{stats}{mean_squared_err} = $mean_squared_err;
    $sdb->{genome}{stats}{mean_squared_err_score} = $mean_squared_err_score;
    $sdb->{genome}{stats}{complexity} = $complexity;
    $sdb->{genome}{stats}{complexity_score} = $complexity_score;
    printn "mean_squared_err_score = $mean_squared_err_score, mean_squared_err = $mean_squared_err, max_mean_squared_err = $max_mean_squared_err, steady_state_slope = $steady_state_slope, steady_state_score = $steady_state_score, sim_score = $sim_score, correlation=$correlation, variance=$variance, stdev=$stdev, magnitude=$magnitude, complexity=$complexity, complexity_score=$complexity_score";

    my $final_score = 0;

    $final_score += ($sim_score > 0) ? 0.1 : 0.1 * (1-exp(-$connectivity_score/300));  	# get full points for connectivity if sim_score > 0
    $final_score += 0.8 * $sim_score;
    $final_score += ($sim_score > 0) ? 0.1 * $complexity_score : 0.0;   # complexity counts whenever sim works

    $final_score = ($final_score < 0) ? 0 : $final_score;  # prevent neg've scores

    if (@problems) {
	printn "Warning: matlab had some problems, tossing out this solution....";
	for (my $i=0; $i < @problems; $i++) {
	    printn "Warning: (from matlab_wait_on) ==> $problems[$i]";
	    last if $i > 10;
	}
	$final_score = -100;
    }

    $sdb->{genome}{score} = ($sdb->{genalg}{score} * ($sdb->{genome}{score_count} - 1) + $final_score) / $sdb->{genome}{score_count};
#    $sdb->{genome}{score} = $final_score;

    printn "final_score = $final_score cumulative_score = $sdb->{genome}{score} count = $sdb->{genome}{score_count}";
}

#######################################################################################
# Package initialization
#######################################################################################
printn "Package Ev1in loaded";

# Package BEGIN must return true value
return 1;
