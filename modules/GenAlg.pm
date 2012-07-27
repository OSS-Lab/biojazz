#-#####################################################################################
#- File:     GenAlg.pm
#- Synopsys: 
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package GenAlg;
use Class::Std;
use base qw();
{
    use Carp;

    use Utils;

    use Globals qw ($verbosity $TAG);

    use ScorCluster;
    use Generation;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %config_ref_of :ATTR(get => 'config_ref', set => 'config_ref', init_arg => 'config_ref');
    my %cluster_ref_of :ATTR(get => 'cluster_ref', set => 'cluster_ref');

    my %current_generation_ref_of :ATTR(get => 'current_generation_ref', set => 'current_generation_ref');
    my %current_generation_number_of :ATTR(get => 'current_generation_number', set => 'current_generation_number');

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
    # Function: START
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub START {
        my ($self, $obj_ID, $arg_ref) = @_;

	my $config_ref = $config_ref_of{$obj_ID};

	# append scoring_class to local_dir to prevent collisions
	my $local_dir = (defined $config_ref->{local_dir} ?
			 "$config_ref->{local_dir}/$config_ref->{scoring_class}" :
			 undef);

	my $host_list_ref = ref $config_ref->{host_list} ? $config_ref->{host_list} : [$config_ref->{host_list}];

	my $cluster_ref = $cluster_ref_of{$obj_ID} = ScorCluster->new({
	    config_file => $config_ref->{config_file},
	    cluster_type => $config_ref->{cluster_type},
	    cluster_size => $config_ref->{cluster_size},
	    host_list => $host_list_ref,
	    nice => $config_ref->{nice},
	    work_dir => $config_ref->{work_dir},
	    local_dir => $local_dir,
	    scoring_class => $config_ref->{scoring_class},
	});

	$cluster_ref->spawn_rrobin();

	# check initializers
	# ...
    }

    #--------------------------------------------------------------------------------------
    # Function: load_current_generation
    # Synopsys: Load generation (i) from disk into current.
    #--------------------------------------------------------------------------------------
    sub load_current_generation {
	my $self = shift; my $obj_ID = ident $self;
	my $number = shift;

	my $config_ref = $config_ref_of{$obj_ID};
	my $current_generation_ref = $current_generation_ref_of{$obj_ID};

	$current_generation_ref->clear_genomes();
	$current_generation_ref->load_generation(
	    dir => "$config_ref->{work_dir}/$TAG/obj",
	    number => $number,
	   );
	$current_generation_number_of{$obj_ID} = $number;
    }

    #--------------------------------------------------------------------------------------
    # Function: save_current_generation
    # Synopsys: Save the current generation to disk.
    #--------------------------------------------------------------------------------------
    sub save_current_generation {
	my $self = shift; my $obj_ID = ident $self;

	my $config_ref = $config_ref_of{$obj_ID};
	my $current_generation_ref = $current_generation_ref_of{$obj_ID};

	$current_generation_ref->save_generation(
	    dir => "$config_ref->{work_dir}/$TAG/obj",
	    number => $current_generation_number_of{$obj_ID},
	   );
    }

    #--------------------------------------------------------------------------------------
    # Function: create_initial_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub create_initial_generation {
	my $self = shift; my $obj_ID = ident $self;

	printn "create_initial_generation: starting...";

	my $config_ref = $config_ref_of{$obj_ID};

	my $current_generation_ref = $current_generation_ref_of{$obj_ID} = Generation->new({});

	if ($config_ref->{initial_genome} =~ /load\s+(\S+)/) {
	    my $file_glob = $1;
	    my @files = glob $file_glob;

	    my $save_tag = $TAG;
	    if ($file_glob !~ /$save_tag/) {
		printn "\nWARNING: you appear to be loading from a different run of BioJazz";
		printn "Continue (y/n)? ";
		my $answer = <>; chomp($answer);
		if ($answer ne "y") {
		    printn "no, then fix the configuration file";
		    exit(1);
		}
	    }
	    if (($file_glob =~ /G(...)_I/) && ($1 != $config_ref->{first_generation})) {
		printn "\nWARNING: Configuration parameter first_generation appears not to be set correctly";
		printn "Continue (y/n)? ";
		my $answer = <>; chomp($answer);
		if ($answer ne "y") {
		    printn "no, then fix the configuration file";
		    exit(1);
		}
	    }		
	    if (!@files) {
		printn "ERROR: nothing to load... ";
		exit(1);
	    }

	    # need to sort on individual number
#	    @files = sort {$a=~/_I(\d+)/; $a_i=$1; $b=~/_I(\d+)/; $b_i=$1; return $a_i <=> $b_i} @files;
	    printn "create_initial_generation: loading initial generation from disk";
	    printn join "\n", @files;
	    $current_generation_ref->retrieve_genomes(
		files => \@files,
		history_flag => 1,
	       );
	} else {
	    printn "create_initial_generation: creating $config_ref->{inum_genomes} individuals";
	    $current_generation_ref->create_random_genomes($config_ref);
	}
	$current_generation_ref->refresh_individual_names($config_ref->{first_generation});
	$current_generation_number_of{$obj_ID} = $config_ref->{first_generation};
	printn "create_initial_generation: done";
    }

    #--------------------------------------------------------------------------------------
    # Function: create_next_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub create_next_generation {
	my $self = shift; my $obj_ID = ident $self;

	my $config_ref = $config_ref_of{$obj_ID};

	my $current_generation_ref = $current_generation_ref_of{$obj_ID};
	my $current_generation_number = $current_generation_number_of{$obj_ID};
	my $current_generation_size = $current_generation_ref->get_num_elements();

	my $next_generation_ref = Generation->new({});
	my $next_generation_number = $current_generation_number + 1;

	printn "create_next_generation: creating generation $next_generation_number";

	# check scores to make sure defined and positive
	my @scores = map {$current_generation_ref->get_element($_)->get_score()} (0..$current_generation_size-1);
	if (grep {!defined $_} @scores) {
	    printn "ERROR: not all scores are defined (if this is first generation, check value of score_initial_generation in config file)";
	    exit(1);
	}
	if (grep {$_ < 0} @scores) {
	    printn "ERROR: all scores must be non-negative";
	    exit(1);
	}

	# rank the current population by sorting on score
	my $ranked_genomes_ref = $current_generation_ref->get_ranked_genomes();

	my @ranked_indices = @{$ranked_genomes_ref->{ranked_indices}};
	my @ranked_genomes = @{$ranked_genomes_ref->{ranked_genomes}};
	my @ranking = @{$ranked_genomes_ref->{ranking}};

	printn "ranked indices:  ".join ",", @ranked_indices;
	printn "ranking:         ".join ",", @ranking;

	# propagate the top genomes unchanged to next generation (elitist strategy)
	my $num_elite = (@ranked_genomes >= $config_ref->{elite_pool_size}) ? $config_ref->{elite_pool_size} : @ranked_genomes;

	for (my $i = 0; $i < $num_elite; $i++) {
	    my $elite_ref = $ranked_genomes[$i];
	    my $elite_name = $elite_ref->get_name();
	    printn "create_next_generation: creating child from elite genome $elite_name";
	    my $duplicate_ref = $elite_ref->duplicate();
	    $duplicate_ref->add_history(sprintf("ELITE REPLICATION: $elite_name -> G%03d_I%02d", $next_generation_number, $i));
	    $duplicate_ref->set_elite_flag(1);
	    $next_generation_ref->add_element($duplicate_ref);
	}

	my $total_children = $config_ref->{max_population} - $num_elite;
	printn "create_next_generation: total_children = $total_children";

	# compute rank-based scores
	my @ranking_scores = ();
	if (!defined $config_ref->{score_by_rank_flag} || $config_ref->{score_by_rank_flag}) {
	    printn "create_next_generation: computing parent ranking_scores";
	    my $ranking_nonviable = $config_ref->{ranking_nonviable};
	    my $num_nonviable = int $ranking_nonviable * $current_generation_size;
	    my $rank_nonviable = $current_generation_size - $num_nonviable; # inclusive
	    my $p1 = ($config_ref->{ranking_pressure})**(1/5);  # per-centile fitness fold-change
	    my $pr = $p1**(100/$current_generation_size);       # per-rank fitness fold-change

	    for (my $i = 0; $i < $current_generation_size; $i++) {
		if ($ranking[$i] >= $rank_nonviable) {
		    $ranking_scores[$i] = 0;
		} else {
		    $ranking_scores[$i] = $pr ** ($rank_nonviable - $ranking[$i] - 1);
		}
	    }
	}

	my @fitness_scores = ($config_ref->{score_by_rank_flag}) ? @ranking_scores : @scores;
	for (my $i = 0; $i < $current_generation_size; $i++) {
	    printn "fitness_score($i) = $fitness_scores[$i]";
	}

	my $total_score = 0;
	for (my $i = 0; $i < $current_generation_size; $i++) {
	    $total_score += $fitness_scores[$i];
	}

	printn "create_next_generation: total_score = $total_score";

	# compute scale_factor and rescale scores
	my $scale_factor;
	$scale_factor = $total_children / $total_score;
	my @scaled_fitness_scores = map {$_ * $scale_factor} @fitness_scores;

	# guiding principle for the following is that no. children is proportional to fitness score
	# (use the roulette wheel sampling approach)
	my $children_created = 0;
	my $parent_pointer = 0;
	my $child_pointer = rand(1);
	for (my $i = 0; $i < $current_generation_size; $i++) {
	    my $num_children = 0;

	    $parent_pointer += $scaled_fitness_scores[$i];
	    while ($parent_pointer > $child_pointer) {
		$child_pointer++;
		$num_children++;
	    }

	    my $parent_ref = $current_generation_ref->get_element($i);

	    printn "create_next_generation: individual $i, ranking=$ranking[$i], fitness_score = $fitness_scores[$i], scaled_fitness_score=$scaled_fitness_scores[$i], num_children=$num_children";
	    for (my $j = 0; $j < $num_children; $j++) {
		my $parent_name = $parent_ref->get_name();
		printn "create_next_generation: creating child $j of individual $parent_name";
		my $child_ref = $parent_ref->duplicate();
		$child_ref->add_history(sprintf("REPLICATION: $parent_name -> G%03d_I%02d", $next_generation_number, $num_elite+$children_created+$j));
		$child_ref->set_elite_flag(0);
		$child_ref->mutate(
		    prob_mutate_params => $config_ref->{prob_mutate_params},
		    prob_mutate_global => $config_ref->{prob_mutate_global},
		    prob_recombination => $config_ref->{prob_recombination},
		    prob_duplicate => $config_ref->{prob_duplicate},
		    prob_delete => $config_ref->{prob_delete},
		    mutation_rate => $config_ref->{mutation_rate},
		   );
		$child_ref->set_score(undef);
		$child_ref->clear_stats();
		$next_generation_ref->add_element($child_ref);
	    }
	    $children_created += $num_children;
	}

	# number of children created may be different from target due to rounding
	printn "create_next_generation: children_created = $children_created";

	if ($children_created != $total_children) {
	    printn "ERROR: create_next_generation -- did not create correct number of children ($children_created != $total_children)";
	    exit(1);
	}

	# change to next generation
	$current_generation_number_of{$obj_ID} = $current_generation_number = $next_generation_number;
	$current_generation_ref = $current_generation_ref_of{$obj_ID} = $next_generation_ref;
	$current_generation_ref->refresh_individual_names($current_generation_number);
    }

    #--------------------------------------------------------------------------------------
    # Function: score_current_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub score_current_generation {
	my $self = shift; my $obj_ID = ident $self;

	my $config_ref = $config_ref_of{$obj_ID};
	my $cluster_ref = $cluster_ref_of{$obj_ID};

	my $current_generation_number = $current_generation_number_of{$obj_ID};

	my $file_glob = Generation->get_generation_glob(
	    dir => "$config_ref->{work_dir}/$TAG/obj",
	    number => $current_generation_number,
	   );

	my @genome_files = (glob $file_glob);

	my %used_nodes = ();
	printn "score_current_generation: scoring generation $current_generation_number ....";

	for (my $i=0; $i < @genome_files; $i++) {
	    my $genome_file = $genome_files[$i];
	    printn "score_current_generation: scoring file $genome_file";

	    my $node_ref = $cluster_ref->get_free_node();
	    # ensure reproducibility independent of node scoring if there is element of randomness
	    # by deriving node scoring seed from main random generator
	    my $seed = int 1_000_000_000 * rand;  # don't make seed bigger or you lose randomness
	    $node_ref->node_print("srand($seed); \$genome_ref = retrieve(\"$genome_file\"); \$scoring_ref->score_genome(\$genome_ref); store(\$genome_ref, \"$genome_file\");\n");
	    $node_ref->node_expect(undef, 'PERL_SHELL');
	    $node_ref->node_print("NODE_READY");
	    $used_nodes{$node_ref->get_node_ID()} = 1;  # mark this node as one we must wait on
	}

	# now wait for scoring to finish
	while (1) {
	    my @used_list = keys %used_nodes;
	    my @busy_list = $cluster_ref->get_busy_node_id_list();
	    my @wait_list = intersection(\@busy_list, \@used_list);
	    printn "score_current_generation: waiting on nodes... @wait_list";
	    last if (@wait_list == 0);
	    sleep 5;		# poll again in 5 seconds....
	}
	printn"score_current_generation: done waiting on nodes...";
    }

    #--------------------------------------------------------------------------------------
    # Function: report_current_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub report_current_generation {
	my $self = shift; my $obj_ID = ident $self;
	my $current_generation_number = $current_generation_number_of{$obj_ID};

	my @genomes = $self->get_current_generation_ref()->get_elements();

	printn "report_current_generation: generation $current_generation_number";
	for (my $i=0; $i < @genomes; $i++) {
	    my $genome_ref = $genomes[$i];
	    printn "individual $i  : ".$genome_ref->sprint_stats();
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: evolve
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub evolve {
	my $self = shift; my $obj_ID = ident $self;
	my $config_ref = $config_ref_of{$obj_ID};

	$self->create_initial_generation();

	if ($config_ref->{score_initial_generation}) {
	    # reset all score/stats for first generation
	    foreach my $genome_model_ref ($current_generation_ref_of{$obj_ID}->get_elements()) {
		$genome_model_ref->clear_stats(preserve => []);
		$genome_model_ref->set_score(undef);
		$genome_model_ref->set_elite_flag(0);
	    }
	}

      GEN_ALG: while (1) {
	    my $current_generation_number = $current_generation_number_of{$obj_ID};
	    printn "GEN_ALG: generation $current_generation_number";
	    if (($current_generation_number > $config_ref->{first_generation}) || $config_ref->{score_initial_generation}) {
		$self->save_current_generation();
		$self->score_current_generation();
		$self->load_current_generation($current_generation_number_of{$obj_ID});
	    }
	    $self->report_current_generation();

	    if ($current_generation_number_of{$obj_ID} + 1 < $config_ref->{num_generations}) {
		$self->create_next_generation();   # this increments current generation counter
	    } else {
		last GEN_ALG;
	    }
	}
    }
}


sub run_testcases {

    $TAG = "GenAlg";
    system("rm -f test/modules/GenAlg/*.log");  # remove old files

    my $config_file = <<END;
#----------------------------------------
# CPU AND CLUSTER SETTINGS
#----------------------------------------
nice = 15
vmem = 750000
local_dir = scoring/localdir

#----------------------------------------
# WORKSPACE AND CUSTOM SCORING MODULE
#----------------------------------------
scoring_class = Scoring
work_dir = scoring

#----------------------------------------
# GENOME PARAMS
#----------------------------------------
# Genome class
radius = 0
kf_max = 1e5
kf_min = 0.1
kb_max = 10
kb_min = 0.001
kp_max = 1000
kp_min = 1

# Gene class
regulated_concentration_width = 4
gene_unused_width = 32
regulated_concentration_max = 1e-3
regulated_concentration_min = 1e-12

# Domain class
RT_transition_rate_width = 5
TR_transition_rate_width = 6
RT_phi_width = 7
domain_unused_width = 8
RT_transition_rate_max = 1e2
RT_transition_rate_min = 1e-2
TR_transition_rate_max = 1e4
TR_transition_rate_min = 1e-4
RT_phi_max = 0.99
RT_phi_min = 0.01

# ProtoDomain class
binding_profile_width = 3
kf_profile_width = 4
kb_profile_width = 5
kp_profile_width = 6
Keq_profile_width = 7
protodomain_unused_width = 10
Keq_ratio_max = 1e2
Keq_ratio_min = 1e-2

#----------------------------------------
# EVOLUTION PARAMS
#----------------------------------------
inum_genomes = 2           # initial number of genomes when generated randomly
max_population = 4         # population size at each generation (except possibly the first generation)
elite_pool_size = 1        # how many individuals to keep unchanged under the elite strategy
num_generations = 2        # total number of generations before evolution stops

first_generation = 0           # first generation number
remove_old_files = 0           # clean out files from a prior run (be careful with this)
score_initial_generation = 1   # whether to score the initial generation (set to 0 if the generation has already been scored)
rescore_elite = 0              # whether to re-score elite individuals (set to 1 if the score has a random component)

initial_genome = random        # generate initial generation randomly

ranking_nonviable = 0.2      # proportion of non-viable low-ranking individuals
ranking_pressure = 1.2       # per-5-centile fold-change in fitness

# mutate everything
#------------------
prob_mutate_params = 0.6
prob_mutate_global = 0.2
prob_recombination = 0.05
prob_duplicate = 0.05
prob_delete = 0.1

# mutate params only
#------------------
#prob_mutate_params = 1.0
#prob_mutate_global = 0.0
#prob_recombination = 0.0
#prob_duplicate = 0.0
#prob_delete = 0.0

mutation_rate = 0.01

#----------------------------------------
# ANC PARAMS
#----------------------------------------
max_external_iterations = 1
max_internal_iterations = 0
max_complex_size = 4
max_species = 40
max_csite_bound_to_msite_number = 2
default_steric_factor = 1e-3
export_graphviz = nothing
#export_graphviz = network,collapse_states,collapse_complexes

#----------------------------------------
# SIMULATION/SCORING PARAMS
#----------------------------------------
plot_input = 0
plot_output = 0
plot_species = 1
plot_phase = 0
plot_min = -1

solver = ode23s

sampling_interval = 0.1
t_final = 100

# MATLAB odeset params
InitialStep = 1e-8
AbsTol = 1e-9
RelTol = 1e-3
MaxStep = 10.0
END

    burp_file("test/modules/GenAlg.cfg", $config_file);

    my $ref = GenAlg->new({
	config_ref => {
	    work_dir => "test/modules",
	    local_dir => "test/modules/localdir",
	    host_list => ["localhost"],
	    cluster_type => "LOCAL",
	    cluster_size => 2,
	    config_file => "test/modules/GenAlg.cfg",
	    scoring_class => "Scoring",
	    nice => 0,
	},
    });

    while ($ref->get_cluster_ref()->get_busy_node_id_list()) {
	sleep 1;
    }

    printn $ref->_DUMP();

    $ref = undef;
    printn;
    printn "LOGFILE test/modules/GenAlg/ScorNode.1.*.log:";

    system("cat test/modules/GenAlg/ScorNode.1.*.log");
}


# Package BEGIN must return true value
return 1;

