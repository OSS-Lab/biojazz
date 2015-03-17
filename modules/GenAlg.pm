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
    use Storable qw(store retrieve);
    use Text::CSV;

    use Utils;

    use Globals qw ($verbosity $TAG);

    use ScorCluster;
    use Generation;

    use GenomeModel;

    use Scoring;



    use FindBin qw($Bin);  # need application path

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
    
    my %score_array_ref_of      :ATTR(get => 'score_array_ref', set => 'score_array_ref');
    my %index_array_ref_of      :ATTR(get => 'index_array_ref', set => 'index_array_ref');

    my %reach_target_flag_of    :ATTR(get => 'reach_target_flag', set => 'reach_target_flag', default => 0);
    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

        # INIT
        $score_array_ref_of{$obj_ID} = [];
        $index_array_ref_of{$obj_ID} = [];

        $reach_target_flag_of{$obj_ID} = 0;
    }

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

        # if the host_list is not asigned value put the host_list 
        # ref return a reference or NULL value
        my $host_list_ref = ref $config_ref->{host_list} ? $config_ref->{host_list} : [$config_ref->{host_list}];

        # assign values to cluster reference
        if ($config_ref->{selection_method} eq "population_based_selection") {
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
        }
        # check initializers
        # ...
    }

    #--------------------------------------------------------------------------------------
    # Function: load_current_generation
    # Synopsys: Load generation (i) from disk into current.
    #--------------------------------------------------------------------------------------
    sub load_current_generation {
        # The ident() utility that produces this unique key is provided by the Class::Std module and is identical in effect to the refaddr() function in the standard Scalar::Util module.
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
                #printn "Continue (y/n)? ";
                #my $answer = <>; chomp($answer);
                #if ($answer ne "y") {
                #    printn "no, then fix the configuration file";
                #    exit(1);
                #}
            }
            if (($file_glob =~ /G(...)_I/) && ($1 != $config_ref->{first_generation})) {
                printn "\nWARNING: Configuration parameter first_generation appears not to be set correctly (not 0 or 1)";
                #printn "Continue (y/n)? ";
                #my $answer = <>; chomp($answer);
                #if ($answer ne "y") {
                #    printn "no, then fix the configuration file";
                #    exit(1);
                #}
            }		
            if (!@files) {
                printn "ERROR: nothing to load... ";
                exit(1);
            }

            # need to sort on individual number
            #@files = sort {$a=~/_I(\d+)/; $a_i=$1; $b=~/_I(\d+)/; $b_i=$1; return $a_i <=> $b_i} @files;
            printn "create_initial_generation: loading initial generation from disk";
            printn join "\n", @files;
            $current_generation_ref->retrieve_genomes(
                files => \@files,
                history_flag => 1,
            );

            my @genome_model_refs = $current_generation_ref_of{$obj_ID}->get_elements();
            if ($config_ref->{score_initial_generation}) {
                my $local_dir = (defined $config_ref->{local_dir} ?
                    "$config_ref->{local_dir}/$config_ref->{scoring_class}" :
                    undef);
                my $defined_local_dir = (defined $local_dir) ? "$local_dir/$TAG" : undef;

                eval("use $config_ref->{scoring_class};");
                if ($@) {print $@; return;}

                my $scoring_ref = $config_ref->{scoring_class}->new({
                        config_file => $config_ref->{config_file},
                        node_ID => 999,
                        work_dir => "$config_ref->{work_dir}/$TAG",
                        local_dir => $defined_local_dir,
                        matlab_startup_options => "-nodesktop -nosplash",  # need jvm
                    });

                # reset all score/stats for first generation
                foreach my $genome_model_ref (@genome_model_refs) {
                    $genome_model_ref->clear_stats(preserve => []);
                    $genome_model_ref->set_score(undef);
                    $genome_model_ref->set_elite_flag(0);
                    if ($config_ref->{selection_method} eq "population_based_selection") {
                        if (!$genome_model_ref->get_number()) {
                            $genome_model_ref->set_number(1);
                        }
                    }
                    $scoring_ref->score_genome($genome_model_ref);
                    $genome_model_ref->static_analyse($config_ref->{rescore_elite});
                    $genome_model_ref->set_elite_flag(1);
                }
            }
            my $loaded_genome_num = 0;
            if ($config_ref->{selection_method} eq "kimura_selection") {
                foreach my $genome_model_ref (@genome_model_refs) {
                    if (!$config_ref->{continue_sim}) {
                        $genome_model_ref->set_number(0);
                        $genome_model_ref->set_stepwise_mutations(0);
                        $genome_model_ref->set_stepwise_point_mutations(0);
                        $genome_model_ref->set_accum_mutations(0);
                        $genome_model_ref->set_accum_point_mutations(0);
                    }
                }
            } elsif ($config_ref->{selection_method} eq "population_based_selection") {
                my $i = 0;
                foreach my $genome_model_ref (@genome_model_refs) {
                    $genome_model_ref->set_mutation_index(undef);
                    if (!$genome_model_ref->get_number()) {
                        $genome_model_ref->set_number(1);
                    }
                    if (!$config_ref->{continue_sim}) {
                        $genome_model_ref->set_stepwise_mutations(0);
                        $genome_model_ref->set_stepwise_point_mutations(0);
                        $genome_model_ref->set_accum_mutations(0);
                        $genome_model_ref->set_accum_point_mutations(0);
                    }
                    my $number = $genome_model_ref->get_number();
                    $loaded_genome_num += $number;
                    my $score = $genome_model_ref->get_score();
                    push @{$score_array_ref_of{$obj_ID}}, ($score) x $number;
                    push @{$index_array_ref_of{$obj_ID}}, ($i) x $number;
                }
                $i++;
            } else {
                confess "The selection method is not set appropriately!";
            }

            if ($config_ref->{selection_method} eq 'population_based_selection') {
                my $inum = $config_ref->{evolve_population};

                if ($inum > $loaded_genome_num) { 
                    for (my $i = 0; $i < ($inum - $loaded_genome_num); $i++) {
                        my $index = int(rand($loaded_genome_num));
                        my $genome_model_ref = $genome_model_refs[$index];
                        my $current_number = $genome_model_ref->get_number();
                        $genome_model_ref->set_number($current_number + 1);

                        if ($config_ref->{selection_method} eq 'population_based_selection') {
                            my $score = $genome_model_ref->get_score();
                            push @{$score_array_ref_of{$obj_ID}}, $score;
                            push @{$index_array_ref_of{$obj_ID}}, $index;
                        }
                    }
                } elsif ($inum < $loaded_genome_num) {
                    confess "The inum is smaller than number of genome loaded for initial generation!";
                }
                # check if those two array have the same size as well as with numbers.
                if (scalar @{$score_array_ref_of{$obj_ID}} != scalar @{$index_array_ref_of{$obj_ID}} 
                    || scalar @{$score_array_ref_of{$obj_ID}} != $config_ref->{evolve_population}) {
                    confess "The evolve_population is not equal to the size of score array and index array for selection.";
                }

           }
        } else {
            printn "create_initial_generation: creating $config_ref->{inum_genomes} individuals";
            $current_generation_ref->create_random_genomes($config_ref);

            my @genome_model_refs = $current_generation_ref_of{$obj_ID}->get_elements();
            if ($config_ref->{score_initial_generation}) {
                my $local_dir = (defined $config_ref->{local_dir} ?
                    "$config_ref->{local_dir}/$config_ref->{scoring_class}" :
                    undef);
                my $defined_local_dir = (defined $local_dir) ? "$local_dir/$TAG" : undef;

                eval("use $config_ref->{scoring_class};");
                if ($@) {print $@; return;}

                my $scoring_ref = $config_ref->{scoring_class}->new({
                        config_file => $config_ref->{config_file},
                        node_ID => 999,
                        work_dir => "$config_ref->{work_dir}/$TAG",
                        local_dir => $defined_local_dir,
                        matlab_startup_options => "-nodesktop -nosplash",  # need jvm
                    });

                # reset all score/stats for first generation
                foreach my $genome_model_ref (@genome_model_refs) {
                    $genome_model_ref->clear_stats(preserve => []);
                    $genome_model_ref->set_score(undef);
                    $genome_model_ref->set_elite_flag(0);

                    $scoring_ref->score_genome($genome_model_ref);
                    $genome_model_ref->static_analyse($config_ref->{rescore_elite});
                    $genome_model_ref->set_mutation_index(undef);
                    $genome_model_ref->set_stepwise_mutations(0);
                    $genome_model_ref->set_stepwise_point_mutations(0);
                    $genome_model_ref->set_accum_mutations(0);
                    $genome_model_ref->set_accum_point_mutations(0);
                    $genome_model_ref->set_elite_flag(1);
                }
            }

            my $population = $config_ref->{evolve_population};
            if ($config_ref->{selection_method} eq "kimura_selection") {
                $population = 1;
            }

            my $loaded_genome_num = 0;
            my $i = 0;
            foreach my $genome_model_ref (@genome_model_refs) {
                if (!defined $genome_model_ref->get_number()) {
                    $genome_model_ref->set_number(1);
                }
                my $number = $genome_model_ref->get_number();
                $loaded_genome_num += $number;

                if ($config_ref->{selection_method} eq 'population_based_selection') {
                    my $score = $genome_model_ref->get_score();
                    push @{$score_array_ref_of{$obj_ID}}, ($score) x $number;
                    push @{$index_array_ref_of{$obj_ID}}, ($i) x $number;
                }
                $i++;
            }

            if ($config_ref->{selection_method} eq 'population_based_selection') {
                if ($population > $loaded_genome_num) {
                    for (my $i = 0; $i < ($population - $loaded_genome_num); $i++) {
                        my $index = int(rand($loaded_genome_num));
                        my $genome_model_ref = $genome_model_refs[$index];
                        my $current_number = $genome_model_ref->get_number();
                        $genome_model_ref->set_number($current_number + 1);

                        my $score = $genome_model_ref->get_score();
                        push @{$score_array_ref_of{$obj_ID}}, $score;
                        push @{$index_array_ref_of{$obj_ID}}, $index;
                    }
                }
                # check if those two array have the same size as well as with numbers.
                if (scalar @{$score_array_ref_of{$obj_ID}} != scalar @{$index_array_ref_of{$obj_ID}} 
                    || scalar @{$score_array_ref_of{$obj_ID}} != $config_ref->{evolve_population}) {
                    confess "The evolve_population is not equal to the size of score array and index array for selection.";
                }
            } elsif ($population < $loaded_genome_num) {
               confess "The inum is smaller than number of genome loaded for initial generation!";
            }

        }
        if ($config_ref->{continue_sim} == 1) {
            my $continue_gen = $config_ref->{continue_init};
            $current_generation_ref->refresh_individual_names($continue_gen);
            $current_generation_number_of{$obj_ID} = $continue_gen;
            printn "create_initial_generation for continue simulation: done";
        } else {
            $current_generation_ref->refresh_individual_names($config_ref->{first_generation});
            $current_generation_number_of{$obj_ID} = $config_ref->{first_generation};
            printn "create_initial_generation: done";
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: random_walk_selection
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub random_walk_selection {
        my $self = shift; my $obj_ID = ident $self;

        my $config_ref = $config_ref_of{$obj_ID};
        my $max_mutate_attempts = $config_ref->{max_mutate_attempts};

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

        printn "@scores" if $verbosity > 1;
        my $local_dir = (defined $config_ref->{local_dir} ?
            "$config_ref->{local_dir}/$config_ref->{scoring_class}" :
            undef);
        my $defined_local_dir = (defined $local_dir) ? "$local_dir/$TAG" : undef;

        eval("use $config_ref->{scoring_class};");
        if ($@) {print $@; return;}

        my $scoring_ref = $config_ref->{scoring_class}->new({
                config_file => $config_ref->{config_file},
                node_ID => 999,
                work_dir => "$config_ref->{work_dir}/$TAG",
                local_dir => $defined_local_dir,
                matlab_startup_options => "-nodesktop -nosplash",  # need jvm
            });

        # start to generate new generation 
        my $effective_population_size = $config_ref->{effective_population_size};
        my $amplifier_alpha = $config_ref->{amplifier_alpha};
        my @mutation_step_nums;

        for (my $i = 0; $i < $current_generation_size; $i++) {
            my $parent_ref = $current_generation_ref->get_element($i); 
            my $parent_name = $parent_ref->get_name();

            # start the kimura selection (random walk)
            my $fixation_p = -1;
            my $mutated_score;
            my $mutation_step_num = 0;
            while ($fixation_p < rand) {

                $current_generation_ref->clear_genomes();
                $current_generation_ref->load_generation(
                    dir => "$config_ref->{work_dir}/$TAG/obj",
                    number => $current_generation_number,
                );
                
                $parent_ref = $current_generation_ref->get_element($i);
                printn "discard previous mutation and reloaded parent genome to mutate and score\n" if $verbosity > 1;
                $parent_ref->set_score(undef);
                $parent_ref->clear_stats();
                $parent_ref->set_elite_flag(0);
                $parent_ref->set_stepwise_mutations(0);
                $parent_ref->set_stepwise_point_mutations(0);
                $parent_ref->mutate(
                    mutation_rate_params => $config_ref->{mutation_rate_params},
                    mutation_rate_global => $config_ref->{mutation_rate_global},
                    gene_duplication_rate => $config_ref->{gene_duplication_rate},
                    gene_deletion_rate => $config_ref->{gene_deletion_rate},
                    domain_duplication_rate => $config_ref->{domain_duplication_rate},
                    domain_deletion_rate => $config_ref->{domain_deletion_rate},
                    recombination_rate => $config_ref->{recombination_rate},
                );
                $scoring_ref->score_genome($parent_ref);
                $parent_ref->static_analyse($config_ref->{rescore_elite});
                $parent_ref->set_elite_flag(1);


                printn "the parent's score is: $scores[$i]";
                my $child_score = $parent_ref->get_score();

                printn "the child's score is; $child_score";

                if ($scores[$i] != 0) {
                    $mutated_score = ($child_score - $scores[$i]) / $scores[$i];
                } else {
                    $mutated_score = 1;
                }

                if ($mutated_score == 0.0) {
                    $fixation_p = 1 / (2 * $effective_population_size);  # prevent divide by zero
                } else {
                    $fixation_p = (1 - exp(-2 * $mutated_score)) / (1 - exp(-4 * $effective_population_size * $mutated_score));
                }

                $fixation_p *= $amplifier_alpha;
                $mutation_step_num++;
                if ($max_mutate_attempts > 0 && $mutation_step_num > $max_mutate_attempts) {
                    my $child_ref = $parent_ref->duplicate();
                    $child_ref->set_number($mutation_step_num);
                    $next_generation_ref->add_element($child_ref);

                    # change to next generation
                    $current_generation_ref->clear_genomes();
                    $current_generation_ref = $current_generation_ref_of{$obj_ID} = $next_generation_ref;
                    $current_generation_number = $current_generation_number_of{$obj_ID} = $next_generation_number;
                    $current_generation_ref->refresh_individual_names($current_generation_number);
                    
                    $self->print_attribute_names();
                    $self->report_current_generation();
                    exit(1);
                }
            }

            # after fix the mutation
            my $child_ref = $parent_ref->duplicate();
            $child_ref->set_number($mutation_step_num);
            $next_generation_ref->add_element($child_ref);
        }


        # change to next generation
        $current_generation_ref->clear_genomes();
        $current_generation_ref = $current_generation_ref_of{$obj_ID} = $next_generation_ref;
        $current_generation_number = $current_generation_number_of{$obj_ID} = $next_generation_number;
        $current_generation_ref->refresh_individual_names($current_generation_number);
        if (defined $config_ref->{target_score}) {
            my $current_generation_size = $current_generation_ref->get_num_elements();
            my @new_scores = map {$current_generation_ref->get_element($_)->get_score()} (0..$current_generation_size-1);
            my $max_score = $new_scores[0];
            for (@new_scores) {$max_score = $_ if $_ > $max_score;}
            if ($max_score >= $config_ref->{target_score}) {
                $self->set_reach_target_flag(1);
            }
        }
    }


    #--------------------------------------------------------------------------------------
    # Function: mutate_current_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------

    sub mutate_current_generation {
        my $self = shift; my $obj_ID = ident $self;
        my $config_ref = $config_ref_of{$obj_ID};

        my $current_generation_ref = $current_generation_ref_of{$obj_ID};
        my $current_generation_number = $current_generation_number_of{$obj_ID};
        my $current_generation_size = $current_generation_ref->get_num_elements();

        my $mutation_rate = $config_ref->{mutation_rate};
        my @genome_model_refs = $current_generation_ref->get_elements();
        printn "MUTATION: mutating the $current_generation_number th generation.";

        my $genotype_num = scalar @genome_model_refs;
        my $population = $config_ref->{evolve_population};
        for (my $i = 0; $i < $population; $i++) {
            my $genome_ref = $genome_model_refs[$index_array_ref_of{$obj_ID}->[$i]];
            my $parent_name = $genome_ref->get_name();
            confess "ERROR: The score of $parent_name is UNDEFINED!" if !defined $genome_ref->get_score();
            my $mutation_count = 0;
            my $general_mutation = rand(1);
            my $hgt_mutation = rand(1);
            if ($general_mutation < $mutation_rate || $hgt_mutation < $config_ref->{hgt_rate}) {
                my $child_ref = $genome_ref->duplicate();
                if ($general_mutation < $mutation_rate) {
                    $child_ref->set_stepwise_mutations(0);
                    $child_ref->set_stepwise_point_mutations(0);
                    printn "MUTATION: mutating genome $parent_name.";
                    $mutation_count = $child_ref->mutate(
                        mutation_rate_params => $config_ref->{mutation_rate_params},
                        mutation_rate_global => $config_ref->{mutation_rate_global},
                        gene_duplication_rate => $config_ref->{gene_duplication_rate},
                        gene_deletion_rate => $config_ref->{gene_deletion_rate},
                        domain_duplication_rate => $config_ref->{domain_duplication_rate},
                        domain_deletion_rate => $config_ref->{domain_deletion_rate},
                        recombination_rate => $config_ref->{recombination_rate},
                    );

                } 
                if ($hgt_mutation < $config_ref->{hgt_rate}) {
                    # now to implement the horizontal gene transfer
                    my $donor_ref = $genome_model_refs[$index_array_ref_of{$obj_ID}->[rand($population)]];
                    my $hgt_sequence = $donor_ref->get_chunk();
                    $child_ref->insert_chunk($hgt_sequence);
                    $mutation_count += length($hgt_sequence);
                }

                if ($mutation_count) {
                    $child_ref->set_score(undef);
                    $child_ref->clear_stats(preserve => []);
                    $child_ref->set_elite_flag(0);
                    $child_ref->set_number(1);
                    $child_ref->set_mutation_index($i);
                    $child_ref->add_history(sprintf("MUTATION: $parent_name -> G%03d_I%02d", $current_generation_number + 1, $genotype_num));
                    $current_generation_ref->add_element($child_ref);
                    my $number = $genome_ref->get_number();
                    $genome_ref->set_number($number - 1);
                    $score_array_ref_of{$obj_ID}->[$i] = undef;
                    $index_array_ref_of{$obj_ID}->[$i] = $genotype_num;
                    $genotype_num++;
                } else {
                    $child_ref->DEMOLISH();
                }
            }
        }
        my $next_generation_number = $current_generation_number + 1;
        $current_generation_number = $current_generation_number_of{$obj_ID} = $next_generation_number;
        $current_generation_ref->refresh_individual_names($next_generation_number);

    }

    #--------------------------------------------------------------------------------------
    # Function: score_mutated_genomes
    # Synopsys: 
    #--------------------------------------------------------------------------------------

    sub score_mutated_genomes {
        my $self = shift; my $obj_ID = ident $self;

        my $config_ref = $config_ref_of{$obj_ID};
        my $cluster_ref = $cluster_ref_of{$obj_ID};

        my $current_generation_ref = $current_generation_ref_of{$obj_ID};
        my $current_generation_number = $current_generation_number_of{$obj_ID};
        my $current_generation_size = $current_generation_ref->get_num_elements();

        my $rescore_elite = $config_ref->{rescore_elite}; die "rescore_elite is not specified in config file!" if !defined $rescore_elite;

        my $file_glob = Generation->get_generation_glob(
            dir => "$config_ref->{work_dir}/$TAG/obj",
            number => $current_generation_number,
        );

        my @genome_files = (glob $file_glob);

        my %used_nodes = ();
        printn "score_current_generation: scoring generation $current_generation_number ....";

        for (my $i=0; $i < $current_generation_size; $i++) {
            my $genome_file = $genome_files[$i];
            printn "score_current_generation: scoring file $genome_file";

            my $node_ref = $cluster_ref->get_free_node();
            # ensure reproducibility independent of node scoring if there is element of randomness
            # by deriving node scoring seed from main random generator
            my $seed = int 1_000_000_000 * rand;  # don't make seed bigger or you lose randomness
            $node_ref->node_print("srand($seed); \$genome_ref = retrieve(\"$genome_file\"); \$scoring_ref->score_genome(\$genome_ref); \$genome_ref->static_analyse($rescore_elite); \$genome_ref->set_elite_flag(1); store(\$genome_ref, \"$genome_file\");\n");
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
            sleep 10;		# poll again in 10 seconds....
        }
        printn"score_current_generation: done waiting on nodes!";
    }



    #--------------------------------------------------------------------------------------
    # Function: population_based_selection
    # Synopsys: 
    #--------------------------------------------------------------------------------------

    sub population_based_selection {
        my $self = shift; my $obj_ID = ident $self;

        my $config_ref = $config_ref_of{$obj_ID};
        my $cluster_ref = $cluster_ref_of{$obj_ID};

        my $current_generation_ref = $current_generation_ref_of{$obj_ID};
        my $current_generation_number = $current_generation_number_of{$obj_ID};
        my $current_generation_size = $current_generation_ref->get_num_elements();
        my $population_size = $config_ref->{evolve_population};

        my $temp_generation_ref = Generation->new({});

        printn "create_next_generation: selecting from generation $current_generation_number";
        # update the scores of mutated genotypes
        my @genome_model_refs = $current_generation_ref->get_elements();
        my $i = 0;
        foreach my $genome_model_ref (@genome_model_refs) {
            if (defined $genome_model_ref->get_mutation_index()) {
                my $mutated_score = $genome_model_ref->get_score();
                my $mutation_index = $genome_model_ref->get_mutation_index();
                $score_array_ref_of{$obj_ID}->[$mutation_index] = $mutated_score;
                $genome_model_ref->set_mutation_index(undef);
                if ($index_array_ref_of{$obj_ID}->[$mutation_index] != $i) {
                    die "There is something wrong that mutation index recorded in genome is not same as the genotype index recoding the corresponding mutated genome";
                }
            } else {
                $genome_model_ref->set_stepwise_mutations(0);
                $genome_model_ref->set_stepwise_point_mutations(0);
 
            }
            $i++;
        }
        my @scores = @{$score_array_ref_of{$obj_ID}};
        my @indice = @{$index_array_ref_of{$obj_ID}};

        # check if @score and @indice are in the same size
        if (scalar @scores != scalar @indice || scalar @scores != $population_size) {
            die "The scores and indice arrays are not in same size or not equal to evolve_population size";
        }

        # check scores to make sure defined and positive
        if (grep {!defined $_} @scores || grep {!defined $_} @indice) {
            printn "ERROR: not all scores and population numbers are defined (if this is first generation, check value of score_initial_generation in config file)";
            exit(1);
        }
        if (grep {$_ < 0} @scores || grep {$_ < 0} @indice) {
            printn "ERROR: all scores and indice must be non-negative";
            exit(1);
        }

        # select the next generation
        my $total_score = 0;
        my @cumulative_scores;
        for (my $i = 0; $i < $population_size; $i++) {
            $total_score += $scores[$i];
            push(@cumulative_scores, $total_score);
        }

        my @selection_count = (0) x $current_generation_size;
        for (my $i = 0; $i < $population_size; $i++) {
            my $threshold = rand() * $total_score;
            my $left = 0;
            my $right = $population_size - 1;
            my $index = 0;

            while ($left < $right) {
                $index = int(($left + $right) / 2);
                if ($cumulative_scores[$index] < $threshold) {
                    $left = $index + 1;
                } else {
                    $right = $index;
                }
            }
            $score_array_ref_of{$obj_ID}->[$i] = $scores[$left];
            $index_array_ref_of{$obj_ID}->[$i] = $indice[$left];
            $selection_count[$indice[$left]]++;
        }

        my $accumulative_count = 0;
        my $not_selected_num = 0;
        for (my $i = 0; $i < $current_generation_size; $i++) {
            if ($selection_count[$i]) {
                my $parent_ref = $current_generation_ref->get_element($i);
                my $parent_name = $parent_ref->get_name();

                my $child_ref = $parent_ref->duplicate();
                $child_ref->add_history(sprintf("REPLICATION: $parent_name has $selection_count[$i] descendants in G%03d_I%02d", $current_generation_number, $i));
                $child_ref->set_number($selection_count[$i]);
                $temp_generation_ref->add_element($child_ref);
                $accumulative_count += $selection_count[$i];
                for (my $j = 0; $j < $population_size; $j++) {
                    if ($index_array_ref_of{$obj_ID}->[$j] == $i) {
                        $index_array_ref_of{$obj_ID}->[$j] = $i - $not_selected_num;
                    }
                }
            } else {
                $not_selected_num++;
            }
        }

        if ($accumulative_count != $population_size) {
            die "The accumulative count is not equal to population size, which means the selection might be wrong!";
        }

        $current_generation_ref->clear_genomes();
        $current_generation_ref = $current_generation_ref_of{$obj_ID} = $temp_generation_ref;
        $current_generation_size = $current_generation_ref->get_num_elements();
        $current_generation_number = $current_generation_number_of{$obj_ID};
        $current_generation_ref->refresh_individual_names($current_generation_number);

        if (grep {$_ >= $current_generation_size} @{$index_array_ref_of{$obj_ID}}) {
            die "The indice array is not set properly after selection: some index are larger than genotype number in generation $current_generation_number";
        }

        if (defined $config_ref->{target_score}) {
            my @new_scores = map {$current_generation_ref->get_element($_)->get_score()} (0..$current_generation_size-1);
            my $max_score = $new_scores[0];
            for (@new_scores) {$max_score = $_ if $_ > $max_score;}
            if ($max_score >= $config_ref->{target_score}) {
                $self->set_reach_target_flag(1);
            }
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: print_attribute_names
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub print_attribute_names {
        my $self = shift; my $obj_ID = ident $self;
        my $current_generation_number = $current_generation_number_of{$obj_ID};
        my $current_generation_ref = $current_generation_ref_of{$obj_ID};
        my $config_ref = $config_ref_of{$obj_ID};

        my @genome_attribute_names = @{$config_ref->{genome_attribute_names}};
        confess "genome_attribute_names is not specified!" if (!@genome_attribute_names);

        if ($genome_attribute_names[0] eq "all") {
            undef @genome_attribute_names;
            @genome_attribute_names = $current_generation_ref->get_attribute_names();
        }

        my $data_dir = "$config_ref->{work_dir}/$TAG/stats";
        my $file_name = sprintf("$data_dir/Generation%03d.csv", $current_generation_number);
        open my $data_file, ">> $file_name" or die "$file_name: $!";
        my $csv = Text::CSV->new({binary => 1, eol => "\n"});

        my $second_attribute = 'Population/MutationSteps';
        if ($config_ref->{selection_method} eq "kimura_selection") {
            $second_attribute = 'Mutation_attempts';
        } elsif ($config_ref->{selection_method} eq "population_based_selection") {
            $second_attribute = 'Population_per_mutant';
        }
 
        my @attribute_names_new = ('Name', $second_attribute, 'Accum_mutations', 
            'Accum_point_mutations', 'Stepwise_mutations', 'Stepwise_point_mutations', @genome_attribute_names);
        $csv->print($data_file, \@attribute_names_new);

        close($data_file) || warn "close failed: $!";

        return 1;
    }

    #--------------------------------------------------------------------------------------
    # Function: report_current_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub report_current_generation {
        my $self = shift; my $obj_ID = ident $self;
        my $current_generation_number = $current_generation_number_of{$obj_ID};
        my $current_generation_ref = $current_generation_ref_of{$obj_ID};
        my $config_ref = $config_ref_of{$obj_ID};

        my @genomes = $self->get_current_generation_ref()->get_elements();

        my @genome_attribute_names = @{$config_ref->{genome_attribute_names}};

        confess "genome_attribute_names is not specified!" if (!@genome_attribute_names);
        if ($genome_attribute_names[0] eq "all") {
            undef @genome_attribute_names;
            @genome_attribute_names = $current_generation_ref->get_attribute_names();
        }

        printn "report_current_generation: generation $current_generation_number";

        my $data_dir = "$config_ref->{work_dir}/$TAG/stats";
        my $file_name = sprintf("$data_dir/Generation%03d.csv", $current_generation_number);
        open my $data_file, ">> $file_name" or die "$file_name: $!";
        my $csv = Text::CSV->new({binary => 1, eol => "\n"});

        my @attributes = ();
        for (my $i=0; $i < @genomes; $i++) {
            my $genome_ref = $genomes[$i];
            push(@attributes, $genome_ref->get_name());
            push(@attributes, $genome_ref->get_number());
            push(@attributes, $genome_ref->get_accum_mutations());
            push(@attributes, $genome_ref->get_accum_point_mutations());
            push(@attributes, $genome_ref->get_stepwise_mutations());
            push(@attributes, $genome_ref->get_stepwise_point_mutations());
            # Here, we output each genome stats into a line 
            # of CSV file
            for (my $j = 0; $j < scalar @genome_attribute_names; $j++) {
                my $attribute_value = $genome_ref->get_stats_ref()->{$genome_attribute_names[$j]};
                push(@attributes, $attribute_value);
            }
            # process the attributes and output
            $csv->print($data_file, \@attributes);

            # destroy @attributes
            undef @attributes;

        }
        close($data_file) || warn "close failed: $!";

        return 1;
    }

    #--------------------------------------------------------------------------------------
    # Function: report_selection
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub report_selection {
        my $self = shift; my $obj_ID = ident $self;
        my $current_generation_number = $current_generation_number_of{$obj_ID};
        my $current_generation_ref = $current_generation_ref_of{$obj_ID};
        my $config_ref = $config_ref_of{$obj_ID};
        
        my $scores_ref = $score_array_ref_of{$obj_ID};
        my $indice_ref = $index_array_ref_of{$obj_ID};

        my $population = $config_ref->{evolve_population};

        if ($population != scalar @{$scores_ref} || $population != scalar @{$indice_ref}) {
            confess "The size of scores array and indice array is not equal to population size!";
        }

        my $csv = Text::CSV->new({binary => 1, eol => "\n"});

        my $data_dir = "$config_ref->{work_dir}/$TAG/report";
        my $scores_file = "$data_dir/selection_scores.csv";
        open my $scores_data, ">> $scores_file" or die "$scores_file: $!";
        $csv->print($scores_data, $scores_ref);
        close($scores_data) || warn "close failed: $!";

        my $indice_file = "$data_dir/selection_indice.csv";
        open my $indice_data, ">> $indice_file" or die "$indice_file: $!";
        $csv->print($indice_data, $indice_ref);
        close($indice_data) || warn "close failed: $!";

        return 1;
    }

    #---  FUNCTION  ----------------------------------------------------------------
    #         NAME: clear_objs
    #   PARAMETERS: ????
    #      RETURNS: ????
    #  DESCRIPTION: ????
    #       THROWS: no exceptions
    #     COMMENTS: none
    #     SEE ALSO: n/a
    #-------------------------------------------------------------------------------

    sub clear_objs {
        my $self = shift; my $obj_ID = ident $self;
        my $generation_number = shift;
        die "No generation number are specified!" if !defined $generation_number;
        my $config_ref = $config_ref_of{$obj_ID};
        my $obj_dir = "$config_ref->{work_dir}/$TAG/obj";
        my $removal_files = sprintf("$obj_dir/G%03d_I*.obj", $generation_number);

        `echo $removal_files    |  xargs rm -f`;
        printn "Removing genome files of the $generation_number th generation" if $verbosity > 1;

        return 1;
    } ## --- end sub clear_pre_gen

    #--------------------------------------------------------------------------------------
    # Function: evolve
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub evolve {
        my $self = shift; my $obj_ID = ident $self;
        my $config_ref = $config_ref_of{$obj_ID};
        my $first_generation = $config_ref->{first_generation};
        my $fossil_epoch = $config_ref->{fossil_epoch};

        $self->create_initial_generation();
        if ($config_ref->{continue_sim} != 1) {
            if ($config_ref->{report_on_fly} == 1) {
                $self->print_attribute_names();
                $self->report_current_generation();
            }
            if ($config_ref->{selection_method} eq "population_based_selection" && $config_ref->{report_selection} == 1) {
                $self->report_selection();
            }
        }
        $self->save_current_generation();

        GEN_ALG: while (1) {
            my $current_generation_number = $current_generation_number_of{$obj_ID};
            printn "GEN_ALG: generation $current_generation_number";

            $self->load_current_generation($current_generation_number_of{$obj_ID});


            if ($current_generation_number_of{$obj_ID} <= $config_ref->{num_generations}
                && $reach_target_flag_of{$obj_ID} != 1) {
                if ($config_ref->{selection_method} eq "kimura_selection") {
                    $self->random_walk_selection();
                    if ($config_ref->{report_on_fly} == 1) {
                        $self->print_attribute_names();
                        $self->report_current_generation();
                    }
                    $self->save_current_generation();
                } elsif ($config_ref->{selection_method} eq "population_based_selection") {
                    $self->mutate_current_generation();
                    $self->save_current_generation();
                    $self->score_mutated_genomes();
                    $self->load_current_generation($current_generation_number_of{$obj_ID});
                    $self->population_based_selection();
                    if ($config_ref->{report_on_fly} == 1) {
                        $self->print_attribute_names();
                        $self->report_current_generation();
                    }
                    if ($config_ref->{report_selection} == 1) {
                        $self->report_selection();
                    }
                    $self->clear_objs($current_generation_number_of{$obj_ID});
                    $self->save_current_generation();
                    # clear genome files in order to relife the storage burdon
                    if (defined $fossil_epoch) {
                        my $current_generation_number = $current_generation_number_of{$obj_ID};
                        if ($current_generation_number != 1 &&
                            (($current_generation_number - $first_generation) % $fossil_epoch) 
                            != 1) {
                            $self->clear_objs($current_generation_number - 1);
                        }
                    }
                } elsif (!$config_ref->{selection_method}) {
                    printn "the selection method is not specified";
                    exit(1);
                } else {
                    printn "Selection method is not recognizible or made-up";
                    exit(1);
                }
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
vmem = 2000000
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

