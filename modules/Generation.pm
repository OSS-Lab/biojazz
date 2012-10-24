#-#####################################################################################
#- File:     Generation.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Generation;
use Class::Std::Storable;
use base qw(Set ClassData);
{
    use Carp;
    use Storable qw(store retrieve);

    use Utils;

    use GenomeModel;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    Generation->set_class_data("ELEMENT_CLASS", "GenomeModel");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: get_generation_glob
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_generation_glob {
	my $class = shift;

	my %args = (
	    dir => undef,
	    number => undef,
	    @_,
	   );
	check_args(\%args, 2);

	my $number = $args{number};
	my $dir = $args{dir};

	my $file_glob = sprintf("$dir/G%03d_I*.obj", $number);

	return $file_glob;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_generation_glob
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_generation_temp {
	my $class = shift;

	my %args = (
	    dir => undef,
	    cluster_size => undef,
	    @_,
	   );
	check_args(\%args, 2);

	my $number = $args{cluster_size};
	my $dir = $args{dir};
	
	my @temp_genome_files = undef;
	for (my $i = 0; $i < $number; $i++) {
	    my $temp_file = sprintf("$dir/Gtemp_I%03d.obj", $i);
	    push(@temp_genome_files, $temp_file);
	}

	return @temp_genome_files;
    }

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
#     #--------------------------------------------------------------------------------------
#     # Function: START
#     # Synopsys: 
#     #--------------------------------------------------------------------------------------
# #    sub BUILD {
#     sub START {
#         my ($self, $obj_ID, $arg_ref) = @_;

# 	# check initializers
# 	# ...
#     }

    #--------------------------------------------------------------------------------------
    # Function: retrieve_genomes
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub retrieve_genomes {
	my $self = shift;
	my %args = (
	    files => undef,
	    history_flag => 0,
	    @_,
	   );
	check_args(\%args, 2);

	my @files = @{$args{files}};
	my $history_flag = $args{history_flag};

	my $genome_class = $self->get_class_data("ELEMENT_CLASS");

	foreach my $file (@files) {
	    my $genome_ref = retrieve("$file");  # retrieve() provided by Storable package
	    $genome_ref->add_history("loaded from $file") if $history_flag;
	    $self->add_element($genome_ref);
	}
    }

    
    sub retrieve_largest_temp_score {
	my $self = shift;
	my %args = (
	    files => undef,
	    @_,
	    );
	check_args(\%args, 1);
	my @files = @{$args{files}};
	my @tempScores = undef;
	foreach my $file (@files) {
	    my $genome_ref = retrieve("$file");
	    push(@tempScores, ($genome_ref->get_score()));
	}
	
	my $index = 0;
	
	$tempScores[$index] > $tempScores[$_] or $index = $_ for 1 .. $#tempScores;
	
	return {
	    score => $tempScores[$index],
	    index => $index,
	};
    }

    #--------------------------------------------------------------------------------------
    # Function: store_genomes
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub store_genomes {
	my $self = shift;
	my %args = (
	    files => undef,
	    history_flag => 0,
	    @_,
	   );
	check_args(\%args, 2);

	my @files = @{$args{files}};
	my $history_flag = $args{history_flag};

	my @genomes = $self->get_elements();

	confess "ERROR: no. of filenames doesn't match no. of genomes" if @files != @genomes;

	for (my $i=0; $i < @genomes; $i++) {
	    my $genome_ref = $genomes[$i];
	    my $file = $files[$i];
	    $genome_ref->add_history("saved to $file") if $history_flag;
	    store($genome_ref, "$file"); # store() provided by Storable package
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: create_random_genomes
    # Synopsys: here is the subroutine that creating the genome.
    #           we could use the similar way to generate genomes of the first generation
    #           Besides the sequence there are some other information in the genome file
    #--------------------------------------------------------------------------------------
    sub create_random_genomes {
	my $self = shift; my $obj_ID = ident $self;
	my $config_ref = shift;

	my $inum = $config_ref->{inum_genomes};

	my $genome_class = $self->get_class_data("ELEMENT_CLASS");
	for (my $i = 0; $i < $inum; $i++) {
	    my $genome_ref  = $genome_class->new({
		name => "UNDEF",
		Genome => {
		    %$config_ref,
		    Gene => {
			%$config_ref,
			Domain => {
			    %$config_ref,
			    ProtoDomain => {
				%$config_ref,
			    },
			},
		    },
		},
	    });
	    $genome_ref->get_sequence_ref()->generate_random_sequence(5000);
	    $genome_ref->add_history("random sequence created");
	    $self->add_element($genome_ref);
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: load_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub load_generation {
	my $self = shift; my $obj_ID = ident $self;
	my %args = (
	    dir => undef,
	    number => undef,
	    @_,
	   );
	check_args(\%args, 2);

	my $number = $args{number};
	my $dir = $args{dir};

	my $class = ref $self;
	my $file_glob = $class->get_generation_glob(
	    dir => $dir,
	    number => $number,
	   );

	my @files = glob $file_glob;
	# need to sort on individual number
	@files = sort {$a=~/_I(\d+)/; my $a_i=$1; $b=~/_I(\d+)/; my $b_i=$1; return $a_i <=> $b_i} @files;
	my $num_genomes = @files;
	printn "load_generation: loading generation $number, with $num_genomes individuals";

	$self->retrieve_genomes(files => \@files);
    }

    #--------------------------------------------------------------------------------------
    # Function: save_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub save_generation {
	my $self = shift; my $obj_ID = ident $self;
	my %args = (
	    dir => undef,
	    number => undef,  # this generation's number by default
	    @_,
	   );
	check_args(\%args, 2);

	my $number = $args{number};
	my $dir = $args{dir};

	my $current_population_size = $self->get_num_elements();

	printn "save_generation: saving generation $number, with $current_population_size individuals";

	my @files = map {sprintf("$dir/G%03d_I%02d.obj", $number, $_)} (0..$current_population_size-1);

	$self->store_genomes(files => \@files);
    }

    #--------------------------------------------------------------------------------------
    # Function: refresh_individual_names
    # Synopsys: Refresh individual names based on given generation number and index.
    #--------------------------------------------------------------------------------------
    sub refresh_individual_names {
	my $self = shift; my $obj_ID = ident $self;
	my $number = shift;

	my @refs = $self->get_elements();
	for (my $i=0; $i < @refs; $i++) {
	    my $ref = $refs[$i];
	    my $name = sprintf("G%03d_I%02d", $number, $i);
	    $ref->set_name($name);
#	    $ref->add_history("name refreshed to $name");
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: clear_genomes
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub clear_genomes {
	my $self = shift;

	$self->clear_elements();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_ranked_genomes
    # Synopsys: Returns ranked_xxx, which is xxx ordered by rank, and ranking, which
    #           is the rank of each genome.
    #--------------------------------------------------------------------------------------
    sub get_ranked_genomes {
	my $self = shift;

	my @genomes = $self->get_elements();
	my @scores = map {$_->get_score()} @genomes;
	my @indexes = (0..$#genomes);

	if (grep {!defined $_} @scores) {
	    printn "ERROR: not all scores are defined (if this is first generation, check value of score_initial_generation in config file)";
	    exit(1);
	}

	my @ranked = map {[$genomes[$_], $scores[$_], $indexes[$_]]} @indexes;
	@ranked = sort {$b->[1] <=> $a->[1]} @ranked;

	my @ranked_genomes = map {$_->[0]} @ranked;
	my @ranked_scores = map {$_->[1]} @ranked;
	my @ranked_indices = map {$_->[2]} @ranked;

	my @ranking = ();
	for (my $i = 0; $i < @genomes; $i++) {
	    $ranking[$ranked_indices[$i]] = $i;
	}

	return {
	    ranking =>        \@ranking,  # rank of each genome
	    ranked =>         \@ranked,   # everything ordered by rank
	    ranked_genomes => \@ranked_genomes,  # genomes ordered by rank
	    ranked_scores =>  \@ranked_scores,   # scores ordered by rank
	    ranked_indices => \@ranked_indices   # ordered indices
	   };
    }


    #--------------------------------------------------------------------------------------
    # Function: get_progeny
    # Synopsys: Get the progeny of each genome in the given generation (normally from G+1)
    #--------------------------------------------------------------------------------------
    sub get_progeny {
	my $self = shift;

	my %args = (
	    dir => undef,
	    child_generation_number => undef,  # this should be G+1
	    @_,
	   );
	check_args(\%args, 2);

	my $child_generation_number = $args{child_generation_number};
	my $dir = $args{dir};

	my $child_gen_ref = Generation->new({});
	$child_gen_ref->load_generation(
	    dir => $dir,
	    number => $child_generation_number,
	   );
	return undef if $child_gen_ref->get_num_elements() == 0;

 	my @child_genomes = $child_gen_ref->get_elements();
 	my @replication_history = map {$_->search_history("REPLICATION")} @child_genomes;

 	my @progeny = map {[]} (0..$self->get_last_element_index());
 	my @progeny_indices = map {[]} (0..$self->get_last_element_index());
 	for (my $i=0; $i < @replication_history; $i++) {
 	    my $line = $replication_history[$i];
 	    $line =~ /(\S+)\s*->\s*(\S+)/;
 	    my $parent_name = $1;
 	    #	    my $child_name = $2;
 	    if (defined $parent_name) {
 		$parent_name =~ /G(\S+)_I(\S+)/;
 		my $parent_index = $2;
 		push @{$progeny[$parent_index]}, $child_genomes[$i];
 		push @{$progeny_indices[$parent_index]}, $i;
 	    }
 	}

 	return {
	    child_gen_ref => $child_gen_ref,
	    progeny => \@progeny,
	    progeny_indices => \@progeny_indices,
	};
    }

    #--------------------------------------------------------------------------------------
    # Function: get_wt_sequence
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_wt_sequence {
	my $self = shift;

	my @wt_sequence = ();
	my @genomes = $self->get_elements();
	foreach my $genome_model_ref (@genomes) {
	    my @sequence = split(//,$genome_model_ref->get_sequence_ref()->sprint());
#	    printn join " ", @sequence[0..30];
	    my $L = @sequence > @wt_sequence ? @sequence : @wt_sequence;
	    for (my $i = 0; $i < $L; $i++) {
		$wt_sequence[$i] = ($wt_sequence[$i] || 0) + $sequence[$i];
	    }
	}
	my $N = @genomes;
	@wt_sequence = map {$_/$N < 0.5 ? 0 : 1} @wt_sequence;
#	printn join " ", @wt_sequence[0..30];
	return join "", @wt_sequence;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_HD_matrix
    # Synopsys: Create a matrix that bins genomes according to HD from WT.
    #--------------------------------------------------------------------------------------
    sub get_HD_matrix {
	my $self = shift;

	my $wt_sequence = $self->get_wt_sequence();
	my @HD_matrix = ();
	my @genomes = $self->get_elements();
	foreach my $genome_model_ref (@genomes) {
	    my $sequence = $genome_model_ref->get_sequence_ref()->sprint();
	    my $HD = BitVector->hamming_distance($wt_sequence, $sequence);
	    push @{$HD_matrix[$HD]}, $genome_model_ref;
	}
	@HD_matrix = map {defined $_ ? $_ : []} @HD_matrix;
	return @HD_matrix;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_pop_vs_HD
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_pop_vs_HD {
	my $self = shift;
	my @HD_matrix = $self->get_HD_matrix();
	my @pop_vs_HD = ();
	for (my $i=0; $i < @HD_matrix; $i++) {
	    $pop_vs_HD[$i] =  scalar(@{$HD_matrix[$i]});
	}
	return @pop_vs_HD;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_attribute_vs_HD
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_attribute_vs_HD {
	my $self = shift;
	my $attribute = shift;

	my @HD_matrix = $self->get_HD_matrix();
	my @attribute_vs_HD = ();
	for (my $i=0; $i < @HD_matrix; $i++) {
	    my @attributes = map {$_->get_stats_ref()->{$attribute}} @{$HD_matrix[$i]};
	    $attribute_vs_HD[$i] = average(@attributes);
	}
	return @attribute_vs_HD;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_attributes
    # Synopsys: This function gets a specific attribute's values for all individuals.
    #--------------------------------------------------------------------------------------
    sub get_attributes {
	my $self = shift;
	my $attribute = shift;

	my @genomes = $self->get_elements();
	my @attributes = map {$_->get_stats_ref()->{$attribute}} @genomes;
	return @attributes;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_attribute_names
    # Synopsys: This function collects defined keys from the stats_ref of
    #           each individual in the generation.
    #--------------------------------------------------------------------------------------
    sub get_attribute_names {
	my $self = shift;

	my @genomes = $self->get_elements();
	my @attribute_names = map {keys %{$_->get_stats_ref()}} @genomes;
	@attribute_names = union(\@attribute_names, []);  # this removes duplicates
	return @attribute_names;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_attribute_stats
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_attribute_stats {
	my $self = shift;
	my $attribute = shift;

	my @genomes = $self->get_elements();
	my @attributes = map {$_->get_stats_ref()->{$attribute}} @genomes;
	
	my $mu = average(@attributes);
	my $sigma = sqrt(variance(@attributes));
	return ($mu, $sigma);
    }
}


sub run_testcases {

    my $ref = Generation->new({

       });
    printn $ref->_DUMP();
}


# Package BEGIN must return true value
return 1;

