#-#####################################################################################
#- File:     GenomeModel.pm
#- Synopsys: A GenomeModel consists of a Sequence and Parsers required to interpret the
#-           information contained in the Sequence.
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package GenomeModel;
use Class::Std::Storable;
use base qw(Model);
{
    use Carp;
    use Storable qw(dclone);

    use Utils;
    use Globals;

    use Genome;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    GenomeModel->set_class_data("PARSER_CLASS", "Genome");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    # simulation parameters
    my %history_ref_of :ATTR(get => 'history_ref', set => 'history_ref');
    my %stats_ref_of   :ATTR(get => 'stats_ref', set => 'stats_ref');
    my %score_of       :ATTR(get => 'score', set => 'score');
    my %elite_flag_of  :ATTR(get => 'elite_flag', set => 'elite_flag', default => 0);

#    my %cell_volume_of :ATTR(get => 'cell_volume', set => 'cell_volume');

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
#    #--------------------------------------------------------------------------------------
#    # Function: XXX
#    # Synopsys: 
#    #--------------------------------------------------------------------------------------
#    sub XXX {
#	my $class = shift;
#    }

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: START
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

	# INIT
	$history_ref_of{$obj_ID} = [];
	$stats_ref_of{$obj_ID} = {};

#	$cell_volume_of{$obj_ID} = $arg_ref->{cell_volume} if exists $arg_ref->{cell_volume};
    }

    #--------------------------------------------------------------------------------------
    # Function: get_stat
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_stat {
	my $self = shift;
	my $stat = shift;
	return $stats_ref_of{ident $self}->{$stat};
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint_stats
    # Synopsys: Return
    #--------------------------------------------------------------------------------------
    sub sprint_stats {
	my $self = shift;

	my $stats_ref = $self->get_stats_ref();

	my $str = "score=".(defined $self->get_score() ? $self->get_score() : "UNDEF");
	foreach my $key (sort keys %{$stats_ref}) {
	    my $value = $stats_ref->{$key};
	    $value = defined $value ? $value : "UNDEF";
	    $str .= ", $key=".$value;
	}
	return $str;
    }

    #--------------------------------------------------------------------------------------
    # Function: set_stat
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub set_stat {
	my $self = shift;
	my $stat = shift;
	my $value = shift;
	$stats_ref_of{ident $self}->{$stat} = $value;
    }

    #--------------------------------------------------------------------------------------
    # Function: clear_stats
    # Synopsys: Clear stats except those indicated, without changing stat_ref.
    #--------------------------------------------------------------------------------------
    sub clear_stats {
	my $self = shift; my $obj_ID = ident $self;
	my %args = (
	    preserve => [],
	    @_,
	   );
	check_args(\%args, 1);
	my $preserve_ref = $args{preserve};
	my $stats_ref = $stats_ref_of{$obj_ID};
	my @preserve = ();
	foreach my $key (@$preserve_ref) {
	    if (exists $stats_ref->{$key}) {
		push @preserve, ($key, $stats_ref->{$key});
	    }
	}
	%$stats_ref = @preserve;
    }

    #--------------------------------------------------------------------------------------
    # Function: duplicate
    # Synopsys: Creates a new genome identical to current one.
    #--------------------------------------------------------------------------------------
    sub duplicate {
	my $self = shift;

	my $child_ref = dclone($self);

	return $child_ref;
    }

    #--------------------------------------------------------------------------------------
    # Function: add_history
    # Synopsys: Add history lines, only preserve last 50 lines by default.
    #--------------------------------------------------------------------------------------
    sub add_history {
	my $self = shift; my $obj_ID = ident $self;
	my $event = shift;
	my $max_lines = shift;

	$max_lines = defined $max_lines ? $max_lines : 50;

	my $history_ref = $history_ref_of{$obj_ID};

	my $first = ($max_lines == -1) ? 0 : $#{$history_ref}-$max_lines+1;
	$first = 0 if $first < 0;

	push @{$history_ref}, $event;
	@{$history_ref} = @{$history_ref}[$first..$#{$history_ref}];
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint_history
    # Synopsys: Print history lines, only last 10 lines by default.
    #--------------------------------------------------------------------------------------
    sub sprint_history {
	my $self = shift; my $obj_ID = ident $self;
	my $max_lines = shift;

	$max_lines = defined $max_lines ? $max_lines : 50;

	my $history_ref = $history_ref_of{$obj_ID};

	my $first = ($max_lines == -1) ? 0 : $#{$history_ref}-$max_lines+1;
	$first = 0 if $first < 0;

	my $history = join "\n",@{$history_ref}[$first..$#{$history_ref}];
	return "HISTORY($max_lines):\n$history\n";
    }

    #--------------------------------------------------------------------------------------
    # Function: search_history
    # Synopsys: Search history for most recent match
    #--------------------------------------------------------------------------------------
    sub search_history {
	my $self = shift; my $obj_ID = ident $self;
	my $pattern = shift;

	my $history_ref = $history_ref_of{$obj_ID};
	for (my $i=$#{$history_ref}; $i >= 0; $i--) {
	    return $history_ref->[$i] if ($history_ref->[$i] =~ /$pattern/);
	}
	return undef;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_genome_parser_ref
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_genome_parser_ref {
	my $self = shift;
	return $self->get_parser_ref();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_gene_parser_ref
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_gene_parser_ref {
	my $self = shift;
	return $self->get_genome_parser_ref()->get_gene_parser_ref();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_domain_parser_ref
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_domain_parser_ref {
	my $self = shift;
	return $self->get_gene_parser_ref()->get_domain_parser_ref();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_protodomain_parser_ref
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_protodomain_parser_ref {
	my $self = shift;
	return $self->get_domain_parser_ref()->get_protodomain_parser_ref();
    }

    #--------------------------------------------------------------------------------------
    # Function: generate_random_genome
    # Synopsys: Generate a genome of a given length in bits.
    #--------------------------------------------------------------------------------------
    sub generate_random_genome {
	my $self = shift;
	my $genome_length = shift;

	printn "generate_random_genome: generating genome of length $genome_length" if ($verbosity >= 1);
	$self->get_sequence_ref()->generate_random_sequence($genome_length);
    }

    #--------------------------------------------------------------------------------------
    # Function: load_genome_sequence
    # Synopsys: Load a genome sequence from a file and clear parser.
    #--------------------------------------------------------------------------------------
    sub load_genome_sequence {
	my $self = shift;

	$self->get_sequence_ref()->load_sequence(@_);
	$self->get_parser_ref()->clear_parser_instances();
    }

    #--------------------------------------------------------------------------------------
    # Function: save_genome
    # Synopsys: Save a genome to a file
    #--------------------------------------------------------------------------------------
    sub save_genome_sequence {
	my $self = shift;
	$self->get_sequence_ref()->save_sequence(@_);
    }

    #--------------------------------------------------------------------------------------
    # Function: check
    # Synopsys: Check if the genome was parsed successfully.
    #--------------------------------------------------------------------------------------
    sub check {
	my $self = shift;
	my $obj_ID = ident $self;

	my @parser_instances = $self->get_parser_ref()->get_object_instances();

	if (@parser_instances == 0) {
	    printn "WARNING: no sequences parsed successfully";
	    return 0;
	} elsif (ident($parser_instances[0]->get_sequence_ref()) != ident($self->get_sequence_ref())) {
	    printn "WARNING: primary sequence of GenomeModel was not parsed successfully";
	    return 0;
	} else {
	    return 1;
	}
    }

    #--------------------------------------------------------------------------------------
    # These functions return ParserInstance(S)
    #--------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------------
    # Function: get_genome
    # Synopsys: Get genome ParserInstance.
    #--------------------------------------------------------------------------------------
    sub get_genome {
	my $self = shift;
	my $obj_ID = ident $self;

	if ($self->get_genome_parser_ref()->get_object_instance_count() != 1) {
	    return undef;
	} else {
	    return $self->get_genome_parser_ref()->get_object_instance_by_index(0);
	}
    }


    #--------------------------------------------------------------------------------------
    # Function: get_genes
    # Synopsys: Return list of gene ParserInstances.
    #--------------------------------------------------------------------------------------
    sub get_genes {
	my $self = shift;

	return $self->get_gene_parser_ref()->get_object_instances();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_gene_names
    # Synopsys: Return list of parsed gene names.
    #--------------------------------------------------------------------------------------
    sub get_gene_names {
	my $self = shift;

	return $self->get_gene_parser_ref()->get_object_instance_names();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_num_genes
    # Synopsys: Return number of parsed genes.
    #--------------------------------------------------------------------------------------
    sub get_num_genes {
	my $self = shift;

	return $self->get_gene_parser_ref()->get_object_instance_count();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_gene_info_by_index
    # Synopsys: Find gene ParserInstance in genome by its index and return info
    #--------------------------------------------------------------------------------------
    sub get_gene_by_index {
	my $self = shift;

	return $self->get_gene_parser_ref()->get_object_instance_by_index(@_);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_gene_by_name
    # Synopsys: Find gene ParserInstance in genome by its name and return info.
    #--------------------------------------------------------------------------------------
    sub get_gene_by_name {
	my $self = shift;
	my $name = shift;

	return $self->get_gene_parser_ref()->lookup_object_instance_by_name($name);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_last_gene_index
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_last_gene_index {
	my $self = shift;

	return $self->get_gene_parser_ref()->get_object_instances() - 1;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_last_gene
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_last_gene {
	my $self = shift;

	return $self->get_gene_by_index($self->get_last_gene_index());
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint_genes
    # Synopsys: Print all genes that have been parsed from genome
    #--------------------------------------------------------------------------------------
    sub sprint_genes {
	my $self = shift;

	my @genes = $self->get_genes();
	return join "\n", (map {$_->sprint(@_)} @genes);
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint_gene_by_name
    # Synopsys: Print info for named gene
    #--------------------------------------------------------------------------------------
    sub sprint_gene_by_name {
	my $self = shift;
	my $gene_name = shift;

	my $gene_ref = $self->get_gene_by_name($gene_name);
	return $gene_ref->sprint(@_);
    }

    #--------------------------------------------------------------------------------------
    # Function: pick_random_gene
    # Synopsys: Pick a gene at random (returns gene index).
    #--------------------------------------------------------------------------------------
    sub pick_random_gene {
	my $self = shift;

	return int rand $self->get_num_genes();
    }


    #--------------------------------------------------------------------------------------
    # Function: get_domain_sequence
    # Synopsys: Given a gene and domain number, extract the sequence from genome.
    #--------------------------------------------------------------------------------------
    sub get_domain_sequence {
	my $self = shift;
	my $gene_index = shift;
	my $domain_index = shift;

	return $self->get_gene_by_index($gene_index)->get_field_sequence(["domains", $domain_index]);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_domains
    # Synopsys: Given a gene index, return the list of domain ParserInstances.
    #--------------------------------------------------------------------------------------
    sub get_domains {
	my $self = shift;
	my $gene_index = shift;

	return @{$self->get_gene_by_index($gene_index)->get_field_ref(["domains"])};
    }

    #--------------------------------------------------------------------------------------
    # These functions are genome mutators
    #--------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------------
    # Function: terminate_genome
    # Synopsys: Terminate tail end of genome by forcing last gene's stop code to
    #           all one's and also the  of genome.  This makes sure that there
    #           are no partial genes which will interfere with duplication.
    #           Requires that genome has been parsed.
    #--------------------------------------------------------------------------------------
    sub terminate_genome {
	my $self = shift;
	my $obj_ID = ident $self;

	confess "ERROR: genome has not been parsed" if (!defined $self->get_genome());

	my $last_gene_ref = $self->get_last_gene();
	my $last_stop_locus = $last_gene_ref->get_field_locus("STOP_CODE");
	my $last_stop_length = $last_gene_ref->get_field_length("STOP_CODE");

	my $genome_ref = $self->get_genome();
	my $post_junk_locus = $genome_ref->get_field_locus("POST_JUNK");
	my $post_junk_length = $genome_ref->get_field_length("POST_JUNK");
	my $post_junk_stop = $genome_ref->get_field_stop("POST_JUNK");

	printn "terminate_genome: setting all bits from position $last_stop_locus to $post_junk_stop";
	my $last_stop = "1" x $last_stop_length;
	$last_gene_ref->set_field("STOP_CODE", $last_stop);
	my $post_junk = "1" x $post_junk_length;
	$genome_ref->set_field("POST_JUNK", $post_junk);
    }

    #--------------------------------------------------------------------------------------
    # Function: set_field_mutation_rates
    # Synopsys: Select the rate at which to mutate fields affecting biochemical params only,
    #           versus fields affecting structure and regulation of network.
    #--------------------------------------------------------------------------------------
    sub set_field_mutation_rates {
	my $self = shift;
	my %args = (
	    mutate_params_rate => 1.0,    # mutation rate for parameter fields
	    mutate_network_rate => 1.0,   # mutation rate for network structure fields
	    @_,
	   );
	check_args(\%args,2);
	my $mutate_params_rate = $args{mutate_params_rate};
	my $mutate_network_rate = $args{mutate_network_rate};

	my $protodomain_parser_ref = $self->get_protodomain_parser_ref();
	$protodomain_parser_ref->set_field_mutation_rate("type", $mutate_network_rate);
	$protodomain_parser_ref->set_field_mutation_rate("substrate_polarity", $mutate_network_rate);
	$protodomain_parser_ref->set_field_mutation_rate("binding_profile", $mutate_network_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kf_profile", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kb_profile", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kp_profile", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("Keq_ratio", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kf_polarity_mask", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kb_polarity_mask", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kf_conformation_mask", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kb_conformation_mask", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("kp_conformation_mask", $mutate_params_rate);
	$protodomain_parser_ref->set_field_mutation_rate("UNUSED", $mutate_network_rate);

	my $domain_parser_ref = $self->get_domain_parser_ref();
	$domain_parser_ref->set_field_mutation_rate("allosteric_flag", $mutate_network_rate);
	$domain_parser_ref->set_field_mutation_rate("RT_transition_rate", $mutate_params_rate);
	$domain_parser_ref->set_field_mutation_rate("TR_transition_rate", $mutate_params_rate);
	$domain_parser_ref->set_field_mutation_rate("RT_phi", $mutate_params_rate);
	$domain_parser_ref->set_field_mutation_rate("UNUSED", $mutate_network_rate);
	$domain_parser_ref->set_field_mutation_rate("protodomains", 1.0);
	$domain_parser_ref->set_field_mutation_rate("protodomains_linker", $mutate_network_rate);

	my $gene_parser_ref = $self->get_gene_parser_ref();
	$gene_parser_ref->set_field_mutation_rate("START_CODE", $mutate_network_rate);
	$gene_parser_ref->set_field_mutation_rate("regulated_concentration", $mutate_params_rate);
	$gene_parser_ref->set_field_mutation_rate("UNUSED", $mutate_network_rate);
	$gene_parser_ref->set_field_mutation_rate("domains", 1.0);
	$gene_parser_ref->set_field_mutation_rate("domains_linker", $mutate_network_rate);
	$gene_parser_ref->set_field_mutation_rate("STOP_CODE", $mutate_network_rate);

	my $genome_parser_ref = $self->get_genome_parser_ref();
	$genome_parser_ref->set_field_mutation_rate("PRE_JUNK", $mutate_network_rate);
	$genome_parser_ref->set_field_mutation_rate("genes", 1.0);
	$genome_parser_ref->set_field_mutation_rate("genes_linker", $mutate_network_rate);
	$genome_parser_ref->set_field_mutation_rate("POST_JUNK", $mutate_network_rate);
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate_genome/mutate_genes
    # Synopsys: Mutate all fields in genome and genes that have been parsed.  This includes
    #           unused fields within and between genes and linkers.
    #           Returns list of mutated gene names and number of bits mutated in each.
    #--------------------------------------------------------------------------------------
    sub mutate_genes {
	my $self = shift;
	my $mutation_rate = shift;
	confess "ERROR: specify a mutation rate" if !defined $mutation_rate;

	my @genes = $self->get_genes();
	my $num_bits;

	return map(@$_, grep($_->[1] != 0, map([$_->get_name(), $_->mutate($mutation_rate)], @genes)));
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate_genome
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub mutate_genome {
	my $self = shift;
	my $mutation_rate = shift;
	confess "ERROR: specify a mutation rate" if !defined $mutation_rate;

	my @genes = $self->get_genes();
	my $num_bits;

	my $genome_parser_ref = $self->get_genome_parser_ref();
	# first mutate only stuff between genes and junk
	my $genes_rate = $genome_parser_ref->get_field_mutation_rate("genes");
	$genome_parser_ref->set_field_mutation_rate("genes", 0.0);
	my $mutated_junk_bits = $self->get_genome->mutate($mutation_rate);
	printn "mutated $mutated_junk_bits bits in junk dna";
	# now mutate genes themselves, restoring saved rate
	$genome_parser_ref->set_field_mutation_rate("genes", $genes_rate);
	return $self->mutate_genes($mutation_rate);
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate_gene_by_name
    # Synopsys: Given name of gene, mutate bits within it.  Returns number of bits mutated.
    #--------------------------------------------------------------------------------------
    sub mutate_gene_by_name {
	my $self = shift;
	my $gene_name = shift;

	my $gene_ref = $self->get_gene_by_name($gene_name);
	confess "ERROR: can't find gene $gene_name" if !defined $gene_ref;
	return $gene_ref->mutate(@_);
    }

    #--------------------------------------------------------------------------------------
    # Function: duplicate_gene
    # Synopsys: Pick a gene at random and duplicate it (copy to end of genome) with
    #           probability given by argument (or always if no argument).
    #           Returns name of chosen gene and locus of duplicate.
    #--------------------------------------------------------------------------------------
    sub duplicate_gene {
	my $self = shift; my $obj_ID = ident $self;
	my $probability = shift;
	$probability = !defined $probability ? 1.0 : 0.0;

	if ((rand 1) >= $probability) {
	    printn "duplicate_gene: no duplication occurred" if $verbosity >= 1;
	    return;
	}
	
	
	my $sequence_ref = $self->get_sequence_ref();
	my $duplicate_locus = $sequence_ref->get_length();

	my $gene_ref = $self->get_gene_by_index($self->pick_random_gene());
	my $gene_name = $gene_ref->get_name();
	my $gene_sequence = $gene_ref->get_sequence();

	printn "duplicate_gene: duplicating gene $gene_name" if $verbosity >= 1;

	# append the gene
	$sequence_ref->splice_subseq($gene_sequence);
	# terminate it
	$sequence_ref->splice_subseq($gene_ref->get_STOP_linker_code(),
				     $sequence_ref->get_length() - length($gene_ref->get_STOP_linker_code()),
				     length($gene_ref->get_STOP_linker_code()));

	return ($gene_name, $duplicate_locus);
    }

    #--------------------------------------------------------------------------------------
    # Function: recombine_genes
    # Synopsys: Recombine two given genes as follows
    #           1D + 1D   -> fuse to given 2D gene
    #           1D + 2D   -> swap LHS with RHS domain chosen at random
    #           2D + 2D   -> swap randomly chosen LHS domain w/ randomly chosen RHS domain
    #           xD + xD   -> general case is to choose a single domain from each gene
    #--------------------------------------------------------------------------------------
    # TBD: Right now only small 2-domain genes are created.  Need to change general case to
    #      choose m<=M and n<=N sequential domains from each gene where M and N are 
    #      the no. domains from each gene.
    #--------------------------------------------------------------------------------------
    sub recombine_genes {
	my $self = shift;
	my $obj_ID = ident $self;
	my ($gene1_index, $gene2_index) = (shift, shift);

	my $gene1_ref = $self->get_gene_by_index($gene1_index);
	my $gene2_ref = $self->get_gene_by_index($gene2_index);
	my $gene1_name = $gene1_ref->get_name();
	my $gene2_name = $gene2_ref->get_name();
	my @gene1_domains = $self->get_domains($gene1_index);
	my @gene2_domains = $self->get_domains($gene2_index);
	my $gene1_num_domains = @gene1_domains;
	my $gene2_num_domains = @gene2_domains;

	printn "recombine_genes: recombining $gene1_name($gene1_index, $gene1_num_domains domains) and $gene2_name($gene2_index, $gene2_num_domains domains)" if $verbosity >= 1;
	printn "recombine_genes: $gene1_name sequence = " . $gene1_ref->get_sequence() if $verbosity >= 2;
	printn "recombine_genes: $gene2_name sequence = " . $gene2_ref->get_sequence() if $verbosity >= 2;

	my $domain1_index = int rand $gene1_num_domains;
	my $domain2_index = int rand $gene2_num_domains;

	printn "recombine_genes: using domains $gene1_name($domain1_index) and $gene2_name($domain2_index)" if $verbosity >= 1;
	my $domain1_sequence = $gene1_domains[$domain1_index]->get_sequence();
	my $domain2_sequence = $gene2_domains[$domain2_index]->get_sequence();
	printn "recombine_genes: domain1_sequence = " . $domain1_sequence if $verbosity >= 2;
	printn "recombine_genes: domain2_sequence = " . $domain2_sequence if $verbosity >= 2;

	my $gene_parser_ref = $self->get_gene_parser_ref();

	my $new_gene_sequence = ($gene_parser_ref->get_gene_start_code() .
				 Sequence->new()->generate_random_sequence($gene_parser_ref->get_regulated_concentration_width()) .
				 Sequence->new()->generate_random_sequence($gene_parser_ref->get_unused_width()) .
				 $domain1_sequence .
				 $gene_parser_ref->get_soft_linker_code() .
				 $domain2_sequence .
				 $gene_parser_ref->get_STOP_linker_code()
				);
	if (@{$gene_parser_ref->get_structure_ref()} != 5) {
	    confess "ERROR: internal error, structure_ref of gene appears to have changed, so new gene sequence may not be built correctly";
	}

	my $new_gene_start = $self->get_sequence_ref()->get_length();

	$self->get_sequence_ref()->splice_subseq("00000000".$new_gene_sequence);
	printn "recombine_genes: new gene sequence = " . $new_gene_sequence if $verbosity >= 2;

	return $new_gene_start;
    }

    #--------------------------------------------------------------------------------------
    # Function: delete_gene (by ref or by name)
    # Synopsys: Delete genes from genome by zeroing out the entire gene!
    #--------------------------------------------------------------------------------------
    sub delete_gene {
	my $self = shift;
	my @genes = @_;

	my @deleted_genes;
	foreach my $gene_ref (@genes) {
	    $gene_ref = ref $gene_ref ? $gene_ref : $self->get_gene_by_name($gene_ref);
	    my $gene_name = $gene_ref->get_name();
	    my $gene_start = $gene_ref->get_locus();
	    my $gene_stop = $gene_ref->get_stop_locus();
	    my $gene_length = $gene_ref->get_length();
	    if (defined $gene_ref) {
		printn "delete_gene: deleting $gene_name start=$gene_start stop=$gene_stop" if $verbosity >= 2;
	    } else {
		printn "ERROR: delete_gene -- no such gene $gene_name";
		exit(1);
	    }
	    $self->get_sequence_ref()->splice_subseq("0" x $gene_length, $gene_start, $gene_length);
	    push @deleted_genes, $gene_ref;
	}
	return @deleted_genes;
    }

    #--------------------------------------------------------------------------------------
    # Function: delete_random_gene
    # Synopsys: Pick a random gene and delete it
    #--------------------------------------------------------------------------------------
    sub delete_random_gene {
	my $self = shift;
	my @gene_refs = $self->get_genes();
	my $gene_ref = $gene_refs[int rand @gene_refs];
	printn "delete_random_gene: deleting gene ".$gene_ref->get_name() if $verbosity >= 1;
	my @deleted_genes = $self->delete_gene($gene_ref);
	confess "ERROR: internal error -- could not delete a gene" if (@deleted_genes != 1) || ($deleted_genes[0] != $gene_ref);
	return $gene_ref;
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub mutate {
	my $self = shift; my $obj_ID = ident $self;
	my %args = (
		    mutation_rate_params => undef,
		    mutation_rate_global => undef,
		    duplication_rate => undef,
		    deletion_rate => undef,
		    recombination_rate => undef,
		    @_,
		   );
	check_args(\%args,5);
	
	my $mutation_rate_params = $args{mutation_rate_params};
	my $mutation_rate_global = $args{mutation_rate_global};
	my $duplication_rate = $args{duplication_rate};
	my $deletion_rate = $args{deletion_rate};
	my $recombination_rate = $args{recombination_rate};

	my $sequence_ref = $self->get_sequence_ref();

	# pre-mutation parsing
	$self->parse();

	if ($mutation_rate_params > 0.0 && $mutation_rate_params <= 1.0) { # mutate parameters
	  printn "mutate: POINT MUTATION (PARAMS)" if $verbosity >= 1;
	  $self->set_field_mutation_rates(
					  mutate_params_rate => 1.0,
					  mutate_network_rate => 0.0,
					 );
	  my @mutated_list = $self->mutate_genes($mutation_rate_params);
	  my $total_bits = 0; for (my $i = 0; $i < @mutated_list; $i++) {$total_bits += $mutated_list[$i] if $i % 2};
	  my $history = "POINT MUTATION (PARAMS) $total_bits bits in genes ". join ",", @mutated_list;
	  printn $history if $verbosity >= 1;
	  $self->add_history($history);
	}
	elsif ($mutation_rate_params != 0.0) {
	  printn "ERROR: mutation_rate_params is not set in proper range";
	  exit;
	}

	if ($mutation_rate_global > 0.0 && $mutation_rate_global <= 1.0) { # mutate the whole network
	  printn "mutate: POINT MUTATION (GLOBAL)" if $verbosity >= 1;
	  $self->set_field_mutation_rates(
					  mutate_params_rate => 1.0,
					  mutate_network_rate => 1.0,
					 );
	  my @mutated_list = $self->mutate_genome($mutation_rate_global);
	  my $total_bits = 0; for (my $i = 0; $i < @mutated_list; $i++) {$total_bits += $mutated_list[$i] if $i % 2};
	  my $history = "POINT MUTATION (GLOBAL) $total_bits bits in genes ". join ",", @mutated_list;
	  printn $history if $verbosity >= 1;
	  $self->add_history($history);
	}
	elsif ($mutation_rate_global != 0.0) {
	  printn "ERROR: mutation_rate_global is not set in proper range";
	  exit;
	}

	if ($duplication_rate > 0.0 && $duplication_rate <= 1.0) {	# duplicate
	  printn "mutate: DUPLICATION" if $verbosity >= 1;
	  my ($duplicated_gene, $duplicate_start) = $self->duplicate_gene();
	  my $duplicate_name = sprintf("G%04d",$duplicate_start);
	  $self->get_gene_parser_ref()->parse(sequence_ref => $sequence_ref, start_pos => $duplicate_start, dont_clear_flag => 1);  # N.B. this only updates gene_parser_ref not genome_parser_ref
	  my $num_bits = $self->mutate_gene_by_name($duplicate_name, $mutation_rate_global); # mutate duplicated gene
	  my $history = "DUPLICATION of gene $duplicated_gene to $duplicate_name, with $num_bits mutations";
	  printn $history if $verbosity >= 1;
	  $self->add_history($history);
	}
	elsif ($duplication_rate != 0.0) {
	  printn "ERROR: duplication_rate is not set in proper range";
	  exit;
	}

	if ($deletion_rate > 0.0 && $deletion_rate <= 1.0) { # delete a gene
	  printn "mutate: DELETION" if $verbosity >= 1;
	  
	  my $deleted_gene_ref = $self->delete_random_gene();
	  my $deleted_gene_name = $deleted_gene_ref->get_name();
	  my $history = "DELETION of gene $deleted_gene_name";
	  printn $history if $verbosity >= 1;
	  $self->add_history($history);
	}
	elsif ($deletion_rate != 0.0) {
	  printn "ERROR: deletion_rate is not set in proper range";
	  exit;
	}

	if ($recombination_rate > 0.0 && $recombination_rate <= 1.0) { # recombine genes
	  printn "mutate: RECOMBINATION" if $verbosity >= 1;
	  my $gene1_index = $self->pick_random_gene();
	  my $gene2_index = $self->pick_random_gene();
	  my $gene1_name = $self->get_gene_by_index($gene1_index)->get_name();
	  my $gene2_name = $self->get_gene_by_index($gene2_index)->get_name();
	  my $recombinatory_start = $self->recombine_genes($gene1_index, $gene2_index);
	  #$self->get_gene_parser_ref()->parse(sequence_ref => $sequence_ref, start_pos => $recombinatory_start, dont_clear_fag => 1); # N.B. this only updates gene_parser_ref not genome_parser_ref
	  my $history = "RECOMBINATION ($gene1_name, $gene2_name) to G$recombinatory_start";
	  printn $history if $verbosity >= 1;
	  $self->add_history($history);
	}
	elsif ($recombination_rate != 0.0) {
	  printn "ERROR: recombination_rate is not set in proper range";
	  exit;
	}

	# post-mutation parsing
	$self->parse();
      }


    #--------------------------------------------------------------------------------------
    # Function: xxx
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub xxx {
	my $self = shift;
    }
}

sub run_testcases {
    $verbosity = 3;

    srand(333);

    printn "PARSING TESTS";
    my $genome_model_ref = GenomeModel->new({
	name => "GenomeModel",
	Genome => {
	    radius => 1,
	    kf_max => 1e5,
	    kf_min => 0.1,
	    kb_max => 10,
	    kb_min => 0.001,
	    kp_max => 1000,
	    kp_min => 1,
	}
    });
    $genome_model_ref->generate_random_genome(1000);
    $genome_model_ref->parse();
    $genome_model_ref->check();

    printn $genome_model_ref->_DUMP();
    printn $genome_model_ref->get_genome_parser_ref()->_DUMP();
    printn $genome_model_ref->get_gene_parser_ref()->_DUMP();
    printn $genome_model_ref->get_domain_parser_ref()->_DUMP();
    printn $genome_model_ref->get_protodomain_parser_ref()->_DUMP();

    printn $genome_model_ref->get_genome->_DUMP();

    printn $genome_model_ref->get_genome()->sprint(colour_flag => 0, hierarchical_flag => 1);
    printn $genome_model_ref->get_genome()->sprint(hierarchical_flag => 0);

    printn "TERMINATION TEST";
    $genome_model_ref->terminate_genome();

    # parse() will auto-clear, but explicit, redundant clear lets us see DEMOLISH
    $genome_model_ref->clear_parser_instances();

    $verbosity = 0; $genome_model_ref->parse(); $verbosity = 3;
    printn $genome_model_ref->get_genome()->sprint(colour_flag => 0, hierarchical_flag => 1);
    printn $genome_model_ref->get_genome()->sprint(hierarchical_flag => 0);

    printn "DIRECT MUTATION TEST";
    $genome_model_ref->get_sequence_ref()->mutate_subseq(1,0,10);
    printn $genome_model_ref->get_sequence_ref()->get_sequence();
    $genome_model_ref->get_sequence_ref()->mutate_subseq(1,0,10);  # put things back as before

    printn "GENE BY NAME TESTS";
    my @gene_names;
    printn join " ", @gene_names = $genome_model_ref->get_gene_names();
    printn $genome_model_ref->sprint_gene_by_name($gene_names[1], hierarchical_flag => 0);
    printn $genome_model_ref->sprint_gene_by_name($gene_names[1], colour_flag => 0);

    printn "MUTATION TESTS";
    $genome_model_ref->set_field_mutation_rates(
	mutate_params_rate => 1.0,    # mutation rate for parameter fields
	mutate_network_rate => 0.0,   # mutation rate for network structure fields
       );
    printn join " ", @gene_names = $genome_model_ref->get_gene_names();
    printn $genome_model_ref->sprint_genes(hierarchical_flag => 0);
    $verbosity = 0; printn join ",", $genome_model_ref->mutate_genes(0.1); $genome_model_ref->parse(); $verbosity = 3;
    printn join " ", @gene_names = $genome_model_ref->get_gene_names();
    printn $genome_model_ref->sprint_genes(hierarchical_flag => 0);

    printn "DUPLICATION TESTS";
    srand(434653644);
    printn join ",",$genome_model_ref->duplicate_gene();
    printn join ",",$genome_model_ref->duplicate_gene();
#!!! test parsing of only duplicates???
    $verbosity = 0; $genome_model_ref->parse(); $verbosity = 3;
    printn $genome_model_ref->get_genome()->sprint(hierarchical_flag => 0);
    printn $genome_model_ref->get_sequence_ref()->sprint();
    printn $genome_model_ref->sprint_genes(hierarchical_flag => 0);

    printn "RECOMBINATION TESTS";
    printn $genome_model_ref->get_domain_sequence(1,0);
    printn $genome_model_ref->get_domain_sequence(0,0);
    $genome_model_ref->recombine_genes(1,0);
    $verbosity = 0; $genome_model_ref->parse(); $verbosity = 3;
    printn $genome_model_ref->sprint_genes(hierarchical_flag => 0);
    printn $genome_model_ref->get_genome()->sprint(hierarchical_flag => 0);
    printn $genome_model_ref->get_sequence_ref()->sprint();

    printn "STORABLE TEST";
    use Storable;
    my $ice_ref = Storable::freeze($genome_model_ref);
    my $water_ref = Storable::thaw($ice_ref);
    printn $genome_model_ref->_DUMP();
    printn $water_ref->_DUMP();
}

# Package BEGIN must return true value
return 1;

    #--------------------------------------------------------------------------------------
## Function: add_gene
## Synopsys: Append a gene to end of genome, given its sequence
    #--------------------------------------------------------------------------------------
#sub add_gene {
#    my $gene_sequence = $_[0];
#    my $new_gene_start = length($sdb->{genome}{sequence});
#    $sdb->{genome}{sequence} = $sdb->{genome}{sequence} . $gene_sequence;
#    printn "add_gene: created gene P$new_gene_start";
#    return $new_gene_start;
#}

    #--------------------------------------------------------------------------------------
## Function: decimate_genes
## Synopsys: Delete parsed genes with a given probability
    #--------------------------------------------------------------------------------------
#sub decimate_genes {
#    my $probability = $_[0];

#    my @return_list = ();

#    if (!defined $probability) {
#	printn "ERROR: delete_parsed_genes -- probability of deletion required";
#	exit;
#    }

#    my @gene_list = get_genes();
#    my $deleted = 0;
#    my $gene;
#    foreach $gene (@gene_list) {
#	if (rand() < $probability) {
#	    printn "delete_parsed_genes: deleting gene $gene";
#	    push @return_list, $gene;
#	    delete_gene($gene);
#	    $deleted++;
#	}
#    }
#    printn "delete_parsed_genes: deleted $deleted genes (".join ",",@return_list.")";
#    return @return_list;
#}


    #--------------------------------------------------------------------------------------
## Function: replace_protodomain_sequence
## Synopsys: Swap out old protodomain sequence with a new one
    #--------------------------------------------------------------------------------------
#sub replace_protodomain_sequence {
#    my $protodomain = shift;
#    my $sequence = shift;

#    my $locus = $sdb->{protodomain_table}{$protodomain}{locus};
#    substr($sdb->{genome}{sequence}, $locus, $protodomain_sequence_length) = $sequence;
#}

