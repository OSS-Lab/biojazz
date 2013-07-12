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
    my %number_of       :ATTR(get => 'number', set => 'number', default => 0);
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
        $number_of{$obj_ID} = 1;
        $score_of{$obj_ID} = 0;
        $elite_flag_of{$obj_ID} = 0;

#	$cell_volume_of{$obj_ID} = $arg_ref->{cell_volume} if exists $arg_ref->{cell_volume};
    }
    
    sub DEMOLISH {
        my $self = shift;

        printn "Genome is detroyed.\n" if $verbosity > 3;
    }

    #===  FUNCTION  ================================================================
    #         NAME: static_analyse
    #      PURPOSE: anaylse the static information of genome and network
    #   PARAMETERS: number of genes, domains, protodomains, rules etc.
    #      RETURNS: 
    #  DESCRIPTION: ????
    #       THROWS: no exceptions
    #     COMMENTS: none
    #     SEE ALSO: n/a
    #===============================================================================
    sub static_analyse {
        my $self = shift; my $obj_ID = ident $self;
        my $rescore_elite = shift;

        my $elite_flag = $self->get_elite_flag();
        if ($elite_flag) {
            printn "This genome has been static analysed before.";
            return if ($rescore_elite == 0);
            printn "Now re-static analyse this genome.";
        }

        my $genome_ref = $self->get_genome();
        my @gene_refs = $self->get_genes();
        my @domain_refs = $self->get_domain_parser_ref()->get_object_instances();
        my @protodomain_refs = $self->get_protodomain_parser_ref()->get_object_instances();

        $stats_ref_of{$obj_ID}->{num_genes} = scalar @gene_refs;
        $stats_ref_of{$obj_ID}->{num_domains} = scalar @domain_refs;
        $stats_ref_of{$obj_ID}->{num_protodomains} = scalar @protodomain_refs;
        # analyse protodomains
        my $num_MDs = 0;
        my $num_CDs = 0;
        my $num_BDs = 0;
        my $num_KDs = 0;
        my $num_PDs = 0;
        my $num_allosteric_protodomains = 0;
        foreach my $protodomain_ref (@protodomain_refs) {
            my $protodomain_translation_ref = $protodomain_ref->get_translation_ref();
            if ($protodomain_translation_ref->{type} eq "msite" ) {
                $num_MDs++;
            } elsif ($protodomain_translation_ref->{type} eq "csite") {
                $num_CDs++;
                if ($protodomain_translation_ref->{substrate_polarity} == 0) {
                    $num_KDs++;
                } elsif ($protodomain_translation_ref->{substrate_polarity} == 1) {
                    $num_PDs++;
                }
            } elsif ($protodomain_translation_ref->{type} eq "bsite") {
                $num_BDs++;
            }

            if ($protodomain_ref->get_allosteric_flag() == 1) {
                $num_allosteric_protodomains++;
            }
        }

        # analyse domains
        my $num_allosteric_domains = 0;
        foreach my $domain_ref (@domain_refs) {
            my $domain_translation_ref = $domain_ref->get_translation_ref();
            if ($domain_translation_ref->{allosteric_flag} == 1) {
                $num_allosteric_domains++;
            }
        }


        @{$stats_ref_of{$obj_ID}}{ "num_phosphorylation_protodomains", 
        "num_catalytic_protodomains", "num_kinase_protodomains", "num_phosphatase_protodomains",
        "num_binding_protodomains", "num_allosteric_protodomains", "num_allosteric_domains",
        } = ($num_MDs, $num_CDs, $num_KDs, $num_PDs, $num_BDs, $num_allosteric_protodomains, 
            $num_allosteric_domains);

        # number of rules
        # N.B.: currently has some problem need to revise
        # $stats_ref_of{$obj_ID}->{num_rules} = scalar @{$genome_ref->get_rules()};

    } ## --- end sub static_analyse

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
        $str .= ", number=".(defined $self->get_number() ? $self->get_number() : "UNDEF");
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
            printn "The genome is not parsered successfully, no genome_parser_ref returned!" if $verbosity > 2;
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


    #---  FUNCTION  ----------------------------------------------------------------
    #         NAME: get_domain_by_index
    #   PARAMETERS: gene_index, domain_index
    #      RETURNS: domain reference
    #     COMMENTS: none
    #-------------------------------------------------------------------------------
    sub get_domain_by_index {
        my	( $self, $gene_index, $domain_index )	= @_;
        
        return $self->get_gene_by_index($gene_index)->get_field_ref(["domains", $domain_index]);
    } ## --- end sub get_domain_by_index

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
    # Function: get_protodomains
    # Synopsys: Given a gene index, return the list of protodomain ParserInstances.
    #--------------------------------------------------------------------------------------
    sub get_protodomains {
        my $self = shift;
        my $gene_index = shift;
        my @domain_refs = $self->get_domains($gene_index);
        my @protodomain_refs;
        foreach my $domain_ref (@domain_refs) {
            push (@protodomain_refs, @{$domain_ref->get_field_ref(["protodomains"])});
        }

        return @protodomain_refs;
    }


    #---  FUNCTION  ----------------------------------------------------------------
    #         NAME: get_protodomain_by_index
    #   PARAMETERS: gene_index, domain_index, protodomain_index
    #      RETURNS: protodomain reference
    #     COMMENTS: none
    #-------------------------------------------------------------------------------
    sub get_protodomain_by_index {
        my	( $self, $gene_index, $domain_index, $protodomain_index )	= @_;
        
        return $self->get_domain_by_index($gene_index, $domain_index)->get_field_ref(["protodomains", $protodomain_index]);
    } ## --- end sub get_protodomain_by_index


    #--------------------------------------------------------------------------------------
    # Function: get_protodomain_sequence
    # Synopsys: Given a gene, domain and protodomain number, extract the sequence from genome.
    #--------------------------------------------------------------------------------------
    sub get_protodomain_sequence {
        my $self = shift;
        my $gene_index = shift;
        my $domain_index = shift;
        my $protodomain_index = shift;

        return $self->get_domain_by_index($gene_index, $domain_index)->get_field_sequence(["protodomains", $protodomain_index]);
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
        my $mutated_junk_bits = $self->get_genome()->mutate($mutation_rate);
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
    # Function: duplicate_random_gene
    # Synopsys: Pick a gene at random and duplicate it (copy to end of genome) with
    #           probability given by argument (or always if no argument).
    #           Returns name of chosen gene and locus of duplicate.
    #--------------------------------------------------------------------------------------
    sub duplicate_random_gene {
        my $self = shift; my $obj_ID = ident $self;
        my $probability = shift;
        $probability = 1.0 if (!defined $probability);

        if ((rand 1) > $probability) {
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
    # Function: duplicate_gene
    # Synopsys: duplicate a gene by its name passed from arguments, if there is no
    #           specified gene name, then it will randomly pick a gene to duplicate
    #--------------------------------------------------------------------------------------
    sub duplicate_gene {
        my $self = shift; my $obj_ID = ident $self;
        my $gene_name = shift;
        my $gene_ref = !defined $gene_name ? ($self->get_gene_by_index($self->pick_random_gene())) : ($self->get_gene_by_name($gene_name));

        my $sequence_ref = $self->get_sequence_ref();
        my $duplicate_locus = $sequence_ref->get_length();

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
                printn "delete_gene: deleting $gene_name start=$gene_start stop=$gene_stop" if $verbosity > 2;
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
    # Function: duplicate_domain
    # Synopsys: duplication of domain is a bit more difficult to implement. Based on the
    #           passed to, try to iterate each domain to decide whether duplicate it or not
    #           and either append it at the end or put it in tandem. 
    #           NOTE: it could easily messup the whole gene, so pay attention
    #--------------------------------------------------------------------------------------
    sub duplicate_domain {
        my $self = shift; my $obj_ID = ident $self;
        my $gene_index = shift;

        if (!defined $gene_index) {
            confess "There is no gene index specified, have no idea how to do domain deletion!";
        }
        my $gene_ref = $self->get_gene_by_index($gene_index);
        my $sequence_ref = $self->get_sequence_ref();
        
        my $gene_name = $gene_ref->get_name();
        my @gene_protodomains = $self->get_protodomains($gene_index);
        my $gene_num_protodomains = scalar @gene_protodomains;
        my $gene_pd1 = int rand ($gene_num_protodomains + 1);
        my $gene_pd2 = int rand ($gene_num_protodomains + 1);
        my $duplicate_num = abs($gene_pd1 - $gene_pd2);

        my $gene_site1; my $gene_site2;
        if ($gene_pd1 == $gene_num_protodomains && $gene_pd2 == $gene_num_protodomains) {
            return $duplicate_num;
        } elsif ($gene_pd1 == $gene_num_protodomains && $gene_pd2 != $gene_num_protodomains) {
            $gene_site1 = $gene_protodomains[$gene_pd1 - 1]->get_locus() + $gene_protodomains[$gene_pd1 - 1]->get_length();
            $gene_site2 = $gene_protodomains[$gene_pd2]->get_locus() + $gene_protodomains[$gene_pd2]->get_length();
        } elsif ($gene_pd1 != $gene_num_protodomains && $gene_pd2 == $gene_num_protodomains) {
            $gene_site1 = $gene_protodomains[$gene_pd1]->get_locus() + $gene_protodomains[$gene_pd1]->get_length();
            $gene_site2 = $gene_protodomains[$gene_pd2 - 1]->get_locus() + $gene_protodomains[$gene_pd2 - 1]->get_length();
        } else {
            $gene_site1 = $gene_protodomains[$gene_pd1]->get_locus();
            $gene_site2 = $gene_protodomains[$gene_pd2]->get_locus();
        }

        my @gene_locus_unsorted = ($gene_site1, $gene_site2);
        my @gene_locus = sort {$a <=> $b} @gene_locus_unsorted;
        my $middle_sequence = $sequence_ref->get_subseq($gene_locus[0], ($gene_locus[1] - $gene_locus[0]));
        $sequence_ref->splice_subseq($middle_sequence, $gene_locus[1]);
 
        return $duplicate_num;
    }


    #--------------------------------------------------------------------------------------
    # Function: delete_domain
    # Synopsys: delete domain basicly based on the doamin deletion rate and gene index
    # 	    return number of protodomains that be deleted.
    #           NOTE: it could easily messup the whole gene, so pay attention
    #--------------------------------------------------------------------------------------
    sub delete_domain {
        my $self = shift; my $obj_ID = ident $self;
        my $gene_index = shift;

        if (!defined $gene_index) {
            confess "There is no gene index specified, have no idea how to do domain deletion!";
        }
        my $gene_ref = $self->get_gene_by_index($gene_index);
        my $sequence_ref = $self->get_sequence_ref();
        
        my $gene_name = $gene_ref->get_name();
        my @gene_protodomains = $self->get_protodomains($gene_index);
        my $gene_num_protodomains = scalar @gene_protodomains;
        my $gene_pd1 = int rand ($gene_num_protodomains + 1);
        my $gene_pd2 = int rand ($gene_num_protodomains + 1);
        my $delete_num = abs($gene_pd1 - $gene_pd2);

        if ($delete_num == $gene_num_protodomains) {
            return (- $delete_num);
        }

        my $gene_site1; my $gene_site2;
        if ($gene_pd1 == $gene_num_protodomains && $gene_pd2 == $gene_num_protodomains) {
            return $delete_num;
        } elsif ($gene_pd1 == $gene_num_protodomains && $gene_pd2 != $gene_num_protodomains) {
            $gene_site1 = $gene_protodomains[$gene_pd1 - 1]->get_locus() + $gene_protodomains[$gene_pd1 - 1]->get_length();
            $gene_site2 = $gene_protodomains[$gene_pd2]->get_locus() + $gene_protodomains[$gene_pd2]->get_length();
        } elsif ($gene_pd1 != $gene_num_protodomains && $gene_pd2 == $gene_num_protodomains) {
            $gene_site1 = $gene_protodomains[$gene_pd1]->get_locus() + $gene_protodomains[$gene_pd1]->get_length();
            $gene_site2 = $gene_protodomains[$gene_pd2 - 1]->get_locus() + $gene_protodomains[$gene_pd2 - 1]->get_length();
        } else {
            $gene_site1 = $gene_protodomains[$gene_pd1]->get_locus();
            $gene_site2 = $gene_protodomains[$gene_pd2]->get_locus();
        }
        
        my @gene_locus_unsorted = ($gene_site1, $gene_site2);
        my @gene_locus = sort {$a <=> $b} @gene_locus_unsorted;
        $sequence_ref->splice_subseq("", $gene_locus[0], ($gene_locus[1] - $gene_locus[0]));
 
        return $delete_num;
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
    # IMPL: Now we try this method. Try to find the sequence from one gene from beginning
    #       to a certain protodomain, then find another sequence from another gene from one
    #       certain domain to the end of it. try to find at least one of the protodomains 
    #       in each gene.
    #--------------------------------------------------------------------------------------
    sub recombine_genes {
        my $self = shift;
        my $obj_ID = ident $self;
        my ($gene1_index, $gene2_index) = (shift, shift);

        my $sequence_ref = $self->get_sequence_ref();
        
        my $gene1_ref = $self->get_gene_by_index($gene1_index);
        my $gene2_ref = $self->get_gene_by_index($gene2_index);
        my $gene1_name = $gene1_ref->get_name();
        my $gene2_name = $gene2_ref->get_name();

        my @gene1_protodomains = $self->get_protodomains($gene1_index);
        my @gene2_protodomains = $self->get_protodomains($gene2_index);
        my $gene1_num_protodomains = scalar @gene1_protodomains;
        my $gene2_num_protodomains = scalar @gene2_protodomains;

        printn "recombine_genes: recombining $gene1_name($gene1_index, $gene1_num_protodomains protodomains) and $gene2_name($gene2_index, $gene2_num_protodomains protodomains)" if $verbosity > 1;
        printn "recombine_genes: $gene1_name sequence = " . $gene1_ref->get_sequence() if $verbosity > 3;
        printn "recombine_genes: $gene2_name sequence = " . $gene2_ref->get_sequence() if $verbosity > 3;

        my $gene1_pd1 = int rand ($gene1_num_protodomains);
        my $gene1_pd2 = int rand ($gene1_num_protodomains);
        my $gene2_pd1 = int rand ($gene2_num_protodomains);
        my $gene2_pd2 = int rand ($gene2_num_protodomains);

        my $gene1_site1;
        my $gene1_site2;
        my $gene2_site1;
        my $gene2_site2;
        
        if ($gene1_pd1 == ($gene1_num_protodomains - 1) || $gene1_pd2 == ($gene1_num_protodomains - 1)
            || $gene2_pd1 == ($gene2_num_protodomains -1) || $gene2_pd2 == ($gene2_num_protodomains -1)) {
            $gene1_site1 = $gene1_protodomains[$gene1_pd1]->get_locus() + $gene1_protodomains[$gene1_pd1]->get_length();
            $gene1_site2 = $gene1_protodomains[$gene1_pd2]->get_locus() + $gene1_protodomains[$gene1_pd2]->get_length();
            $gene2_site1 = $gene2_protodomains[$gene2_pd1]->get_locus() + $gene2_protodomains[$gene2_pd1]->get_length();
            $gene2_site2 = $gene2_protodomains[$gene2_pd2]->get_locus() + $gene2_protodomains[$gene2_pd2]->get_length();

        } else {
            if ($gene1_pd1 == 0 || $gene1_pd2 == 0 || $gene2_pd1 == 0 || $gene2_pd2 == 0) {
                $gene1_site1 = $gene1_protodomains[$gene1_pd1]->get_locus();
                $gene1_site2 = $gene1_protodomains[$gene1_pd2]->get_locus();
                $gene2_site1 = $gene2_protodomains[$gene2_pd1]->get_locus();
                $gene2_site2 = $gene2_protodomains[$gene2_pd2]->get_locus();
            } else {
                if (rand() < 0.5) {
                    $gene1_site1 = $gene1_protodomains[$gene1_pd1]->get_locus();
                    $gene1_site2 = $gene1_protodomains[$gene1_pd2]->get_locus();
                    $gene2_site1 = $gene2_protodomains[$gene2_pd1]->get_locus();
                    $gene2_site2 = $gene2_protodomains[$gene2_pd2]->get_locus();
                } else {
                    $gene1_site1 = $gene1_protodomains[$gene1_pd1]->get_locus() + $gene1_protodomains[$gene1_pd1]->get_length();
                    $gene1_site2 = $gene1_protodomains[$gene1_pd2]->get_locus() + $gene1_protodomains[$gene1_pd2]->get_length();
                    $gene2_site1 = $gene2_protodomains[$gene2_pd1]->get_locus() + $gene2_protodomains[$gene2_pd1]->get_length();
                    $gene2_site2 = $gene2_protodomains[$gene2_pd2]->get_locus() + $gene2_protodomains[$gene2_pd2]->get_length();

                }
            }
        }

        my @gene1_locus_unsorted = ($gene1_site1, $gene1_site2);
        my @gene2_locus_unsorted = ($gene2_site1, $gene2_site2);
        my @gene1_locus = sort {$a <=> $b} @gene1_locus_unsorted;
        my @gene2_locus = sort {$a <=> $b} @gene2_locus_unsorted;

        my $middle_sequence = $sequence_ref->get_subseq($gene1_locus[0], ($gene1_locus[1]-$gene1_locus[0]));
        my $gene2_orig = $gene2_ref->get_locus();
        my $gene2_end = $gene2_orig + $gene2_ref->get_length();
        my $left_sequence = $sequence_ref->get_subseq($gene2_orig, ($gene2_locus[0] - $gene2_orig));
        my $right_sequence = $sequence_ref->get_subseq($gene2_locus[1], ($gene2_end - $gene2_locus[1]));

        my $gene_parser_ref = $self->get_gene_parser_ref();
        my $new_gene_sequence = ($left_sequence .
            $middle_sequence .
            $right_sequence . 
            $gene_parser_ref->get_STOP_linker_code());

        my $new_gene_start = $sequence_ref->get_length();

        $self->get_sequence_ref()->splice_subseq("00000000".$new_gene_sequence);
        printn "recombine_genes: new gene sequence = " . $new_gene_sequence if $verbosity > 1;

        return $new_gene_start;
    }

    #--------------------------------------------------------------------------------------
    # Function: mutate
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub mutate 
    {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            mutation_rate_params => undef,
            mutation_rate_global => undef,
            gene_duplication_rate => undef,
            gene_deletion_rate => undef,
            domain_duplication_rate => undef,
            domain_deletion_rate => undef,
            recombination_rate => undef,
            @_,
        );
        check_args(\%args,7);

        my $mutation_rate_params = $args{mutation_rate_params};
        my $mutation_rate_global = $args{mutation_rate_global};
        my $gene_duplication_rate = $args{gene_duplication_rate};
        my $gene_deletion_rate = $args{gene_deletion_rate};
        my $domain_duplication_rate = $args{domain_duplication_rate};
        my $domain_deletion_rate = $args{domain_deletion_rate};
        my $recombination_rate = $args{recombination_rate};


        my $mutation_count = 0;

        # pre-mutation parsing
        $self->parse();
        $self->check();
        my $num_genes = $self->get_num_genes();
        ###########################
        # gene duplication
        # by appending genes based on gene_duplication_rate
        if ($gene_duplication_rate > 0.0 && $gene_duplication_rate <= 1.0) 	# duplicate genes
        {
            my @gene_refs = $self->get_genes();
            my @gene_names = map $_->get_name(), @gene_refs;

            my $gene_duplication_count = 0;
            foreach my $gene_name (@gene_names) {
                if ( rand() < $gene_duplication_rate ) {
                    printn "mutate: GENE_DUPLICATION" if $verbosity >= 1;
                    my ($duplicated_gene, $duplicate_start) = $self->duplicate_gene($gene_name);

                    my $duplicate_name = sprintf("G%04d",$duplicate_start);
                    my $history = "DUPLICATION of gene $duplicated_gene"; 
                    printn $history if $verbosity >= 1;
                    $self->add_history($history);

                    # post-mutation parsing
                    $self->parse();
                    $gene_duplication_count++;
                }
            }
            $mutation_count += $gene_duplication_count;
        }
        elsif ($gene_duplication_rate != 0.0) {
            printn "ERROR: gene_duplication_rate is not set in proper range";
            exit;
        }
        $self->check();
 
        #######################
        # Domain shuffling
        if ($recombination_rate > 0.0 && $recombination_rate <= 1.0) { # recombine genes
            for (my $i = 0; $i < $num_genes; $i ++) {
                if (rand() < $recombination_rate) {
                    printn "mutate: RECOMBINATION" if $verbosity >= 1;
                    my $gene1_index = int rand $num_genes;
                    my $gene2_index = int rand $num_genes;
                    my $gene1_name = $self->get_gene_by_index($gene1_index)->get_name();
                    my $gene2_name = $self->get_gene_by_index($gene2_index)->get_name();
                    my $recombinatory_start = $self->recombine_genes($gene1_index, $gene2_index);
                    my $history = "RECOMBINATION ($gene1_name, $gene2_name) to G$recombinatory_start";
                    printn $history if $verbosity >= 1;
                    $self->add_history($history);

                    # post-mutation parsing
                    $self->parse();
                    $mutation_count++;
                }
            }
        }
        elsif ($recombination_rate != 0.0) {
            printn "ERROR: recombination_rate is not set in proper range";
            exit;
        }
        $self->check();
 
        #######################
        # DOMAIN DUPLICATION

        # choose a gene to duplicate it, and duplicate domains in it, then 
        # append the gene to the end and delete the orginal one. 
        if ($domain_duplication_rate > 0.0 && $domain_duplication_rate <= 1.0) 	# duplicate domains
        {

            my @gene_refs = $self->get_genes();
            my $domain_duplication_count = 0;

            for (my $i = 0; $i < $num_genes; $i ++) {
                my $gene_ref = $self->get_gene_by_index($i);
                my $gene_name = $gene_ref->get_name();
                if ($gene_name ne ($gene_refs[$i]->get_name())) {
                    confess "The index of gene instances are not consistent with index of gene_refs array";
                }
                if (rand() < $domain_duplication_rate) {
                    my $duplicate_num = $self->duplicate_domain($i);
                    printn "Mutate: duplicated $duplicate_num domains in gene $gene_name" if $verbosity > 1;

                    # update the parser instances and gene_refs after reparsing
                    if ($duplicate_num) {
                        # post domain_duplication parsing
                        $self->parse();
                        undef @gene_refs;
                        @gene_refs = $self->get_genes();
                        $domain_duplication_count += $duplicate_num;
                    }
                }
            }

            $mutation_count += $domain_duplication_count;

        }
        elsif ($domain_duplication_rate != 0.0) 
        {
            printn "ERROR: domain_duplication_rate is not set in proper range";
            exit;
        }
        $self->check();
 
        #######################
        # DOMAIN DELETION

        my @gene_need_delete_indice;
        if ($domain_deletion_rate > 0.0 && $domain_deletion_rate <= 1.0)  # delete domains
        {
            my @gene_refs = $self->get_genes();
            for (my $i = 0; $i < $num_genes; $i ++) {
                my $gene_ref = $self->get_gene_by_index($i);
                my $gene_name = $gene_ref->get_name();
                if ($gene_name ne $gene_refs[$i]->get_name()) {
                    confess "The index of gene instances are not consistent with index of gene_refs array";
                }
                if (rand() < $domain_deletion_rate) {
                    my $deleted_num = $self->delete_domain($i);
                    printn "Mutate: deleted $deleted_num domains in gene $gene_name" if $verbosity > 1;

                    # update the parser instances and gene_refs after reparsing
                    if ($deleted_num < 0) {
                        push(@gene_need_delete_indice, $i);
                    }
                    elsif ($deleted_num > 0) {
                        # post domain_duplication parsing
                        $self->parse();
                        undef @gene_refs;
                        @gene_refs = $self->get_genes();
                        $mutation_count += $deleted_num;
                    }
                }
            }
        }
        elsif ($domain_deletion_rate != 0.0) 
        {
            printn "ERROR: domain_deletion_rate is not set in proper range";
            exit;
        }
        $self->check();
 
        ########################
        # GENE_DELETION
        if ($gene_deletion_rate > 0.0 && $gene_deletion_rate <= 1.0)  # delete genes
        {
            printn "mutate: GENE_DELETION" if $verbosity > 1;

            my $num_genes = $self->get_num_genes();
            my $deleted_gene_num = 0;

            if (scalar @gene_need_delete_indice > 0 && scalar @gene_need_delete_indice < $num_genes) {
                my @gene_indice = sort {$b <=> $a} @gene_need_delete_indice;
                foreach my $gene_index (@gene_indice) {
                    my $gene_need_delete = $self->get_gene_by_index($gene_index - $deleted_gene_num);
                    $self->delete_gene($gene_need_delete);
                    my $deleted_gene_name = $gene_need_delete->get_name();
                    my $history = "DELETION of gene $deleted_gene_name";
                    printn $history if $verbosity >= 1;
                    $self->add_history($history);

                    # post-mutation parsing
                    $self->parse();

                    $deleted_gene_num++;
                }


                my $num_genes_left = $num_genes - scalar @gene_need_delete_indice;
                my $current_num_genes = $self->get_num_genes();
                confess "The number of genes after deletion is not consistant with of left ones." if ($current_num_genes != $num_genes_left);
                my $deleted_gene_num2 = 0;
                GENE_DELETION: {
                    for (my $i = 0; $i < $current_num_genes; $i++) {
                        if ($current_num_genes - $deleted_gene_num2 == 1) {
                            last GENE_DELETION;
                        }
                        elsif (rand() < $gene_deletion_rate) {
                            my $deleted_gene_ref = $self->delete_random_gene();
                            my $deleted_gene_name = $deleted_gene_ref->get_name();
                            my $history = "DELETION of gene $deleted_gene_name";
                            printn $history if $verbosity >= 1;
                            $self->add_history($history);

                            # post-mutation parsing
                            $self->parse();

                            # conuting the deleted gene number to calculate number of rest genes
                            $deleted_gene_num2++;
                        }
                    }
                }

                $mutation_count += $deleted_gene_num;
                $mutation_count += $deleted_gene_num2;


            } elsif (@gene_need_delete_indice == $num_genes) {
                DELETE_GENES: {
                    for (my $i = 0; $i < $num_genes; $i++) {
                        if ($num_genes - $deleted_gene_num == 1) {
                            last DELETE_GENES;
                        }
                        elsif (rand() < $gene_deletion_rate) {
                            my $deleted_gene_ref = $self->delete_random_gene();
                            my $deleted_gene_name = $deleted_gene_ref->get_name();
                            my $history = "DELETION of gene $deleted_gene_name";
                            printn $history if $verbosity >= 1;
                            $self->add_history($history);

                            # post-mutation parsing
                            $self->parse();

                            # conuting the deleted gene number to calculate number of rest genes
                            $deleted_gene_num++;
                        }
                    }
                }

                $mutation_count += $deleted_gene_num;
            }
        }
        elsif ($gene_deletion_rate != 0.0) {
            printn "ERROR: gene_deletion_rate is not set in proper range";
            exit;
        }
        $self->check();
 
 
        ########################
        if ($mutation_rate_params > 0.0 && $mutation_rate_params <= 1.0)  # mutate parameters
        {
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

            # post-mutation parsing
            $self->parse();
            $mutation_count += $total_bits;
        }
        elsif ($mutation_rate_params != 0.0) 
        {
            printn "ERROR: mutation_rate_params is not set in proper range";
            exit;
        }
        $self->check();
 
        ########################
        if ($mutation_rate_global > 0.0 && $mutation_rate_global <= 1.0)  # mutate the whole network
        {
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

            # post-mutation parsing
            $self->parse();
            $mutation_count += $total_bits;
        }
        elsif ($mutation_rate_global != 0.0) 
        {
            printn "ERROR: mutation_rate_global is not set in proper range";
            exit;
        }
        $self->check();
 

        return $mutation_count;
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

