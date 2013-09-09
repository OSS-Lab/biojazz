#-#####################################################################################
#- File:     Network.pm
#- Synopsys: Package of routines for interaction network analysis and searching.
#-           Uses adjacency matrices where element i,j is 1 if nodes are
#-           connected, 0 otherwise.
#-
#- REF: http://www.math.ucdavis.edu/~daddel/linear_algebra_appl/Applications/GraphTheory/GraphTheory_9_17/node9.html
#-
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Network;
use Class::Std::Storable;
use base qw();
{
    use Carp;
    use Storable qw(dclone);

    use Utils;
    use Globals;

    use Matrix;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %adjacency_matrix_ref_of :ATTR(get => 'adjacency_matrix_ref', set => 'adjacency_matrix_ref');
    my %adjacency_matrix_node_refs_of :ATTR(get => 'adjacency_matrix_node_refs', set => 'adjacency_matrix_node_refs');
    my %adjacency_matrix_node_states_of :ATTR(get => 'adjacency_matrix_node_states', set => 'adjacency_matrix_node_states');

    my %connectivity_matrix_ref_of :ATTR(get => 'connectivity_matrix_ref', set => 'connectivity_matrix_ref');

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
    # Function: reset_network
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub reset_network {
        my $self = shift; my $obj_ID = ident $self;

        $adjacency_matrix_ref_of{$obj_ID}{protodomains} = [Matrix->new({})];
        $adjacency_matrix_node_refs_of{$obj_ID}{protodomains} = [];
        $adjacency_matrix_node_states_of{$obj_ID}{protodomains} = [];

        $adjacency_matrix_ref_of{$obj_ID}{domains} = [Matrix->new({})];
        $adjacency_matrix_node_refs_of{$obj_ID}{domains} = [];
        $adjacency_matrix_node_states_of{$obj_ID}{domains} = [];

        $adjacency_matrix_ref_of{$obj_ID}{genes} = [Matrix->new({})];
        $adjacency_matrix_node_refs_of{$obj_ID}{genes} = [];
        $adjacency_matrix_node_states_of{$obj_ID}{genes} = [];
    }

    #--------------------------------------------------------------------------------------
    # Function: STORABLE_freeze_pre
    # Synopsys: hook provided by Class::Std::Storable
    #--------------------------------------------------------------------------------------
    sub STORABLE_freeze_pre: CUMULATIVE {
        my ($self, $clone_flag) = @_;

        # we reset the network because it contains refs to ParserInstances which
        # get cleared when Parser object is stored
        $self->reset_network();
    }

    #--------------------------------------------------------------------------------------
    # Function: build_network
    # Synopsys: Build adjacency matrix for protodomains, domains and genes.
    #--------------------------------------------------------------------------------------
    sub build_network {
        my $self = shift; my $obj_ID = ident $self;

        $self->reset_network();

        # any msites must be included twice, i.e. once for each state
        my @temp_refs = $self->get_gene_parser_ref()->get_domain_parser_ref()->get_protodomain_parser_ref()->get_object_instances();
        @temp_refs = map {$_->get_translation_ref()->{type} eq "msite" ? ([$_,0],[$_,1]) : [$_,0]} @temp_refs;

        my @protodomain_refs = map {$_->[0]} @temp_refs;
        my @protodomain_states = map {$_->[1]} @temp_refs;

        my @is_substrate_flag = map { 0 } @protodomain_refs;
        my @is_product_flag   = map { 0 } @protodomain_refs;

        my %rules_cache = ();

        # INIT AND BUILD PROTODOMAIN ADJACENCY MATRIX
        my $pd_size = @protodomain_refs;
        $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($pd_size-1,$pd_size-1, 0);
        for (my $i=0; $i < @protodomain_refs; $i++) {
            my $x_ref = $protodomain_refs[$i];
            my $x_state = $protodomain_states[$i];
            my $x_translation_ref = $x_ref->get_translation_ref();
            for (my $j=$i; $j < @protodomain_refs; $j++) {
                my $y_ref = $protodomain_refs[$j];
                my $y_state = $protodomain_states[$j];
                my $y_translation_ref = $y_ref->get_translation_ref();
                my @rules = ProtoDomainInstance->create_canbindrules(
                    x_ref => $x_ref,
                    y_ref => $y_ref,
                    radius => $self->get_radius(),
                    kf_max => $self->get_kf_max(),
                    kf_min => $self->get_kf_min(),
                    kb_max => $self->get_kb_max(),
                    kb_min => $self->get_kb_min(),
                    kp_max => $self->get_kp_max(),
                    kp_min => $self->get_kp_min(),
                );
                $rules_cache{$x_ref}{$y_ref} = \@rules;

                if (@rules) {
                    my $catalytic_flag = ((($x_translation_ref->{type} eq "csite") && ($y_translation_ref->{type} eq "msite")) ||
                        (($y_translation_ref->{type} eq "csite") && ($x_translation_ref->{type} eq "msite"))) ? 1 : 0;

                    if (!$catalytic_flag) {
                        # set symmetric elements if x binds to y
                        $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($i,$j, 1);
                        $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($j,$i, 1);
                    } else {
                        if ($x_translation_ref->{type} eq "csite") {
                            confess "ERROR: unexpected condition" if $y_translation_ref->{type} ne "msite";
                            if ($x_translation_ref->{substrate_polarity} == $y_state) {
                                # set symmetric elements if x binds to y
                                $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($i,$j, 1);
                                $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($j,$i, 1);
                                # set uni-arrow from substrate to product
                                my $jj = ($y_state == 1) ? $j-1 : $j+1;
                                $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($j,$jj, 1);
                                $is_substrate_flag[$j] = 1;
                                $is_product_flag[$jj] = 1;
                            }
                        } elsif ($y_translation_ref->{type} eq "csite") {
                            confess "ERROR: unexpected condition" if $x_translation_ref->{type} ne "msite";
                            if ($y_translation_ref->{substrate_polarity} == $x_state) {
                                # set symmetric elements if x binds to y
                                $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($i,$j, 1);
                                $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($j,$i, 1);
                                # set uni-arrow from substrate to product
                                my $ii = ($x_state == 1) ? $i-1 : $i+1;
                                $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($i,$ii, 1);
                                $is_substrate_flag[$i] = 1;
                                $is_product_flag[$ii] = 1;
                            }
                        } else {
                            confess "ERROR: unexpected condition";
                        }
                    }
                }

                my $x_domain_ref = $protodomain_refs[$i]->get_upper_ref();
                my $y_domain_ref = $protodomain_refs[$j]->get_upper_ref();

                # set symmetric elements if x is in same domain as y and this domain is allosteric
                if ($i != $j) {
                    if (($x_domain_ref == $y_domain_ref) && ($x_domain_ref->get_translation_ref->{allosteric_flag})) {
                        $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($i,$j, 1);
                        $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->set_element($j,$i, 1);
                    }
                }
            }
        }

        # any msite is potentially prunable because we don't know
        #   i) only msite states of 0 have initial value
        #  ii) it takes a kinase to create msite_state of 1
        # iii) if they will exist at steady-state (e.g. if only the kinase exists, they msite_state of 0 will go to nil concentration)
        my @protodomain_prunable_flags = map {$_->get_translation_ref()->{type} eq "msite" ? 1 : 0} @protodomain_refs;
        for (my $i=0; $i < @protodomain_prunable_flags; $i++) {
            if ($protodomain_prunable_flags[$i]) {
                if ($protodomain_states[$i] == 1) {
                    $protodomain_prunable_flags[$i] = 0 if $is_product_flag[$i];
                } else {
                    $protodomain_prunable_flags[$i] = 0 if $is_product_flag[$i] || !$is_substrate_flag[$i];
                }
            }
        }

        # now we prune any msite nodes whose concentration will be zero at steady-state
        my @prunable_msite_indices = ();
        map {push @prunable_msite_indices, $_ if $protodomain_prunable_flags[$_]} (0..$#protodomain_prunable_flags);
        my $num_pruned = 0;
        for (my $i=0; $i < @prunable_msite_indices; $i++) {
            my $prune_index = $prunable_msite_indices[$i] - $num_pruned;
            splice @protodomain_refs, $prune_index, 1;
            splice @protodomain_states, $prune_index, 1;
            $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->delete_row($prune_index);
            $adjacency_matrix_ref_of{$obj_ID}{protodomains}->[0]->delete_column($prune_index);
            $num_pruned++;
        }

        # assign ref/state lists to attributes
        $adjacency_matrix_node_refs_of{$obj_ID}{protodomains} = [@protodomain_refs];
        $adjacency_matrix_node_states_of{$obj_ID}{protodomains} = [@protodomain_states];

        # INIT AND BUILD DOMAIN AND GENE ADJACENCY MATRICES
        my %domain_matrix_hash;  # quick store of connected domains
        my %gene_matrix_hash;    # quick store of connected genes
        for (my $i=0; $i < @protodomain_refs; $i++) {
            my $x_ref = $protodomain_refs[$i];
            for (my $j=$i; $j < @protodomain_refs; $j++) {
                my $y_ref = $protodomain_refs[$j];
                my @rules = @{$rules_cache{$x_ref}{$y_ref}};
                if (@rules) {
                    my $x_domain_ref = $protodomain_refs[$i]->get_upper_ref();
                    my $y_domain_ref = $protodomain_refs[$j]->get_upper_ref();
                    my $x_gene_ref = $x_domain_ref->get_upper_ref();
                    my $y_gene_ref = $y_domain_ref->get_upper_ref();
                    $domain_matrix_hash{$x_domain_ref->get_name()}{$y_domain_ref->get_name()} = 1;
                    $domain_matrix_hash{$y_domain_ref->get_name()}{$x_domain_ref->get_name()} = 1;
                    $gene_matrix_hash{$x_gene_ref->get_name()}{$y_gene_ref->get_name()} = 1;
                    $gene_matrix_hash{$y_gene_ref->get_name()}{$x_gene_ref->get_name()} = 1;
                }
            }
        }

        my @domains = $self->get_gene_parser_ref()->get_domain_parser_ref()->get_object_instances();
        $adjacency_matrix_node_refs_of{$obj_ID}{domains} = [@domains];
        $adjacency_matrix_ref_of{$obj_ID}{domains}->[0]->set_element(@domains-1,@domains-1, 0);   # init matrix
        for (my $i=0; $i < @domains; $i++) {
            for (my $j=0; $j < @domains; $j++) {
                my $x_domain_ref = $domains[$i];
                my $y_domain_ref = $domains[$j];
                if ($domain_matrix_hash{$x_domain_ref->get_name()}{$y_domain_ref->get_name()}) {
                    $adjacency_matrix_ref_of{$obj_ID}{domains}->[0]->set_element($i,$j, 1);
                } else {
                    $adjacency_matrix_ref_of{$obj_ID}{domains}->[0]->set_element($i,$j, 0);
                }
            }
        }

        my @genes = $self->get_gene_parser_ref()->get_object_instances();
        $adjacency_matrix_node_refs_of{$obj_ID}{genes} = [@genes];
        $adjacency_matrix_ref_of{$obj_ID}{genes}->[0]->set_element(@genes-1,@genes-1, 0);   # init matrix
        for (my $i=0; $i < @genes; $i++) {
            for (my $j=0; $j < @genes; $j++) {
                my $x_gene_ref = $genes[$i];
                my $y_gene_ref = $genes[$j];
                if ($gene_matrix_hash{$x_gene_ref->get_name()}{$y_gene_ref->get_name()}) {
                    $adjacency_matrix_ref_of{$obj_ID}{genes}->[0]->set_element($i,$j, 1);
                } else {
                    $adjacency_matrix_ref_of{$obj_ID}{genes}->[0]->set_element($i,$j, 0);
                }
            }
        }

        # now build the high-order adjacency matrices and connectivity matrix
        $self->compute_high_order_adjacency_matrix(key => "protodomains", order => 12);
        $self->compute_connectivity_matrix(key => "protodomains");
        $self->compute_high_order_adjacency_matrix(key => "genes", order => 8);
        $self->compute_connectivity_matrix(key => "genes");

    }

    #--------------------------------------------------------------------------------------
    # Function: compute_high_order_adjacency_matrix
    # Synopsys: Computes n-th order adjacency matrix, wherein element i,j gives the
    #           number of n-step connections between nodes i and j.
    #--------------------------------------------------------------------------------------
    sub compute_high_order_adjacency_matrix {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            key => undef,
            order => undef,
            @_,
        );
        check_args(\%args, 2);

        my $key = $args{key};
        my $order = $args{order};

        confess "ERROR: no matrix defined for $key" if !ref $adjacency_matrix_ref_of{$obj_ID}{$key};

        printn "compute_high_order_adjacency_matrix: computing to order $order" if $verbosity >= 1;

        for (my $i = 1; $i < $order; $i++) {  # if order is 1, nothing to do...
            $adjacency_matrix_ref_of{$obj_ID}{$key}->[$i] = Matrix->mmult(
                $adjacency_matrix_ref_of{$obj_ID}{$key}->[$i-1],
                $adjacency_matrix_ref_of{$obj_ID}{$key}->[0]
            );
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: compute_connectivity_matrix
    # Synopsys: Use adjacency matrices to determine if nodes i,j are connected
    #           by a path <= n, where n is the highest order adjacency matrix that was
    #           previously computed. Element (i,j) is set to  1 iff connected.
    #--------------------------------------------------------------------------------------
    sub compute_connectivity_matrix {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            key => undef,
            @_,
        );
        check_args(\%args, 1);
        my $key = $args{key};

        printn "compute_connectivity_matrix: computing" if $verbosity >= 1;

        $connectivity_matrix_ref_of{$obj_ID}{$key} = dclone($adjacency_matrix_ref_of{$obj_ID}{$key}->[0]);
        if (@{$adjacency_matrix_ref_of{$obj_ID}{$key}} > 1) {
            for (my $i = 1; $i < @{$adjacency_matrix_ref_of{$obj_ID}{$key}}; $i++) {
                $connectivity_matrix_ref_of{$obj_ID}{$key} = Matrix->madd(
                    $connectivity_matrix_ref_of{$obj_ID}{$key},
                    $adjacency_matrix_ref_of{$obj_ID}{$key}->[$i],
                );
            }
        }
        $connectivity_matrix_ref_of{$obj_ID}{$key} = Matrix->mnonzero($connectivity_matrix_ref_of{$obj_ID}{$key});
    }

    #--------------------------------------------------------------------------------------
    # Function: get_adjacent (!!! change to get_adjacent_fanout)
    # Synopsys: Using adjacency matrix, return list of nodes immediately connected to argument.
    #--------------------------------------------------------------------------------------
    sub get_adjacent {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            key => undef,
            ref => undef,
            state => "UNDEF",
            @_,
        );
        check_args(\%args, 3);
        my $key = $args{key};
        my $ref = $args{ref};
        my $state = $args{state};

        # find index corresponding to $ref
        my $index;
        for ($index = 0; $index < @{$adjacency_matrix_node_refs_of{$obj_ID}{$key}}; $index++) {
            last if (($adjacency_matrix_node_refs_of{$obj_ID}{$key}->[$index] == $ref) &&
                ($state eq "UNDEF" ||
                    (defined $adjacency_matrix_node_states_of{$obj_ID}{$key}->[$index] &&
                        $adjacency_matrix_node_states_of{$obj_ID}{$key}->[$index] == $state))) ;
        }

        if ($index == @{$adjacency_matrix_node_refs_of{$obj_ID}{$key}}) {
            return ();  # i.e. ref not in adjacency matrix, so not connected at all -- return empty
        }

        my @return_list = ();
        my $matrix_ref = $adjacency_matrix_ref_of{$obj_ID}{$key}->[0]->get_matrix_ref();
        my $matrix_row_ref = $matrix_ref->[$index];

        for (my $i=0; $i < @{$matrix_row_ref}; $i++) {
            my $elem = $matrix_row_ref->[$i];
            if ($elem != 0) {
                my $node_ref = $adjacency_matrix_node_refs_of{$obj_ID}{$key}->[$i];
                my $node_state = $adjacency_matrix_node_states_of{$obj_ID}{$key}->[$i];
                my $sprint = $node_ref->get_name().(defined $node_state ? "/$node_state" : "");
                push @return_list, [$node_ref, $node_state, $sprint];
            }
        }
        return @return_list;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_connected (!!!should really be called get_fanout???)
    # Synopsys: Using connectivity matrix, return list of nodes
    #           that are connected (same subnet) to argument.
    #--------------------------------------------------------------------------------------
    sub get_connected {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            key => undef,
            ref => undef,
            state => "UNDEF",
            @_,
        );
        check_args(\%args, 3);
        my $key = $args{key};
        my $ref = $args{ref};
        my $state = $args{state};

        # find index corresponding to $ref
        my $index;
        for ($index = 0; $index < @{$adjacency_matrix_node_refs_of{$obj_ID}{$key}}; $index++) {
            last if (($adjacency_matrix_node_refs_of{$obj_ID}{$key}->[$index] == $ref) &&
                ($state eq "UNDEF" ||
                    (defined $adjacency_matrix_node_states_of{$obj_ID}{$key}->[$index] &&
                        $adjacency_matrix_node_states_of{$obj_ID}{$key}->[$index] == $state))) ;
        }

        if ($index == @{$adjacency_matrix_node_refs_of{$obj_ID}{$key}}) {
            return ();  # i.e. ref not in adjacency matrix, so not connected at all -- return empty
        }

        my @return_list = ();
        my $matrix_ref = $connectivity_matrix_ref_of{$obj_ID}{$key}->get_matrix_ref();
        my $matrix_row_ref = $matrix_ref->[$index];

        for (my $i=0; $i < @{$matrix_row_ref}; $i++) {
            my $elem = $matrix_row_ref->[$i];
            if ($elem != 0) {
                my $node_ref = $adjacency_matrix_node_refs_of{$obj_ID}{$key}->[$i];
                my $node_state = $adjacency_matrix_node_states_of{$obj_ID}{$key}->[$i];
                my $sprint = $node_ref->get_name().(defined $node_state ? "/$node_state" : "");
                push @return_list, [$node_ref, $node_state, $sprint];
            }
        }
        return @return_list;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_isolated_genes
    # Synopsys: Find genes that don't interact with anything, based on computed
    #           adjacency matrix.
    #--------------------------------------------------------------------------------------
    sub get_isolated_genes {
        my $self = shift; my $obj_ID = ident $self;

        printn "get_isolated_genes: finding non-interacting genes..." if $verbosity >= 2;

        my @genes;
        for (my $i=0; $i < @{$adjacency_matrix_node_refs_of{$obj_ID}{genes}}; $i++) {
            my $gene_ref = $adjacency_matrix_node_refs_of{$obj_ID}{genes}->[$i];
            my $gene_name = $gene_ref->get_name();
            if (!grep($_==1, @{$adjacency_matrix_ref_of{$obj_ID}{genes}->[0]->get_matrix_ref()->[$i]})) {
                printn "get_isolated_genes: found isolated gene $gene_name" if $verbosity >= 3;
                push @genes, $gene_ref;
            } else {
                printn "get_isolated_genes: gene $gene_name is not isolated" if $verbosity >= 3;
            }
        }
        printn "get_isolated_genes: found isolated genes: ". join ",", map($_->get_name(), @genes) if $verbosity >= 2;
        return @genes
    }

    #--------------------------------------------------------------------------------------
    # Function: prune_isolated_genes
    # Synopsys: Find isolated genes and reset their export_flag.
    #--------------------------------------------------------------------------------------
    sub prune_isolated_genes {
        my $self = shift; my $obj_ID = ident $self;

        my @isolated_genes = $self->get_isolated_genes();
        map {$_->set_export_flag(0)} @isolated_genes;
        printn "prune_isolated_genes: pruned genes: ". join ",", map($_->get_name(), @isolated_genes) if $verbosity >= 1;
        return @isolated_genes;
    }

    #--------------------------------------------------------------------------------------
    # Function: find_adjacent_csites
    # Synopsys: Returns a list of csites adjacent to given protodomain.
    #--------------------------------------------------------------------------------------
    sub find_adjacent_csites {
        my $self = shift; my $obj_ID = ident $self;
        my $protodomain_ref = shift;
        my $polarity = shift;	# 0 for kinase, 1 for phosphatase

        confess "ERROR: protodomain_ref not defined" if (!defined $protodomain_ref);

        printn "find_adjacent_enzymes: finding csites adjacent to protodomain ".$protodomain_ref->get_name() if $verbosity >= 2;

        my @csite_list = ();
        foreach my $get_adjacent_ref ($self->get_adjacent(key => "protodomains", ref => $protodomain_ref, state => $polarity)) {
            my $ligand_protodomain_ref = $get_adjacent_ref->[0];
            next if ($ligand_protodomain_ref->get_translation_ref()->{type} ne "csite");
            next if ($ligand_protodomain_ref->get_translation_ref()->{substrate_polarity} != $polarity);
            printn "find_adjacent_enzymes: ".$ligand_protodomain_ref->get_name()." modifies (polarity=$polarity) ".$protodomain_ref->get_name();
            push @csite_list, $ligand_protodomain_ref;
        }
        return @csite_list;
    }

}


sub run_testcases {
    printn "NO TESTCASES!!!";
}


# Package BEGIN must return true value
return 1;


# #######################################################################################
# # Function: find_input_ligand
# # Synopsys: Returns list of simple proteins that can be used as inputs to system.
# #           The compile_ligand_lists() routine must be run first.
# #######################################################################################
# sub find_input_ligand {
#     printn "find_input_ligand: searching for input ligand candidates" if ($verbosity >= 1);

#     my @return_list = ();
#     my ($protein, $domain, $protodomain);
#     foreach $protein (sort keys %{$sdb->{protein_table}}) {
# 	if (@{$sdb->{protein_table}{$protein}{domain_list}} != 1) {next;}      # one domain
# 	if (@{$sdb->{protein_table}{$protein}{ligand_list}[0]} != 1) {next;}   # one protodomain
# 	if (@{$sdb->{protein_table}{$protein}{ligand_list}[0][0]} == 0) {next;}   # at least one interaction
# 	$domain = $sdb->{protein_table}{$protein}{domain_list}[0];
# 	$protodomain = $sdb->{domain_table}{$domain}{protodomain_list}[0];
# 	if ($sdb->{protodomain_table}{$protodomain}{protodomain_type} eq "msite") {next;}   # not an msite
# 	if ($sdb->{protodomain_table}{$protodomain}{protodomain_type} eq "csite") {next;}   # not a  csite
# 	if ($sdb->{protodomain_table}{$protodomain}{bsite_default_state} ne "on") {next;}   # must be bindable :)
# 	if (protodomain_self_binding($protodomain, $bsite_profile_width)) {next;} # must bind to something else than itself
# 	push @return_list, $protein;
#     };
#     if (@return_list >= 1) {
# 	printn "find_input_ligand: found input ligands (@return_list)" if ($verbosity >= 1);
#     } else {
# 	printn "find_input_ligand: found no input ligands" if ($verbosity >= 1);
#     }
#     return @return_list;
# }

# #######################################################################################
# # Function: find_output_tfs
# # Synopsys: Returns list of simple proteins that can be used as outputs to system
# #######################################################################################
# # Detailed Description:
# # ---------------------
# # We are trying to find candidate transcription factors whose binding site
# # will turn on once phosphorylated.  Look for a single protodomain protein with a
# # modification site and whose binding site is on by default.
# #######################################################################################
# sub find_output_tfs {
#     my @search_list;
#     if (defined $_[0]) {
# 	printn "find_output_tfs: searching for output tf candidates in given protein list" if ($verbosity >= 1);
# 	printn "find_output_tfs: searching in " . join ",", @_ if ($verbosity >= 1);
# 	@search_list = @_;
#     } else {
# 	printn "find_output_tfs: searching for output tf candidates in all proteins" if ($verbosity >= 1);
# 	@search_list = (sort keys %{$sdb->{protein_table}});
#     }

#     my @return_list = ();
#     my ($protein, $domain, $protodomain, $ligand_protein, $ligand_domain, $ligand_protodomain, $kinase_count);
#   PROTEIN: foreach $protein (@search_list) {
# 	printn "find_output_tfs: checking protein $protein";
# 	if (get_protodomains_in_protein($protein) != 1) {next;}         # one protodomain
# 	if (network_get_adjacent_protein($protein) == 0) {next;}        # at least one interaction
# 	$domain = (get_domains_in_protein($protein))[0];
# 	$protodomain = (get_protodomains_in_protein($protein))[0];
# 	if ($sdb->{protodomain_table}{$protodomain}{protodomain_type} ne "msite") {next;}   # phosphorylatable
# 	if ($sdb->{protodomain_table}{$protodomain}{bsite_default_state} ne "on") {next;}       # must be bindable :)
# 	if (protodomain_self_binding($protodomain, $bsite_profile_width)) {next;}             # must bind to something else than itself
# 	$kinase_count = 0;
#       LIGAND: foreach $ligand_protodomain (network_get_adjacent_protodomain($protodomain)) {
# 	    $ligand_domain = (in_domain($ligand_protodomain))[0];
# 	    $ligand_protein = (in_protein($ligand_domain))[0];
# 	    if ($sdb->{protodomain_table}{$ligand_protodomain}{protodomain_type} ne "csite") {next;}
# #	    if ($sdb->{protodomain_table}{$ligand_protodomain}{bsite_default_state} ne "on") {next;}
# #	    if ($sdb->{protodomain_table}{$ligand_protodomain}{csite_default_state} ne "on") {next;}
# 	    printn "find_output_tfs: kinase that phosphorylates $protein:$domain:$protodomain is $ligand_protein:$ligand_domain:$ligand_protodomain";
# 	    $kinase_count++;
# 	}
# 	if ($kinase_count == 0) {next;}  # must have a connected kinase
# 	push @return_list, $protein;
#     };
#     if (@return_list >= 1) {
# 	printn "find_output_tfs: found outputs (@return_list)" if ($verbosity >= 1);
#     } else {
# 	printn "find_output_tfs: found no outputs" if ($verbosity >= 1);
#     }
#     return @return_list;
# }

