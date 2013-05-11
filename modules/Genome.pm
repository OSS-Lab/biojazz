#-#####################################################################################
#- File:     Genome.pm
#- Synopsys: 
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Genome;
use Class::Std::Storable;
use base qw(Parser Network);
{
    use Carp;

## For manipulation of bit vectors
#use Bit::Vector;

    use Utils;
    use Globals qw($verbosity $VERSION $RELEASE_DATE);

#use BitString;

    use GenomeInstance;

    use Gene;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    Genome->set_class_data("ELEMENT_CLASS", "Gene");
    Genome->set_class_data("INSTANCE_CLASS", "GenomeInstance");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    # interaction rule parameters
    my %radius_of :ATTR(get => 'radius', set => 'radius', init_arg => 'radius');
    my %kf_max_of :ATTR(get => 'kf_max', set => 'kf_max', init_arg => 'kf_max');
    my %kf_min_of :ATTR(get => 'kf_min', set => 'kf_min', init_arg => 'kf_min');
    my %kb_max_of :ATTR(get => 'kb_max', set => 'kb_max', init_arg => 'kb_max');
    my %kb_min_of :ATTR(get => 'kb_min', set => 'kb_min', init_arg => 'kb_min');
    my %kp_max_of :ATTR(get => 'kp_max', set => 'kp_max', init_arg => 'kp_max');
    my %kp_min_of :ATTR(get => 'kp_min', set => 'kp_min', init_arg => 'kp_min');

    # translated rules
    my %rules_of :ATTR(get => 'rules', set => 'rules');

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

        # DEFAULTS

        # INIT
        $rules_of{$obj_ID} = [];

        # STRUCTURE
        my $gene_parser_ref = Gene->new({
                %{$arg_ref->{Gene} || {}},
                name => "Gene",
            });
        $self->add_element($gene_parser_ref);

        my $gene_start_code = $gene_parser_ref->get_gene_start_code();
        my $gene_min_length = $gene_parser_ref->get_min_length();

        # set structure/mutatation rates
        $self->set_structure_ref([
                ["PRE_JUNK", "^.*?(?=$gene_start_code)"],
                ["genes",
                    $gene_parser_ref,
                    "\\G.*?(?=$gene_start_code.{".($gene_min_length-length($gene_start_code))."})"],
                ["POST_JUNK", "\\G.*"],
            ]);

        $self->set_linker_ref({
                genes => "0000", # necessary for untranscribe routine
            });

        $self->set_mutation_rate_ref({  # don't mutate JUNK (to preserve number of genes)
                "PRE_JUNK" => 0.0,
                "genes" => 1.0,
                "genes_linker" => 0.0,
                "POST_JUNK" => 0.0,
            });
    }

    #--------------------------------------------------------------------------------------
    # Function: get_gene_parser_ref
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_gene_parser_ref {
        my $self = shift;
        return $self->get_element(0);
    }

    #--------------------------------------------------------------------------------------
    # Function: untranslate_field
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub untranslate_field {
        my $self = shift; my $obj_ID = ident $self;
        my $field_name = shift;
        my $field_value = shift;

        confess "ERROR: unexpected reference" if ref $field_name;

        my $field_sequence = "";

        SWITCH: {
            if ($field_name =~ /JUNK/) {
                # no untranslation necessary
                $field_sequence .= defined $field_value ? $field_value : "";
                last SWITCH;
            }
            confess "ERROR: translate_field() -- unknown field $field_name";
        }			# SWITCH
        return $field_sequence;
    }

    #--------------------------------------------------------------------------------------
    # Function: export_anc
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub export_anc {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            max_external_iterations => -1,
            max_internal_iterations => -1,
            max_complex_size => -1,
            max_species => 50,
            max_csite_bound_to_msite_number => 1,
            default_max_count => -1,
            default_steric_factor => 1e-3,
            export_graphviz => "nothing",
            equations => [],
            compartment_volume => "UNDEF",
            ode_event_times => "UNDEF",
            SS_timescale => "UNDEF",
            t_final => "UNDEF",
            t_vector => "UNDEF",
            matlab_ode_solver => "UNDEF",
            matlab_odeset_options => "UNDEF",
            @_,
        );
        check_args(\%args, 16);

        my $date = `date`; chomp($date);
        my $str = <<END;
##########################################################################
# Allosteric Network Compiler (ANC) model file
# Created by BioJazz version $VERSION released on $RELEASE_DATE
# $date
##########################################################################

##########################################################################
# PARAMETERS
##########################################################################
END

        # APPEND PARAMETERS
        $str .= "\$max_external_iterations = $args{max_external_iterations}\n" if defined $args{max_external_iterations};
        $str .= "\$max_internal_iterations = $args{max_internal_iterations}\n" if defined $args{max_internal_iterations};
        $str .= "\$max_complex_size = $args{max_complex_size}\n" if defined $args{max_complex_size};
        $str .= "\$max_species = $args{max_species}\n" if defined $args{max_species};
        $str .= "\$max_csite_bound_to_msite_number = $args{max_csite_bound_to_msite_number}\n" if defined $args{max_csite_bound_to_msite_number};
        $str .= "\$default_steric_factor = $args{default_steric_factor}\n" if defined $args{default_steric_factor};
        $str .= "\$export_graphviz = \"$args{export_graphviz}\"\n" if defined $args{export_graphviz};

        $str .= <<END;
        \n
##########################################################################
# OBJECTS:
##########################################################################
END

        # APPEND PROTEINS, DOMAINS, PROTODOMAINS
        my @genomes = $self->get_object_instances();
        foreach my $genome_instance_ref (@genomes) {
            $str .= $genome_instance_ref->export_anc($args{default_max_count});
        }

        # APPEND RULES
        my @genes = map {$_->[0]} (map {@{$_->get_transcription_ref->{genes}}} @genomes);
        @genes = grep {$_->get_export_flag()} @genes;
        my @domains = map {$_->[0]} (map(@{$_->get_transcription_ref->{domains}}, @genes));
        my @protodomains = map {$_->[0]} (map(@{$_->get_transcription_ref->{protodomains}}, @domains));

        $str .= <<END;
        \n
##########################################################################
# RULES:
##########################################################################
MODEL:

END
        $rules_of{$obj_ID} = [];
#	my @protodomains = $self->get_gene_parser_ref()->get_domain_parser_ref()->get_protodomain_parser_ref()->get_object_instances();
        for (my $i=0; $i < @protodomains; $i++) {
            for (my $j=$i; $j < @protodomains; $j++) {
                my @rules = ProtoDomainInstance->create_canbindrules(
                    x_ref => $protodomains[$i],
                    y_ref => $protodomains[$j],
                    radius => $self->get_radius(),
                    kf_max => $self->get_kf_max(),
                    kf_min => $self->get_kf_min(),
                    kb_max => $self->get_kb_max(),
                    kb_min => $self->get_kb_min(),
                    kp_max => $self->get_kp_max(),
                    kp_min => $self->get_kp_min(),
                );
                push @{$rules_of{$obj_ID}}, @rules;
            }
        }
        printn "export_anc: created ".scalar(@{$rules_of{$obj_ID}})." rules" if $verbosity >= 1;
        $str .= join "\n", @{$rules_of{$obj_ID}};

        $str .= <<END;
        \n
##########################################################################
# EQUATIONS:
##########################################################################
EQN:

END
        $str .= join "\n", @{$args{equations}};


        $str .= <<END;
        \n
##########################################################################
# CONFIG:
##########################################################################
CONFIG:

END
        $str .= "compartment_volume = $args{compartment_volume};\n" if $args{compartment_volume} ne "UNDEF";
        $str .= "ode_event_times = $args{ode_event_times};\n" if $args{ode_event_times} ne "UNDEF";
        $str .= "SS_timescale = $args{SS_timescale};\n" if $args{SS_timescale} ne "UNDEF";
        $str .= "t_final = $args{t_final};\n" if $args{t_final} ne "UNDEF";
        $str .= "t_vector = $args{t_vector};\n" if $args{t_vector} ne "UNDEF";
        $str .= "matlab_ode_solver = $args{matlab_ode_solver};\n" if $args{matlab_ode_solver} ne "UNDEF";
        $str .= "matlab_odeset_options = $args{matlab_odeset_options};\n" if $args{matlab_odeset_options} ne "UNDEF";

        return $str;
    }
}

sub run_testcases {

    $verbosity = 2;

    printn "CREATION";
    srand(98275);
    my $seq_ref = Sequence->new({});
    $seq_ref->generate_random_sequence(2000);

    my $ref = Genome->new({
            name => "Genome",
            radius => 1,
            kf_max => 1e5,
            kf_min => 0.1,
            kb_max => 10,
            kb_min => 0.001,
            kp_max => 1000,
            kp_min => 1,
        });

    printn $ref->_DUMP();
    my $iref = $ref->parse(sequence_ref => $seq_ref);
    $iref->translate();
    printn $iref->_DUMP();

    printn "ANC EXPORT";
    my $anc_model = $ref->export_anc(
        max_species => 200,
    );
    printn $anc_model;

    burp_file("test/modules/Genome.mod", $anc_model);

    printn "TEST ADJACENCY/CONNECTIVITY MATRIX";
    $verbosity = 1;
    $ref->build_network();
    printn join(",", map {$_->get_name()} @{$ref->get_adjacency_matrix_node_refs()->{protodomains}});
    printn $ref->get_adjacency_matrix_ref()->{protodomains}->[0]->sprint_matrix();
    printn join(",", map {$_->get_name()} @{$ref->get_adjacency_matrix_node_refs()->{domains}});
    printn $ref->get_adjacency_matrix_ref()->{domains}->[0]->sprint_matrix();
    printn join(",", map {$_->get_name()} @{$ref->get_adjacency_matrix_node_refs()->{genes}});
    printn $ref->get_adjacency_matrix_ref()->{genes}->[0]->sprint_matrix();
    $ref->compute_high_order_adjacency_matrix(key => "genes", order => 4);
    printn $ref->get_adjacency_matrix_ref()->{genes}->[1]->sprint_matrix();
    printn $ref->get_adjacency_matrix_ref()->{genes}->[2]->sprint_matrix();
    printn $ref->get_adjacency_matrix_ref()->{genes}->[3]->sprint_matrix();
    $ref->compute_connectivity_matrix(key => "genes");
    printn $ref->get_connectivity_matrix_ref()->{genes}->sprint_matrix();

    printn "TEST GET_ADJACENT_*()";
    my ($pdref) = ($ref->get_gene_parser_ref->get_domain_parser_ref->get_protodomain_parser_ref->get_object_instances())[2];
    printn $pdref->get_name();
    printn join ",", (map {$_->[2]} $ref->get_adjacent(key => "protodomains", ref => $pdref));

    my ($gref) = ($ref->get_gene_parser_ref->get_object_instances())[2];
    printn $gref->get_name();
    printn join ",", (map {$_->[2]} $ref->get_adjacent(key => "genes", ref => $gref));
    printn join ",", (map {$_->[2]} $ref->get_connected(key => "genes", ref => $gref));

    printn "TEST ISOLATED GENE FIND/EXPORT";
    $ref->prune_isolated_genes();
    my $pruned_anc_model = $ref->export_anc(
        max_species => 200,
    );
    printn $pruned_anc_model;
    burp_file("test/modules/Genome.pruned.mod", $pruned_anc_model);

    system("tkdiff test/modules/Genome.mod test/modules/Genome.pruned.mod");

    printn "STORABLE TEST";
    $ref->set_rules(undef);  # the rules don't get DUMPed properly
    use Storable;
    my $ice_ref = Storable::freeze($ref);
    my $water_ref = Storable::thaw($ice_ref);
    printn $ref->_DUMP();
    printn $water_ref->_DUMP();
}

# Package BEGIN must return true value
return 1;

