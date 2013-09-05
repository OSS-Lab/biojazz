#-#####################################################################################
#- File:     Gene.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Gene;
use Class::Std::Storable;
use base qw(Parser);
{
    use Carp;
    use Utils;

    use BitString;
    use GeneInstance;
    use Domain;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    Gene->set_class_data("ELEMENT_CLASS", "Domain");
    Gene->set_class_data("INSTANCE_CLASS", "GeneInstance");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    # parsing parameters
    my %gene_start_code_of :ATTR(get => 'gene_start_code', set => 'gene_start_code');
    my %soft_linker_code_of :ATTR(get => 'soft_linker_code', set => 'soft_linker_code');
    my %regulated_concentration_width_of :ATTR(get => 'regulated_concentration_width', set => 'regulated_concentration_width');
    my %unused_width_of :ATTR(get => 'unused_width', set => 'unused_width');

    # parsing parameters (derived)
    my %gene_start_code_length_of :ATTR(get => 'gene_start_code_length', set => 'gene_start_code_length');
    my %stop_linker_code_of :ATTR(get => 'stop_linker_code', set => 'stop_linker_code');  # regexp for !hard && !soft
    my %STOP_linker_code_of :ATTR(get => 'STOP_linker_code', set => 'STOP_linker_code');  # all 1s pattern

    # scaling parameters
    my %regulated_concentration_max_of :ATTR(get => 'regulated_concentration_max', set => 'regulated_concentration_max');
    my %regulated_concentration_min_of :ATTR(get => 'regulated_concentration_min', set => 'regulated_concentration_min');

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
    # Function: 
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

        # DEFAULTS
        $gene_start_code_of{$obj_ID} = "01111110";
        $soft_linker_code_of{$obj_ID} = "001";
        $regulated_concentration_width_of{$obj_ID} = 10;
        $unused_width_of{$obj_ID} = 4;

        $regulated_concentration_max_of{$obj_ID} = 1e-3;
        $regulated_concentration_min_of{$obj_ID} = 1e-3;

        # INIT
        $gene_start_code_of{$obj_ID} = $arg_ref->{gene_start_code} if exists $arg_ref->{gene_start_code};
        $soft_linker_code_of{$obj_ID} = $arg_ref->{soft_linker_code} if exists $arg_ref->{soft_linker_code};
        $regulated_concentration_width_of{$obj_ID} = $arg_ref->{regulated_concentration_width} if exists $arg_ref->{regulated_concentration_width};
        $unused_width_of{$obj_ID} = $arg_ref->{unused_width} if exists $arg_ref->{unused_width};

        $regulated_concentration_max_of{$obj_ID} = $arg_ref->{regulated_concentration_max} if exists $arg_ref->{regulated_concentration_max};
        $regulated_concentration_min_of{$obj_ID} = $arg_ref->{regulated_concentration_min} if exists $arg_ref->{regulated_concentration_min};

        # STRUCTURE
        my $domain_parser_ref = Domain->new({
                %{$arg_ref->{Domain} || {}},
                name => "Domain",
            });
        $self->add_element($domain_parser_ref);

        my $domain_min_length = $domain_parser_ref->get_min_length();

        # COMPUTE DERIVED ATTRIBUTES
        my $gene_start_code_length = $gene_start_code_length_of{$obj_ID} = length($gene_start_code_of{$obj_ID});
        my $hard_linker_code = $domain_parser_ref->get_hard_linker_code();

        my $soft_linker_code = $soft_linker_code_of{$obj_ID};
        my $linker_code_length = length($hard_linker_code);
        if ($linker_code_length != length($soft_linker_code)) {
            confess "ERROR: START -- linker codes must be same length";
        };
        # stop_linker_code regexp is whatever is on tail end of gene
        # (n.b. the pattern may be a hard/soft code if there were not
        # enough bits left for a protodomain/domain)
        my $stop_linker_code = $stop_linker_code_of{$obj_ID} = ".{0,$linker_code_length}";
        # this "hard" STOP code is a CONSTANT consisting of the all 1's pattern
        my $STOP_linker_code = $STOP_linker_code_of{$obj_ID} = "1" x $linker_code_length;

        $self->set_structure_ref([
                ["START_CODE", "$gene_start_code_of{$obj_ID}"],
                ["regulated_concentration", "\\G.{$regulated_concentration_width_of{$obj_ID}}"],
                ["UNUSED", "\\G.{$unused_width_of{$obj_ID}}"],
                ["domains", $domain_parser_ref, "\\G$soft_linker_code_of{$obj_ID}(?=.{$domain_min_length})"],
                ["STOP_CODE", "\\G$stop_linker_code"],
            ]);

        $self->set_linker_ref({
                domains => $soft_linker_code_of{$obj_ID}, # necessary for untranscribe routine
            });

        $self->set_mutation_rate_ref({   # don't mutate START/STOP (to preserve number of genes)
                "START_CODE" => 0.0,
                "regulated_concentration" => 1.0,
                "domains" => 1.0,
                "domains_linker" => 1.0,
                "STOP_CODE" => 0.0,
            });
    }

    #--------------------------------------------------------------------------------------
    # Function: get_domain_parser_ref
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_domain_parser_ref {
        my $self = shift;
        return $self->get_element(0);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_min_length
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_min_length {
        my $self = shift; my $obj_ID = ident $self;

        my $min_length = 0;
        $min_length += length($gene_start_code_of{$obj_ID});
        $min_length += $regulated_concentration_width_of{$obj_ID};
        $min_length += $unused_width_of{$obj_ID};
        $min_length += $self->get_domain_parser_ref()->get_min_length();

        return $min_length;
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
            if ($field_name eq "UNUSED") {
                # no untranslation necessary, just 0-padding or truncation
                my $lc_field_name = lc $field_name;  # for UNUSED field
                my $field_width = eval("\$self->get_${lc_field_name}_width()");
                confess "ERROR: internal error -- field_width not defined (field $field_name)" if !defined $field_width;
                $field_sequence .= "0"x($field_width-length($field_value)).$field_value;  # 0-pad
                $field_sequence = substr $field_sequence, -$field_width; # truncate
                last SWITCH;
            }
            if ($field_name eq "START_CODE") {
                $field_sequence .= $gene_start_code_of{$obj_ID};
                last SWITCH;
            }
            if ($field_name eq "STOP_CODE") {
                $field_sequence .= $STOP_linker_code_of{$obj_ID};
                last SWITCH;
            }
            if ($field_name eq "regulated_concentration") {
                my $field_width = eval("\$self->get_${field_name}_width()");
                confess "ERROR: internal error -- field_width not defined" if !defined $field_width;
                my $max = eval("\$self->get_${field_name}_max()");
                my $min = eval("\$self->get_${field_name}_min()");
                confess "ERROR: internal error -- max not defined" if !defined $max;
                confess "ERROR: internal error -- min not defined" if !defined $min;
                confess "ERROR: given $field_name=$field_value gt max=$max" if $field_value > $max;
                confess "ERROR: given $field_name=$field_value lt min=$min" if $field_value < $min;
                my $field_int_value = round2int(loglinear_inv($min, $max, (2**$field_width) - 1, $field_value));
                my $field_bin_value = dec2bin($field_int_value);
                $field_sequence .= "0"x($field_width-length($field_bin_value)).$field_bin_value;  # 0-pad
                $field_sequence = substr $field_sequence, -$field_width; # truncate
                last SWITCH;
            }
            confess "ERROR: translate_field() -- unknown field $field_name";
        }			# SWITCH
        return $field_sequence;
    }

}


sub run_testcases {
    use Data::Dumper;

    srand(365345);
    my $seq_ref = Sequence->new({});
    $seq_ref->generate_random_sequence(1000);

    my $ref = Gene->new({
            name => "Gene",
        });

    printn $ref->_DUMP();

    printn "Gene min length is: ".$ref->get_min_length();

    my $iref = $ref->parse(sequence_ref => $seq_ref);
    $iref->translate();
    printn $iref->_DUMP();
    printn $iref->export_anc();

    printn "STORABLE TEST";
    use Storable;
    my $ice_ref = Storable::freeze($ref);
    my $water_ref = Storable::thaw($ice_ref);
    printn $ref->_DUMP();
    printn $water_ref->_DUMP();

    printn "UNTRANSLATE TEST";
    my $untranslated_ref = $ref->untranslate({
            START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
            regulated_concentration => 7e-4,
            UNUSED => "1010",
            domains => [
                {
                    allosteric_flag => 1,
                    RT_transition_rate => 1e1,
                    TR_transition_rate => 1e-1,
                    RT_phi => 0.45,
                    protodomains => [
                        {
                            type => "csite",
                            substrate_polarity => 1,
                            binding_profile => "010",
                            kf_profile => "0011",
                            kb_profile => "111011",
                            kp_profile => "100011",
                            Keq_ratio => 1.0e-3,
                            kf_polarity_mask => "0000",
                            kb_polarity_mask => "110011",
                            kf_conformation_mask => "1111",
                            kb_conformation_mask => "001100",
                            kp_conformation_mask => "010101",
                            UNUSED => "1000000001",
                        },
                        {
                            type => "msite",
                            substrate_polarity => 1,
                            binding_profile => "110",
                            kf_profile => "0011",
                            kb_profile => "111011",
                            kp_profile => "100011",
                            Keq_ratio => 1.0e-2,
                            kf_polarity_mask => "0000",
                            kb_polarity_mask => "110011",
                            kf_conformation_mask => "1111",
                            kb_conformation_mask => "001100",
                            kp_conformation_mask => "010101",
                            UNUSED => "1000000011",
                        },
                        {
                            type => "bsite",
                            substrate_polarity => 0,
                            binding_profile => "111",
                            kf_profile => "0011",
                            kb_profile => "111011",
                            kp_profile => "100011",
                            Keq_ratio => 1.0e-1,
                            kf_polarity_mask => "0000",
                            kb_polarity_mask => "110011",
                            kf_conformation_mask => "1111",
                            kb_conformation_mask => "001100",
                            kp_conformation_mask => "010101",
                            UNUSED => "1000000111",
                        },
                    ],
                    UNUSED => "100",
                },
                {
                    allosteric_flag => 0,
                    RT_transition_rate => 1e2,
                    TR_transition_rate => 1e-2,
                    RT_phi => 0.55,
                    protodomains => [
                        {
                            type => "msite",
                            substrate_polarity => 1,
                            binding_profile => "1111",
                            kf_profile => "0011",
                            kb_profile => "111011",
                            kp_profile => "100011",
                            Keq_ratio => 1.0e-3,
                            kf_polarity_mask => "0000",
                            kb_polarity_mask => "110011",
                            kf_conformation_mask => "1111",
                            kb_conformation_mask => "001100",
                            kp_conformation_mask => "010101",
                            UNUSED => "1100000001",
                        },
                    ],
                    UNUSED => "101",
                },
            ],
        });
    printn Dumper($untranslated_ref);
    my $untranscribed_ref = $ref->untranscribe($untranslated_ref);
    printn $untranscribed_ref->get_sequence();
    my $irref = $ref->parse(sequence_ref => $untranscribed_ref);
    $irref->translate();
    printn $irref->_DUMP();
}


# Package BEGIN must return true value
return 1;

