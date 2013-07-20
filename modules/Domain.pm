#-#####################################################################################
#- File:     Domain.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Domain;
use Class::Std::Storable;
use base qw(Parser);
{
    use Carp;
    use Utils;

    use BitString;
    use DomainInstance;
    use ProtoDomain;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    Domain->set_class_data("ELEMENT_CLASS", "ProtoDomain");
    Domain->set_class_data("INSTANCE_CLASS", "DomainInstance");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    # parsing parameters
    my %hard_linker_code_of :ATTR(get => 'hard_linker_code', set => 'hard_linker_code');
    my %allosteric_flag_width_of :ATTR(get => 'allosteric_flag_width', set => 'allosteric_flag_width');
    my %RT_transition_rate_width_of :ATTR(get => 'RT_transition_rate_width', set => 'RT_transition_rate_width');
    my %TR_transition_rate_width_of :ATTR(get => 'TR_transition_rate_width', set => 'TR_transition_rate_width');
    my %RT_phi_width_of :ATTR(get => 'RT_phi_width', set => 'RT_phi_width');
    my %unused_width_of :ATTR(get => 'unused_width', set => 'unused_width');

    # translation scaling factors
    my %RT_transition_rate_max_of :ATTR(get => 'RT_transition_rate_max', set => 'RT_transition_rate_max');
    my %RT_transition_rate_min_of :ATTR(get => 'RT_transition_rate_min', set => 'RT_transition_rate_min');
    my %TR_transition_rate_max_of :ATTR(get => 'TR_transition_rate_max', set => 'TR_transition_rate_max');
    my %TR_transition_rate_min_of :ATTR(get => 'TR_transition_rate_min', set => 'TR_transition_rate_min');
    my %RT_phi_max_of :ATTR(get => 'RT_phi_max', set => 'RT_phi_max');
    my %RT_phi_min_of :ATTR(get => 'RT_phi_min', set => 'RT_phi_min');

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
    # Function: BUILD
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

        # DEFAULTS
        $hard_linker_code_of{$obj_ID} = "000";
        $allosteric_flag_width_of{$obj_ID} = 1;
        $RT_transition_rate_width_of{$obj_ID} = 10;
        $TR_transition_rate_width_of{$obj_ID} = 10;
        $RT_phi_width_of{$obj_ID} = 10;
        $unused_width_of{$obj_ID} = 4;

        $RT_transition_rate_max_of{$obj_ID} = 1e3;
        $RT_transition_rate_min_of{$obj_ID} = 1e-3;
        $TR_transition_rate_max_of{$obj_ID} = 1e3;
        $TR_transition_rate_min_of{$obj_ID} = 1e-3;
        $RT_phi_max_of{$obj_ID} = 1.0;
        $RT_phi_min_of{$obj_ID} = 0.0;

        # INIT
        $hard_linker_code_of{$obj_ID} = $arg_ref->{hard_linker_code} if exists $arg_ref->{hard_linker_code};
        $allosteric_flag_width_of{$obj_ID} = $arg_ref->{allosteric_flag_width} if exists $arg_ref->{allosteric_flag_width};
        $RT_transition_rate_width_of{$obj_ID} = $arg_ref->{RT_transition_rate_width} if exists $arg_ref->{RT_transition_rate_width};
        $TR_transition_rate_width_of{$obj_ID} = $arg_ref->{TR_transition_rate_width} if exists $arg_ref->{TR_transition_rate_width};
        $RT_phi_width_of{$obj_ID} = $arg_ref->{RT_phi_width} if exists $arg_ref->{RT_phi_width};
        $unused_width_of{$obj_ID} = $arg_ref->{unused_width} if exists $arg_ref->{unused_width};

        $RT_transition_rate_max_of{$obj_ID} = $arg_ref->{RT_transition_rate_max} if exists $arg_ref->{RT_transition_rate_max};
        $RT_transition_rate_min_of{$obj_ID} = $arg_ref->{RT_transition_rate_min} if exists $arg_ref->{RT_transition_rate_min};
        $TR_transition_rate_max_of{$obj_ID} = $arg_ref->{TR_transition_rate_max} if exists $arg_ref->{TR_transition_rate_max};
        $TR_transition_rate_min_of{$obj_ID} = $arg_ref->{TR_transition_rate_min} if exists $arg_ref->{TR_transition_rate_min};
        $RT_phi_max_of{$obj_ID} = $arg_ref->{RT_phi_max} if exists $arg_ref->{RT_phi_max};
        $RT_phi_min_of{$obj_ID} = $arg_ref->{RT_phi_min} if exists $arg_ref->{RT_phi_min};

        # DEFAULT/INIT

        # STRUCTURE
        my $protodomain_parser_ref = ProtoDomain->new({
                %{$arg_ref->{ProtoDomain} || {}},
                name => "ProtoDomain",
            });
        $self->add_element($protodomain_parser_ref);

        my $protodomain_sequence_length = $protodomain_parser_ref->get_length();

        $self->set_structure_ref([
                ["allosteric_flag",              "\\G.{1}"],
                ["RT_transition_rate",           "\\G.{$RT_transition_rate_width_of{$obj_ID}}"],
                ["TR_transition_rate",           "\\G.{$TR_transition_rate_width_of{$obj_ID}}"],
                ["RT_phi",                       "\\G.{$RT_phi_width_of{$obj_ID}}"],
                ["UNUSED",                       "\\G.{$unused_width_of{$obj_ID}}"],
                ["protodomains", $protodomain_parser_ref, "\\G$hard_linker_code_of{$obj_ID}(?=.{$protodomain_sequence_length})"],
            ]);

        $self->set_linker_ref({
                protodomains => $hard_linker_code_of{$obj_ID}, # necessary for untranscribe routine
            });
    }

    #--------------------------------------------------------------------------------------
    # Function: get_protodomain_parser_ref
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_protodomain_parser_ref {
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
        $min_length += 1;  # allosteric_flag
        $min_length += $RT_transition_rate_width_of{$obj_ID};
        $min_length += $TR_transition_rate_width_of{$obj_ID};
        $min_length += $RT_phi_width_of{$obj_ID};
        $min_length += $unused_width_of{$obj_ID};
        $min_length += $self->get_protodomain_parser_ref()->get_length();

        return $min_length;
    }

    #--------------------------------------------------------------------------------------
    # Function: untranslate_field
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub untranslate_field {
        my $self = shift;
        my $field_name = shift;
        my $field_value = shift;

        confess "ERROR: unexpected reference" if ref $field_name;

        my $field_sequence = "";

        SWITCH: {
            if (($field_name eq "allosteric_flag") ||
                ($field_name eq "UNUSED")
            ) {
                # no untranslation necessary, just 0-padding or truncation
                my $lc_field_name = lc $field_name;  # for UNUSED field
                my $field_width = eval("\$self->get_${lc_field_name}_width()");
                confess "ERROR: internal error -- field_width not defined (field $field_name)" if !defined $field_width;
                $field_sequence .= "0"x($field_width-length($field_value)).$field_value;  # 0-pad
                $field_sequence = substr $field_sequence, -$field_width; # truncate
                last SWITCH;
            }
            if (($field_name eq "RT_transition_rate") ||
                ($field_name eq "TR_transition_rate")
            ) {
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
            if ($field_name eq "RT_phi") {
                my $field_width = eval("\$self->get_${field_name}_width()");
                confess "ERROR: internal error -- field_width not defined" if !defined $field_width;
                my $max = eval("\$self->get_${field_name}_max()");
                my $min = eval("\$self->get_${field_name}_min()");
                confess "ERROR: internal error -- max not defined" if !defined $max;
                confess "ERROR: internal error -- min not defined" if !defined $min;
                confess "ERROR: given $field_name=$field_value gt max=$max" if $field_value > $max;
                confess "ERROR: given $field_name=$field_value lt min=$min" if $field_value < $min;
                my $field_int_value = round2int(linear_inv($min, $max, (2**$field_width) - 1, $field_value));
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

    srand(36141471);  # this seed makes it such that we get 2 protodomains

    my $seq_ref = Sequence->new({});
    $seq_ref->generate_random_sequence(1000);

    my $ref = Domain->new({
            name => "Domain",
        });

    printn $ref->_DUMP();
    printn "Domain min length is: ".$ref->get_min_length();

    my $iref = $ref->parse(sequence_ref => $seq_ref);
#    printn $iref->sprint(hierarchical_flag => 0, context_flag => 1);
#    printn $seq_ref->sprint();
    $iref->translate();
    printn $iref->_DUMP();
    map {printn $_->_DUMP()} @{$iref->get_protodomain_parser_ref->get_instances()};
    printn $iref->export_anc();

    printn "STORABLE TEST";
    use Storable;
    my $ice_ref = Storable::freeze($ref);
    my $water_ref = Storable::thaw($ice_ref);
    printn $ref->_DUMP();
    printn $water_ref->_DUMP();

    printn "UNTRANSLATE TEST";
    my $untranslated_ref = $ref->untranslate({
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

