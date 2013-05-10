#-#####################################################################################
#- File:     ProtoDomain.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package ProtoDomain;
use Class::Std::Storable;
use base qw(Parser);
{
    use Carp;
    use Globals;
    use Utils;

    use BitString;
    use ProtoDomainInstance;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    ProtoDomain->set_class_data("ELEMENT_CLASS", "NONE");
    ProtoDomain->set_class_data("INSTANCE_CLASS", "ProtoDomainInstance");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    # parsing
    my %type_width_of :ATTR(get => 'type_width', set => 'type_width');
    my %substrate_polarity_width_of :ATTR(get => 'substrate_polarity_width', set => 'substrate_polarity_width');
    my %binding_profile_width_of :ATTR(get => 'binding_profile_width', set => 'binding_profile_width', init_arg => 'binding_profile_width');
    my %kf_profile_width_of :ATTR(get => 'kf_profile_width', set => 'kf_profile_width', init_arg => 'kf_profile_width');
    my %kb_profile_width_of :ATTR(get => 'kb_profile_width', set => 'kb_profile_width', init_arg => 'kb_profile_width');
    my %kp_profile_width_of :ATTR(get => 'kp_profile_width', set => 'kp_profile_width', init_arg => 'kp_profile_width');
    my %Keq_ratio_width_of :ATTR(get => 'Keq_ratio_width', set => 'Keq_ratio_width', init_arg => 'Keq_ratio_width');
    my %unused_width_of :ATTR(get => 'unused_width', set => 'unused_width');

    # translation scaling factors
    my %Keq_ratio_max_of :ATTR(get => 'Keq_ratio_max', set => 'Keq_ratio_max');
    my %Keq_ratio_min_of :ATTR(get => 'Keq_ratio_min', set => 'Keq_ratio_min');

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
    # Function: BUILD
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

        # DEFAULTS
        $type_width_of{$obj_ID} = 2;
        $substrate_polarity_width_of{$obj_ID} = 1;
        $binding_profile_width_of{$obj_ID} = 4;
        $kf_profile_width_of{$obj_ID} = 6;
        $kb_profile_width_of{$obj_ID} = 4;
        $kp_profile_width_of{$obj_ID} = 4;
        $Keq_ratio_width_of{$obj_ID} = 10;
        $unused_width_of{$obj_ID} = 23;

        $Keq_ratio_max_of{$obj_ID} = 1e3;
        $Keq_ratio_min_of{$obj_ID} = 1e-3;

        # INIT
        $type_width_of{$obj_ID} = $arg_ref->{type_width} if exists $arg_ref->{type_width};
        $substrate_polarity_width_of{$obj_ID} = $arg_ref->{substrate_polarity_width} if exists $arg_ref->{substrate_polarity_width};
        $binding_profile_width_of{$obj_ID} = $arg_ref->{binding_profile_width} if exists $arg_ref->{binding_profile_width};
        $kf_profile_width_of{$obj_ID} = $arg_ref->{kf_profile_width} if exists $arg_ref->{kf_profile_width};
        $kb_profile_width_of{$obj_ID} = $arg_ref->{kb_profile_width} if exists $arg_ref->{kb_profile_width};
        $kp_profile_width_of{$obj_ID} = $arg_ref->{kp_profile_width} if exists $arg_ref->{kp_profile_width};
        $Keq_ratio_width_of{$obj_ID} = $arg_ref->{Keq_ratio_width} if exists $arg_ref->{Keq_ratio_width};
        $unused_width_of{$obj_ID} = $arg_ref->{unused_width} if exists $arg_ref->{unused_width};

        $Keq_ratio_max_of{$obj_ID} = $arg_ref->{Keq_ratio_max} if exists $arg_ref->{Keq_ratio_max};
        $Keq_ratio_min_of{$obj_ID} = $arg_ref->{Keq_ratio_min} if exists $arg_ref->{Keq_ratio_min};

        # STRUCTURE
        $self->set_structure_ref([
                ["type",                   "\\G.{$type_width_of{$obj_ID}}"],
                ["substrate_polarity",     "\\G.{$substrate_polarity_width_of{$obj_ID}}"],
                ["binding_profile",        "\\G.{$binding_profile_width_of{$obj_ID}}"],
                ["kf_profile",             "\\G.{$kf_profile_width_of{$obj_ID}}"],
                ["kb_profile",             "\\G.{$kb_profile_width_of{$obj_ID}}"],
                ["kp_profile",             "\\G.{$kp_profile_width_of{$obj_ID}}"],
                ["Keq_ratio",              "\\G.{$Keq_ratio_width_of{$obj_ID}}"],
                ["kf_polarity_mask",       "\\G.{$kf_profile_width_of{$obj_ID}}"],
                ["kb_polarity_mask",       "\\G.{$kb_profile_width_of{$obj_ID}}"],
                # there is no kp_polarity_mask
                ["kf_conformation_mask",   "\\G.{$kf_profile_width_of{$obj_ID}}"],
                ["kb_conformation_mask",   "\\G.{$kb_profile_width_of{$obj_ID}}"],
                ["kp_conformation_mask",   "\\G.{$kp_profile_width_of{$obj_ID}}"],
                ["UNUSED",                 "\\G.{$unused_width_of{$obj_ID}}"],
            ]);
    }

    #--------------------------------------------------------------------------------------
    # Function: get_length
    # Synopsys: Compute the sum of lengths of fixed-length fields.  Excludes repeats
    #           and sub-sequences.  Includes {##}-type and bit-string patterns.
    #--------------------------------------------------------------------------------------
    sub get_length {
        my $self = shift;
        my $obj_ID = ident $self;

        my $length = 0;
        my $structure_ref = $self->get_structure_ref();

        foreach my $field_ref (@$structure_ref) {
            if ($field_ref->[1] =~ /\{(\S+)\}/) {
                my $field_length = $1;
                $length += $field_length;
            } else {
                confess "ERROR: can't compute length";
            }
        }
        return $length;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_*_width
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_kf_polarity_mask_width {my $self = shift; return $self->get_kf_profile_width()}
    sub get_kf_conformation_mask_width {my $self = shift; return $self->get_kf_profile_width()}
    sub get_kb_polarity_mask_width {my $self = shift; return $self->get_kb_profile_width()}
    sub get_kb_conformation_mask_width {my $self = shift; return $self->get_kb_profile_width()}
    sub get_kp_conformation_mask_width {my $self = shift; return $self->get_kp_profile_width()}

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
            if (($field_name eq "substrate_polarity") ||
                ($field_name eq "binding_profile") ||
                ($field_name eq "kf_profile") ||
                ($field_name eq "kb_profile") ||
                ($field_name eq "kp_profile") ||
                ($field_name =~ /polarity_mask/) ||
                ($field_name =~ /conformation_mask/) ||
                ($field_name eq "UNUSED")
            ) {
                # no untranslation necessary, just 0-padding or truncation
                my $lc_field_name = lc $field_name;  # for UNUSED field
                my $field_width = eval("\$self->get_${lc_field_name}_width()");
                confess "ERROR: internal error -- field_width not defined" if !defined $field_width;
                $field_sequence .= "0"x($field_width-length($field_value)).$field_value;  # 0-pad
                $field_sequence = substr $field_sequence, -$field_width; # truncate
                last SWITCH;
            }
            if (($field_name eq "Keq_ratio")
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
            if ($field_name eq "type") {
                if ($field_value eq "csite") {
                    $field_sequence = "10";
                    last SWITCH;
                }
                if ($field_value eq "bsite") {
                    $field_sequence = "00";
                    last SWITCH;
                }
                if ($field_value eq "msite") {
                    $field_sequence = "01";
                    last SWITCH;
                }
                confess "ERROR: unknown $field_name $field_value";
            }
            confess "ERROR: translate_field() -- unknown field $field_name";
        }			# SWITCH

        return $field_sequence;
    }
}

sub run_testcases {

    use Data::Dumper;

    $verbosity = 2;

    srand(1001);
    my $seq_ref = Sequence->new({});
    $seq_ref->generate_random_sequence(1000);
    my $ref = ProtoDomain->new({
            name => "ProtoDomain",
        });
    printn $ref->_DUMP();
    printn "Protodomain length is: ".$ref->get_length();
    my $iref = $ref->parse(sequence_ref => $seq_ref);
    $iref->translate();
    printn $iref->_DUMP();
    printn $iref->export_anc();

    # start at different position, don't clear
    my $i2ref = $ref->parse(sequence_ref => $seq_ref, start_pos => 200, dont_clear_flag => 1);
    $i2ref->translate();
    printn $i2ref->_DUMP();
    printn $i2ref->export_anc();

    $iref->set_allosteric_flag(1);
    $i2ref->set_allosteric_flag(1);
    printn join "\n", (ProtoDomainInstance->create_canbindrules(
            x_ref => $iref,
            y_ref => $i2ref,
            radius => 100,
            kf_max => 1e5,
            kf_min => 0.1,
            kb_max => 10,
            kb_min => 0.001,
            kp_max => 1e4,
            kp_min => 0.1,
        ));

    printn "STORABLE TEST";
    use Storable;
    my $ice_ref = Storable::freeze($ref);
    my $water_ref = Storable::thaw($ice_ref);
    printn $ref->_DUMP();
    printn $water_ref->_DUMP();
    map {printn $_->_DUMP()} ($water_ref->get_object_instances());

    printn "UNTRANSLATE TEST";
    my $untranslated_ref = $ref->untranslate({
            type => "csite",
            substrate_polarity => 1,
            binding_profile => "110",
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

