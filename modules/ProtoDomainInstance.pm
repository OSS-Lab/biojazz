#-#####################################################################################
#- File:     ProtoDomainInstance.pm
#- Synopsys: Instance of ProtoDomain
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package ProtoDomainInstance;
use Class::Std::Storable;
use base qw(ParserInstance);
{
    use Carp;

    use Utils;
    use Globals;
    use BitString;

    use BindingProfile;

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %allosteric_flag_of                    :ATTR(get => 'allosteric_flag', set => 'allosteric_flag', default => 0);

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: create_canbindrules
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub create_canbindrules {
        my $class = shift;
        my %args = (
            x_ref => undef,
            y_ref => undef,
            radius => undef,
            kf_max => undef,
            kf_min => undef,
            kb_max => undef,
            kb_min => undef,
            kp_max => undef,
            kp_min => undef,
            @_,
        );

        check_args(\%args, 9);

        my $x_ref = $args{x_ref};
        my $y_ref = $args{y_ref};

        confess "ERROR: x_ref sequence not translated" if (!$x_ref->get_translation_valid_flag());
        confess "ERROR: y_ref sequence not translated" if (!$y_ref->get_translation_valid_flag());

        confess "ERROR: x_ref and y_ref were not created by the same parser" if $x_ref->get_parent_ref() != $y_ref->get_parent_ref();

        # ascending SORT in order csite, msite, then bsite
        my %type_rank = (
            csite => 0,
            msite => 1,
            bsite => 2,
        );
        ($x_ref, $y_ref) = sort {$type_rank{$a->get_translation_ref()->{type}} <=> $type_rank{$b->get_translation_ref()->{type}}} ($x_ref, $y_ref);
        my $xy_swapped = ($x_ref == $args{x_ref}) ? 0 : 1;

        my $x_name = $x_ref->get_name();
        my $y_name = $y_ref->get_name();
        my $x_translation_ref = $x_ref->get_translation_ref();
        my $y_translation_ref = $y_ref->get_translation_ref();
        my $catalytic_flag = (($x_translation_ref->{type} eq "csite") && ($y_translation_ref->{type} eq "msite")) ? 1 : 0;

        my @state_table = (
            []
        );
        my $x_msite_flag = ($x_translation_ref->{type} eq "msite");
        my $x_allosteric_flag = ($x_ref->get_allosteric_flag());
        my $y_msite_flag = ($y_translation_ref->{type} eq "msite");
        my $y_allosteric_flag = ($y_ref->get_allosteric_flag());

        @state_table = map {$x_msite_flag ? ([@$_, "0"], [@$_, "1"]) : [@$_, "0"]} @state_table;
        @state_table = map {$x_allosteric_flag ? ([@$_, "R"], [@$_, "T"]) : [@$_, "."]} @state_table;
        @state_table = map {$y_msite_flag ? ([@$_, "0"], [@$_, "1"]) : [@$_, "0"]} @state_table;
        @state_table = map {$y_allosteric_flag ? ([@$_, "R"], [@$_, "T"]) : [@$_, "."]} @state_table;

        my @rules;

        my $x_binding_profile = $x_translation_ref->{binding_profile};
        my $y_binding_profile = $y_translation_ref->{binding_profile};
        my $binding_profile_width = $x_ref->get_parent_ref()->get_binding_profile_width();
        confess "ERROR: inconsistent binding_profile_width" if $binding_profile_width != $y_ref->get_parent_ref()->get_binding_profile_width();
        my $binding_profile_mismatch = BindingProfile->mismatch($x_binding_profile,
            $y_binding_profile,
        );
        printn "create_canbindrule: binding_profile_mismatch=$binding_profile_mismatch for $x_name($x_binding_profile) and $y_name($y_binding_profile)" if $verbosity >= 2;

        if ($binding_profile_mismatch <= $args{radius}) {
            if ($verbosity >= 2) {
                if (!$xy_swapped) {
                    printn "create_canbindrules: $x_name($x_binding_profile) and $y_name($y_binding_profile) interact (binding_profile_mismatch=$binding_profile_mismatch)";
                } else {
                    printn "create_canbindrules: $y_name($y_binding_profile) and $x_name($x_binding_profile) interact (binding_profile_mismatch=$binding_profile_mismatch)";
                }
            }

            foreach my $state_ref (@state_table) {
                my ($x_msite_state, $x_allosteric_state, $y_msite_state, $y_allosteric_state) = @$state_ref;

                next if ($catalytic_flag && ($y_msite_state != $x_translation_ref->{substrate_polarity}));  # skip if incorrect substrate

                my $x_kf_profile = $x_translation_ref->{kf_profile};
                my $y_kf_profile = $y_translation_ref->{kf_profile};
                my $kf_profile_width = $x_ref->get_parent_ref()->get_kf_profile_width();
                confess "ERROR: inconsistent kf_profile_width" if $kf_profile_width != $y_ref->get_parent_ref()->get_kf_profile_width();
                my $x_kf_polarity_mask = $x_msite_state ? $x_translation_ref->{kf_polarity_mask} : "0" x $kf_profile_width;
                my $y_kf_polarity_mask = $y_msite_state ? $y_translation_ref->{kf_polarity_mask} : "0" x $kf_profile_width;
                my $x_kf_conformation_mask = $x_allosteric_state eq "T" ? $x_translation_ref->{kf_conformation_mask} : "0" x $kf_profile_width;
                my $y_kf_conformation_mask = $y_allosteric_state eq "T" ? $y_translation_ref->{kf_conformation_mask} : "0" x $kf_profile_width;
                my $kf_profile_mismatch = BindingProfile->mismatch($x_kf_profile,
                    $y_kf_profile,
                    $x_kf_polarity_mask,
                    $y_kf_polarity_mask,
                    $x_kf_conformation_mask,
                    $y_kf_conformation_mask,
                );
                # kf/kb mismatch represents fold-change from starting from kf_max/kb_min to kf_min/kb_max
                my $kf = loglinear($args{kf_max}, $args{kf_min}, $kf_profile_width, $kf_profile_mismatch);
                printn "create_canbindrule: state=(@$state_ref) kf=$kf kf_profile_mismatch=$kf_profile_mismatch for $x_name($x_kf_profile, $x_kf_polarity_mask, $x_kf_conformation_mask) and $y_name($y_kf_profile, $y_kf_polarity_mask, $y_kf_conformation_mask)" if $verbosity >= 2;

                my $x_kb_profile = $x_translation_ref->{kb_profile};
                my $y_kb_profile = $y_translation_ref->{kb_profile};
                my $kb_width = $x_ref->get_parent_ref()->get_kb_profile_width();
                confess "ERROR: inconsistent kb_width" if $kb_width != $y_ref->get_parent_ref()->get_kb_profile_width();
                my $x_kb_polarity_mask = $x_msite_state ? $x_translation_ref->{kb_polarity_mask} : "0" x $kb_width;
                my $y_kb_polarity_mask = $y_msite_state ? $y_translation_ref->{kb_polarity_mask} : "0" x $kb_width;
                my $x_kb_conformation_mask = $x_allosteric_state eq "T" ? $x_translation_ref->{kb_conformation_mask} : "0" x $kb_width;
                my $y_kb_conformation_mask = $y_allosteric_state eq "T" ? $y_translation_ref->{kb_conformation_mask} : "0" x $kb_width;
                my $kb_mismatch = BindingProfile->mismatch($x_kb_profile,
                    $y_kb_profile,
                    $x_kb_polarity_mask,
                    $y_kb_polarity_mask,
                    $x_kb_conformation_mask,
                    $y_kb_conformation_mask,
                );
                # kf/kb mismatch represents fold-change from starting from kf_max/kb_min to kf_min/kb_max
                my $kb = loglinear($args{kb_min}, $args{kb_max}, $kb_width, $kb_mismatch);
                printn "create_canbindrule: state=(@$state_ref) kb=$kb kb_mismatch=$kb_mismatch for $x_name($x_kb_profile, $x_kb_polarity_mask, $x_kb_conformation_mask) and $y_name($y_kb_profile, $y_kb_polarity_mask, $y_kb_conformation_mask)" if $verbosity >= 2;

                my $kp_line = "";
                if ($catalytic_flag) {
                    my $kp_profile = $x_translation_ref->{kp_profile};
                    my $kp_profile_width = $x_ref->get_parent_ref()->get_kp_profile_width();
                    my $kp_conformation_mask = $x_allosteric_state eq "T" ? $x_translation_ref->{kp_conformation_mask} : "0" x $kp_profile_width;
                    my $kp_masked = BindingProfile->XOR($kp_profile, $kp_conformation_mask)->sprint();
                    my $kp = loglinear($args{kp_min}, $args{kp_max}, (2**$kp_profile_width)-1, bin2dec($kp_masked));
                    printn "create_canbindrule: state=(@$state_ref) kp_final = $kp kp_masked=$kp_masked for $x_name($kp_profile, $kp_conformation_mask)" if $verbosity >= 2;
                    $kp_line = "  kp => $kp,\n";
                }

                my $str = <<END;
CanBindRule : {
  name => "${x_name} ${y_name} (@$state_ref)",
  ligand_names => ["$x_name", "$y_name"],
  ligand_msite_states => ["$x_msite_state", "$y_msite_state"],
  ligand_allosteric_labels => ["$x_allosteric_state", "$y_allosteric_state"],
  kf => $kf,
  kb => $kb,
$kp_line}
END

                push @rules, $str;
            }
        }
        return @rules;
    }

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: translate_field
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub translate_field {
        my $self = shift;
        my $field_name = shift;

        confess "ERROR: unexpected repeat sequence" if ref $field_name;

        my $transcription_ref = $self->get_transcription_ref();

        my $field_ref = $transcription_ref->{$field_name};
        my $field_sequence = $field_ref->[0];
        my $field_int_value = bin2dec($field_sequence);

        my $new_value;
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
                # no translation necessary
                $new_value = $field_sequence;
                last SWITCH;
            }
            if (($field_name eq "Keq_ratio")
            ) {
                my ($max, $min);
                eval "\$max = \$self->get_parent_ref->get_${field_name}_max";
                eval "\$min = \$self->get_parent_ref->get_${field_name}_min";
                confess "ERROR: internal error -- max not defined" if !defined $max;
                confess "ERROR: internal error -- min not defined" if !defined $min;
                $new_value = loglinear($min, $max, (2**length($field_sequence)) - 1, $field_int_value);
                last SWITCH;
            }
            if ($field_name eq "type") {
                if (($field_sequence eq "10") || ($field_sequence eq "11")) {
                    $new_value = "csite";
                    last SWITCH;
                }
                if ($field_sequence eq "00") {
                    $new_value = "bsite";
                    last SWITCH;
                }
                if ($field_sequence eq "01") {
                    $new_value = "msite";
                    last SWITCH;
                }
                confess "ERROR: unknown type $field_sequence";
            }
            confess "ERROR: translate_field() -- unknown field $field_name";
        }			# SWITCH

        return $new_value;
    }

    #--------------------------------------------------------------------------------------
    # Function: export_anc
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub export_anc {
        my $self = shift;

        confess "ERROR: sequence not translated" if (!$self->get_translation_valid_flag());


        my $name = $self->get_name();
        my $type = $self->get_translation_ref->{type};
        my $Keq_ratio = $self->get_translation_ref->{Keq_ratio};

        my $str = <<END;
ReactionSite : {
    name => "$name",
    type => "$type",
    Keq_ratio => $Keq_ratio,
}
END

        return $str;
    }
}


sub run_testcases {

}


# Package BEGIN must return true value
return 1;

