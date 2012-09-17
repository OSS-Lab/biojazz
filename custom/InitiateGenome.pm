#-#####################################################################################
#- File:      InitiateGenome.pm
#- Synopsis: Initate the genome model and save as genome boject for later retrieving
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package InitiateGenome;
use Class::Std::Storable;
use base qw();
{
    
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG $config_ref);

    use GenomeModel;
    use Storable qw(store retrieve);

    
    sub initiate_genome {

	my $genome_model_ref = GenomeModel -> new({
	name => "Ultrasensitive",
	Genome => {
	    radius => $config_ref->{radius},
	    kf_max => $config_ref->{kf_max},
	    kf_min => $config_ref->{kf_min},
	    kb_max => $config_ref->{kb_max},
	    kb_min => $config_ref->{kb_min},
	    kp_max => $config_ref->{kp_max},
	    kp_min => $config_ref->{kp_min},
	    Gene => {
		regulated_concentration_width => $config_ref->{regulated_concentration_width},
		unused_width => $config_ref->{gene_unused_width},
		regulated_concentration_max => $config_ref->{regulated_concentration_max},
		regulated_concentration_min => $config_ref->{regulated_concentration_min},
		Domain => {
		    RT_transition_rate_width => $config_ref->{RT_transition_rate_width},
		    TR_transition_rate_width => $config_ref->{TR_transition_rate_width},
		    RT_phi_width => $config_ref->{RT_phi_width},
		    unused_width => $config_ref->{domain_unused_width},
		    RT_transition_rate_max => $config_ref->{RT_transition_rate_max},
		    RT_transition_rate_min => $config_ref->{RT_transition_rate_min},
		    TR_transition_rate_max => $config_ref->{TR_transition_rate_max},
		    TR_transition_rate_min => $config_ref->{TR_transition_rate_min},
		    RT_phi_max => $config_ref->{RT_phi_max},
		    RT_phi_min => $config_ref->{RT_phi_min},
		    ProtoDomain => {
			binding_profile_width => $config_ref->{binding_profile_width},
			kf_profile_width => $config_ref->{kf_profile_width},
			kb_profile_width => $config_ref->{kb_profile_width},
			kp_profile_width => $config_ref->{kp_profile_width},
			Keq_profile_width => $config_ref->{Keq_profile_width},
			unused_width => $config_ref->{protodomain_unused_width},
			Keq_ratio_max => $config_ref->{Keq_ratio_max},
			Keq_ratio_min => $config_ref->{Keq_ratio_min},
		    },
		},
	    },
	},
    });

    # CONFIGURE/CREATE GENOME
    my $lg_binding_profile = $config_ref->{lg_binding_profile};
    my $tg_binding_profile = $config_ref->{tg_binding_profile};
    my $sequence_ref = $genome_model_ref->get_genome_parser_ref()->create_sequence({
	PRE_JUNK => undef,   # undef == ""
	POST_JUNK => "0000",
	genes => [
	    {
		START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
		regulated_concentration => 1.0, # uM
		UNUSED => "0000",
		domains => [
		    {
			allosteric_flag => 1,
			RT_transition_rate => 0.01,
			TR_transition_rate => 1.00,
			RT_phi => 1.0,
			protodomains => [
			    {
				type => "bsite",
				substrate_polarity => 0,
				binding_profile => BindingProfile->binding_complement($lg_binding_profile)->sprint(),
				kf_profile => "00000000000000000000",
				kb_profile => "11000000110000001000",
				kp_profile => "00011111000111110011",
				Keq_ratio => 1.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "11111100111111001110",
				kb_conformation_mask => "0",
				kp_conformation_mask => "0",
				UNUSED => "0",
			    },
			    {
				type => "csite",
				substrate_polarity => 0,
				binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
				kf_profile => "11111100111111001110",
				kb_profile => "11000000110000001000",
				kp_profile => "11111100111111001110",
				Keq_ratio => 1.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "11111100111111001110",
				kb_conformation_mask => "0",
				kp_conformation_mask => "0",
				UNUSED => "0",
			    },
			   ],
			UNUSED => "0",
		    },
		   ],
	    },
	    {
		START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
		regulated_concentration => 0.1, # uM
		UNUSED => "0000",
		domains => [
		    {
			allosteric_flag => 0,
			RT_transition_rate => 1.0,
			TR_transition_rate => 1.0,
			RT_phi => 0.0,
			protodomains => [
			    {
				type => "csite",
				substrate_polarity => 1,
				binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
				kf_profile => "11111100111111001110",
				kb_profile => "11000000110000001000",
				kp_profile => "11111100111111001110",
				Keq_ratio => 2.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "0",
				kb_conformation_mask => "0",
				kp_conformation_mask => "0",
				UNUSED => "0",
			    },
			    ],
			UNUSED => "0",
		    },
		    ],
	    },
	    ],

										   });

    printn "sequence=".$sequence_ref->get_sequence();
    $genome_model_ref->set_sequence_ref($sequence_ref);
	
	printn "$config_ref->{initial_genome} is specified\n";

	my $initial_genome_file = "test/custom/Ultrasensitive.obj";


	if ($config_ref->{initial_genome} =~ /load\s+(\S+)/) {
	    $initial_genome_file = $1;
	    
	} else {
	    printn "file to save the initiated genome is not specified\n";
	}
	
	store($genome_model_ref, "$initial_genome_file");
	
	return 1;
	
    }



}
