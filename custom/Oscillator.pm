#-#####################################################################################
#- File:     Oscillator.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Oscillator;
use Class::Std::Storable;
use base qw(Ultrasensitive);
{
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG);

    use Stimulus;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
#    my %A_of :ATTR(get => 'A', set => 'A', init_arg => 'A');  # constructor must supply initialization
#    my %B_of :ATTR(get => 'B', set => 'B', default => 'yyy'); # default value is yyy

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
    # Function: score_genome
    # Synopsys: Score for ultrasensitivity.
    #--------------------------------------------------------------------------------------
#    sub score_genome {
#	my $self = shift;
#	my $genome_model_ref = shift;
#
#    }

}


sub run_testcases {
    $verbosity = 1;

    $TAG = "test";
    srand(113);

    my $config_file = <<END;
#----------------------------------------
# CPU AND CLUSTER SETTINGS
#----------------------------------------
nice = 15
vmem = 450000

#----------------------------------------
# WORKSPACE AND CUSTOM SCORING MODULE
#----------------------------------------
scoring_class = Oscillator
work_dir = oscillator

#----------------------------------------
# GENOME PARAMS
#----------------------------------------

# Scaling: all concentrations in uM, all 2nd-order rates in uM^-1 s^-1

# Genome class
radius = 1
kf_max = 1e3    # uM^-1 s^-1
kf_min = 1e-3
kb_max = 1e3
kb_min = 1e-3
kp_max = 1e3
kp_min = 1e-3

# Gene class
regulated_concentration_width = 8
gene_unused_width = 8
regulated_concentration_max = 1e3    # 1mM
regulated_concentration_min = 1e-3   # 1nM  ~ 1 molecule in prokaryote

# Domain class
RT_transition_rate_width = 8
TR_transition_rate_width = 8
RT_phi_width = 8
domain_unused_width = 8
RT_transition_rate_max = 1e2
RT_transition_rate_min = 1e-2
TR_transition_rate_max = 1e2
TR_transition_rate_min = 1e-2
RT_phi_max = 1.0
RT_phi_min = 0.0

# ProtoDomain class
binding_profile_width = 8
kf_profile_width = 8
kb_profile_width = 8
kp_profile_width = 8
Keq_profile_width = 8
protodomain_unused_width = 8
Keq_ratio_max = 1e2
Keq_ratio_min = 1e-2

#----------------------------------------
# ANC PARAMS
#----------------------------------------
max_external_iterations = -1
max_internal_iterations = -1
max_complex_size = 4
max_species = 160
max_csite_bound_to_msite_number = 1
default_max_count = 2          # this prevents polymerization (see ANC manual)
default_steric_factor = 1e3
#export_graphviz = nothing
export_graphviz = network,collapse_states,collapse_complexes

#----------------------------------------
# GENERIC SIMULATION/SCORING PARAMETERS
#----------------------------------------
plot_input = 1
plot_output = 1
plot_species = 0
plot_phase = 1
plot_min = -1

sprint_transcript = 1

# MATLAB solver
solver = ode23s

# MATLAB odeset params
InitialStep = 1e-8
AbsTol = 1e-9
RelTol = 1e-3
MaxStep = 1.0

# MATLAB sampling interval
sampling_interval = 0.1

#----------------------------------------
# APPLICATION-SPECIFIC PARAMETERS
#----------------------------------------
lg_binding_profile = 01000011
tg_binding_profile = 10100011

steady_state_threshold = 0.01   # 1% relative change
steady_state_delta_t = 20

iss_check_time = 200.0          # initial steady-state check
iss_check_score_threshold = 0.99

delta_threshold = 0.02  # relative

LG_level = 10
LG_period = 4000       # in s
LG_strength = 0.5     # in Hz
LG_duty = 62.5
LG_rftime = 1000
LG_delay = 1000
LG_steps = 3

TG_init = 200  # uM

hill_n = 8
hill_k = 5

stimulus = staircase_equation
#stimulus = ramp_equation


END

    burp_file("test/modules/Oscillator.cfg", $config_file);

    my $scoring_ref = Oscillator->new({
	node_ID => 98,
	config_file => "test/modules/Oscillator.cfg",
	work_dir    => "test/modules",
	matlab_startup_options => "-nodesktop -nosplash",
    });

    printn $scoring_ref->_DUMP();

    my $config_ref = {};
    read_config($config_ref, "test/modules/Oscillator.cfg");

    use GenomeModel;
    my $genome_model_ref = GenomeModel->new({
	name => "Oscillator",
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
	PRE_JUNK => undef, POST_JUNK => "0000",
	genes => [
	    {
		START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
		regulated_concentration => 2.0, # uM
		UNUSED => "0000",
		domains => [
		    {
			allosteric_flag => 1,
			RT_transition_rate => 0.01,
			TR_transition_rate => 1.0,
			RT_phi => 1.0,
			protodomains => [
			    {
				type => "bsite",
				substrate_polarity => 0,
				binding_profile => BindingProfile->binding_complement($lg_binding_profile)->sprint(),
				kf_profile => "00000000",
				kb_profile => "11000000",
				kp_profile => "11111",
				Keq_ratio => 1.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "11111111",
				kb_conformation_mask => "00110000",
				kp_conformation_mask => "0",
				UNUSED => "0",
			    },
			    {
				type => "csite",
				substrate_polarity => 1,
				binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
				kf_profile => "11111111",
				kb_profile => "11000000",
				kp_profile => "11100000",
				Keq_ratio => 1.0,
				kf_polarity_mask => "0",
				kb_polarity_mask => "0",
				kf_conformation_mask => "11111100",
				kb_conformation_mask => "11000000",
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
		regulated_concentration => 1.0, # uM
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
				substrate_polarity => 0,
				binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
				kf_profile => "11111111",  # gives kf=1000, kb=31.6, kp=186 hence Km=~0.2
				kb_profile => "11000000",
				kp_profile => "11100000",
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

    # save the genome object
    use Storable qw(store retrieve);
    store($genome_model_ref, "test/custom/Oscillator.obj");

    $scoring_ref->score_genome($genome_model_ref);
    printn $genome_model_ref->_DUMP();
    sleep 20;
}


# Package BEGIN must return true value
return 1;

