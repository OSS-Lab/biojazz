##########################################################################
# Allosteric Network Compiler (ANC) model file
# Created by BioJazz version UNKNOWN released on UNKNOWN
# Thu Jun  5 14:55:27 EDT 2008
##########################################################################

##########################################################################
# PARAMETERS
##########################################################################
$max_external_iterations = -1
$max_internal_iterations = -1
$max_complex_size = 4
$max_species = 160
$max_csite_bound_to_msite_number = 1
$default_steric_factor = 1e3
$export_graphviz = "network,collapse_states,collapse_complexes"


##########################################################################
# OBJECTS:
##########################################################################
#-------------------------
MODEL:  # G0000
#-------------------------
ProtoDomain : {
    name => "PD0057",
    type => "bsite",
    Keq_ratio => 1.00451178020473,
}

ProtoDomain : {
    name => "PD0155",
    type => "csite",
    Keq_ratio => 1.00451178020473,
}
Domain : {
    name => "D0024",
    allosteric_flag => 1,
    RT_transition_rate => 0.01,
    TR_transition_rate => 1.01822354964927,
    RT_phi => 1,
    elements => [PD0057,PD0155],
}
Protein : {  # IC of G0000 = 1.96841944728661
    name => "G0000",
    elements => [D0024],
    max_count => 2,
}


ProtoDomain : {
    name => "PD0314",
    type => "csite",
    Keq_ratio => 1.9912245340072,
}
Domain : {
    name => "D0281",
    allosteric_flag => 0,
    RT_transition_rate => 1.01822354964927,
    TR_transition_rate => 1.01822354964927,
    RT_phi => 0,
    elements => [PD0314],
}
Protein : {  # IC of G0257 = 1.02745948544618
    name => "G0257",
    elements => [D0281],
    max_count => 2,
}


#-------------------------
INIT:  # G0000
#-------------------------
G0000 = 1.96841944728661
G0257 = 1.02745948544618

#-------------------------
MODEL:  # LG0000
#-------------------------
ProtoDomain : {
    name => "LPD0057",
    type => "bsite",
    Keq_ratio => 0.01,
}
Domain : {
    name => "LD0024",
    allosteric_flag => 0,
    RT_transition_rate => 0.01,
    TR_transition_rate => 0.01,
    RT_phi => 0,
    elements => [LPD0057],
}
Protein : {  # IC of LG0000 = 0.001
    name => "LG0000",
    elements => [LD0024],
    max_count => 2,
}


#-------------------------
INIT:  # LG0000
#-------------------------
LG0000 = 0.001

#-------------------------
MODEL:  # TG0000
#-------------------------
ProtoDomain : {
    name => "TPD0057",
    type => "msite",
    Keq_ratio => 0.01,
}
Domain : {
    name => "TD0024",
    allosteric_flag => 0,
    RT_transition_rate => 0.01,
    TR_transition_rate => 0.01,
    RT_phi => 0,
    elements => [TPD0057],
}
Protein : {  # IC of TG0000 = 196.841944728661
    name => "TG0000",
    elements => [TD0024],
    max_count => 2,
}


#-------------------------
INIT:  # TG0000
#-------------------------
TG0000 = 196.841944728661



##########################################################################
# RULES:
##########################################################################
MODEL:

CanBindRule : {
  name => "PD0057 LPD0057 (0 R 0 .)",
  ligand_class => ProtoDomain,
  ligand_names => ["PD0057", "LPD0057"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_states => ["R", "."],
  kf => 0.001,
  kb => 31.6227766016838,
}

CanBindRule : {
  name => "PD0057 LPD0057 (0 T 0 .)",
  ligand_class => ProtoDomain,
  ligand_names => ["PD0057", "LPD0057"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_states => ["T", "."],
  kf => 1000,
  kb => 1,
}

CanBindRule : {
  name => "PD0155 TPD0057 (0 R 1 .)",
  ligand_class => ProtoDomain,
  ligand_names => ["PD0155", "TPD0057"],
  ligand_msite_states => ["0", "1"],
  ligand_allosteric_states => ["R", "."],
  kf => 1000,
  kb => 31.6227766016838,
  kp => 186.461097142696,
}

CanBindRule : {
  name => "PD0155 TPD0057 (0 T 1 .)",
  ligand_class => ProtoDomain,
  ligand_names => ["PD0155", "TPD0057"],
  ligand_msite_states => ["0", "1"],
  ligand_allosteric_states => ["T", "."],
  kf => 0.0316227766016838,
  kb => 1000,
  kp => 186.461097142696,
}

CanBindRule : {
  name => "PD0314 TPD0057 (0 . 0 .)",
  ligand_class => ProtoDomain,
  ligand_names => ["PD0314", "TPD0057"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_states => [".", "."],
  kf => 1000,
  kb => 31.6227766016838,
  kp => 186.461097142696,
}


##########################################################################
# EQUATIONS:
##########################################################################
EQN:

null -> LG0000_x; clamp_source_LG0000_x="(0.5*10*0.5*(1 + square((t-1000)/4000*2*pi, 37.5)) + 0.5*10*0.5*(1 + square((t-1500)/4000*2*pi, 37.5)) + 0.5*10*0.5*(1 + square((t-2000)/4000*2*pi, 37.5)))/3"
LG0000_x -> null; clamp_sink_LG0000_x=0.5

##########################################################################
# CONFIG:
##########################################################################
CONFIG:

event_times = 1000 1500 2000 2500 3000 3500;
t_final = 4000;
t_vector = [t0:0.1:tf];
matlab_ode_solver = ode23s;
matlab_odeset_options = odeset('InitialStep', 1e-8, 'AbsTol', 1e-9, 'RelTol', 1e-3, 'MaxStep', 1.0);
