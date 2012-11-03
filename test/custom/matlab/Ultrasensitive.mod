##########################################################################
# Allosteric Network Compiler (ANC) model file
# Created by BioJazz version UNKNOWN released on UNKNOWN
# Sat  3 Nov 2012 14:12:33 GMT
##########################################################################

##########################################################################
# PARAMETERS
##########################################################################
$max_external_iterations = -1
$max_internal_iterations = -1
$max_complex_size = 4
$max_species = 80
$max_csite_bound_to_msite_number = 1
$default_steric_factor = 1e3
$export_graphviz = "network,collapse_states,collapse_complexes"


##########################################################################
# OBJECTS:
##########################################################################
#-------------------------
MODEL:  # G0000
#-------------------------
ReactionSite : {
    name => "PD0057",
    type => "bsite",
    Keq_ratio => 1.00451178020473,
}

ReactionSite : {
    name => "PD0227",
    type => "csite",
    Keq_ratio => 1.00451178020473,
}
AllostericStructure : {
    name => "D0022",
    allosteric_flag => 1,
    allosteric_transition_rates => [0.01,1.00451178020473],
    Phi => 1,
    elements => [PD0057,PD0227],
}
Structure : {  # IC of G0000 = 1.00677529813686
    name => "G0000",
    elements => [D0022],
    max_count => 2,
}


ReactionSite : {
    name => "PD0458",
    type => "csite",
    Keq_ratio => 1.9912245340072,
}
AllostericStructure : {
    name => "D0423",
    allosteric_flag => 0,
    allosteric_transition_rates => [1.00451178020473,1.00451178020473],
    Phi => 0,
    elements => [PD0458],
}
Structure : {  # IC of G0401 = 0.1
    name => "G0401",
    elements => [D0423],
    max_count => 2,
}


# Initial concentrations
Init : {
	structure => G0000,
	IC => 1.00677529813686,
}
Init : {
	structure => G0401,
	IC => 0.1,
}

#-------------------------
MODEL:  # LG0000
#-------------------------
ReactionSite : {
    name => "LPD0057",
    type => "bsite",
    Keq_ratio => 0.01,
}
AllostericStructure : {
    name => "LD0022",
    allosteric_flag => 0,
    allosteric_transition_rates => [0.01,0.01],
    Phi => 0,
    elements => [LPD0057],
}
Structure : {  # IC of LG0000 = 0.001
    name => "LG0000",
    elements => [LD0022],
    max_count => 2,
}


# Initial concentrations
Init : {
	structure => LG0000,
	IC => 0.001,
}

#-------------------------
MODEL:  # TG0000
#-------------------------
ReactionSite : {
    name => "TPD0057",
    type => "msite",
    Keq_ratio => 0.01,
}
AllostericStructure : {
    name => "TD0022",
    allosteric_flag => 0,
    allosteric_transition_rates => [0.01,0.01],
    Phi => 0,
    elements => [TPD0057],
}
Structure : {  # IC of TG0000 = 1.00677529813686
    name => "TG0000",
    elements => [TD0022],
    max_count => 2,
}


# Initial concentrations
Init : {
	structure => TG0000,
	IC => 1.00677529813686,
}



##########################################################################
# RULES:
##########################################################################
MODEL:

CanBindRule : {
  name => "PD0057 LPD0057 (0 R 0 .)",
  ligand_names => ["PD0057", "LPD0057"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "."],
  kf => 0.001,
  kb => 31.6227766016838,
}

CanBindRule : {
  name => "PD0057 LPD0057 (0 T 0 .)",
  ligand_names => ["PD0057", "LPD0057"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "."],
  kf => 31.6227766016838,
  kb => 31.6227766016838,
}

CanBindRule : {
  name => "PD0227 TPD0057 (0 R 0 .)",
  ligand_names => ["PD0227", "TPD0057"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "."],
  kf => 31.6227766016838,
  kb => 31.6227766016838,
  kp => 515.952796467086,
}

CanBindRule : {
  name => "PD0227 TPD0057 (0 T 0 .)",
  ligand_names => ["PD0227", "TPD0057"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "."],
  kf => 0.001,
  kb => 31.6227766016838,
  kp => 515.952796467086,
}

CanBindRule : {
  name => "PD0458 TPD0057 (0 . 1 .)",
  ligand_names => ["PD0458", "TPD0057"],
  ligand_msite_states => ["0", "1"],
  ligand_allosteric_labels => [".", "."],
  kf => 31.6227766016838,
  kb => 31.6227766016838,
  kp => 515.952796467086,
}


##########################################################################
# EQUATIONS:
##########################################################################
EQN:

null -> LG0000; clamp_source_LG0000="(+(event_flags(1) && ~event_flags(4))*min((t-event_times(1))/1000, 1)*3.33333333333333*4.0+event_flags(4)*max(1-(t-event_times(4))/1000, 0)*3.33333333333333*4.0+(event_flags(2) && ~event_flags(5))*min((t-event_times(2))/1000, 1)*3.33333333333333*4.0+event_flags(5)*max(1-(t-event_times(5))/1000, 0)*3.33333333333333*4.0+(event_flags(3) && ~event_flags(6))*min((t-event_times(3))/1000, 1)*3.33333333333333*4.0+event_flags(6)*max(1-(t-event_times(6))/1000, 0)*3.33333333333333*4.0)"
LG0000 -> null; clamp_sink_LG0000=4.0

##########################################################################
# CONFIG:
##########################################################################
CONFIG:

ode_event_times = ~ ~ ~ ~ ~ ~ ~;
SS_timescale = 500.0;
t_final = 20000;
t_vector = [t0:1:tf];
matlab_ode_solver = ode23s;
matlab_odeset_options = odeset('InitialStep', 1e-8, 'AbsTol', 1e-9, 'RelTol', 1e-3, 'MaxStep', 500.0);
