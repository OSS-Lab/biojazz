##########################################################################
# Allosteric Network Compiler (ANC) model file
# Created by BioJazz version UNKNOWN released on UNKNOWN
# Mon 17 Sep 2012 14:12:36 BST
##########################################################################

##########################################################################
# PARAMETERS
##########################################################################
$max_external_iterations = 1
$max_internal_iterations = 0
$max_complex_size = 4
$max_species = 40
$max_csite_bound_to_msite_number = 2
$default_steric_factor = 1e-3
$export_graphviz = "network,collapse_states,collapse_complexes"


##########################################################################
# OBJECTS:
##########################################################################
#-------------------------
MODEL:  # G0000
#-------------------------
ReactionSite : {
    name => "PD0155",
    type => "msite",
    Keq_ratio => 0.0220844211988575,
}
AllostericStructure : {
    name => "D0128",
    allosteric_flag => 1,
    allosteric_transition_rates => [0.144974067037263,537.228111832403],
    Phi => 0.743070866141732,
    elements => [PD0155],
}
Structure : {  # IC of G0084 = 6.30957344480193e-11
    name => "G0084",
    elements => [D0128],
    max_count => -1,
}


ReactionSite : {
    name => "PD0345",
    type => "csite",
    Keq_ratio => 0.0120812379151953,
}
AllostericStructure : {
    name => "D0318",
    allosteric_flag => 1,
    allosteric_transition_rates => [0.475794431400941,401.027913949521],
    Phi => 0.99,
    elements => [PD0345],
}
Structure : {  # IC of G0274 = 1e-06
    name => "G0274",
    elements => [D0318],
    max_count => -1,
}


ReactionSite : {
    name => "PD1368",
    type => "bsite",
    Keq_ratio => 0.317297226293716,
}
AllostericStructure : {
    name => "D1341",
    allosteric_flag => 0,
    allosteric_transition_rates => [22.6380340952145,1291.54966501488],
    Phi => 0.372677165354331,
    elements => [PD1368],
}
Structure : {  # IC of G1297 = 1.58489319246111e-05
    name => "G1297",
    elements => [D1341],
    max_count => -1,
}


ReactionSite : {
    name => "PD1917",
    type => "csite",
    Keq_ratio => 0.084466834159386,
}
AllostericStructure : {
    name => "D1890",
    allosteric_flag => 0,
    allosteric_transition_rates => [2.82886943462597,166.810053720006],
    Phi => 0.28007874015748,
    elements => [PD1917],
}
Structure : {  # IC of G1846 = 0.000251188643150957
    name => "G1846",
    elements => [D1890],
    max_count => -1,
}


# Initial concentrations
Init : {
	structure => G0084,
	IC => 6.30957344480193e-11,
}
Init : {
	structure => G0274,
	IC => 1e-06,
}
Init : {
	structure => G1297,
	IC => 1.58489319246111e-05,
}
Init : {
	structure => G1846,
	IC => 0.000251188643150957,
}

#-------------------------
MODEL:  # LG0000
#-------------------------
ReactionSite : {
    name => "LPD0071",
    type => "bsite",
    Keq_ratio => 0.100225335130631,
}
AllostericStructure : {
    name => "LD0044",
    allosteric_flag => 0,
    allosteric_transition_rates => [100,0.0107583589854218],
    Phi => 0.550157480314961,
    elements => [LPD0071],
}
Structure : {  # IC of LG0000 = 1e-12
    name => "LG0000",
    elements => [LD0044],
    max_count => -1,
}


# Initial concentrations
Init : {
	structure => LG0000,
	IC => 1e-12,
}



##########################################################################
# RULES:
##########################################################################
MODEL:

CanBindRule : {
  name => "PD0345 PD0155 (0 R 0 R)",
  ligand_names => ["PD0345", "PD0155"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "R"],
  kf => 3162.27766016838,
  kb => 0.00630957344480193,
  kp => 8.96150501946605,
}

CanBindRule : {
  name => "PD0345 PD0155 (0 R 0 T)",
  ligand_names => ["PD0345", "PD0155"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "T"],
  kf => 100,
  kb => 1.58489319246112,
  kp => 8.96150501946605,
}

CanBindRule : {
  name => "PD0345 PD0155 (0 T 0 R)",
  ligand_names => ["PD0345", "PD0155"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "R"],
  kf => 100,
  kb => 0.0398107170553498,
  kp => 51.7947467923121,
}

CanBindRule : {
  name => "PD0345 PD0155 (0 T 0 T)",
  ligand_names => ["PD0345", "PD0155"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "T"],
  kf => 3.16227766016838,
  kb => 0.251188643150958,
  kp => 51.7947467923121,
}

CanBindRule : {
  name => "PD0155 LPD0071 (0 R 0 .)",
  ligand_names => ["PD0155", "LPD0071"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "."],
  kf => 3.16227766016838,
  kb => 1.58489319246112,
}

CanBindRule : {
  name => "PD0155 LPD0071 (0 T 0 .)",
  ligand_names => ["PD0155", "LPD0071"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "."],
  kf => 100,
  kb => 0.00630957344480193,
}

CanBindRule : {
  name => "PD0155 LPD0071 (1 R 0 .)",
  ligand_names => ["PD0155", "LPD0071"],
  ligand_msite_states => ["1", "0"],
  ligand_allosteric_labels => ["R", "."],
  kf => 3162.27766016838,
  kb => 0.00630957344480193,
}

CanBindRule : {
  name => "PD0155 LPD0071 (1 T 0 .)",
  ligand_names => ["PD0155", "LPD0071"],
  ligand_msite_states => ["1", "0"],
  ligand_allosteric_labels => ["T", "."],
  kf => 0.1,
  kb => 1.58489319246112,
}

CanBindRule : {
  name => "PD1917 PD1368 (0 . 0 .)",
  ligand_names => ["PD1917", "PD1368"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => [".", "."],
  kf => 3.16227766016838,
  kb => 10,
}


##########################################################################
# EQUATIONS:
##########################################################################
EQN:

null -> LG_0000_x; clamp_source_LG_0000_x="(0.5*0.001*1*(1 + square((t-10)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-15)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-20)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-25)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-30)/100*2*pi, 50)))/5"
LG_0000_x -> null; clamp_sink_LG_0000_x=1

##########################################################################
# CONFIG:
##########################################################################
CONFIG:

t_final = 100;
t_vector = [t0:0.1:tf];
matlab_ode_solver = ode23s;
matlab_odeset_options = odeset('InitialStep', 1e-8, 'AbsTol', 1e-9, 'RelTol', 1e-3, 'MaxStep', 10.0);
