##########################################################################
# Allosteric Network Compiler (ANC) model file
# Created by BioJazz version UNKNOWN released on UNKNOWN
# Mon 17 Sep 2012 14:10:19 BST
##########################################################################

##########################################################################
# PARAMETERS
##########################################################################
$max_external_iterations = -1
$max_internal_iterations = -1
$max_complex_size = -1
$max_species = 200
$max_csite_bound_to_msite_number = 1
$default_steric_factor = 0.001
$export_graphviz = "nothing"


##########################################################################
# OBJECTS:
##########################################################################
#-------------------------
MODEL:  # G0000
#-------------------------
ReactionSite : {
    name => "PD0217",
    type => "csite",
    Keq_ratio => 3.04698957090351,
}
AllostericStructure : {
    name => "D0154",
    allosteric_flag => 0,
    allosteric_transition_rates => [0.00145956309563553,131.89690478391],
    Phi => 0.84652981427175,
    elements => [PD0217],
}
Structure : {  # IC of G0110 = 6.30957344480193e-08
    name => "G0110",
    elements => [D0154],
    max_count => -1,
}


ReactionSite : {
    name => "PD0496",
    type => "bsite",
    Keq_ratio => 27.5348489224465,
}
AllostericStructure : {
    name => "D0433",
    allosteric_flag => 0,
    allosteric_transition_rates => [0.260868888272655,232.578151753603],
    Phi => 0.00782013685239492,
    elements => [PD0496],
}
Structure : {  # IC of G0389 = 1.58489319246111e-11
    name => "G0389",
    elements => [D0433],
    max_count => -1,
}


ReactionSite : {
    name => "PD0689",
    type => "csite",
    Keq_ratio => 0.424195427069203,
}
AllostericStructure : {
    name => "D0626",
    allosteric_flag => 1,
    allosteric_transition_rates => [0.653505529858691,0.00800250227816105],
    Phi => 0.63049853372434,
    elements => [PD0689],
}
Structure : {  # IC of G0582 = 1.58489319246111e-05
    name => "G0582",
    elements => [D0626],
    max_count => -1,
}


ReactionSite : {
    name => "PD0988",
    type => "bsite",
    Keq_ratio => 0.00407352770587525,
}
AllostericStructure : {
    name => "D0925",
    allosteric_flag => 1,
    allosteric_transition_rates => [0.778924925037525,139.217879549033],
    Phi => 0.161290322580645,
    elements => [PD0988],
}
Structure : {  # IC of G0881 = 1e-12
    name => "G0881",
    elements => [D0925],
    max_count => -1,
}


ReactionSite : {
    name => "PD1758",
    type => "csite",
    Keq_ratio => 2.38945645393108,
}
AllostericStructure : {
    name => "D1695",
    allosteric_flag => 1,
    allosteric_transition_rates => [0.708663202813242,348.756245878594],
    Phi => 0.127077223851417,
    elements => [PD1758],
}
Structure : {  # IC of G1651 = 1.58489319246111e-11
    name => "G1651",
    elements => [D1695],
    max_count => -1,
}


# Initial concentrations
Init : {
	structure => G0110,
	IC => 6.30957344480193e-08,
}
Init : {
	structure => G0389,
	IC => 1.58489319246111e-11,
}
Init : {
	structure => G0582,
	IC => 1.58489319246111e-05,
}
Init : {
	structure => G0881,
	IC => 1e-12,
}
Init : {
	structure => G1651,
	IC => 1.58489319246111e-11,
}



##########################################################################
# RULES:
##########################################################################
MODEL:

CanBindRule : {
  name => "PD0217 PD0217 (0 . 0 .)",
  ligand_names => ["PD0217", "PD0217"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => [".", "."],
  kf => 100000,
  kb => 10,
}

CanBindRule : {
  name => "PD0689 PD0988 (0 R 0 R)",
  ligand_names => ["PD0689", "PD0988"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "R"],
  kf => 1000,
  kb => 1,
}

CanBindRule : {
  name => "PD0689 PD0988 (0 R 0 T)",
  ligand_names => ["PD0689", "PD0988"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "T"],
  kf => 100,
  kb => 10,
}

CanBindRule : {
  name => "PD0689 PD0988 (0 T 0 R)",
  ligand_names => ["PD0689", "PD0988"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "R"],
  kf => 10,
  kb => 0.1,
}

CanBindRule : {
  name => "PD0689 PD0988 (0 T 0 T)",
  ligand_names => ["PD0689", "PD0988"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "T"],
  kf => 1,
  kb => 1,
}

CanBindRule : {
  name => "PD0689 PD1758 (0 R 0 R)",
  ligand_names => ["PD0689", "PD1758"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "R"],
  kf => 1000,
  kb => 0.1,
}

CanBindRule : {
  name => "PD0689 PD1758 (0 R 0 T)",
  ligand_names => ["PD0689", "PD1758"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["R", "T"],
  kf => 10000,
  kb => 0.1,
}

CanBindRule : {
  name => "PD0689 PD1758 (0 T 0 R)",
  ligand_names => ["PD0689", "PD1758"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "R"],
  kf => 10,
  kb => 0.01,
}

CanBindRule : {
  name => "PD0689 PD1758 (0 T 0 T)",
  ligand_names => ["PD0689", "PD1758"],
  ligand_msite_states => ["0", "0"],
  ligand_allosteric_labels => ["T", "T"],
  kf => 100,
  kb => 1,
}


##########################################################################
# EQUATIONS:
##########################################################################
EQN:



##########################################################################
# CONFIG:
##########################################################################
CONFIG:

