# Scaffolded MAPK cascade patterned after yeast
# mating pathway (see Lodish et al.)
#
# Ad-hoc regulatory logic
#

###################################
MODEL:
###################################

#-------------------------------------------------
# COMPILE PARAMETERS
#-------------------------------------------------
#$max_csite_bound_to_msite_number = 1;

#$kf_1st_order_rate_cutoff = 1e-3;
#$kf_2nd_order_rate_cutoff = 1e-3;    # ensures scaffold-only phosphorylations

#-------------------------------------------------
# MODEL PARAMETERS
#-------------------------------------------------
Parameter : {
	name => Keq_S,
	value => 1e3,
}
Parameter : {
	name => kf_S,
	value => 1,
}
Parameter : {
	name => kb_S,
	value => "kf_S/Keq_S",
}

Parameter : {
	name => kf_KNX,
	value => 0.1,
}
Parameter : {
	name => kb_KNX,
	value => 0.5,
}
Parameter : {
	name => kp_KNX,
	value => 10.0,
}

Parameter : {
	name => kf_MAPKKK,
	value => 0.5e-3,
}
Parameter : {
	name => kb_MAPKKK,
	value => 0.5,
}
Parameter : {
	name => kp_MAPKKK,
	value => 10.0,
}

Parameter : {
	name => kf_MAPKK,
	value => 0.5e-3,
}
Parameter : {
	name => kb_MAPKK,
	value => 0.5,
}
Parameter : {
	name => kp_MAPKK,
	value => 10.0,
}

Parameter : {
	name => kf_MAPK,
	value => 1,
}
Parameter : {
	name => kb_MAPK,
	value => 0.5,
}
Parameter : {
	name => kp_MAPK,
	value => 10.0,
}

Parameter : {
	name => kf_P,
	value => 1,
}
Parameter : {
	name => kb_P,
	value => 0.5,
}
Parameter : {
	name => kp_P,
	value => 10.0,
}

Parameter : {
	name => kf_PY,
	value => 0.1,
}
Parameter : {
	name => kb_PY,
	value => 0.5,
}
Parameter : {
	name => kp_PY,
	value => 10.0,
}

#-------------------------------------------------
# INPUT
#-------------------------------------------------
ReactionSite: {
	name => KNX,
	type => "csite",
}
Structure: {name => X, elements => [KNX], max_count => 1}

#-------------------------------------------------
# OUTPUT
#-------------------------------------------------
ReactionSite: {
	name => MY,
	type => "msite",
}
Structure: {name => Y, elements => [MY], max_count => 1}

#-------------------------------------------------
# MAPKKK
#-------------------------------------------------
ReactionSite: {
	name => M1,
	type => "msite",
}
ReactionSite: {
	name => KN1,
	type => "csite",
}
ReactionSite: {
	name => S1n,
	type => "bsite",
}
# Protein
Structure: {name => KKK, elements => [S1n,M1,KN1], max_count => 1}

#-------------------------------------------------
# MAPKK
#-------------------------------------------------
ReactionSite: {
	name => M2a,
	type => "msite",
}
ReactionSite: {
	name => M2b,
	type => "msite",
}
ReactionSite: {
	name => KN2,
	type => "csite",
}
ReactionSite: {
	name => S2n,
	type => "bsite",
}
Structure: {name => KK, elements => [S2n,M2a,M2b,KN2], max_count => 1}

#-------------------------------------------------
# MAPK
#-------------------------------------------------
ReactionSite: {
	name => M3a,
	type => "msite",
}
ReactionSite: {
	name => M3b,
	type => "msite",
}
ReactionSite: {
	name => KN3,
	type => "csite",
}
ReactionSite: {
	name => S3n,
	type => "bsite",
}
# Protein
Structure: {name => K, elements => [S3n,M3a,M3b,KN3], max_count => 1}

#-------------------------------------------------
# PHOSPHATASE
#-------------------------------------------------
ReactionSite: {
	name => P1,
	type => "csite",
}
Structure: {name => P1, elements => [P1], max_count => 1}
ReactionSite: {
	name => P2,
	type => "csite",
}
Structure: {name => P2, elements => [P2], max_count => 1}
ReactionSite: {
	name => P3,
	type => "csite",
}
Structure: {name => P3, elements => [P3], max_count => 1}

ReactionSite: {
	name => PY,
	type => "csite",
}
Structure: {name => PY, elements => [PY], max_count => 1}


#-------------------------------------------------
# MAPK SCAFFOLD
#-------------------------------------------------
ReactionSite: {
	name => S1,
	type => "bsite",
}
ReactionSite: {
	name => S2,
	type => "bsite",
}
ReactionSite: {
	name => S3,
	type => "bsite",
}
Structure: {name => SS, elements => [S1, S2, S3], max_count => 1}

#-------------------------------------------------
# RULES
#-------------------------------------------------
# Scaffold association and dissociation
CanBindRule : {
	ligand_names => ['S1', 'S1n'], 
	kf => kf_S,
	kb => kb_S,
	steric_factor => 1e6,
	# constrain to free MAPKKK
	association_constraints => [
		'$R->get_in_object()->get_element(1)->is_bound() == 0',
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
}
CanBindRule : {
	ligand_names => ['S2', 'S2n'], 
	kf => kf_S,
	kb => kb_S,
	steric_factor => 1e6,
	# constrain to free MAPKK
	association_constraints => [
		'$R->get_in_object()->get_element(1)->is_bound() == 0',
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
}
CanBindRule : {
	ligand_names => ['S3', 'S3n'], 
	kf => kf_S,
	kb => kb_S,
	steric_factor => 1e6,
	# constrain to free MAPK
	association_constraints => [
		'$R->get_in_object()->get_element(1)->is_bound() == 0',
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
	],
}

# KNX (input) phosphorylates M1 (MAPKKK)
CanBindRule : {
	ligand_names => ['KNX', 'M1'], 
	ligand_msite_states => ['.', 0],
	kf => kf_KNX,
	kb => kb_KNX,
	kp => kp_KNX,
	association_constraints => [
		# not bound to substrate
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
}

# KN1 (MAPKKK) phosphorylates M2 (MAPKK) if M1 is phosphorylated
CanBindRule : {
	ligand_names => ['KN1', 'M2a'], 
	ligand_msite_states => ['.', 0],
	kf => kf_MAPKKK,
	kb => kb_MAPKKK,
	kp => kp_MAPKKK,
	association_constraints => [
		'$L->get_in_object()->get_element(1)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(1)->is_bound() == 0',
		# MAPKK not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPKK not bound to phosphatase
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
	steric_factor => 1e6,
}
CanBindRule : {
	ligand_names => ['KN1', 'M2b'], 
	ligand_msite_states => ['.', 0],
	kf => kf_MAPKKK,
	kb => kb_MAPKKK,
	kp => kp_MAPKKK,
	association_constraints => [
		'$L->get_in_object()->get_element(1)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(1)->is_bound() == 0',
		# MAPKK not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPKK not bound to phosphatase
		'$R->get_in_object()->get_element(1)->is_bound() == 0',
	],
	steric_factor => 1e6,
}

# KN2 (MAPKK) phosphorylates M3 (MAPK) if both M2a and M2b are phosphorylated
CanBindRule : {
	ligand_names => ['KN2', 'M3a'], 
	ligand_msite_states => ['.', 0],
	kf => kf_MAPKK,
	kb => kb_MAPKK,
	kp => kp_MAPKK,
	steric_factor => 1e6,
	association_constraints => [
		'$L->get_in_object()->get_element(1)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(1)->is_bound() == 0',
		'$L->get_in_object()->get_element(2)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(2)->is_bound() == 0',
		# MAPK not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPK not bound to phosphatase
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
}
CanBindRule : {
	ligand_names => ['KN2', 'M3b'], 
	ligand_msite_states => ['.', 0],
	kf => kf_MAPKK,
	kb => kb_MAPKK,
	kp => kp_MAPKK,
	steric_factor => 1e6,
	association_constraints => [
		'$L->get_in_object()->get_element(1)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(1)->is_bound() == 0',
		'$L->get_in_object()->get_element(2)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(2)->is_bound() == 0',
		# MAPK not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPK not bound to phosphatase
		'$R->get_in_object()->get_element(1)->is_bound() == 0',
	],
}

# KN3 (MAPK) phosphorylates MY (output) if both M3 sites are phosphorylated
CanBindRule : {
	ligand_names => ['KN3', 'MY'], 
	ligand_msite_states => ['.', 0],
	kf => kf_MAPK,
	kb => kb_MAPK,
	kp => kp_MAPK,
	steric_factor => 1e6,
	association_constraints => [
		'$L->get_in_object()->get_element(1)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(1)->is_bound() == 0',
		'$L->get_in_object()->get_element(2)->get_msite_state() == 1',
		'$L->get_in_object()->get_element(2)->is_bound() == 0',
	],
}

# MAPK phosphatases
CanBindRule : {
	ligand_names => ['P1', 'M1'],
	ligand_msite_states => ['.', 1],
	kf => kf_P,
	kb => kb_P,
	kp => kp_P,
	association_constraints => [
		## MAPKKK is not on scaffold
		#'$R->get_in_object()->get_element(0)->is_bound() == 0',
		# MAPKKK is not bound to substrate
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
}
CanBindRule : {
	ligand_names => ['P2', 'M2a'],
	ligand_msite_states => ['.', 1],
	kf => kf_P,
	kb => kb_P,
	kp => kp_P,
	association_constraints => [
		## MAPKK is not on scaffold
		#'$R->get_in_object()->get_element(0)->is_bound() == 0',
		# MAPKK is not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPKK not bound to kinase
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
}
CanBindRule : {
	ligand_names => ['P2', 'M2b'],
	ligand_msite_states => ['.', 1],
	kf => kf_P,
	kb => kb_P,
	kp => kp_P,
	association_constraints => [
		## MAPKK is not on scaffold
		#'$R->get_in_object()->get_element(0)->is_bound() == 0',
		# MAPKK is not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPKK not bound to kinase
		'$R->get_in_object()->get_element(1)->is_bound() == 0',
	],
}
CanBindRule : {
	ligand_names => ['P3', 'M3a'],
	ligand_msite_states => ['.', 1],
	kf => kf_P,
	kb => kb_P,
	kp => kp_P,
	association_constraints => [
		## MAPK is not on scaffold
		#'$R->get_in_object()->get_element(0)->is_bound() == 0',
		# MAPK is not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPK not bound to kinase
		'$R->get_in_object()->get_element(2)->is_bound() == 0',
	],
}
CanBindRule : {
	ligand_names => ['P3', 'M3b'],
	ligand_msite_states => ['.', 1],
	kf => kf_P,
	kb => kb_P,
	kp => kp_P,
	association_constraints => [
		## MAPK is not on scaffold
		#'$R->get_in_object()->get_element(0)->is_bound() == 0',
		# MAPK is not bound to substrate
		'$R->get_in_object()->get_element(3)->is_bound() == 0',
		# MAPK not bound to kinase
		'$R->get_in_object()->get_element(1)->is_bound() == 0',
	],
}

CanBindRule : {
	ligand_names => ['PY', 'MY'],
	ligand_msite_states => ['.', 1],
	kf => kf_PY,
	kb => kb_PY,
	kp => kp_PY,
}

#-----------------------------------------------------
# STIMULUS
#-----------------------------------------------------
Stimulus : {
	structure => 'X',
	type => "dose_response",
	strength => 1000,
	range => [1e-3,1e3],
	steps => 60,
	log_steps => 1,
}

#Stimulus : {
#	structure => 'Y',
#	state => "[,0]",
#	type => "clamp",
#	strength => 0.1,
#	concentration => 1,
#}

#Stimulus : {
#	structure => 'Y',
#	state => "[,1]",
#	type => "sink",
#	strength => 0.1,
#}

#-------------------------------------------------
# PROBE
#-------------------------------------------------
Probe : {
	name => X,
	structure => X,
}
Probe : {
	name => Y0,
	structure => Y,
	state => "[,0]",
}
Probe : {
	name => Y1,
	structure => Y,
	state => "[,1]",
}

Probe : {
	name => p_KKK1,
	classes => ComplexInstance,
	filters => [
 		'$_->get_exported_name() =~ /KKK1/',
        ],
}
Probe : {
	name => p_KK11,
	classes => ComplexInstance,
	filters => [
 		'$_->get_exported_name() =~ /(_|^)KK11/',
        ],
}
Probe : {
	name => p_K11,
	classes => ComplexInstance,
	filters => [
 		'$_->get_exported_name() =~ /(_|^)K11/',
        ],
}


#-------------------------------------------------
# INIT
#-------------------------------------------------
Init : {
	structure => X,
	IC => 0,
}
Init : {
	structure => Y,
	IC => 1,
}
Init : {
	structure => SS,
	IC => 1,
}
Init : {
	structure => KKK,
	IC => 1,
}
Init : {
	structure => KK,
	IC => 1,
}
Init : {
	structure => K,
	IC => 1,
}
Init : {
	structure => P1,
	IC => 1,
}
Init : {
	structure => P2,
	IC => 1,
}
Init : {
	structure => P3,
	IC => 1,
}
Init : {
	structure => PY,
	IC => 1,
}

################################
CONFIG:
################################
t_final = 100000
t_vector = [0:1:tf]

matlab_ode_solver = ode15s

