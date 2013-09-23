#####################################################################################
# File:     Makefile
# Synopsys: BioJazz project. Runs tests and examples.
# Author:   Julien Ollivier
######################################################################################

######################################################################################
# SETUP
######################################################################################
ifeq ($(strip $(BIOJAZZ_HOME)),)
	#text-if-empty
	BIOJAZZ_HOME = .
endif

ANC_CMD = $(ANC_HOME)/anc.pl
FACILE_CMD = $(FACILE_HOME)/facile.pl

DIFF_CMD = git diff HEAD --

######################################################################################
# TESTING
######################################################################################
main_modules = \
			   test/modules/BitVector.log \
			   test/modules/BindingProfile.log \
			   test/modules/Sequence.log \
			   test/modules/Parser.log \
			   test/modules/Model.log \
			   test/modules/ProtoDomain.log \
			   test/modules/Domain.log \
			   test/modules/Gene.log \
			   test/modules/Genome.log \
			   test/modules/GenomeModel.log \
			   test/modules/Generation.log \
			   test/modules/MatlabDriver.log \
			   test/modules/Network.log \
			   test/modules/Scoring.log \
			   test/modules/ScorNode.log \
			   test/modules/ScorCluster.log \
			   test/modules/GenAlg.log \
			   test/modules/History.log \

custom_modules = \
				 test/custom/Template.log \
				 test/custom/Ultrasensitive.log \
				 test/custom/Oscillator.log \
				 test/custom/Linear.log \

test : test_all_modules test_main

test_main : FORCE
	-chmod -f +w test/biojazz.log
	-$(BIOJAZZ_HOME)/biojazz.pl --config=config/template.cfg --tag=test_main --cluster_type="LOCAL" --cluster_size=2 --seed=1092946010 2>&1 | tee -i test/biojazz.log
	-$(DIFF_CMD) test/biojazz.log

test_all_modules : test_main_modules test_custom_modules

test_main_modules : $(main_modules)
	@echo "Done running module tests."

test_custom_modules : $(custom_modules)
	@echo "Done running applicatin-specific custom module tests."

test/modules/%.log : FORCE
	@echo "Running $* testcase..."
	-mkdir -p test/modules
	-chmod -f +w test/modules/$*.log
	-perl -I$(ANC_HOME)/base -I$(BIOJAZZ_HOME)/modules -M$* -e '$*::run_testcases()' 2>&1 | cat > test/modules/$*.log
	-$(DIFF_CMD) test/modules/$*.log

test/custom/%.log : FORCE
	@echo "Running $* testcase..."
	-mkdir -p test/custom
	-chmod -f +w test/custom/$*.log
	-chmod -f +w test/custom/matlab/$*.mod
	-chmod -f +w test/custom/matlab/$*.eqn
	-perl -I$(ANC_HOME)/base -I$(BIOJAZZ_HOME)/modules -Icustom -M$* -e '$*::run_testcases()' 2>&1 | cat > test/custom/$*.log
	-$(DIFF_CMD) test/custom/$*.log
	-$(DIFF_CMD) test/custom/matlab/$*.mod
	-$(DIFF_CMD) test/custom/matlab/$*.eqn

test_custom : FORCE
	@echo "Running $(module) testcase..."
	-mkdir -p test/custom
	-chmod -f +w test/custom/$(module).log
	-chmod -f +w test/custom/matlab/$(module).mod
	-chmod -f +w test/custom/matlab/$(module).eqn
	-perl -I$(ANC_HOME)/base -I$(BIOJAZZ_HOME)/modules -Icustom -M$(module) -e '$(module)::run_testcases()' 2>&1 | tee test/custom/$(module).log

######################################################################################
# CODE ARCHITECTURE
######################################################################################

biojazz_class_hierarchy : FORCE
	-perl -I$(ANC_HOME)/base -Imodules -MUtils -e 'print Utils::sprint_class_hierarchy("biojazz_class_hierarchy.png", GenomeModel, Genome, Gene, Domain, ProtoDomain, GenomeInstance, GeneInstance, DomainInstance, ProtoDomainInstance, BindingProfile, BioJazz, BitString, BitVector, Generation, Model, Network, Sequence, Stimulus, Scoring, History, ScorNode, ScorCluster, Parser, ParserInstance)'

######################################################################################
# EVOLUTION, SHELL, ETC.
######################################################################################
cluster_type = LOCAL
cluster_size = 1
seed = -1
tag = test
config=ultrasensitive.cfg

evolve : FORCE
	@echo "ANC_HOME =" $(ANC_HOME)
	@echo "FACILE_HOME =" $(FACILE_HOME)
	mkdir -p $(basename $(config))
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --seed=$(seed) 2>&1 | tee -a -i $(basename $(config))/$(basename $(config)).$(tag).log --rescore

score : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --genome=$(genome) --score --store --shell 

score_generation : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --generation=$(generation) --shell

rescore : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --rescore --post_evolution

shell : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --seed=$(seed) --shell

load : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --genome=$(genome) --shell

collect : collect_from_networks

collect_from_logfile : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --command='collect_history_from_logfile("$(basename $(config))/$(basename $(config)).$(tag).log");'

collect_from_genomes : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --command='collect_history_from_genomes();'

collect_from_networks : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --command='collect_info_from_networks();'

######################################################################################
# ADMIN/UTILS
######################################################################################
retag : FORCE
	@mv $(basename $(config))/$(old_tag) $(basename $(config))/$(new_tag)
	@mv $(basename $(config))/$(basename $(config)).$(old_tag).log $(basename $(config))/$(basename $(config)).$(new_tag).log
	@echo "" | tee -a $(basename $(config))/$(basename $(config)).$(new_tag).log
	@echo "# ***************************************" | tee -a $(basename $(config))/$(basename $(config)).$(new_tag).log
	@echo "# RE-TAGGED FROM $(old_tag) to $(new_tag)" | tee -a $(basename $(config))/$(basename $(config)).$(new_tag).log
	@echo "# ***************************************" | tee -a $(basename $(config))/$(basename $(config)).$(new_tag).log
	@echo "" | tee -a $(basename $(config))/$(basename $(config)).$(new_tag).log
	@mv $(basename $(config))/$(basename $(config)).$(old_tag).xls $(basename $(config))/$(basename $(config)).$(new_tag).xls

ipckill : FORCE
	gcc -o ipckill_c ipckill.c
	chmod -f +x ipckill_c

clean: FORCE
	@rm -rf $(basename $(config))/$(tag)
	@rm -rf $(basename $(config))/$(basename $(config)).$(tag).log

clean_all: clean
	@rm -rf $(basename $(config))/scratch/$(tag)

archive:
	@mv $(basename $(config))/$(basename $(config)).$(tag).log ~/Archive/ultrasensitive
	@gzip -r $(basename $(config))/$(tag)
	@mv $(basename $(config))/$(tag) ~/Archive/ultrasensitive


######################################################################################
# DUMMY TARGET
######################################################################################

FORCE:

