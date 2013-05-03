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
	test/modules/ExcelSheet.log \
	test/modules/ExcelBook.log \
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
	-perl -I$(ANC_HOME)/base -Imodules -MUtils -e 'print Utils::sprint_class_hierarchy("biojazz_class_hierarchy.png", GenomeModel, Genome, Gene, Domain, ProtoDomain, GenomeInstance, GeneInstance, DomainInstance, ProtoDomainInstance)'

######################################################################################
# EVOLUTION, SHELL, ETC.
######################################################################################
cluster_type = SSH
cluster_size = 5
seed = -1
tag = test
config=ultrasensitive.cfg

evolve : FORCE
	@echo "ANC_HOME =" $(ANC_HOME)
	@echo "FACILE_HOME =" $(FACILE_HOME)
	mkdir -p $(basename $(config))
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --seed=$(seed) 2>&1 | tee -a -i $(basename $(config))/$(basename $(config)).$(tag).log

score : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --seed=$(seed) --shell --genome=$(genome) --score

shell : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --seed=$(seed) --shell

load : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --cluster_type=$(cluster_type) --cluster_size=$(cluster_size) --seed=$(seed) --shell --genome=$(genome)

collect : collect_from_genomes

collect_from_logfile : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --command='collect_history_from_logfile("$(basename $(config))/$(basename $(config)).$(tag).log"); save_history();'
collect_from_genomes : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --command='collect_history_from_genomes(); save_history();'
analyze : FORCE
	$(BIOJAZZ_HOME)/biojazz.pl --config=config/$(config) --tag=$(tag) --command='load_history(); export_history();'

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
	@rm -rf $(basename $(config))/history.$(tag).xls

clean_all: clean
	@rm -rf $(basename $(config))/scratch/$(tag)

archive:
	@mv $(basename $(config))/$(basename $(config)).$(tag).log ~/Archive/ultrasensitive
	@mv $(basename $(config))/history.$(tag).xls ~/Archive/ultrasensitive
	@gzip -r $(basename $(config))/$(tag)
	@mv $(basename $(config))/$(tag) ~/Archive/ultrasensitive


######################################################################################
# DUMMY TARGET
######################################################################################

FORCE:

