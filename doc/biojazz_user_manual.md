### Introduction

BioJazz is a tool for evolving and designing biochemical reaction networks using genetic algorithm (GA).  Typically, a BioJazz user wishes to evolve or design a small network or motif which accomplishes a specific function, such as a switch or an oscillator module.  The network comprises a set of proteins whose attributes are encoded in a network's ''genome''.  The ''genome'' is a binary string which  contains all the information necessary to determine how many proteins are present in the network, their structure, which proteins interact and the biochemical parameters of their interaction. 

BioJazz implements a genetic algorithm which through a random process of replication, mutation, and selection (with either population-based or random walk style), attempts to incrementally improve how well those "genomes" perform a user-specified function. By encoding the network in a fashion that mimics the way nature does, BioJazz can use a larger variety of mutational operators than do traditional GAs (which use point mutations and crossover), such as gene duplications, gene deletions, and domain shuffling.  Thus, BioJazz has the ability to change and evolve networks with respect to both topology and biochemical parameters, by starting from a designed network ''_de novo_'', or a partially or completely functional seed network.

Much of BioJazz's ability to design realistic networks comes from the accompanying Allosteric Network Compiler (ANC).  ANC is a stand-alone, rule-based compiler which has the ability to turn a high-level description of allosteric proteins into the corresponding set of biochemical equations.  The proteins can exhibit many of the behaviours observed in nature, such as co-localization, allosteric transitions, binding and catalytic reactions.  BioJazz is likely the first tool to couple a rule-based compiler with an evolutionary algorithm.  

The principal user input consists of a scoring function, which evaluates a particular genome against the desired functionality, and returns a score reflecting the network's performance. This score is compared against the score of other networks to determine whether the network survives to the next generation and replicates.  Generally speaking, this involves applying a stimulus to the network and evaluating it's response.  Simulation of the network is accomplished by integrating a set of ordinary differential equations (ODEs) in ''Matlab''.  The required Matlab files are automatically generated from ANC's output using a tool called Facile.

BioJazz is highly configurable.  For example, the user can specify evolutionary parameters such as mutation rates.  Also, the user may restrict BioJazz to changing a subset of the network's attributes.  This is useful to "freeze" the network topology, with the effect that only the network's biochemical parameters and not its structure are allowed to evolve.  

While the genetic algorithm itself is not very tasking, scoring each individual of a population of genomes may require a lot of processing power. Therefore, BioJazz has an integrated capability to use workstation clusters to speed the computation.

The main features of BioJazz are:

* evolves both network topology and connection weights
* designs a network ''_de novo_'', or starting from user-specified seed network
* uses workstation clusters to speed up the design
* produces a human-readable model of the winning network
* highly configurable and parameterized

###Getting Started

####Installation

BioJazz requires the ANC and Facile tools.  You can tell BioJazz where to get them by setting the ANC_HOME and FACILE_HOME environment variables to point to the appropriate directories.  It is recommended to add the following lines to your ''~/.bashrc'' file:
```sh
 export ANC_HOME    = ~/workspace/anc
 export FACILE_HOME = ~/workspace/facile
 alias  anc='$ANC_HOME/anc.pl'
 alias  facile='$FACILE_HOME/facile.pl'
 
 export BIOJAZZ_HOME = ~/workspace/biojazz
 alias  biojazz='$BIOJAZZ_HOME/biojazz.pl'
```

BioJazz requires Matlab to be installed on all nodes used for computation, and assumes matlab can be started with the command ''matlab''. Here is an example of configuration in "~/.bashrc" file (on Mac OS X):
```sh
 export MATLAB_HOME = /Applications/MATLAB_R2011b.app/bin
 alias  matlab='$MATLAB_HOME/matlab'
 export PATH = $MATLAB_HOME:$PATH

 DYLD_LIBRARY_PATH = /Applications/MATLAB_R2011b.app/bin/maci64:/Applications/MATLAB_R2011b.app/sys/os/maci64:/Applications/MATLAB_R2011b.app/runtime/maci64:$DYLD_LIBRARY_PATH
 export DYLD_LIBRARY_PATH
```
Note that if you decide to use a cluster of workstations, these installation instructions apply to all workstations used.

#####CPAN modules
CPAN is an internet database of Perl modules.  BioJazz/ANC/Facile use several of them and they must be installed prior to use.  You will need system administrator priviledges to install these modules (or see  for instructions on how to install them in your home directory).  You or your sysadmin will typically need to run the following commands on each system used (use sudo as prefix if available, if you don't have a admin privilege [here is a solution](http://twiki.org/cgi-bin/view/TWiki/HowToInstallCpanModules#Install_CPAN_modules_into_your_l) let you install perl modules in your user directory):
```sh
 cpan -i Class::Std
 cpan -i Class::Std::Storable
 cpan -i String::CRC32
 cpan -i Expect
 cpan -i Carp
 cpan -i WeakRef
 cpan -i IPC::Shareable
 cpan -i Linux::Pid
 cpan -i Text::CSV
```

Test your installation by running Facile, ANC and BioJazz without any arguments:
```sh
 $FACILE_HOME/facile.pl
 $ANC_HOME/anc.pl
 $BIOJAZZ_HOME/biojazz.pl
```

An error will be reported if any of the required modules are still missing.  Simply run CPAN again to install the missing module.

#####GraphViz

If you would like ANC to generate diagrams of the reaction network and species, you will the ''dot'' application and the following CPAN module:
```sh
 cpan -i GraphViz
```

#####SSH Agent

####Workspace Creation
Depending on your specific application, BioJazz will require some customized configuration and scoring functions.  Also, during the course of a single design run, BioJazz will generate a large number of files.  For this reason, the user must create a properly configured workspace which will contain the appropriate configuration files, scoring functions, and design files.

To facilitate this, BioJazz can create the workspace for and populate it with the required directories and with template files to get you started.  To do this, run the following command:
```sh
 biojazz --command='create_workspace("bjazz")'
```

This will create the directory `bjazz` and various sub-directories including `config` and `custom`.  Your configuration files go in the `config` directory, while your custom scoring functions go in the `custom` directory.  At this point, the user should familiarize him/herself with some the template files that are provided, and try to run BioJazz.  

The template file will try to design a network which contains a high concentration of dimers, and demonstrates how to use some of the functions available to the user.
```sh
 cd bjazz
 less config/template.cfg   # template configuration file
 less config/Template.pm    # template application-specific scoring function
```

####Running BioJazz

After taking a look around, try to run BioJazz.  The `cluster_type` and `cluster_size` arguments override the specification contained in the configuration file, and will launch both slave nodes of the cluster on your machine.
```sh
 biojazz --config=config/template.cfg --tag=first_try --cluster_type="LOCAL" --cluster_size=2
```

This will evolve the network for only a couple generations.  The `tag` argument is very important.  In BioJazz, each design attempt is associated with a specific, user-specified tag.  BioJazz will create a directory in your workspace containing all the results and other files generated during the optimization.  This allows the user to attempt several optimizations simultaneously without fear of accidental loss of files.  The name of the design's working directory is `work_dir/tag`.  The `work_dir` parameter is specified in your configuration file (and has a value of `template` in this example).  Thus the results of the above run are contained in the directory `template/first_try`.
```sh
 [user@host bjazz]$ ls -la template/test_main/
 total 168
 drwx------ 5 user group  4096 2008-06-03 14:53 .
 drwx------ 3 user group  4096 2008-06-03 14:51 ..
 drwx------ 2 user group  4096 2008-06-03 14:53 matlab
 drwx------ 2 user group  4096 2008-06-03 14:53 obj
 -rw------- 1 user group 71904 2008-06-03 14:53 ScorNode.0.20080603_145159.log
 -rw------- 1 user group 67215 2008-06-03 14:53 ScorNode.1.20080603_145159.log
 drwx------ 2 user group  4096 2008-06-03 14:51 source_2008-06-03-14:51:58
```

The `obj` directory contains all the genomes generated in a machine-readable form.  The `matlab` contains the models generated by ANC, and the `matlab` scripts generated by Facile.  The `ScorNode*` files are a log of the activity of each node in the cluster as work to compile and score genomes.  The `source*` directory is a snapshot of the source code used for that run such as your configuration and custom scoring files.

Now try modifying the configuration file to use other available workstations and run BioJazz again...

####Using a Makefile


####Scoring a specific genome

####Collecting and Analyzing Statistics

####BioJazz Shell


####Writing and Testing Application Specific Functions

#####Template Files

#####Testing

#####Debugging when things go wrong


####Workspace Directory Structure

```sh
 bjazz                               # workspace home
   config                            # configuration files
   custom                            # application-specific modules and functions (incl. scoring function)
   test/custom                       # recommended location for test results of custom modules
   test/modules                      # BioJazz module test results
   mydesign                          # application-specific directory
     mydesign.08jun01.log            # master node logfile
     08jun01                         # results directory for run with TAG=08jun01
       ScorNode.i.timestamp.log      # log file for slave node i
       matlab                        # ANC genome models, eqn files, and matlab files
       obj                           # genome objects in binary form
```

###Evolutionary Algorithm

BioJazz uses a genetic algorithm (GA) to search for a functional network.  The idea behind a genetic algorithm is to encode each candidate network as a ''genome'', which is just a sequence of bits having the value 0 or 1.  Given a population of individual network genomes, the GA calls a user-specified scoring function on each of the individuals in the population to determine it's ''fitness score''.  After computing the fitness scores, the GA creates a new generation of individual networks by cloning and mutating the previous generation in such a way that fitter individuals tend to have more descendants than less fit individuals.  Thus, depending on its relative fitness score, the fate of a particular individual network may be to:

  1. pass unchanged to the next generation (elite duplication)

  2. generate one or more mutated clones (asexual duplication)

  3. generate one or more children

  4. have no descendants at all (elimination)

Networks with a higher score has a better chance of fates i), ii) or iii) (i.e. surviving) than networks with a relatively lower score.

The number of creation created is such that the population size is kept constant from generation to generation.  The population size is specified via the parameter ''max_population''.

The following sections detail each step in the genetic algorithm.

####Initial Generation

The initial generation can be either generated randomly or loaded from disk, as specified by the `initial_generation` parameter of the configuration file.  In the random case, the user can also specify the number of individuals to create (parameter `inum_genomes`) and the genome length (parameter -- currently fixed at 5000).  Loading from disk is useful to resume work on a partially completed design starting from the last generation created, or to load hand-crafted seed designs.  The following shows some examples for each case:
```
    initial_genome = random                                    # random generation
    initial_genome = load test/modules/Ultrasensitive.obj      # load a hand-crafted network
    initial_genome = load ultrasensitive/test/obj/G427_I*.obj  # load all individuals of generation 427 of previous run
```

Regardless of how the initial generation is created, each individual is then stored under the following name in the working directory of the design:
```
    $DESIGN_WORK/obj/G`ggg`_I`ii`.obj
```
Where `ggg` is the generation number and `ii` is the individual number.

####Scoring

* the score_genome method of a user-defined scoring class is called to score individual networks
* the scoring class inherits many useful method from ''Scoring''
* describe Scoring class, its use of MatlabDriver, etc.

####Replication

There are two modes of replication in BioJazz.  They are rank-based and fitness-based replication.

#####Rank-based replication

In the rank-based mode, the fitness score returned by the user-defined evaluation routine are not used as is to compute the actual fitness of each network (the actual fitness being the expectation value for the number of children).  This mode is appropriate when small numerical differences in the score do not meaningfully correspond to the number of children an individual should have.  Therefore, BioJazz uses the fitness score of each individual to produce a ranking of the entire generation.

Children are generated according to the following procedure.  First, the top ''elite_pool_size'' individuals are copied as is to the next generation, without any mutations.  Next, an additional (''max_population'' - ''elite_pool_size'') children need to be created to keep the population size constant.  

#####Fitness-based replication

In this mode, the computed fitness dictates the relative proportion of the available pool of children that a genome is to have.

Children are generated according to the following procedure.  First, the top ''elite_pool_size'' individuals are copied as is to the next generation, without any mutations.  Next, an additional (''max_population'' - ''elite_pool_size'') children need to be created to keep the population size constant.  Each individual in the current generation is a potential parent and its expected number of children is proportional to the fitness computed from the ranking as described above.  

To decide which individuals become parents a sampling procedure called the "roulette wheel" is used.  Imagine a roulette wheel where each parent has a slot whose size is proportional to the fitness....


####Mutation Operators
#####Point Mutation
#####Gene Duplication
#####Gene Deletion
#####Domain Duplication
#####Domain Deletion
#####Domain Shuffling

###Network Encoding Model

####Complementarity Encoding Scheme

The complementarity encoding scheme is currently implemented in BioJazz.

| Field Name | Length (L) | Allowed Values RegExp | Scaling | Dynamic Range | Resolution | Description |
|-------------|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
|**Genome** |  |  |  |  |  |  |
|PRE_JUNK | Any | [01]* | | | |Untranslated sequence preceding first gene|
|*genes*| 1 or more genes | Any legal gene | | | | One or more genes separated by untranslated sub-sequences|
|POST_JUNK | Any | [01]* | | | | Untranslated sequence following last gene|
|**Genw** |  |  |  |  |  |  |
|START_CODE | 8 | 01111110 | | | | Untranslated sequence preceding first gene|
|concentration | 10 | [01]{L} | Log-linear | 1e3/1e-3 | 1.36% | Encodes the initial concentration of the protein|
|UNUSED | 4 | [01]{L} | | | | Reserved field|
|*domains* | 1 or more domains | Any legal domain | | | | One or more domains separated by a soft linker pattern ("001")|
|STOP_CODE | 3 | 111 | | | | Terminates the gene|
|**Domain** |  |  |  |  |  |  |
|allosteric flag | 1 | [01]{1} | | | | Maps directly to ANC domain attribute of the same name|
|R <--> T rates | 10 | [01]{L} | Loglinear | 1e+2/1e-2 | 0.90% | Embedded in ANC model|
|&Phi; | 10 | [01]{10} | Linear | 0.0 - 1.0 | ~1e-3 | Scaled to a number between 0 and 1 (inclusive), embedded in ANC model|
|UNUSED | 4 | [01]{L} | | | | Reserved field|
|*protodomains* | 1 or more protodomains | Any legal protodomain | | | |One or more protodomains separated by a hard linker pattern ("000")|
|**ProtoDomain** |  |  |  |  |  |  |
|type | 2 | [01]{2} | | | | Maps to protodomain type in ANC model, determines enzyme polarity. 00=bsite, 01=msite, 10=csite, 11=csite|
|substrate polarity | 1 | [01]{1} | | | | If the protodomain is a csite, determines whether it modifies (0, kinase) or unmodifies (1, phosphatase) the substrate|
|binding profile | 10 | [01]{L} | | | | Determines ligands pairs. If the binding profiles of two protodomains are sufficiently complementary, then a binding reaction occurs using kinetic rates calculated from the following attributes|
|kf profile | 20 | [01]{L} | Log-linear profile | 1e3/1e-3 | 1.99x | Determines the association kinetics of two protodomains|
|kb profile | 20 | [01]{L} | Loglinear-profile | 1e3/1e-3 | 1.99x | Determines the dissociation kinetics of two protodomains|
|kp profile | 10 | [01]{L} | Loglinear | 1e3/1e-3 | 1.36% | For csite protodomains only, determines rate of product reaction when modifying|
|Keq ratio (&Gamma;) | 10 | [01]{10} | Loglinear | 1e2/1e-2 | 0.90% | Determines allosteric effect of msite modification. Embedded in ANC model|
|kf polarity mask | 20 | [01]{L} | | | | XORed with "kf profile" to determine profile of modified (msite=1) version of protodomain|
|kf conformation mask | 20 | [01]{L} | | | | XORed with "kf profile" to determine new profile of 'T' conformation|
|kb polarity mask | 20 | [01]{L} | | | | XORed with "kb profile" to determine profile of modified (msite=1) version of protodomain|
|kb conformation mask | 20 | [01]{8} | | | | XORed with "kb profile" to determine new profile of 'T' conformation|
|kp conformation mask | 10 | [01]{8} | | | | XORed with "kp profile" to determine new profile of 'T' conformation|
| UNUSED | 4 | [01]{L} | | | | Reserved field|


####Field Scaling

All the above fields correspond to parameter values whose scale and dynamic range is given in a configuration file.  This section deals with the issue of translating an integer (binary encoded) scale into parameter values.  In all sections, x is the integer value corresponding to a binary field, which goes from 0 to x_max=2^L-1 where L is the length of the field.  Also, y_min is the minimum value of the parameter, and y_max is the maximum value.  For example, the regulated_concentration parameter can be made to range from 1e-3 (y_min) to 1e+3 (y_max).

#####Linearly Scaled Fields

Relevant field(s): &Phi;.

Here the value of the parameter is y(x) = y_min + (y_max - y_min)*(x/x_max);

The resolution R(n) for an n-bit field is given by:

R(n) = (y_max-y_min)/(2^n-1)

#####Loglinear Scaled Fields

Relevant field(s): kp_profile, RT/TR_transition_rate, Keq_ratio, regulated_concentration

y(x) = y_min * (y_max/y_min)^(x/x_max), or

log(y) = log(y0) + (x/x_max)*log(y_max/y_min)

The fold-change between two values of x which differ by 1.0 is:

FC = y(x+1)/y(x) = (y_max/y_min)^(1/x_max).

For example, an L=8-bit field associated with a dynamic range of y_max/y_min = 1e6 implies xmin=0, x_max=2^L-1=255 and a fold-change per increment of x of:

FC(8bits) = 1e6^(1/255) = 1.0557

In other words, you can change by 5.57% increments.  For 10 bits,

FC(10bits) = 1e6^(1/1023) = 1.0136.

#####Loglinear Scaled Profiles

Fields: kf_profile, kb_profile

Profiles do not directly encode parameters, but instead are compared to each other after applying a special complement() function to either one of the arguments.  The result of this comparison is used to obtain a Hamming distance of 0 to at most L.  This Hamming distance is then scaled log-linearly.  Thus if two profiles are 8 bits long, x_min=0 and x_max=L=8, so for a dynamic range of 1e6 we have

FC = 1e6^(1/8) = 5.6234

I.e. a 560% change per increment.  More appropriate would be to use 10-30 bits:

FC(10) = 1e6^(1/10) = 3.9811

FC(20) = 1e6^(1/20) = 1.9953

FC(30) = 1e6^(1/30) = 1.5849

I.e. a 2x or 60% change per increment.

######Improvement

Note that there is some inefficiency in this method of coding because we need extremely long fields to get good resolution,  because x_max=L and not 2^L.

Currently for profiles X and Y, we have:

MISMATCH_VECTOR = XOR(X, complement(Y))
MISMATCH = ones(MISMATCH_VECTOR)

Equivalently, 

MISMATCH = HAMMING(X, complement(Y))

An alternative would scale the MISMATCH_VECTOR, which has a range of 2^L, instead of MISMATCH, which has a range of just L.  Caveat is that we must get the same thing if the profiles are commuted, but in fact the MISMATCH_VECTOR gets flipped, so the msb becomes the lsb.  One way around this is to designate the 2 outer bits as msbs, and the 2 most inner bits as lsbs.  You break the MISMATCH_VECTOR in two, flip the 2nd half and add to the first.  This gets us a range of 0 to [2^(L/2)-1]*2.

###Scoring Functions

### Authors and Contributors
The project is firstly developed by Julien Ollivier, then modified and currently maintained by Song Feng (@LifeWorks).

### Support or Contact
If have any problems, please contact @LifeWorks or email lifeworks.song@gmail.com