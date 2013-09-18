### Introduction
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

### Authors and Contributors
The project is firstly developed by Julien Ollivier, then modified and currently maintained by Song Feng (@LifeWorks).

### Support or Contact
If have any problems, please contact @LifeWorks or email lifeworks.song@gmail.com