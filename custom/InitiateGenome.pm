#-#####################################################################################
#- File:      InitiateGenome.pm
#- Synopsis: Initate the genome model and save as genome boject for later retrieving
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package InitiateGenome;
use Class::Std::Storable;
use base qw(GenomeModel);
{
    
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG);

    use GenomeModel;
    use storable qw(store, retrieve);

    
    sub initiate_genome {

	my $initiate_ref = $config_ref

    }






}
