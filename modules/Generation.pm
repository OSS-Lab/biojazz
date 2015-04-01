#-#####################################################################################
#- File:     Generation.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Generation;
use Class::Std::Storable;
use base qw(Set ClassData);
{
    use Carp;
    use Storable qw(store retrieve);

    use Utils;

    use GenomeModel;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################
    Generation->set_class_data("ELEMENT_CLASS", "GenomeModel");

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: get_generation_glob
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_generation_glob {
        my $class = shift;

        my %args = (
            dir => undef,
            number => undef,
            @_,
        );
        check_args(\%args, 2);

        my $number = $args{number};
        my $dir = $args{dir};

        my $file_glob = sprintf("$dir/G%03d_I*.obj", $number);

        return $file_glob;
    }

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: retrieve_genomes
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub retrieve_genomes {
        my $self = shift;
        my %args = (
            files => undef,
            history_flag => 0,
            @_,
        );
        check_args(\%args, 2);

        my @files = @{$args{files}};
        my $history_flag = $args{history_flag};

        my $genome_class = $self->get_class_data("ELEMENT_CLASS");

        foreach my $file (@files) {
            my $genome_ref = retrieve("$file");  # retrieve() provided by Storable package
            $genome_ref->add_history("loaded from $file") if $history_flag;
            $self->add_element($genome_ref);
        }
    }


    #--------------------------------------------------------------------------------------
    # Function: store_genomes
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub store_genomes {
        my $self = shift;
        my %args = (
            files => undef,
            history_flag => 0,
            @_,
        );
        check_args(\%args, 2);

        my @files = @{$args{files}};
        my $history_flag = $args{history_flag};

        my @genomes = $self->get_elements();

        confess "ERROR: no. of filenames doesn't match no. of genomes" if @files != @genomes;

        for (my $i=0; $i < @genomes; $i++) {
            my $genome_ref = $genomes[$i];
            my $file = $files[$i];
            $genome_ref->add_history("saved to $file") if $history_flag;
            store($genome_ref, "$file"); # store() provided by Storable package
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: create_random_genomes
    # Synopsys: here is the subroutine that creating the genome.
    #           we could use the similar way to generate genomes of the first generation
    #           Besides the sequence there are some other information in the genome file
    #--------------------------------------------------------------------------------------
    sub create_random_genomes {
        my $self = shift; my $obj_ID = ident $self;
        my $config_ref = shift;

        my $inum = $config_ref->{inum_genomes};
        if ($config_ref->{selection_method} eq "kimura_selection") {
            $inum = 1;
        }

        my $genome_class = $self->get_class_data("ELEMENT_CLASS");
        for (my $i = 0; $i < $inum; $i++) {
            my $genome_ref  = $genome_class->new({
                    name => "UNDEF",
                    Genome => {
                        %$config_ref,
                        Gene => {
                            %$config_ref,
                            Domain => {
                                %$config_ref,
                                ProtoDomain => {
                                    %$config_ref,
                                },
                            },
                        },
                    },
                });
            my $gene_min_length = $genome_ref->get_gene_parser_ref->get_min_length();
            my $seq_length = int ((rand(16) + 8) * $gene_min_length);
            $genome_ref->get_sequence_ref()->generate_random_sequence($seq_length);
            $genome_ref->set_number(1);
            $genome_ref->add_history("random sequence created");
            $self->add_element($genome_ref);
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: load_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub load_generation {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            dir => undef,
            number => undef,
            @_,
        );
        check_args(\%args, 2);

        my $number = $args{number};
        my $dir = $args{dir};

        my $class = ref $self;
        my $file_glob = $class->get_generation_glob(
            dir => $dir,
            number => $number,
        );

        my @files = glob $file_glob;
        # need to sort on individual number
        @files = sort {$a=~/_I(\d+)/; my $a_i=$1; $b=~/_I(\d+)/; my $b_i=$1; return $a_i <=> $b_i} @files;
        my $num_genomes = @files;
        printn "load_generation: loading generation $number, with $num_genomes individuals";

        $self->retrieve_genomes(files => \@files);
    }

    #--------------------------------------------------------------------------------------
    # Function: save_generation
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub save_generation {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            dir => undef,
            number => undef,  # this generation's number by default
            @_,
        );
        check_args(\%args, 2);

        my $number = $args{number};
        my $dir = $args{dir};

        my $current_population_size = $self->get_num_elements();

        printn "save_generation: saving generation $number, with $current_population_size individuals";

        my @files = map {sprintf("$dir/G%03d_I%02d.obj", $number, $_)} (0..$current_population_size-1);

        $self->store_genomes(files => \@files);
    }

    #--------------------------------------------------------------------------------------
    # Function: refresh_individual_names
    # Synopsys: Refresh individual names based on given generation number and index.
    #--------------------------------------------------------------------------------------
    sub refresh_individual_names {
        my $self = shift; my $obj_ID = ident $self;
        my $number = shift;

        my @refs = $self->get_elements();
        for (my $i=0; $i < @refs; $i++) {
            my $ref = $refs[$i];
            my $name = sprintf("G%03d_I%02d", $number, $i);
            $ref->set_name($name);
            $ref->add_history("name refreshed to $name");
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: clear_genomes
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub clear_genomes {
        my $self = shift;

        $self->clear_elements();
    }

    #--------------------------------------------------------------------------------------
    # Function: get_attributes
    # Synopsys: This function gets a specific attribute's values for all individuals.
    #--------------------------------------------------------------------------------------
    sub get_attributes {
        my $self = shift;
        my $attribute = shift;

        my @genomes = $self->get_elements();
        my @attributes = map {$_->get_stats_ref()->{$attribute}} @genomes;
        return @attributes;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_attribute_names
    # Synopsys: This function collects defined keys from the stats_ref of
    #           each individual in the generation.
    #--------------------------------------------------------------------------------------
    sub get_attribute_names {
        my $self = shift;

        my @genomes = $self->get_elements();
        my @attribute_names = map {keys %{$_->get_stats_ref()}} @genomes;
        @attribute_names = union(\@attribute_names, []);  # this removes duplicates
        return @attribute_names;
    }

    #--------------------------------------------------------------------------------------
    # Function: get_attribute_stats
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub get_attribute_stats {
        my $self = shift;
        my $attribute = shift;

        my @genomes = $self->get_elements();
        my @attributes = map {$_->get_stats_ref()->{$attribute}} @genomes;

        my $mu = average(@attributes);
        my $sigma = sqrt(variance(@attributes));
        return ($mu, $sigma);
    }
}


sub run_testcases {

    my $ref = Generation->new({

        });
    printn $ref->_DUMP();
}


# Package BEGIN must return true value
return 1;

