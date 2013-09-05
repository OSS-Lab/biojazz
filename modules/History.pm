#-#####################################################################################
#- File:     History.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package History;
use Class::Std::Storable;
use base qw();
{
    use Carp;

    use Text::CSV;
    use Utils;
    use Generation;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    # genome history : $ghistory_ref->{attribute}[$generation][$individual]
    my %ghistory_ref_of :ATTR(get => 'ghistory_ref', set => 'ghistory_ref');
    # population history : $phistory_ref->{attribute}[$generation]
    my %phistory_ref_of :ATTR(get => 'phistory_ref', set => 'phistory_ref');
    my %attribute_names_ref_of :ATTR(get => 'attribute_names_ref', set => 'attribute_names_ref');

    my %last_generation_of :ATTR(get => 'last_generation', set => 'last_generation');

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: BUILD
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

        $ghistory_ref_of{$obj_ID} = {};
        $phistory_ref_of{$obj_ID} = {};
    }

    #--------------------------------------------------------------------------------------
    # Function: collect_history_from_logfile
    # Synopsys: Collect history of genomes (i.e. their scores and other attributes).
    #--------------------------------------------------------------------------------------
    sub collect_history_from_logfile {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            logfile => undef,
            @_,
        );
        check_args(\%args,1);

        my $logfile = $args{logfile};

        my $ghistory_ref = $ghistory_ref_of{$obj_ID};

        open (LOGFILE, "< $logfile") or die "Can't open file $logfile\n";

        my $generation;
        my $individual;
        my @stats;
        my %attribute_names = ();
        while (<LOGFILE>) {
            my $line = $_;
            if ($line =~ /report_current_generation:\s+generation\s+(\S+)/) {
                $generation = $1 + 0;
                %attribute_names = ();
            }
            if ($line =~ /individual\s+(\S+)\s+:/) {
                $individual = $1;
                while ($line =~ /(\S+)\s*=\s*(\S*\d)/g) {
                    $ghistory_ref->{$1}[$generation][$individual] = $2;
                    $attribute_names{$1} = 1;
                }
                $attribute_names_ref_of{$obj_ID}[$generation] = [keys %attribute_names];
            }
        }
        $last_generation_of{$obj_ID} = $generation;
        close LOGFILE;
    }

    #--------------------------------------------------------------------------------------
    # Function: collect_history_from_genomes
    # Synopsys: Collect history of genomes (i.e. their scores and other attributes).
    #--------------------------------------------------------------------------------------
    sub collect_history_from_genomes {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            genomes_dir => undef,
            max_generations => -1,
            @_,
        );
        check_args(\%args,2);
        my $max_generations = $args{max_generations};

        my $ghistory_ref = $ghistory_ref_of{$obj_ID};

        my $generation = 0;
        my $progeny_ref = undef;
        while(1) {
            last if ($max_generations != -1) && ($generation >= $max_generations);
            my $gen_ref;
            if (defined $progeny_ref) {
                $gen_ref = $progeny_ref->{child_gen_ref};
            } else {
                $gen_ref = Generation->new({});
                $gen_ref->load_generation(
                    dir => $args{genomes_dir},
                    number => $generation,
                );
            }
            last if ($gen_ref->get_num_elements() == 0);
            $last_generation_of{$obj_ID} = $generation;

            # this function scans for the keys in each stats_ref of
            # each individual in the generation
            my @genome_attribute_names = $gen_ref->get_attribute_names();

            $attribute_names_ref_of{$obj_ID}[$generation] = \@genome_attribute_names;
            foreach (my $i=0; $i < @genome_attribute_names; $i++) {
                my $attribute_name = $genome_attribute_names[$i];
                my @attributes = $gen_ref->get_attributes($attribute_name);
                $ghistory_ref->{$attribute_name}[$generation] = [@attributes];
            }

            # PROGENY
            $progeny_ref = $gen_ref->get_progeny(
                dir => $args{genomes_dir},
                child_generation_number => $generation + 1,
            );
            if (defined $progeny_ref) {
                my @progeny_indices = @{$progeny_ref->{progeny_indices}};
                $ghistory_ref->{progeny_indices}[$generation] = [@progeny_indices];
                $ghistory_ref->{num_progeny}[$generation] = [map {scalar(@$_)} @progeny_indices];
            }
            $generation++;
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: analyze
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub analyze {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            population_attribute_names => [],
            @_,
        );
        check_args(\%args,1);

        my $ghistory_ref = $ghistory_ref_of{$obj_ID};
        my $phistory_ref = $phistory_ref_of{$obj_ID};

        my @population_attribute_names = @{$args{population_attribute_names}};

        my $last_generation = $last_generation_of{$obj_ID};
        for (my $generation = 0; $generation <= $last_generation; $generation++) {
            next if !defined $ghistory_ref->{score}[$generation];
            my @scores = @{$ghistory_ref->{score}[$generation]};
            @scores = map {defined $_ ? $_ : 0.0} @scores;
            my @score_sorted_indices = sort {$scores[$b] <=> $scores[$a]} (0..$#scores);
            my $top_scorer_index = $score_sorted_indices[0];
            $phistory_ref->{"top_score_index"}[$generation] = $top_scorer_index;
            my @attribute_names = @{$attribute_names_ref_of{$obj_ID}[$generation]};
            foreach my $population_attribute_name (@population_attribute_names) {
                $population_attribute_name =~ /^(\S+?)_(\S+)/;
                my $analysis_type = $1;
                my $attribute_name = $2;
                next if $population_attribute_name =~ /_index$/;  # skip
                next if $population_attribute_name =~ /slice.*_children$/;  # skip
                next if $population_attribute_name =~ /slice.*_robustness$/;  # skip
                next if $population_attribute_name =~ /slice.*_parents$/;  # skip
                next if $population_attribute_name =~ /slice.*_max$/;  # skip
                next if $population_attribute_name =~ /slice.*_avg$/;  # skip
                next if $population_attribute_name =~ /slice.*_min$/;  # skip
                next if !defined $ghistory_ref->{$attribute_name}[$generation];
                my @attributes = @{$ghistory_ref->{$attribute_name}[$generation]};
                SWITCH : {
                    if ($analysis_type =~ /^slice(\S+)/) {  # analysis of population within % of max for this attribute
                        my $percent_top = $1;
                        my @attributes = map {defined $_ ? $_ : 0.0} @{$ghistory_ref->{$attribute_name}[$generation]};
                        my @sorted_indices = sort {$attributes[$b] <=> $attributes[$a]} (0..$#attributes);
                        my $range_max = $attributes[$sorted_indices[0]];
                        my $range_min = $range_max*(1-$percent_top/100);

                        my @slice_indices = grep {($range_min <= $attributes[$_]) && ($attributes[$_] <= $range_max)} (0..$#attributes);
                        my @slice_attributes = map {$attributes[$_]} @slice_indices;

                        my @slice_progeny_indices = map {
                        if (defined $ghistory_ref->{progeny_indices}[$generation][$_]) {
                        @{$ghistory_ref->{progeny_indices}[$generation][$_]};
                        } else {
                        ();
                        }
                        } @slice_indices;
                        my $num_progeny = @slice_progeny_indices;
                        my @slice_progeny_attributes = map {$ghistory_ref->{$attribute_name}[$generation+1][$_]} @slice_progeny_indices;
                        my $num_children_in_range = 0;
                        foreach my $progeny_attribute (@slice_progeny_attributes) {
                            next if !defined $progeny_attribute;
                            next if $progeny_attribute < $range_min;
                            next if $progeny_attribute > $range_max;
                            $num_children_in_range++;
                        }
                        $phistory_ref->{"slice${percent_top}_${attribute_name}_children"}[$generation] = $num_progeny;
                        $phistory_ref->{"slice${percent_top}_${attribute_name}_robustness"}[$generation] = $num_children_in_range/$num_progeny if $num_progeny != 0;
                        $phistory_ref->{"slice${percent_top}_${attribute_name}_parents"}[$generation] = scalar(@slice_attributes);
                        $phistory_ref->{"slice${percent_top}_${attribute_name}_max"}[$generation] = $range_max;
                        $phistory_ref->{"slice${percent_top}_${attribute_name}_avg"}[$generation] = average(@slice_attributes);
                        $phistory_ref->{"slice${percent_top}_${attribute_name}_min"}[$generation] = $range_min;
                        last SWITCH;
                    }
                    if ($analysis_type eq "top") {
                        $phistory_ref->{"top_${attribute_name}"}[$generation] = $attributes[$top_scorer_index];
                        last SWITCH;
                    }
                    if ($analysis_type eq "max") {
                        my @sorted_indices = sort {
                        return 0  if !defined $attributes[$b] && !defined $attributes[$a];
                        return -1 if !defined $attributes[$b];
                        return +1 if !defined $attributes[$a];
                        $attributes[$b] <=> $attributes[$a]
                        } (0..$#attributes);
                        $phistory_ref->{"max_${attribute_name}"}[$generation] = $attributes[$sorted_indices[0]];
                        $phistory_ref->{"max_${attribute_name}_index"}[$generation] = $sorted_indices[0];
                        last SWITCH;
                    }
                    if ($analysis_type eq "median") {
                        $phistory_ref->{"median_${attribute_name}"}[$generation] = median(grep {defined $_} @attributes);
                        last SWITCH;
                    }
                    if ($analysis_type eq "avg") {
                        $phistory_ref->{"avg_${attribute_name}"}[$generation] = average(grep {defined $_} @attributes);
                        last SWITCH;
                    }
                    if ($analysis_type eq "sigma") {
                        $phistory_ref->{"sigma_${attribute_name}"}[$generation] = sqrt(variance(grep {defined $_} @attributes));
                        last SWITCH;
                    }
                    # default
                    printn "ERROR: unknown analysis_type $analysis_type";
                }
            }
        }
    }

    #--------------------------------------------------------------------------------------
    # Function: export
    # Synopsys: Export to excel spreadsheet.
    #--------------------------------------------------------------------------------------
    sub export {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            filename => undef,
            genome_attribute_names => [],
            population_attribute_names => [],
            @_,
        );
        check_args(\%args,3);

        # ANALYSIS
        $self->analyze(
            population_attribute_names => $args{population_attribute_names},
        );

        my $ghistory_ref = $ghistory_ref_of{$obj_ID};
        my $phistory_ref = $phistory_ref_of{$obj_ID};

        my $book_ref = ExcelBook->new({
                name => $args{filename},
            });

        # EXPORT OF POPULATION STATISTICS
        my $sheet_ref = $book_ref->add_sheet(
            name => "population",
            title => uc "POPULATION HISTORY",
            column_labels => ["GENERATION"],
        );
        my @population_attribute_names = @{$args{population_attribute_names}};
        foreach (my $k=0; $k < @population_attribute_names; $k++) {
            if (!defined $phistory_ref->{$population_attribute_names[$k]}) {
                printn "No data to export for $population_attribute_names[$k]";
                next;
            }
            my @population_attributes = @{$phistory_ref->{$population_attribute_names[$k]}};
            printn "History::export: exporting ".@population_attributes." generations of $population_attribute_names[$k] data";
            $sheet_ref->set_column($k+1, @population_attributes);
        }

        $sheet_ref->set_column(0, (0..$last_generation_of{$obj_ID}));
        $sheet_ref->set_column_labels([
                "GENERATION",
                @population_attribute_names,
            ]);

        # EXPORT OF PER-INDIVIDUAL HISTORY
        my @genome_attribute_names = @{$args{genome_attribute_names}};
        foreach (my $i=0; $i < @genome_attribute_names; $i++) {
            my $genome_attribute_name = $genome_attribute_names[$i];
            my $genome_attribute_ref = $ghistory_ref->{$genome_attribute_name};

            if (!defined $genome_attribute_ref) {
                printn "No data to export for $genome_attribute_name";
                next;
            }

            my $sheet_ref = $book_ref->add_sheet(
                name => "$genome_attribute_name",
                title => uc "HISTORY OF $genome_attribute_name",
                column_labels => ["GENERATION"],
            );

            my $max_population = 0;
            printn "History::export: exporting ".@{$genome_attribute_ref}." generations of $genome_attribute_name data";
            for (my $j=0; $j < @{$genome_attribute_ref}; $j++) {  # loop over generations
                if (defined $genome_attribute_ref->[$j]) {
                    my @genome_attributes = @{$genome_attribute_ref->[$j]};
                    $max_population = (@genome_attributes > $max_population) ? @genome_attributes : $max_population;
                    $sheet_ref->set_row($j, $j, @genome_attributes);
                } else {
                    printn "WARNING: no data for generation $j, attribute $genome_attribute_name";
                }
            }
            $sheet_ref->set_column_labels([
                    "GENERATION",
                    map {sprintf("I%02d",$_)} (0..$max_population-1),
                ]);
        }

        # WRITE IT OUT....
        $book_ref->write();
    }
}


sub run_testcases {
    printn "NO TESTCASES!!!";
}


# Package BEGIN must return true value
return 1;

