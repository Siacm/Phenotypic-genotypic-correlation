#!usr/bin/perl
use warnings;
use strict;

use Term::ANSIColor qw(:constants);
use Text::CSV_XS;
use less '-R';
use Getopt::Std;
use vars qw($opt_a $opt_b $opt_c $opt_d $opt_e $opt_f $opt_g $opt_h $opt_i $opt_j $opt_k $opt_m $opt_n $opt_o);
getopts('abcdefghijkmno');

#variables
my (%column_index, %mismatched_pheno, %mismatched_pheno2, %mismatched_pheno3, %mismatched_geno, %mismatched_decrease, %matched_geno, %matched_geno2,%matched_geno3, %matched_decrease, %matched_decrease2, %matched_decrease3, %maybe_sens, %maybe_resis, %maybe_decrease, %hash_insuff);
my (%misp_serovar, %misp2_serovar, %misp3_serovar, %misg_serovar, %misd_serovar, %matp_serovar, %matg_serovar, %matg2_serovar, %matg3_serovar, %matd_serovar, %matd2_serovar, %matd3_serovar, %maybesens_serovar, %maybedecrease_serovar, %mayberesis_serovar);
my (%dupmisp_serovar, %dupmisp2_serovar, %dupmisp3_serovar, %dupmisg_serovar, %dupmatp_serovar, %dupmisd_serovar, %dupmatg_serovar, %dupmatg2_serovar, %dupmatg3_serovar, %dupmatd_serovar, %dupmatd2_serovar, %dupmatd3_serovar, %dupmays_serovar, %dupmayd_serovar, %dupmayr_serovar);
my (%matchs_rpattern, %matchr_rpattern, %matchr2_rpattern, %matchr3_rpattern, %matchd_rpattern, %matchd2_rpattern, %matchd3_rpattern, %miss_rpattern, %miss2_rpattern, %miss3_rpattern, %misr_rpattern, %misd_rpattern, %maybes_rpattern, %mayber_rpattern, %maybed_rpattern);
my (%anti_matched_geno, %anti_matched_geno2, %anti_matched_geno3, %anti_matched_decrease, %anti_matched_decrease2, %anti_matched_decrease3, %anti_maybe_resis, %anti_maybe_decrease);
my (%rev_misp, %rev_misp2, %rev_misp3, %rev_misg, %rev_misd, %rev_matp, %rev_matg, %rev_matg2, %rev_matg3, %rev_matd, %rev_matd2, %rev_matd3, %rev_maybes, %rev_mayber, %rev_maybed);
my (@table, @matched_pheno,@i_sensitive,@insufficient);
my $num_pheno_discordant=0;
my $num_pheno2_discordant=0;
my $num_pheno3_discordant=0;
my $num_geno_discordant=0;
my $num_decrease_discordant=0;
my $num_pheno_match=0;
my $num_geno_match=0;
my $num_geno2_match=0;
my $num_geno3_match=0;
my $num_decrease_match=0;
my $num_decrease2_match=0;
my $num_decrease3_match=0;
my $num_geno_maybe=0;
my $num_pheno_maybe=0;
my $num_decrease_maybe=0;
my $num_insuff=0;
my ($input, $output);

if ($opt_h) {
    print "\nINSTRUCTIONS\n\n";
    print "This file can only be used for reading contents of a csv file with quotes. The file must only contain one antibiotic class with the corresponding resistance genes for that specific class.\n";
    print "Sample ID must be in column 1, Sample serovar in column 4, and Sample ST in column 7.\n";
    print "Phenotypic results must be reported as S/R/I/NA, genotypic results must be reported as YES if the resistance gene is present.\n\n";
    print "Two command-line arguments are necessary, except for the help option: \n";
    print "-h                Brings up the help guide\n";
    print "-a <filename.csv> <outfile.tsv>      Prints the number of serovars for each resistance gene profiles i.e mcr-1: ST19 (n=2) ST34(n=3)\n";
    print "-b <filename.csv> <outfile.tsv>      Prints the number of samples for each resistance gene profiles i.e mcr-1 (n=12) \n";
    print "-c <filename.csv> <outfile.tsv>      Prints the count for each serovar type i.e ST19 (n=12) \n";
    print "-d <filename.csv> <outfile.tsv>      Prints sample IDs for each resistance gene profile i.e mcr-1: 2018-0000 2018-0001\n";
    print "-e <filename.csv> <outfile.tsv>      Prints sample IDs with their resistance gene profile i.e 2018-0000: mcr-1, 2018-0001: mcr-1, blaTEM\n";
    print "-f <filename.csv> <outfile.tsv>        Prints the sample IDs for each serovar type for all discordant and concordant pheno/genotypic data i.e ST19: 2018-0000 2018-0001 \n";
    print "-g <filename.csv> <outfile.tsv>        Prints the sample IDs for each resistant antbiotic i.e ampicillin i.e 2018-0000 2018-0001 \n";
    print "-i <filename.csv> <outfile.tsv>      Prints the number of serovars for the overall resistance gene profiles i.e [mcr-1, blaTEM]: ST19(n=2)\n";
    print "-j <filename.csv> <outfile.tsv>      Prints the number of samples for overall resistance gene profiles i.e [mcr-1, blaTEM] (n=1) \n";
    print "-k <filename.csv> <outfile.tsv>        Prints the serovar ST and corresponding Sample IDs and resistance gene profile i.e ST19: 2018-0000 - mcr-1, blaTEM\n";
    print "-m <filename.csv> <outfile.tsv>        Prints the serovar ST and resistance gene profile including the count ST19: [mcr-1, blaTEM] (n=2), [blaTEM] (n=1) \n";
    print "-n <filename.csv> <outfile.tsv>        Prints single resistance genes and serovar ST, i.e mef(B) ST19 (n=1), mef(A) ST19 (n=4)\n";
    print "-o <filename.csv> <outfile.tsv>        Prints the serovar ST and single genes, i.e ST19 mef(B) (n=1), ST34 mef(A) (n=2) \n";
    exit;
}
elsif (@ARGV !=2 ) {                        #if no flag or file is entered
    print "\nNo valid input. Please enter -h for help options\n\n";
    exit;
}

$input = $ARGV[0];    #entered file
$output= $ARGV[1];    #output file

if (defined $input) {
    unless ($input=~ /.csv$/){
        print "\nInvalid file format. Please enter -h for help options\n";
        exit;
    }
}

print "Please wait...\n";
print "Printing results in $output\n\n";


open(my $INFILE, '<:encoding(UTF-8)', $input) or die "Cannot open file $input\n\n";
open(OUTFILE,">>$output") or die "Cannot open file $output\n\n";

#get column headers in hash
my $csv_lines = Text::CSV_XS->new({ binary => 1, sep_char => ','});
my @column_headers = @{$csv_lines->getline($INFILE)};
my $i=0;
for my $header (@column_headers) {
    $column_index{$i}=$header;
    $i++;
}

while (my $rows=<$INFILE>) { #while loop to go through each line of input file
    if ($csv_lines->parse($rows)){  #checks whether the lines can be parsed
    my $line= $csv_lines->string($rows);  #each line is converted to a string and assigned to a variable
    my @elements = $csv_lines->fields(); #each word in each line
        if ($line =~/^(?=.*(\"\bR\b\")+)(?=.*(\"\bYES\b\").*(\"\bYES\b\").*(\"\bYES\b\"))/gi) {    #resistant to abx phenotypically with at least two resistance genes
            push (@{$matchr3_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$matchr3_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $n=0;$n<@elements;$n++){ #loop through each word
                if ($elements[$n] =~/\bR\b/g){ #if word matches uppercase R, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_matched_geno3{$elements[0]}},$n) unless grep{$_ == $n} @{$anti_matched_geno3{$elements[0]}};
                }
                elsif ($elements[$n] =~/\bYES\b/gi){  #if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                    push (@{$matched_geno3{$elements[0]}},$n) unless grep{$_ == $n} @{$matched_geno3{$elements[0]}};
                    next if exists $dupmatg3_serovar{$elements[0]};
                    push (@{$matg3_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$matg3_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmatg3_serovar{$elements[0]}=1;
                }
            }
            $num_geno3_match=keys%matched_geno3; #counts the number of keys -gives the count of each matched geno
        }
        elsif ($line =~/^(?=.*(\"\bR\b\")+)(?=.*(\"\bYES\b\").*(\"\bYES\b\"))/gi) {    #resistant to abx phenotypically with at least two resistance genes
            push (@{$matchr2_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$matchr2_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $n=0;$n<@elements;$n++){ #loop through each word
                if ($elements[$n] =~/\bR\b/g){ #if word matches uppercase R, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_matched_geno2{$elements[0]}},$n) unless grep{$_ == $n} @{$anti_matched_geno2{$elements[0]}};
                }
                elsif ($elements[$n] =~/\bYES\b/gi){  #if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                    push (@{$matched_geno2{$elements[0]}},$n) unless grep{$_ == $n} @{$matched_geno2{$elements[0]}};
                    next if exists $dupmatg2_serovar{$elements[0]};
                    push (@{$matg2_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$matg2_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmatg2_serovar{$elements[0]}=1;
                }
            }
            $num_geno2_match=keys%matched_geno2; #counts the number of keys -gives the count of each matched geno
        }
        elsif ($line =~/^(?=.*(\"\bR\b\")+)(?=.*(\"\bYES\b\"))/gi) {    #resistant to abx phenotypically with resistance genes
            push (@{$matchr_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$matchr_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $n=0;$n<@elements;$n++){ #loop through each word
                if ($elements[$n] =~/\bR\b/g){ #if word matches uppercase R, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_matched_geno{$elements[0]}},$n) unless grep{$_ == $n} @{$anti_matched_geno{$elements[0]}};
                }
                elsif ($elements[$n] =~/\bYES\b/gi){  #if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                    push (@{$matched_geno{$elements[0]}},$n) unless grep{$_ == $n} @{$matched_geno{$elements[0]}};
                    next if exists $dupmatg_serovar{$elements[0]};
                    push (@{$matg_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$matg_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmatg_serovar{$elements[0]}=1;
                }
            }
            $num_geno_match=keys%matched_geno; #counts the number of keys -gives the count of each matched geno
        }
        elsif ($line =~/^(?=.*(\"\bS\b\")+)(?=.*(\"\bYES\b\").*(\"\bYES\b\").*(\"\bYES\b\"))/gi) {#sensitive to abx phenotypically but has at least two resistance genes
            push (@{$miss3_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$miss3_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $l=0;$l<@elements;$l++){ #loop through each word
                next unless $elements[$l] =~/\bYES\b/gi;#if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                push (@{$mismatched_pheno3{$elements[0]}},$l) unless grep{$_ == $l} @{$mismatched_pheno3{$elements[0]}};
                next if exists $dupmisp3_serovar{$elements[0]};
                push (@{$misp3_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$misp3_serovar{$elements[2]." ST".$elements[3]}};
                $dupmisp3_serovar{$elements[0]}=1;
            }
            $num_pheno3_discordant=keys %mismatched_pheno3; #counts the number of keys -gives the count of each mismatched pheno
        }
        elsif ($line =~/^(?=.*(\"\bS\b\")+)(?=.*(\"\bYES\b\").*(\"\bYES\b\"))/gi) {#sensitive to abx phenotypically but has at least two resistance genes
            push (@{$miss2_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$miss2_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $l=0;$l<@elements;$l++){ #loop through each word
                next unless $elements[$l] =~/\bYES\b/gi;#if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                push (@{$mismatched_pheno2{$elements[0]}},$l) unless grep{$_ == $l} @{$mismatched_pheno2{$elements[0]}};
                next if exists $dupmisp2_serovar{$elements[0]};
                push (@{$misp2_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$misp2_serovar{$elements[2]." ST".$elements[3]}};
                $dupmisp2_serovar{$elements[0]}=1;
            }
            $num_pheno2_discordant=keys %mismatched_pheno2; #counts the number of keys -gives the count of each mismatched pheno
        }
        elsif ($line =~/^(?=.*(\"\bS\b\")+)(?=.*(\"\bYES\b\"))/gi) {#sensitive to abx phenotypically but has resistance genes
            push (@{$miss_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$miss_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $l=0;$l<@elements;$l++){ #loop through each word
                next unless $elements[$l] =~/\bYES\b/gi;#if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                push (@{$mismatched_pheno{$elements[0]}},$l) unless grep{$_ == $l} @{$mismatched_pheno{$elements[0]}};
                next if exists $dupmisp_serovar{$elements[0]};
                push (@{$misp_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$misp_serovar{$elements[2]." ST".$elements[3]}};
                $dupmisp_serovar{$elements[0]}=1;
            }
            $num_pheno_discordant=keys %mismatched_pheno; #counts the number of keys -gives the count of each mismatched pheno
        }
        elsif ($line =~/^(?=.*(\"\bD\b\")+)(?=.*(\"\bYES\b\").*(\"\bYES\b\").*(\"\bYES\b\"))/gi) {    #decreased to abx phenotypically with two or more resistance genes
            push (@{$matchd3_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$matchd3_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $n=0;$n<@elements;$n++){ #loop through each word
                if ($elements[$n] =~/\bD\b/g){ #if word matches uppercase D, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_matched_decrease3{$elements[0]}},$n) unless grep{$_ == $n} @{$anti_matched_decrease3{$elements[0]}};
                }
                elsif ($elements[$n] =~/\bYES\b/gi){  #if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                    push (@{$matched_decrease3{$elements[0]}},$n) unless grep{$_ == $n} @{$matched_decrease3{$elements[0]}};
                    next if exists $dupmatd3_serovar{$elements[0]};
                    push (@{$matd3_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$matd3_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmatd3_serovar{$elements[0]}=1;
                }
            }
            $num_decrease3_match=keys%matched_decrease3; #counts the number of keys -gives the count of each matched decreased
        }
        elsif ($line =~/^(?=.*(\"\bD\b\")+)(?=.*(\"\bYES\b\").*(\"\bYES\b\"))/gi) {    #decreased to abx phenotypically with two or more resistance genes
            push (@{$matchd2_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$matchd2_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $n=0;$n<@elements;$n++){ #loop through each word
                if ($elements[$n] =~/\bD\b/g){ #if word matches uppercase D, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_matched_decrease2{$elements[0]}},$n) unless grep{$_ == $n} @{$anti_matched_decrease2{$elements[0]}};
                }
                elsif ($elements[$n] =~/\bYES\b/gi){  #if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                    push (@{$matched_decrease2{$elements[0]}},$n) unless grep{$_ == $n} @{$matched_decrease2{$elements[0]}};
                    next if exists $dupmatd2_serovar{$elements[0]};
                    push (@{$matd2_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$matd2_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmatd2_serovar{$elements[0]}=1;
                }
            }
            $num_decrease2_match=keys%matched_decrease2; #counts the number of keys -gives the count of each matched decreased
        }
        elsif ($line =~/^(?=.*(\"\bD\b\")+)(?=.*(\"\bYES\b\"))/gi) {    #decreased to abx phenotypically with one resistance gene
            push (@{$matchd_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$matchd_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $n=0;$n<@elements;$n++){ #loop through each word
                if ($elements[$n] =~/\bD\b/g){ #if word matches uppercase D, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_matched_decrease{$elements[0]}},$n) unless grep{$_ == $n} @{$anti_matched_decrease{$elements[0]}};
                }
                elsif ($elements[$n] =~/\bYES\b/gi){  #if word matches yes, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                    push (@{$matched_decrease{$elements[0]}},$n) unless grep{$_ == $n} @{$matched_decrease{$elements[0]}};
                    next if exists $dupmatd_serovar{$elements[0]};
                    push (@{$matd_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$matd_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmatd_serovar{$elements[0]}=1;
                }
            }
            $num_decrease_match=keys%matched_decrease; #counts the number of keys -gives the count of each matched decreased
        }
        elsif ($line =~/^(?=.*(\"\bS\b\")+)(?=.*(\"\bMAYBE\b\"))/gi) {    #sensitive to abx phenotypically but has resistance genes
            push (@{$maybes_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$maybes_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $j=0;$j<@elements;$j++){ #loop through each word
                next unless $elements[$j] =~/\bmaybe\b/gi; #if word matches maybe, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistance gene
                push (@{$maybe_sens{$elements[0]}},$j) unless grep{$_ == $j} @{$maybe_sens{$elements[0]}};
                next if exists $dupmays_serovar{$elements[0]};
                push (@{$maybesens_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$maybesens_serovar{$elements[2]." ST".$elements[3]}};
                $dupmays_serovar{$elements[0]}=1;
            }
            $num_pheno_maybe=keys %maybe_sens; #counts the number of keys -gives the count of each maybe pheno
        }
        elsif ($line =~ /^(?=.*(\"\bR\b\")+)(?=.*(\"\bMAYBE\b\"))/gi) {    #resistant to abx phenotypically but has no resistance genes
            push (@{$mayber_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$mayber_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $k=0;$k<@elements;$k++){ #loop through each word
                if ($elements[$k] =~/\bR\b/g){ #if word matches uppercase R, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_maybe_resis{$elements[0]}},$k) unless grep{$_ == $k} @{$anti_maybe_resis{$elements[0]}};
                }
                elsif ($elements[$k] =~/\bmaybe\b/gi){ #if word matches uppercase maybe, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistant gene.
                    push (@{$maybe_resis{$elements[0]}},$k) unless grep{$_ == $k} @{$maybe_resis{$elements[0]}};
                    next if exists $dupmayr_serovar{$elements[0]};
                    push (@{$mayberesis_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$mayberesis_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmayr_serovar{$elements[0]}=1;
                }
            }
            $num_geno_maybe=keys%maybe_resis; #counts the number of keys -gives the count of each maybe geno
        }
        elsif ($line =~ /^(?=.*(\"\bD\b\")+)(?=.*(\"\bMAYBE\b\"))/gi) {    #resistant to abx phenotypically but has no resistance genes
            push (@{$maybed_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$maybed_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $k=0;$k<@elements;$k++){ #loop through each word
                if ($elements[$k] =~/\bD\b/g){ #if word matches uppercase D, store column number in antibiotic hash=sample ID. removed /i in match to avoid matching with lowercase serovar r
                    push (@{$anti_maybe_decrease{$elements[0]}},$k) unless grep{$_ == $k} @{$anti_maybe_decrease{$elements[0]}};
                }
                elsif ($elements[$k] =~/\bmaybe\b/gi){ #if word matches uppercase mayve, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistant gene.
                    push (@{$maybe_decrease{$elements[0]}},$k) unless grep{$_ == $k} @{$maybe_decrease{$elements[0]}};
                    next if exists $dupmayd_serovar{$elements[0]};
                    push (@{$maybedecrease_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$maybedecrease_serovar{$elements[2]." ST".$elements[3]}};
                    $dupmayd_serovar{$elements[0]}=1;
                }
            }
            $num_decrease_maybe=keys%maybe_decrease; #counts the number of keys -gives the count of each maybe decreased
        }
        elsif ($line =~/^(?=.*(\"\bS\b\")+)(?!.*(\"\bYES\b\"))/gi) {#sensitive to abx phenotypically with no resistance genes, hash is serovar(column2&3)=sample ID(column1)
            push (@{$matchs_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$matchs_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            next if exists $dupmatp_serovar{$elements[0]};
            push (@{$matp_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$matp_serovar{$elements[2]." ST".$elements[3]}};
            $dupmatp_serovar{$elements[0]}=1;
        }
        elsif ($line =~ /^(?=.*(\"\bR\b\")+)(?!.*(\"\bYES\b\"))/gi) {    #resistant to abx phenotypically but has no resistance genes
            push (@{$misr_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$misr_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $m=0;$m<@elements;$m++){ #loop through each word
                next unless $elements[$m] =~/\bR\b/g; #if word matches uppercase R, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistant antibiotic. removed /i in match to avoid matching with lowercase serovar r
                push (@{$mismatched_geno{$elements[0]}},$m) unless grep{$_ == $m} @{$mismatched_geno{$elements[0]}};
                next if exists $dupmisp_serovar{$elements[0]};
                push (@{$misg_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$misg_serovar{$elements[2]." ST".$elements[3]}};
                $dupmisg_serovar{$elements[0]}=1;
            }
            $num_geno_discordant=keys%mismatched_geno; #counts the number of keys -gives the count of each mismatched geno
        }
        elsif ($line =~ /^(?=.*(\"\bD\b\")+)(?!.*(\"\bYES\b\"))/gi) {    #resistant to abx phenotypically but has no resistance genes
            push (@{$misd_rpattern{$elements[0]}},$elements[4]) unless grep{$_ eq $elements[4]} @{$misd_rpattern{$elements[0]}};#key=ID, value =overall resistance pattern
            for (my $m=0;$m<@elements;$m++){ #loop through each word
                next unless $elements[$m] =~/\bD\b/g; #if word matches uppercase D, continue. first hash is serovar(column2&3)=sample ID(column1), second hash is sample ID=resistant antibiotic. removed /i in match to avoid matching with lowercase serovar d
                push (@{$mismatched_decrease{$elements[0]}},$m) unless grep{$_ == $m} @{$mismatched_decrease{$elements[0]}};
                next if exists $dupmisd_serovar{$elements[0]};
                push (@{$misd_serovar{$elements[2]." ST".$elements[3]}},$elements[0]) unless grep{$_ eq $elements[0]} @{$misd_serovar{$elements[2]." ST".$elements[3]}};
                $dupmisd_serovar{$elements[0]}=1;
            }
            $num_decrease_discordant=keys%mismatched_decrease; #counts the number of keys -gives the count of each mismatched decreased
        }
        else { #all samples without phenotypic data
            push (@insufficient,$rows);
            foreach my $i_insuff(@insufficient) {
                my @id_insuff= split ",", $i_insuff;
                $hash_insuff{$_}++ foreach $id_insuff[0];
            }
            $num_insuff=keys%hash_insuff;
        }
    }
    else {
        print "cant parse\n";
    }
}

print OUTFILE "SUMMARY: for $input \n\n";

for my $hp_key(keys %matp_serovar){
    $num_pheno_match += @{$matp_serovar{$hp_key}};  #to get the total number of phenotypic matches as not done above
}
print OUTFILE "There are $num_pheno_match samples that are phenotypically sensitive and concurred genotypically.\n";
print OUTFILE "There are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n";
print OUTFILE "There are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n";
print OUTFILE "There are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n";
print OUTFILE "There are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
print OUTFILE "There are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n";
print OUTFILE "There are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n";
print OUTFILE "There are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n";
print OUTFILE "There are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n";
print OUTFILE "There are $num_geno_discordant samples that are  phenotypically resistant but genotypically sensitive.\n\n";
print OUTFILE "There are $num_decrease_match samples that are with reduced susceptibility and concurred genotypically at most once.\n";
print OUTFILE "There are $num_decrease2_match samples that are with reduced susceptibility and concurred genotypically at most twice.\n";
print OUTFILE "There are $num_decrease3_match samples that are with reduced susceptibility and concurred genotypically at least thrice.\n";
print OUTFILE "There are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n";
print OUTFILE "There are $num_decrease_discordant samples with reduced susceptibility but genotypically sensitive.\n\n";

print OUTFILE "#####"x30;

if ($opt_a){ #didnt want to create new hash with resistance genes: serovar because there is no unique id to use for key in hash -hence could have duplicates
    print OUTFILE "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    foreach my $mispheno_key (sort keys %mismatched_pheno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno{$mispheno_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
        }
        my ($rev_mispserovar) = grep { grep { $_ eq $mispheno_key } @{$misp_serovar{$_}} } keys %misp_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp{$reverse_key}},$rev_mispserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_mispkey (sort keys %rev_misp) { #for each serovar in resistance gene profile
        my $count_revmisp=@{$rev_misp{$rev_mispkey}}; #count total
        print OUTFILE "$rev_mispkey\tTotal = $count_revmisp\n";
        #print OUTFILE "$rev_mispkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp{$rev_mispkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    foreach my $mispheno2_key (sort keys %mismatched_pheno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno2{$mispheno2_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
        }
        my ($rev_misp2serovar) = grep { grep { $_ eq $mispheno2_key } @{$misp2_serovar{$_}} } keys %misp2_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp2{$reverse_key}},$rev_misp2serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misp2key (sort keys %rev_misp2) { #for each serovar in resistance gene profile
        my $count_revmisp2=@{$rev_misp2{$rev_misp2key}}; #count total
        print OUTFILE "$rev_misp2key\tTotal = $count_revmisp2\n";
        #print OUTFILE "$rev_misp2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp2{$rev_misp2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    foreach my $mispheno3_key (sort keys %mismatched_pheno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno3{$mispheno3_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
        }
        my ($rev_misp3serovar) = grep { grep { $_ eq $mispheno3_key } @{$misp3_serovar{$_}} } keys %misp3_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp3{$reverse_key}},$rev_misp3serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misp3key (sort keys %rev_misp3) { #for each serovar in resistance gene profile
        my $count_revmisp3=@{$rev_misp3{$rev_misp3key}}; #count total
        print OUTFILE "$rev_misp3key\tTotal = $count_revmisp3\n";
        #print OUTFILE "$rev_misp3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp3{$rev_misp3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    foreach my $matchgeno_key (sort keys %matched_geno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno{$matchgeno_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_matgserovar) = grep { grep { $_ eq $matchgeno_key } @{$matg_serovar{$_}} } keys %matg_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg{$reverse_key}},$rev_matgserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matgkey (sort keys %rev_matg) { #for each serovar in resistance gene profile
        my $count_revmatg=@{$rev_matg{$rev_matgkey}};
        print OUTFILE "$rev_matgkey\tTotal = $count_revmatg\n";
        #print OUTFILE "$rev_matgkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg{$rev_matgkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    foreach my $matchgeno2_key (sort keys %matched_geno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno2{$matchgeno2_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_matg2serovar) = grep { grep { $_ eq $matchgeno2_key } @{$matg2_serovar{$_}} } keys %matg2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg2{$reverse_key}},$rev_matg2serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg2key (sort keys %rev_matg2) { #for each serovar in resistance gene profile
        my $count_revmatg2=@{$rev_matg2{$rev_matg2key}};
        print OUTFILE "$rev_matg2key\tTotal = $count_revmatg2\n";
        #print OUTFILE "$rev_matg2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg2{$rev_matg2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $matchgeno3_key (sort keys %matched_geno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno3{$matchgeno3_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_matg3serovar) = grep { grep { $_ eq $matchgeno3_key } @{$matg3_serovar{$_}} } keys %matg3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg3{$reverse_key}},$rev_matg3serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg3key (sort keys %rev_matg3) { #for each serovar in resistance gene profile
        my $count_revmatg3=@{$rev_matg3{$rev_matg3key}};
        print OUTFILE "$rev_matg3key\tTotal = $count_revmatg3\n";
        #print OUTFILE "$rev_matg3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg3{$rev_matg3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples that are with reduced susceptibility and concurred genotypically at most once.\n\n";
    foreach my $matchdecrease_key (sort keys %matched_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease{$matchdecrease_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_matdserovar) = grep { grep { $_ eq $matchdecrease_key } @{$matd_serovar{$_}} } keys %matd_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd{$reverse_key}},$rev_matdserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matdkey (sort keys %rev_matd) { #for each serovar in resistance gene profile
        my $count_revmatd=@{$rev_matd{$rev_matdkey}};
        print OUTFILE "$rev_matdkey\tTotal = $count_revmatd\n";
        #print OUTFILE "$rev_matdkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd{$rev_matdkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples that are with reduced susceptibility and concurred genotypically at most twice.\n\n";
    foreach my $matchdecrease2_key (sort keys %matched_decrease2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease2{$matchdecrease2_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_matd2serovar) = grep { grep { $_ eq $matchdecrease2_key } @{$matd2_serovar{$_}} } keys %matd2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd2{$reverse_key}},$rev_matd2serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd2key (sort keys %rev_matd2) { #for each serovar in resistance gene profile
        my $count_revmatd2=@{$rev_matd2{$rev_matd2key}};
        print OUTFILE "$rev_matd2key\tTotal = $count_revmatd2\n";
        #print OUTFILE "$rev_matd2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd2{$rev_matd2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples that are with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    foreach my $matchdecrease3_key (sort keys %matched_decrease3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease3{$matchdecrease3_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_matd3serovar) = grep { grep { $_ eq $matchdecrease3_key } @{$matd3_serovar{$_}} } keys %matd3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd3{$reverse_key}},$rev_matd3serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd3key (sort keys %rev_matd3) { #for each serovar in resistance gene profile
        my $count_revmatd3=@{$rev_matd3{$rev_matd3key}};
        print OUTFILE "$rev_matd3key\tTotal = $count_revmatd3\n";
        #print OUTFILE "$rev_matd3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd3{$rev_matd3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    foreach my $maybesens_key (sort keys %maybe_sens) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_sens{$maybesens_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_maybesserovar) = grep { grep { $_ eq $maybesens_key } @{$maybesens_serovar{$_}} } keys %maybesens_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_maybes{$reverse_key}},$rev_maybesserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_maybeskey (sort keys %rev_maybes) { #for each serovar in resistance gene profile
        my $count_revmaybes=@{$rev_maybes{$rev_maybeskey}};
        print OUTFILE "$rev_maybeskey\tTotal = $count_revmaybes\n";
        #print OUTFILE "$rev_maybeskey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybes{$rev_maybeskey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $mayberesis_key (sort keys %maybe_resis) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_resis{$mayberesis_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_mayberserovar) = grep { grep { $_ eq $mayberesis_key } @{$mayberesis_serovar{$_}} } keys %mayberesis_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_mayber{$reverse_key}},$rev_mayberserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_mayberkey (sort keys %rev_mayber) { #for each serovar in resistance gene profile
        my $count_revmayber=@{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "$rev_mayberkey\tTotal = $count_revmayber\n";
        #print OUTFILE "$rev_mayberkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_mayber{$rev_mayberkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    foreach my $maybedecrease_key (sort keys %maybe_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_decrease{$maybedecrease_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        my ($rev_maybedserovar) = grep { grep { $_ eq $maybedecrease_key } @{$maybedecrease_serovar{$_}} } keys %maybedecrease_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_maybed{$reverse_key}},$rev_maybedserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_maybedkey (sort keys %rev_maybed) { #for each serovar in resistance gene profile
        my $count_revmaybed=@{$rev_maybed{$rev_maybedkey}};
        #print OUTFILE "$rev_maybedkey\tTotal = $count_revmaybed\n";
        print OUTFILE "$rev_maybedkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybed{$rev_maybedkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_b){
    print OUTFILE  "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    foreach my $mispheno_key (sort keys %mismatched_pheno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno{$mispheno_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_misp{$reverse_key}},$mispheno_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_mispkey (sort keys %rev_misp) {
        my $count_revmisp=@{$rev_misp{$rev_mispkey}}; #counts the total number in the array for each serovar
        print OUTFILE "$rev_mispkey:\t$count_revmisp\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    foreach my $mispheno2_key (sort keys %mismatched_pheno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno2{$mispheno2_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_misp2{$reverse_key}},$mispheno2_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_misp2key (sort keys %rev_misp2) {
        my $count_revmisp2=@{$rev_misp2{$rev_misp2key}}; #counts the total number in the array for each serovar
        print OUTFILE "$rev_misp2key:\t$count_revmisp2\n";
    }
    print OUTFILE  "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    foreach my $mispheno3_key (sort keys %mismatched_pheno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno3{$mispheno3_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_misp3{$reverse_key}},$mispheno3_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_misp3key (sort keys %rev_misp3) {
        my $count_revmisp3=@{$rev_misp3{$rev_misp3key}}; #counts the total number in the array for each serovar
        print OUTFILE "$rev_misp3key:\t$count_revmisp3\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    foreach my $matchgeno_key (sort keys %matched_geno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno{$matchgeno_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matg{$reverse_key}},$matchgeno_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_matgkey (sort keys %rev_matg) {
        my $count_revmatg=@{$rev_matg{$rev_matgkey}};
        print OUTFILE "$rev_matgkey:\t$count_revmatg\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    foreach my $matchgeno2_key (sort keys %matched_geno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno2{$matchgeno2_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matg2{$reverse_key}},$matchgeno2_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_matg2key (sort keys %rev_matg2) {
        my $count_revmatg2=@{$rev_matg2{$rev_matg2key}};
        print OUTFILE "$rev_matg2key:\t$count_revmatg2\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $matchgeno3_key (sort keys %matched_geno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno3{$matchgeno3_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matg3{$reverse_key}},$matchgeno3_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_matg3key (sort keys %rev_matg3) {
        my $count_revmatg3=@{$rev_matg3{$rev_matg3key}};
        print OUTFILE "$rev_matg3key:\t$count_revmatg3\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptibility and concurred genotypically at most once.\n\n";
    foreach my $matchdecrease_key (sort keys %matched_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease{$matchdecrease_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matd{$reverse_key}},$matchdecrease_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_matdkey (sort keys %rev_matd) {
        my $count_revmatd=@{$rev_matd{$rev_matdkey}};
        print OUTFILE "$rev_matdkey:\t$count_revmatd\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptibility and concurred genotypically at most twice.\n\n";
    foreach my $matchdecrease2_key (sort keys %matched_decrease2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease2{$matchdecrease2_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matd2{$reverse_key}},$matchdecrease2_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_matd2key (sort keys %rev_matd2) {
        my $count_revmatd2=@{$rev_matd2{$rev_matd2key}};
        print OUTFILE "$rev_matd2key:\t$count_revmatd2\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    foreach my $matchdecrease3_key (sort keys %matched_decrease3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease3{$matchdecrease3_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matd3{$reverse_key}},$matchdecrease3_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_matd3key (sort keys %rev_matd3) {
        my $count_revmatd3=@{$rev_matd3{$rev_matd3key}};
        print OUTFILE "$rev_matd3key:\t$count_revmatd3\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    foreach my $maybesens_key (sort keys %maybe_sens) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_sens{$maybesens_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_maybes{$reverse_key}},$maybesens_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_maybeskey (sort keys %rev_maybes) {
        my $count_revmaybes=@{$rev_maybes{$rev_maybeskey}};
        print OUTFILE "$rev_maybeskey:\t$count_revmaybes\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $mayberesis_key (sort keys %maybe_resis) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_resis{$mayberesis_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_mayber{$reverse_key}},$mayberesis_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_mayberkey (sort keys %rev_mayber) {
        my $count_revmayber=@{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "$rev_mayberkey:\t$count_revmayber\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    foreach my $maybedecrease_key (sort keys %maybe_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_decrease{$maybedecrease_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_maybed{$reverse_key}},$maybedecrease_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_maybedkey (sort keys %rev_maybed) {
        my $count_revmaybed=@{$rev_maybed{$rev_maybedkey}};
        print OUTFILE "$rev_maybedkey:\t$count_revmaybed\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_c) {
    print OUTFILE  "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    for my $misp (sort keys %misp_serovar) { # print OUTFILE count for each servoar in mismatch pheno
        my $count_misp=@{$misp_serovar{$misp}}; #counts the total number in the array for each serovar
        print OUTFILE "$misp:\t$count_misp\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    for my $misp2 (sort keys %misp2_serovar) { # print OUTFILE count for each servoar in mismatch pheno
        my $count_misp2=@{$misp2_serovar{$misp2}}; #counts the total number in the array for each serovar
        print OUTFILE "$misp2:\t$count_misp2\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    for my $misp3 (sort keys %misp3_serovar) { # print OUTFILE count for each servoar in mismatch pheno
        my $count_misp3=@{$misp3_serovar{$misp3}}; #counts the total number in the array for each serovar
        print OUTFILE "$misp3:\t$count_misp3\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_discordant samples that are  phenotypically resistant but genotypically sensitive.\n\n";
    for my $misg (sort keys %misg_serovar) { # print OUTFILE count for each servoar in mismatch geno
        my $count_misg=@{$misg_serovar{$misg}};
        print OUTFILE "$misg:\t$count_misg\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_discordant samples with reduced susceptibility but genotypically sensitive.\n\n";
    for my $misd (sort keys %misd_serovar) { # print OUTFILE count for each servoar in mismatch decreased susceptibility
        my $count_misd=@{$misd_serovar{$misd}};
        print OUTFILE "$misd:\t$count_misd\n";
    }
    #for my $hp_key(keys %matp_serovar){
    #    $num_pheno_match += @{$matp_serovar{$hp_key}};  #to get the total number of phenotypic matches as not done above
    #}
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_match samples that are phenotypically sensitive and concurred genotypically.\n\n";
    for my $matp (sort keys %matp_serovar) { # print OUTFILE count for each servoar in match pheno
        my $count_matp=@{$matp_serovar{$matp}};
        print OUTFILE "$matp:\t$count_matp\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    for my $matg (sort keys %matg_serovar) { # print OUTFILE count for each servoar in match geno
        my $count_matg=@{$matg_serovar{$matg}};
        print OUTFILE "$matg:\t$count_matg\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    for my $matg2 (sort keys %matg2_serovar) { # print OUTFILE count for each servoar in match geno
        my $count_matg2=@{$matg2_serovar{$matg2}};
        print OUTFILE "$matg2:\t$count_matg2\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    for my $matg3 (sort keys %matg3_serovar) { # print OUTFILE count for each servoar in match geno
        my $count_matg3=@{$matg3_serovar{$matg3}};
        print OUTFILE "$matg3:\t$count_matg3\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptibility and concurred genotypically at most once.\n\n";
    for my $matd (sort keys %matd_serovar) { # print OUTFILE count for each servoar in match geno
        my $count_matd=@{$matd_serovar{$matd}};
        print OUTFILE "$matd:\t$count_matd\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptibility and concurred genotypically at most twice.\n\n";
    for my $matd2 (sort keys %matd2_serovar) { # print OUTFILE count for each servoar in match geno
        my $count_matd2=@{$matd2_serovar{$matd2}};
        print OUTFILE "$matd2:\t$count_matd2\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    for my $matd3 (sort keys %matd3_serovar) { # print OUTFILE count for each servoar in match geno
        my $count_matd3=@{$matd3_serovar{$matd3}};
        print OUTFILE "$matd3:\t$count_matd3\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    for my $maybes (sort keys %maybesens_serovar) { # print OUTFILE count for each servoar in sensitive pheno, maybe geno
        my $count_maybes=@{$maybesens_serovar{$maybes}};
        print OUTFILE "$maybes:\t$count_maybes\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    for my $mayber (sort keys %mayberesis_serovar) { # print OUTFILE count for each servoar in resistant pheno, maybe geno
        my $count_mayber=@{$mayberesis_serovar{$mayber}};
        print OUTFILE "$mayber:\t$count_mayber\n";
    }
    print OUTFILE "\n";
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    for my $maybed (sort keys %maybedecrease_serovar) { # print OUTFILE count for each servoar in resistant pheno, maybe geno
        my $count_maybed=@{$maybedecrease_serovar{$maybed}};
        print OUTFILE "$maybed:\t$count_maybed\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_d){
    print OUTFILE  "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $mispheno_key (sort keys %mismatched_pheno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno{$mispheno_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_misp{$reverse_key}},$mispheno_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_mispkey (sort keys %rev_misp) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_mispkey\n";
        print OUTFILE "$_\n" foreach @{$rev_misp{$rev_mispkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $mispheno2_key (sort keys %mismatched_pheno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno2{$mispheno2_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_misp2{$reverse_key}},$mispheno2_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_misp2key (sort keys %rev_misp2) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_misp2key\n";
        print OUTFILE "$_\n" foreach @{$rev_misp2{$rev_misp2key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $mispheno3_key (sort keys %mismatched_pheno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno3{$mispheno3_key}}){
            $reverse_key.=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_misp3{$reverse_key}},$mispheno3_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_misp3key (sort keys %rev_misp3) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_misp3key\n";
        print OUTFILE "$_\n" foreach @{$rev_misp3{$rev_misp3key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $matchgeno_key (sort keys %matched_geno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno{$matchgeno_key}}){
            $reverse_key.=$column_index{$_}." "; #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_matg{$reverse_key}},$matchgeno_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_matgkey (sort keys %rev_matg) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_matgkey\n";
        print OUTFILE "$_\n" foreach @{$rev_matg{$rev_matgkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $matchgeno2_key (sort keys %matched_geno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno2{$matchgeno2_key}}){
            $reverse_key.=$column_index{$_}." "; #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_matg2{$reverse_key}},$matchgeno2_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_matg2key (sort keys %rev_matg2) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_matg2key\n";
        print OUTFILE "$_\n" foreach @{$rev_matg2{$rev_matg2key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $matchgeno3_key (sort keys %matched_geno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno3{$matchgeno3_key}}){
            $reverse_key.=$column_index{$_}." "; #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_matg3{$reverse_key}},$matchgeno3_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_matg3key (sort keys %rev_matg3) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_matg3key\n";
        print OUTFILE "$_\n" foreach @{$rev_matg3{$rev_matg3key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptibility and concurred genotypically at most once.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $matchdecrease_key (sort keys %matched_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease{$matchdecrease_key}}){
            $reverse_key.=$column_index{$_}." "; #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_matd{$reverse_key}},$matchdecrease_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_matdkey (sort keys %rev_matd) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_matdkey\n";
        print OUTFILE "$_\n" foreach @{$rev_matd{$rev_matdkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptibility and concurred genotypically at most twice.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $matchdecrease2_key (sort keys %matched_decrease2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease2{$matchdecrease2_key}}){
            $reverse_key.=$column_index{$_}." "; #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_matd2{$reverse_key}},$matchdecrease2_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_matd2key (sort keys %rev_matd2) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_matd2key\n";
        print OUTFILE "$_\n" foreach @{$rev_matd2{$rev_matd2key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $matchdecrease3_key (sort keys %matched_decrease3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease3{$matchdecrease3_key}}){
            $reverse_key.=$column_index{$_}." "; #string to store the resistance genes profile to put as key in reverse hash
        }
        push (@{$rev_matd3{$reverse_key}},$matchdecrease3_key); #reverse hash that has the resistance genes as the key and the sample as value
    }
    foreach my $rev_matd3key (sort keys %rev_matd3) { #print OUTFILEs the reverse hash and the sample IDs
        print OUTFILE "$rev_matd3key\n";
        print OUTFILE "$_\n" foreach @{$rev_matd3{$rev_matd3key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $maybesens_key (sort keys %maybe_sens) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_sens{$maybesens_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_maybes{$reverse_key}},$maybesens_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_maybeskey (sort keys %rev_maybes) {
        print OUTFILE "$rev_maybeskey\n";
        print OUTFILE "$_\n" foreach @{$rev_maybes{$rev_maybeskey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $mayberesis_key (sort keys %maybe_resis) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_resis{$mayberesis_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_mayber{$reverse_key}},$mayberesis_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_mayberkey (sort keys %rev_mayber) {
        print OUTFILE "$rev_mayberkey\n";
        print OUTFILE "$_\n" foreach @{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    print OUTFILE "Resistance genes and corresponding sample IDs\n";
    foreach my $maybedecrease_key (sort keys %maybe_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_decrease{$maybedecrease_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_maybed{$reverse_key}},$maybedecrease_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_maybedkey (sort keys %rev_maybed) {
        print OUTFILE "$rev_maybedkey\n";
        print OUTFILE "$_\n" foreach @{$rev_maybed{$rev_maybedkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_e){
    print OUTFILE  "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $mispheno_key (sort keys %mismatched_pheno) { #print OUTFILE sample ID and its resistance gene profile in mismatched pheno
        print OUTFILE "$mispheno_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$mismatched_pheno{$mispheno_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $mispheno2_key (sort keys %mismatched_pheno2) { #print OUTFILE sample ID and its resistance gene profile in mismatched pheno
        print OUTFILE "$mispheno2_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$mismatched_pheno2{$mispheno2_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $mispheno3_key (sort keys %mismatched_pheno3) { #print OUTFILE sample ID and its resistance gene profile in mismatched pheno
        print OUTFILE "$mispheno3_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$mismatched_pheno3{$mispheno3_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $matchgeno_key (sort keys %matched_geno) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$matchgeno_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$matched_geno{$matchgeno_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $matchgeno2_key (sort keys %matched_geno2) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$matchgeno2_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$matched_geno2{$matchgeno2_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $matchgeno3_key (sort keys %matched_geno3) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$matchgeno3_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$matched_geno3{$matchgeno3_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptibility and concurred genotypically at most once.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $matchdecrease_key (sort keys %matched_decrease) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$matchdecrease_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$matched_decrease{$matchdecrease_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptibility and concurred genotypically at most twice.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $matchdecrease2_key (sort keys %matched_decrease2) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$matchdecrease2_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$matched_decrease2{$matchdecrease2_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $matchdecrease3_key (sort keys %matched_decrease3) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$matchdecrease3_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$matched_decrease3{$matchdecrease3_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $maybesens_key (sort keys %maybe_sens) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$maybesens_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$maybe_sens{$maybesens_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $mayberesis_key (sort keys %maybe_resis) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$mayberesis_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$maybe_resis{$mayberesis_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    print OUTFILE "ID\tDiscordant_genes\n";
    foreach my $maybedecrease_key (sort keys %maybe_decrease) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        print OUTFILE "$maybedecrease_key\t";
        print OUTFILE "$column_index{$_}\t" foreach @{$maybe_decrease{$maybedecrease_key}}; #for each value (column number), looks through keys in %column hash for respective value
        print OUTFILE "\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_f){
    print OUTFILE  "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $misp_key (sort keys %misp_serovar) { #for each serovar in mismatch pheno, print OUTFILE the sample ID
        print OUTFILE  "$misp_key\n";
        print OUTFILE "$_\n" foreach @{$misp_serovar{$misp_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $misp2_key (sort keys %misp2_serovar) { #for each serovar in mismatch pheno, print OUTFILE the sample ID
        print OUTFILE  "$misp2_key\n";
        print OUTFILE "$_\n" foreach @{$misp2_serovar{$misp2_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $misp3_key (sort keys %misp3_serovar) { #for each serovar in mismatch pheno, print OUTFILE the sample ID
        print OUTFILE  "$misp3_key\n";
        print OUTFILE "$_\n" foreach @{$misp3_serovar{$misp3_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_discordant samples that are  phenotypically resistant but genotypically sensitive.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $misg_key (sort keys %misg_serovar) { #for each serovar in mismatch geno, print OUTFILE the sample ID
        print OUTFILE  "$misg_key\n";
        print OUTFILE "$_\n" foreach @{$misg_serovar{$misg_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_discordant samples with reduced susceptibility but genotypically sensitive.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $misd_key (sort keys %misd_serovar) { #for each serovar in mismatch geno, print OUTFILE the sample ID
        print OUTFILE  "$misd_key\n";
        print OUTFILE "$_\n" foreach @{$misd_serovar{$misd_key}};
        print OUTFILE "\n";
    }
    #for my $matp_key(keys %matp_serovar){
    #    $num_pheno_match += @{$matp_serovar{$matp_key}};  #counts the total number in the array for each serovar -gives the count of matched pheno
    #}
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_match samples that are phenotypically sensitive and concurred genotypically.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $matp_serovarkeys (sort keys %matp_serovar){ #for each serovar in match pheno, print OUTFILE the sample ID
        print OUTFILE "$matp_serovarkeys\n";
        print OUTFILE "$_\n" foreach @{$matp_serovar{$matp_serovarkeys}};
        print OUTFILE  "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $matg_key (sort keys %matg_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$matg_key\n";
        print OUTFILE "$_\n" foreach @{$matg_serovar{$matg_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $matg2_key (sort keys %matg2_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$matg2_key\n";
        print OUTFILE "$_\n" foreach @{$matg2_serovar{$matg2_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $matg3_key (sort keys %matg3_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$matg3_key\n";
        print OUTFILE "$_\n" foreach @{$matg3_serovar{$matg3_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptibility and concurred genotypically at most once.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $matd_key (sort keys %matd_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$matd_key\n";
        print OUTFILE "$_\n" foreach @{$matd_serovar{$matd_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptibility and concurred genotypically at most twice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $matd2_key (sort keys %matd2_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$matd2_key\n";
        print OUTFILE "$_\n" foreach @{$matd2_serovar{$matd2_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $matd3_key (sort keys %matd3_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$matd3_key\n";
        print OUTFILE "$_\n" foreach @{$matd3_serovar{$matd3_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are phenotypically sensitive and maybe resistant genotypically.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $maybes_key (sort keys %maybesens_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$maybes_key\n";
        print OUTFILE "$_\n" foreach @{$maybesens_serovar{$maybes_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $mayber_key (sort keys %mayberesis_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$mayber_key\n";
        print OUTFILE "$_\n" foreach @{$mayberesis_serovar{$mayber_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs\n";
    foreach my $maybed_key (sort keys %maybedecrease_serovar) { #for each serovar in match geno, print OUTFILE the sample ID
        print OUTFILE  "$maybed_key\n";
        print OUTFILE "$_\n" foreach @{$maybedecrease_serovar{$maybed_key}};
        print OUTFILE "\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_g){
    print OUTFILE "\nThere are $num_geno_discordant samples that are  phenotypically resistant but genotypically sensitive.\n\n";
    foreach my $misgeno_key (sort keys %mismatched_geno) {
        my $reverse_key=();
        foreach (@{$mismatched_geno{$misgeno_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_misg{$reverse_key}},$misgeno_key);
    }
    foreach my $rev_misgkey (sort keys %rev_misg) {
        my $count_rev_misgkey=@{$rev_misg{$rev_misgkey}};
        print OUTFILE "$rev_misgkey\n";
        print OUTFILE "$_\n" foreach @{$rev_misg{$rev_misgkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    foreach my $anti_matgkey (sort keys %matched_geno) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_matched_geno{$anti_matgkey}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matg{$reverse_key}},$anti_matgkey);
    }
    foreach my $rev_matgkey (sort keys %rev_matg) {
        my $count_rev_matgkey=@{$rev_matg{$rev_matgkey}};
        print OUTFILE "$rev_matgkey\n";
        print OUTFILE "$_\n" foreach @{$rev_matg{$rev_matgkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    foreach my $anti_matg2key (sort keys %matched_geno2) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_matched_geno2{$anti_matg2key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matg2{$reverse_key}},$anti_matg2key);
    }
    foreach my $rev_matg2key (sort keys %rev_matg2) {
        my $count_rev_matg2key=@{$rev_matg2{$rev_matg2key}};
        print OUTFILE "$rev_matg2key\n";
        print OUTFILE "$_\n" foreach @{$rev_matg2{$rev_matg2key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $anti_matg3key (sort keys %matched_geno3) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_matched_geno3{$anti_matg3key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matg3{$reverse_key}},$anti_matg3key);
    }
    foreach my $rev_matg3key (sort keys %rev_matg3) {
        my $count_rev_matg3key=@{$rev_matg3{$rev_matg3key}};
        print OUTFILE "$rev_matg3key\n";
        print OUTFILE "$_\n" foreach @{$rev_matg3{$rev_matg3key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $anti_mayberkey (sort keys %maybe_resis) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_maybe_resis{$anti_mayberkey}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_mayber{$reverse_key}}, $anti_mayberkey);
    }
    foreach my $rev_mayberkey(sort keys %rev_mayber) {
        my $count_rev_mayberkey=@{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "$rev_mayberkey\n";
        print OUTFILE "$_\n" foreach @{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_discordant samples with reduced susceptiblity but genotypically sensitive.\n\n";
    foreach my $misdecrease_key (sort keys %mismatched_decrease) {
        my $reverse_key=();
        foreach (@{$mismatched_decrease{$misdecrease_key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_misd{$reverse_key}},$misdecrease_key);
    }
    foreach my $rev_misdkey (sort keys %rev_misd) {
        my $count_rev_misdkey=@{$rev_misd{$rev_misdkey}};
        print OUTFILE "$rev_misdkey\n";
        print OUTFILE "$_\n" foreach @{$rev_misd{$rev_misdkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptiblity and concurred genotypically at most once.\n\n";
    foreach my $anti_matdkey (sort keys %matched_decrease) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_matched_decrease{$anti_matdkey}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matd{$reverse_key}},$anti_matdkey);
    }
    foreach my $rev_matdkey (sort keys %rev_matd) {
        my $count_rev_matdkey=@{$rev_matd{$rev_matdkey}};
        print OUTFILE "$rev_matdkey\n";
        print OUTFILE "$_\n" foreach @{$rev_matd{$rev_matdkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptiblity and concurred genotypically at most twice.\n\n";
    foreach my $anti_matd2key (sort keys %matched_decrease2) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_matched_decrease2{$anti_matd2key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matd2{$reverse_key}},$anti_matd2key);
    }
    foreach my $rev_matd2key (sort keys %rev_matd2) {
        my $count_rev_matd2key=@{$rev_matd2{$rev_matd2key}};
        print OUTFILE "$rev_matd2key\n";
        print OUTFILE "$_\n" foreach @{$rev_matd2{$rev_matd2key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptiblity and concurred genotypically at least thrice.\n\n";
    foreach my $anti_matd3key (sort keys %matched_decrease3) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_matched_decrease3{$anti_matd3key}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_matd3{$reverse_key}},$anti_matd3key);
    }
    foreach my $rev_matd3key (sort keys %rev_matd3) {
        my $count_rev_matd3key=@{$rev_matd3{$rev_matd3key}};
        print OUTFILE "$rev_matd3key\n";
        print OUTFILE "$_\n" foreach @{$rev_matd3{$rev_matd3key}};
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptiblity and maybe resistant genotypically.\n\n";
    foreach my $anti_maybedkey (sort keys %maybe_decrease) { #print OUTFILE sample ID and its resistance gene profile in matched geno
        my $reverse_key=();
        foreach (@{$anti_maybe_decrease{$anti_maybedkey}}){
            $reverse_key.=$column_index{$_}." ";
        }
        push (@{$rev_maybed{$reverse_key}}, $anti_maybedkey);
    }
    foreach my $rev_maybedkey(sort keys %rev_maybed) {
        my $count_rev_maybedkey=@{$rev_maybed{$rev_maybedkey}};
        print OUTFILE "$rev_maybedkey\n";
        print OUTFILE "$_\n" foreach @{$rev_maybed{$rev_maybedkey}};
        print OUTFILE "\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_i){ #didnt want to create new hash with resistance genes: serovar because there is no unique id to use for key in hash -hence could have duplicates
    print OUTFILE "\nThere are $num_pheno_match samples that are phenotypically sensitive and concurred genotypically.\n\n";
    foreach my $matchs_key (sort keys %matchs_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchs_rpattern{$matchs_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matpserovar) = grep { grep { $_ eq $matchs_key } @{$matp_serovar{$_}} } keys %matp_serovar; # finds the value (sample ID) in matp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matp{$reverse_key}},$rev_matpserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matpkey (sort keys %rev_matp) { #for each serovar in resistance gene profile
        my $count_revmatp=@{$rev_matp{$rev_matpkey}}; #count total
        print OUTFILE "$rev_matpkey\tTotal = $count_revmatp\n";
        #print OUTFILE "$rev_matpkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matp{$rev_matpkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    foreach my $matchr_key (sort keys %matchr_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr_rpattern{$matchr_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matgserovar) = grep { grep { $_ eq $matchr_key } @{$matg_serovar{$_}} } keys %matg_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg{$reverse_key}},$rev_matgserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matgkey (sort keys %rev_matg) { #for each serovar in resistance gene profile
        my $count_revmatg=@{$rev_matg{$rev_matgkey}}; #count total
        print OUTFILE "$rev_matgkey\tTotal = $count_revmatg\n";
        #print OUTFILE "$rev_matgkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg{$rev_matgkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    foreach my $matchr2_key (sort keys %matchr2_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr2_rpattern{$matchr2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matg2serovar) = grep { grep { $_ eq $matchr2_key } @{$matg2_serovar{$_}} } keys %matg2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg2{$reverse_key}},$rev_matg2serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg2key (sort keys %rev_matg2) { #for each serovar in resistance gene profile
        my $count_revmatg2=@{$rev_matg2{$rev_matg2key}}; #count total
        print OUTFILE "$rev_matg2key\tTotal = $count_revmatg2\n";
        #print OUTFILE "$rev_matg2key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg2{$rev_matg2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $matchr3_key (sort keys %matchr3_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr3_rpattern{$matchr3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matg3serovar) = grep { grep { $_ eq $matchr3_key } @{$matg3_serovar{$_}} } keys %matg3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg3{$reverse_key}},$rev_matg3serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg3key (sort keys %rev_matg3) { #for each serovar in resistance gene profile
        my $count_revmatg3=@{$rev_matg3{$rev_matg3key}}; #count total
        print OUTFILE "$rev_matg3key\tTotal = $count_revmatg3\n";
        #print OUTFILE "$rev_matg3key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg3{$rev_matg3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples that are with reduced susceptibility and concurred genotypically at most once.\n\n";
    foreach my $matchd_key (sort keys %matchd_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd_rpattern{$matchd_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matdserovar) = grep { grep { $_ eq $matchd_key } @{$matd_serovar{$_}} } keys %matd_serovar; # finds the value (sample ID) in matd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd{$reverse_key}},$rev_matdserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matdkey (sort keys %rev_matd) { #for each serovar in resistance gene profile
        my $count_revmatd=@{$rev_matd{$rev_matdkey}}; #count total
        print OUTFILE "$rev_matdkey\tTotal = $count_revmatd\n";
        #print OUTFILE "$rev_matdkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd{$rev_matdkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples that are with reduced susceptibility and concurred genotypically at most twice.\n\n";
    foreach my $matchd2_key (sort keys %matchd2_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd2_rpattern{$matchd2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matd2serovar) = grep { grep { $_ eq $matchd2_key } @{$matd2_serovar{$_}} } keys %matd2_serovar; # finds the value (sample ID) in matd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd2{$reverse_key}},$rev_matd2serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd2key (sort keys %rev_matd2) { #for each serovar in resistance gene profile
        my $count_revmatd2=@{$rev_matd2{$rev_matd2key}}; #count total
        print OUTFILE "$rev_matd2key\tTotal = $count_revmatd2\n";
        #print OUTFILE "$rev_matd2key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd2{$rev_matd2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples that are with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    foreach my $matchd3_key (sort keys %matchd3_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd3_rpattern{$matchd3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matd3serovar) = grep { grep { $_ eq $matchd3_key } @{$matd3_serovar{$_}} } keys %matd3_serovar; # finds the value (sample ID) in matd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd3{$reverse_key}},$rev_matd3serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd3key (sort keys %rev_matd3) { #for each serovar in resistance gene profile
        my $count_revmatd3=@{$rev_matd3{$rev_matd3key}}; #count total
        print OUTFILE "$rev_matd3key\tTotal = $count_revmatd3\n";
        #print OUTFILE "$rev_matd3key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd3{$rev_matd3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    foreach my $miss_key (sort keys %miss_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss_rpattern{$miss_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_mispserovar) = grep { grep { $_ eq $miss_key } @{$misp_serovar{$_}} } keys %misp_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp{$reverse_key}},$rev_mispserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_mispkey (sort keys %rev_misp) { #for each serovar in resistance gene profile
        my $count_revmisp=@{$rev_misp{$rev_mispkey}}; #count total
        print OUTFILE "$rev_mispkey\tTotal = $count_revmisp\n";
        #print OUTFILE "$rev_mispkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp{$rev_mispkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    foreach my $miss2_key (sort keys %miss2_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss2_rpattern{$miss2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misp2serovar) = grep { grep { $_ eq $miss2_key } @{$misp2_serovar{$_}} } keys %misp2_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp2{$reverse_key}},$rev_misp2serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misp2key (sort keys %rev_misp2) { #for each serovar in resistance gene profile
        my $count_revmisp2=@{$rev_misp2{$rev_misp2key}}; #count total
        print OUTFILE "$rev_misp2key\tTotal = $count_revmisp2\n";
        #print OUTFILE "$rev_misp2key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp2{$rev_misp2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    foreach my $miss3_key (sort keys %miss3_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss3_rpattern{$miss3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misp3serovar) = grep { grep { $_ eq $miss3_key } @{$misp3_serovar{$_}} } keys %misp3_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp3{$reverse_key}},$rev_misp3serovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misp3key (sort keys %rev_misp3) { #for each serovar in resistance gene profile
        my $count_revmisp3=@{$rev_misp3{$rev_misp3key}}; #count total
        print OUTFILE "$rev_misp3key\tTotal = $count_revmisp3\n";
        #print OUTFILE "$rev_misp3key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp3{$rev_misp3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_discordant samples that are phenotypically resistant but genotypically sensitive.\n\n";
    foreach my $misr_key (sort keys %misr_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$misr_rpattern{$misr_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misgserovar) = grep { grep { $_ eq $misr_key } @{$misg_serovar{$_}} } keys %misg_serovar; # finds the value (sample ID) in misg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misg{$reverse_key}},$rev_misgserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misgkey (sort keys %rev_misg) { #for each serovar in resistance gene profile
        my $count_revmisg=@{$rev_misg{$rev_misgkey}}; #count total
        print OUTFILE "$rev_misgkey\tTotal = $count_revmisg\n";
        #print OUTFILE "$rev_misgkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misg{$rev_misgkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_discordant samples that are phenotypically sensitive but maybe genotypically resistant.\n\n";
    foreach my $misd_key (sort keys %misd_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$misd_rpattern{$misd_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misdserovar) = grep { grep { $_ eq $misd_key } @{$misd_serovar{$_}} } keys %misd_serovar; # finds the value (sample ID) in misd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misd{$reverse_key}},$rev_misdserovar); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misdkey (sort keys %rev_misd) { #for each serovar in resistance gene profile
        my $count_revmisd=@{$rev_misd{$rev_misdkey}}; #count total
        print OUTFILE "$rev_misdkey\tTotal = $count_revmisd\n";
        #print OUTFILE "$rev_misdkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misd{$rev_misdkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    foreach my $maybes_key (sort keys %maybes_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$maybes_rpattern{$maybes_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_maysserovar) = grep { grep { $_ eq $maybes_key } @{$maybesens_serovar{$_}} } keys %maybesens_serovar; # finds the value (sample ID) in maybesens_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_maybes{$reverse_key}},$rev_maysserovar); #push each serovar to the corresponding sample ID for the overall resistance profile
    }
    for my $rev_maybeskey (sort keys %rev_maybes) { #for each serovar in resistance gene profile
        my $count_revmaybes=@{$rev_maybes{$rev_maybeskey}}; #count total
        print OUTFILE "$rev_maybeskey\tTotal = $count_revmaybes\n";
        #print OUTFILE "$rev_maybeskey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybes{$rev_maybeskey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $mayber_key (sort keys %mayber_rpattern) { #goes through the overallresistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$mayber_rpattern{$mayber_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_mayrserovar) = grep { grep { $_ eq $mayber_key } @{$mayberesis_serovar{$_}} } keys %mayberesis_serovar; # finds the value (sample ID) in mayberesis_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_mayber{$reverse_key}},$rev_mayrserovar); #push each serovar to the corresponding sample ID for that overall resistance profile
    }
    for my $rev_mayberkey (sort keys %rev_mayber) { #for each serovar in resistance gene profile
        my $count_revmayber=@{$rev_mayber{$rev_mayberkey}}; #count total
        print OUTFILE "$rev_mayberkey\tTotal = $count_revmayber\n";
        #print OUTFILE "$rev_mayberkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_mayber{$rev_mayberkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    foreach my $maybed_key (sort keys %maybed_rpattern) { #goes through the overallresistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$maybed_rpattern{$maybed_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_maydserovar) = grep { grep { $_ eq $maybed_key } @{$maybedecrease_serovar{$_}} } keys %maybedecrease_serovar; # finds the value (sample ID) in maybedecrease_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_maybed{$reverse_key}},$rev_maydserovar); #push each serovar to the corresponding sample ID for that overall resistance profile
    }
    for my $rev_maybedkey (sort keys %rev_maybed) { #for each serovar in resistance gene profile
        my $count_revmaybed=@{$rev_maybed{$rev_maybedkey}}; #count total
        print OUTFILE "$rev_maybedkey\tTotal = $count_revmaybed\n";
        #print OUTFILE "$rev_maybedkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybed{$rev_maybedkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
}

if ($opt_j){
    print OUTFILE "\nThere are $num_pheno_match samples that are phenotypically sensitive and concurred genotypically.\n\n";
    foreach my $matchs_key (sort keys %matchs_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchs_rpattern{$matchs_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_matp{$reverse_key}},$matchs_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_matpkey (sort keys %rev_matp) {
        my $count_revmatp=@{$rev_matp{$rev_matpkey}}; #counts the total number in the array for each serovar
        print OUTFILE "$rev_matpkey:\t$count_revmatp\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    foreach my $matchr_key (sort keys %matchr_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr_rpattern{$matchr_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_matg{$reverse_key}},$matchr_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matgkey (sort keys %rev_matg) {
        my $count_revmatg=@{$rev_matg{$rev_matgkey}};
        print OUTFILE "$rev_matgkey:\t$count_revmatg\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    foreach my $matchr2_key (sort keys %matchr2_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr2_rpattern{$matchr2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_matg2{$reverse_key}},$matchr2_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg2key (sort keys %rev_matg2) {
        my $count_revmatg2=@{$rev_matg2{$rev_matg2key}};
        print OUTFILE "$rev_matg2key:\t$count_revmatg2\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $matchr3_key (sort keys %matchr3_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr3_rpattern{$matchr3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_matg3{$reverse_key}},$matchr3_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg3key (sort keys %rev_matg3) {
        my $count_revmatg3=@{$rev_matg3{$rev_matg3key}};
        print OUTFILE "$rev_matg3key:\t$count_revmatg3\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptibility and concurred genotypically at most once.\n\n";
    foreach my $matchd_key (sort keys %matchd_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd_rpattern{$matchd_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_matd{$reverse_key}},$matchd_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matdkey (sort keys %rev_matd) {
        my $count_revmatd=@{$rev_matd{$rev_matdkey}};
        print OUTFILE "$rev_matdkey:\t$count_revmatd\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptibility and concurred genotypically at most twice.\n\n";
    foreach my $matchd2_key (sort keys %matchd2_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd2_rpattern{$matchd2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_matd2{$reverse_key}},$matchd2_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd2key (sort keys %rev_matd2) {
        my $count_revmatd2=@{$rev_matd2{$rev_matd2key}};
        print OUTFILE "$rev_matd2key:\t$count_revmatd2\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    foreach my $matchd3_key (sort keys %matchd3_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd3_rpattern{$matchd3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_matd3{$reverse_key}},$matchd3_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd3key (sort keys %rev_matd3) {
        my $count_revmatd3=@{$rev_matd3{$rev_matd3key}};
        print OUTFILE "$rev_matd3key:\t$count_revmatd3\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    foreach my $miss_key (sort keys %miss_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss_rpattern{$miss_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_misp{$reverse_key}},$miss_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_mispkey (sort keys %rev_misp) {
        my $count_revmisp=@{$rev_misp{$rev_mispkey}}; #counts the total number in the array for each serovar
        print OUTFILE "$rev_mispkey:\t$count_revmisp\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    foreach my $miss2_key (sort keys %miss2_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss2_rpattern{$miss2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_misp2{$reverse_key}},$miss2_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_misp2key (sort keys %rev_misp2) {
        my $count_revmisp2=@{$rev_misp2{$rev_misp2key}}; #counts the total number in the array for each serovar
        print OUTFILE "$rev_misp2key:\t$count_revmisp2\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    foreach my $miss3_key (sort keys %miss3_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss3_rpattern{$miss3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_misp3{$reverse_key}},$miss3_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_misp3key (sort keys %rev_misp3) {
        my $count_revmisp3=@{$rev_misp3{$rev_misp3key}}; #counts the total number in the array for each serovar
        print OUTFILE "$rev_misp3key:\t$count_revmisp3\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_discordant samples that are phenotypically resistant but genotypically sensitive.\n\n";
    foreach my $misr_key (sort keys %misr_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$misr_rpattern{$misr_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_misg{$reverse_key}},$misr_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misgkey (sort keys %rev_misg) {
        my $count_revmisg=@{$rev_misg{$rev_misgkey}};
        print OUTFILE "$rev_misgkey:\t$count_revmisg\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_discordant samples with reduced susceptibility and concurred genotypically.\n\n";
    foreach my $misd_key (sort keys %misd_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$misd_rpattern{$misd_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_misd{$reverse_key}},$misd_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misdkey (sort keys %rev_misd) {
        my $count_revmisd=@{$rev_misd{$rev_misdkey}};
        print OUTFILE "$rev_misdkey:\t$count_revmisd\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    foreach my $maybes_key (sort keys %maybes_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$maybes_rpattern{$maybes_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_maybes{$reverse_key}},$maybes_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_maybeskey (sort keys %rev_maybes) {
        my $count_revmaybes=@{$rev_maybes{$rev_maybeskey}};
        print OUTFILE "$rev_maybeskey:\t$count_revmaybes\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $mayber_key (sort keys %mayber_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$mayber_rpattern{$mayber_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_mayber{$reverse_key}},$mayber_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_mayberkey (sort keys %rev_mayber) {
        my $count_revmayber=@{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "$rev_mayberkey:\t$count_revmayber\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    foreach my $maybed_key (sort keys %maybed_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$maybed_rpattern{$maybed_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        push (@{$rev_maybed{$reverse_key}},$maybed_key); #string to store the resistance genes profile to put as key in reverse hash
    }
    for my $rev_maybedkey (sort keys %rev_maybed) {
        my $count_revmaybed=@{$rev_maybed{$rev_maybedkey}};
        print OUTFILE "$rev_maybedkey:\t$count_revmaybed\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_k){
    print OUTFILE  "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $misp_key (sort keys %misp_serovar) { #for each serovar in mismatch pheno where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$misp_key\n";
        foreach my $misp_serovar_key (@{$misp_serovar{$misp_key}}){ #for each sample ID in mismatch pheno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$misp_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$miss_rpattern{$misp_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $misp2_key (sort keys %misp2_serovar) { #for each serovar in mismatch pheno where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$misp2_key\n";
        foreach my $misp2_serovar_key (@{$misp2_serovar{$misp2_key}}){ #for each sample ID in mismatch pheno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$misp2_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$miss2_rpattern{$misp2_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $misp3_key (sort keys %misp3_serovar) { #for each serovar in mismatch pheno where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$misp3_key\n";
        foreach my $misp3_serovar_key (@{$misp3_serovar{$misp3_key}}){ #for each sample ID in mismatch pheno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$misp3_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$miss3_rpattern{$misp3_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_discordant samples that are  phenotypically resistant but genotypically sensitive.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $misg_key (sort keys %misg_serovar) { #for each serovar in mismatch geno where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$misg_key\n";
        foreach my $misg_serovar_key (@{$misg_serovar{$misg_key}}){ #for each sample ID in mismatch geno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$misg_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$misr_rpattern{$misg_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE  "\nThere are $num_decrease_discordant samples that are phenotypically sensitive but genotypically resistant.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $misd_key (sort keys %misd_serovar) { #for each serovar in mismatch decrease where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$misd_key\n";
        foreach my $misd_serovar_key (@{$misd_serovar{$misd_key}}){#for each sample ID in mismatch decrease, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$misd_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$misd_rpattern{$misd_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "There are $num_pheno_match samples that are phenotypically sensitive and concurred genotypically.\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $matp_key (sort keys %matp_serovar) { #for each serovar in match pheno where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$matp_key\n";
        foreach my $matp_serovar_key (@{$matp_serovar{$matp_key}}){#for each sample ID in match pheno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$matp_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$matchs_rpattern{$matp_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $matg_key (sort keys %matg_serovar) { #for each serovar in match geno where serovar(key)=sample ID(value), print serovar
        print OUTFILE "$matg_key\n";
        foreach my $matg_serovar_key (@{$matg_serovar{$matg_key}}){#for each sample ID in match geno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$matg_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$matchr_rpattern{$matg_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $matg2_key (sort keys %matg2_serovar) { #for each serovar in match geno where serovar(key)=sample ID(value), print serovar
        print OUTFILE "$matg2_key\n";
        foreach my $matg2_serovar_key (@{$matg2_serovar{$matg2_key}}){#for each sample ID in match geno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$matg2_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$matchr2_rpattern{$matg2_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $matg3_key (sort keys %matg3_serovar) { #for each serovar in match geno where serovar(key)=sample ID(value), print serovar
        print OUTFILE "$matg3_key\n";
        foreach my $matg3_serovar_key (@{$matg3_serovar{$matg3_key}}){#for each sample ID in match geno, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$matg3_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$matchr3_rpattern{$matg3_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples with reduced susceptibility and concurred genotypically at most once.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $matd_key (sort keys %matd_serovar) { #for each serovar in match decrease where serovar(key)=sample ID(value), print serovar
        print OUTFILE "$matd_key\n";
        foreach my $matd_serovar_key (@{$matd_serovar{$matd_key}}){#for each sample ID in match decrease, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$matd_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$matchd_rpattern{$matd_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples with reduced susceptibility and concurred genotypically at most twice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $matd2_key (sort keys %matd2_serovar) { #for each serovar in match decrease where serovar(key)=sample ID(value), print serovar
        print OUTFILE "$matd2_key\n";
        foreach my $matd2_serovar_key (@{$matd2_serovar{$matd2_key}}){#for each sample ID in match decrease, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$matd2_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$matchd2_rpattern{$matd2_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $matd3_key (sort keys %matd3_serovar) { #for each serovar in match decrease where serovar(key)=sample ID(value), print serovar
        print OUTFILE "$matd3_key\n";
        foreach my $matd3_serovar_key (@{$matd3_serovar{$matd3_key}}){#for each sample ID in match decrease, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$matd3_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$matchd3_rpattern{$matd3_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are phenotypically sensitive and maybe resistant genotypically.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $maybesens_key (sort keys %maybesens_serovar) { #for each serovar in maybe sens where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$maybesens_key\n";
        foreach my $maybesens_serovar_key (@{$maybesens_serovar{$maybesens_key}}){ #for each sample ID in match sens, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$maybesens_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$maybes_rpattern{$maybesens_serovar_key}};
        }
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $mayberesis_key (sort keys %mayberesis_serovar) { #for each serovar in maybe resis where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$mayberesis_key\n";
        foreach my $mayberesis_serovar_key (@{$mayberesis_serovar{$mayberesis_key}}){#for each sample ID in maybe resis, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$mayberesis_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$mayber_rpattern{$mayberesis_serovar_key}};
        }
        print OUTFILE "\n";
    }
    
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    print OUTFILE "Serovar ST and corresponding Sample IDs and resistance gene profile\n";
    foreach my $maybedecrease_key (sort keys %maybedecrease_serovar) { #for each serovar in maybe decrease where serovar(key)=sample ID(value), print serovar
        print OUTFILE  "$maybedecrease_key\n";
        foreach my $maybedecrease_serovar_key (@{$maybedecrease_serovar{$maybedecrease_key}}){ #for each sample ID in maybe decrease, go through the corresponding overall resistance profile i.e SAMPLE_ID(key)=overall resistance profile(value)
            print OUTFILE "$maybedecrease_serovar_key:\t";
            print OUTFILE "$_\n" foreach @{$maybed_rpattern{$maybedecrease_serovar_key}};
        }
        print OUTFILE "\n";
    }
}

if ($opt_m){ #didnt want to create new hash with resistance genes: serovar because there is no unique id to use for key in hash -hence could have duplicates
    print OUTFILE "\nThere are $num_pheno_match samples that are phenotypically sensitive and concurred genotypically.\n\n";
    foreach my $matchs_key (sort keys %matchs_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchs_rpattern{$matchs_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matpserovar) = grep { grep { $_ eq $matchs_key } @{$matp_serovar{$_}} } keys %matp_serovar; # finds the value (sample ID) in matp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matp{$rev_matpserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matpkey (sort keys %rev_matp) { #for each serovar in resistance gene profile
        my $count_revmatp=@{$rev_matp{$rev_matpkey}}; #count total
        print OUTFILE "$rev_matpkey\tTotal = $count_revmatp\n";
        #print OUTFILE "$rev_matpkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matp{$rev_matpkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at most once.\n\n";
    foreach my $matchr_key (sort keys %matchr_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr_rpattern{$matchr_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matgserovar) = grep { grep { $_ eq $matchr_key } @{$matg_serovar{$_}} } keys %matg_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg{$rev_matgserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matgkey (sort keys %rev_matg) { #for each serovar in resistance gene profile
        my $count_revmatg=@{$rev_matg{$rev_matgkey}}; #count total
        print OUTFILE "$rev_matgkey\tTotal = $count_revmatg\n";
        #print OUTFILE "$rev_matgkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg{$rev_matgkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at most twice.\n\n";
    foreach my $matchr2_key (sort keys %matchr2_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr2_rpattern{$matchr2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matg2serovar) = grep { grep { $_ eq $matchr2_key } @{$matg2_serovar{$_}} } keys %matg2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg2{$rev_matg2serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg2key (sort keys %rev_matg2) { #for each serovar in resistance gene profile
        my $count_revmatg2=@{$rev_matg2{$rev_matg2key}}; #count total
        print OUTFILE "$rev_matg2key\tTotal = $count_revmatg2\n";
        #print OUTFILE "$rev_matg2key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg2{$rev_matg2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $matchr3_key (sort keys %matchr3_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchr3_rpattern{$matchr3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matg3serovar) = grep { grep { $_ eq $matchr3_key } @{$matg3_serovar{$_}} } keys %matg3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matg3{$rev_matg3serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matg3key (sort keys %rev_matg3) { #for each serovar in resistance gene profile
        my $count_revmatg3=@{$rev_matg3{$rev_matg3key}}; #count total
        print OUTFILE "$rev_matg3key\tTotal = $count_revmatg3\n";
        #print OUTFILE "$rev_matg3key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg3{$rev_matg3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples that are with reduced susceptibility and concurred genotypically at most once.\n\n";
    foreach my $matchd_key (sort keys %matchd_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd_rpattern{$matchd_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matdserovar) = grep { grep { $_ eq $matchd_key } @{$matd_serovar{$_}} } keys %matd_serovar; # finds the value (sample ID) in matd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd{$rev_matdserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matdkey (sort keys %rev_matd) { #for each serovar in resistance gene profile
        my $count_revmatd=@{$rev_matd{$rev_matdkey}}; #count total
        print OUTFILE "$rev_matdkey\tTotal = $count_revmatd\n";
        #print OUTFILE "$rev_matdkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd{$rev_matdkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples that are with reduced susceptibility and concurred genotypically at most twice.\n\n";
    foreach my $matchd2_key (sort keys %matchd2_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd2_rpattern{$matchd2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matd2serovar) = grep { grep { $_ eq $matchd2_key } @{$matd2_serovar{$_}} } keys %matd2_serovar; # finds the value (sample ID) in matd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd2{$rev_matd2serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd2key (sort keys %rev_matd2) { #for each serovar in resistance gene profile
        my $count_revmatd2=@{$rev_matd2{$rev_matd2key}}; #count total
        print OUTFILE "$rev_matd2key\tTotal = $count_revmatd2\n";
        #print OUTFILE "$rev_matd2key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd2{$rev_matd2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples that are with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    foreach my $matchd3_key (sort keys %matchd3_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$matchd3_rpattern{$matchd3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_matd3serovar) = grep { grep { $_ eq $matchd3_key } @{$matd3_serovar{$_}} } keys %matd3_serovar; # finds the value (sample ID) in matd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_matd3{$rev_matd3serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_matd3key (sort keys %rev_matd3) { #for each serovar in resistance gene profile
        my $count_revmatd3=@{$rev_matd3{$rev_matd3key}}; #count total
        print OUTFILE "$rev_matd3key\tTotal = $count_revmatd3\n";
        #print OUTFILE "$rev_matd3key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd3{$rev_matd3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at most once.\n\n";
    foreach my $miss_key (sort keys %miss_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss_rpattern{$miss_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_mispserovar) = grep { grep { $_ eq $miss_key } @{$misp_serovar{$_}} } keys %misp_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp{$rev_mispserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_mispkey (sort keys %rev_misp) { #for each serovar in resistance gene profile
        my $count_revmisp=@{$rev_misp{$rev_mispkey}}; #count total
        print OUTFILE "$rev_mispkey\tTotal = $count_revmisp\n";
        #print OUTFILE "$rev_mispkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp{$rev_mispkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at most twice.\n\n";
    foreach my $miss2_key (sort keys %miss2_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss2_rpattern{$miss2_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misp2serovar) = grep { grep { $_ eq $miss2_key } @{$misp2_serovar{$_}} } keys %misp2_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp2{$rev_misp2serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misp2key (sort keys %rev_misp2) { #for each serovar in resistance gene profile
        my $count_revmisp2=@{$rev_misp2{$rev_misp2key}}; #count total
        print OUTFILE "$rev_misp2key\tTotal = $count_revmisp2\n";
        #print OUTFILE "$rev_misp2key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp2{$rev_misp2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    foreach my $miss3_key (sort keys %miss3_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$miss3_rpattern{$miss3_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misp3serovar) = grep { grep { $_ eq $miss3_key } @{$misp3_serovar{$_}} } keys %misp3_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misp3{$rev_misp3serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misp3key (sort keys %rev_misp3) { #for each serovar in resistance gene profile
        my $count_revmisp3=@{$rev_misp3{$rev_misp3key}}; #count total
        print OUTFILE "$rev_misp3key\tTotal = $count_revmisp3\n";
        #print OUTFILE "$rev_misp3key:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp3{$rev_misp3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_discordant samples that are phenotypically resistant but genotypically sensitive.\n\n";
    foreach my $misr_key (sort keys %misr_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$misr_rpattern{$misr_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misgserovar) = grep { grep { $_ eq $misr_key } @{$misg_serovar{$_}} } keys %misg_serovar; # finds the value (sample ID) in misg_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misg{$rev_misgserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misgkey (sort keys %rev_misg) { #for each serovar in resistance gene profile
        my $count_revmisg=@{$rev_misg{$rev_misgkey}}; #count total
        print OUTFILE "$rev_misgkey\tTotal = $count_revmisg\n";
        #print OUTFILE "$rev_misgkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misg{$rev_misgkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_discordant samples that are phenotypically sensitive but maybe genotypically resistant.\n\n";
    foreach my $misd_key (sort keys %misd_rpattern) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$misd_rpattern{$misd_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_misdserovar) = grep { grep { $_ eq $misd_key } @{$misd_serovar{$_}} } keys %misd_serovar; # finds the value (sample ID) in misd_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_misd{$rev_misdserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
    }
    for my $rev_misdkey (sort keys %rev_misd) { #for each serovar in resistance gene profile
        my $count_revmisd=@{$rev_misd{$rev_misdkey}}; #count total
        print OUTFILE "$rev_misdkey\tTotal = $count_revmisd\n";
        #print OUTFILE "$rev_misdkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_misd{$rev_misdkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    foreach my $maybes_key (sort keys %maybes_rpattern) { #goes through the overall resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$maybes_rpattern{$maybes_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_maysserovar) = grep { grep { $_ eq $maybes_key } @{$maybesens_serovar{$_}} } keys %maybesens_serovar; # finds the value (sample ID) in maybesens_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_maybes{$rev_maysserovar}},$reverse_key); #push each serovar to the corresponding sample ID for the overall resistance profile
    }
    for my $rev_maybeskey (sort keys %rev_maybes) { #for each serovar in resistance gene profile
        my $count_revmaybes=@{$rev_maybes{$rev_maybeskey}}; #count total
        print OUTFILE "$rev_maybeskey\tTotal = $count_revmaybes\n";
        #print OUTFILE "$rev_maybeskey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybes{$rev_maybeskey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $mayber_key (sort keys %mayber_rpattern) { #goes through the overallresistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$mayber_rpattern{$mayber_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_mayrserovar) = grep { grep { $_ eq $mayber_key } @{$mayberesis_serovar{$_}} } keys %mayberesis_serovar; # finds the value (sample ID) in mayberesis_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_mayber{$rev_mayrserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that overall resistance profile
    }
    for my $rev_mayberkey (sort keys %rev_mayber) { #for each serovar in resistance gene profile
        my $count_revmayber=@{$rev_mayber{$rev_mayberkey}}; #count total
        print OUTFILE "$rev_mayberkey\tTotal = $count_revmayber\n";
        #print OUTFILE "$rev_mayberkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_mayber{$rev_mayberkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    foreach my $maybed_key (sort keys %maybed_rpattern) { #goes through the overallresistance genes profile for each key sample ID
        my $reverse_key=();
        foreach my $values (@{$maybed_rpattern{$maybed_key}}){
            if ($values){
                $reverse_key=$values;  #string to store the overall resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            }
            else {
                $reverse_key="No resistance genes";
            }
        }
        my ($rev_maydserovar) = grep { grep { $_ eq $maybed_key } @{$maybedecrease_serovar{$_}} } keys %maybedecrease_serovar; # finds the value (sample ID) in maybedecrease_servovar hash to print OUTFILE out the key(serovar)
        push (@{$rev_maybed{$rev_maydserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that overall resistance profile
    }
    for my $rev_maybedkey (sort keys %rev_maybed) { #for each serovar in resistance gene profile
        my $count_revmaybed=@{$rev_maybed{$rev_maybedkey}}; #count total
        print OUTFILE "$rev_maybedkey\tTotal = $count_revmaybed\n";
        #print OUTFILE "$rev_maybedkey:";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybed{$rev_maybedkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
}

if ($opt_n){ #didnt want to create new hash with resistance genes: serovar because there is no unique id to use for key in hash -hence could have duplicates
    print OUTFILE "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at least once.\n\n";
    foreach my $mispheno_key (sort keys %mismatched_pheno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno{$mispheno_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            my ($rev_mispserovar) = grep { grep { $_ eq $mispheno_key } @{$misp_serovar{$_}} } keys %misp_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_misp{$reverse_key}},$rev_mispserovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_mispkey (sort keys %rev_misp) { #for each serovar in resistance gene profile
        my $count_revmisp=@{$rev_misp{$rev_mispkey}}; #count total
        print OUTFILE "$rev_mispkey\tTotal = $count_revmisp\n";
        #print OUTFILE "$rev_mispkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp{$rev_mispkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at least twice.\n\n";
    foreach my $mispheno2_key (sort keys %mismatched_pheno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno2{$mispheno2_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            my ($rev_misp2serovar) = grep { grep { $_ eq $mispheno2_key } @{$misp2_serovar{$_}} } keys %misp2_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_misp2{$reverse_key}},$rev_misp2serovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_misp2key (sort keys %rev_misp2) { #for each serovar in resistance gene profile
        my $count_revmisp2=@{$rev_misp2{$rev_misp2key}}; #count total
        print OUTFILE "$rev_misp2key\tTotal = $count_revmisp2\n";
        #print OUTFILE "$rev_misp2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp2{$rev_misp2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    foreach my $mispheno3_key (sort keys %mismatched_pheno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno3{$mispheno3_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr3
            my ($rev_misp3serovar) = grep { grep { $_ eq $mispheno3_key } @{$misp3_serovar{$_}} } keys %misp3_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_misp3{$reverse_key}},$rev_misp3serovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_misp3key (sort keys %rev_misp3) { #for each serovar in resistance gene profile
        my $count_revmisp3=@{$rev_misp3{$rev_misp3key}}; #count total
        print OUTFILE "$rev_misp3key\tTotal = $count_revmisp3\n";
        #print OUTFILE "$rev_misp3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp3{$rev_misp3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at least once.\n\n";
    foreach my $matchgeno_key (sort keys %matched_geno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno{$matchgeno_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matgserovar) = grep { grep { $_ eq $matchgeno_key } @{$matg_serovar{$_}} } keys %matg_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matg{$reverse_key}},$rev_matgserovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matgkey (sort keys %rev_matg) { #for each serovar in resistance gene profile
        my $count_revmatg=@{$rev_matg{$rev_matgkey}};
        print OUTFILE "$rev_matgkey\tTotal = $count_revmatg\n";
        #print OUTFILE "$rev_matgkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg{$rev_matgkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at least twice.\n\n";
    foreach my $matchgeno2_key (sort keys %matched_geno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno2{$matchgeno2_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matg2serovar) = grep { grep { $_ eq $matchgeno2_key } @{$matg2_serovar{$_}} } keys %matg2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matg2{$reverse_key}},$rev_matg2serovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matg2key (sort keys %rev_matg2) { #for each serovar in resistance gene profile
        my $count_revmatg2=@{$rev_matg2{$rev_matg2key}};
        print OUTFILE "$rev_matg2key\tTotal = $count_revmatg2\n";
        #print OUTFILE "$rev_matg2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg2{$rev_matg2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $matchgeno3_key (sort keys %matched_geno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno3{$matchgeno3_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matg3serovar) = grep { grep { $_ eq $matchgeno3_key } @{$matg3_serovar{$_}} } keys %matg3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matg3{$reverse_key}},$rev_matg3serovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matg3key (sort keys %rev_matg3) { #for each serovar in resistance gene profile
        my $count_revmatg3=@{$rev_matg3{$rev_matg3key}};
        print OUTFILE "$rev_matg3key\tTotal = $count_revmatg3\n";
        #print OUTFILE "$rev_matg3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg3{$rev_matg3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples that are with reduced susceptibility and concurred genotypically at least once.\n\n";
    foreach my $matchdecrease_key (sort keys %matched_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease{$matchdecrease_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matdserovar) = grep { grep { $_ eq $matchdecrease_key } @{$matd_serovar{$_}} } keys %matd_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matd{$reverse_key}},$rev_matdserovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matdkey (sort keys %rev_matd) { #for each serovar in resistance gene profile
        my $count_revmatd=@{$rev_matd{$rev_matdkey}};
        print OUTFILE "$rev_matdkey\tTotal = $count_revmatd\n";
        #print OUTFILE "$rev_matdkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd{$rev_matdkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples that are with reduced susceptibility and concurred genotypically at least twice.\n\n";
    foreach my $matchdecrease2_key (sort keys %matched_decrease2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease2{$matchdecrease2_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matd2serovar) = grep { grep { $_ eq $matchdecrease2_key } @{$matd2_serovar{$_}} } keys %matd2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matd2{$reverse_key}},$rev_matd2serovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matd2key (sort keys %rev_matd2) { #for each serovar in resistance gene profile
        my $count_revmatd2=@{$rev_matd2{$rev_matd2key}};
        print OUTFILE "$rev_matd2key\tTotal = $count_revmatd2\n";
        #print OUTFILE "$rev_matd2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd2{$rev_matd2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples that are with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    foreach my $matchdecrease3_key (sort keys %matched_decrease3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease3{$matchdecrease3_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matd3serovar) = grep { grep { $_ eq $matchdecrease3_key } @{$matd3_serovar{$_}} } keys %matd3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matd3{$reverse_key}},$rev_matd3serovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matd3key (sort keys %rev_matd3) { #for each serovar in resistance gene profile
        my $count_revmatd3=@{$rev_matd3{$rev_matd3key}};
        print OUTFILE "$rev_matd3key\tTotal = $count_revmatd3\n";
        #print OUTFILE "$rev_matd3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd3{$rev_matd3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    foreach my $maybesens_key (sort keys %maybe_sens) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_sens{$maybesens_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_maybesserovar) = grep { grep { $_ eq $maybesens_key } @{$maybesens_serovar{$_}} } keys %maybesens_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_maybes{$reverse_key}},$rev_maybesserovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_maybeskey (sort keys %rev_maybes) { #for each serovar in resistance gene profile
        my $count_revmaybes=@{$rev_maybes{$rev_maybeskey}};
        print OUTFILE "$rev_maybeskey\tTotal = $count_revmaybes\n";
        #print OUTFILE "$rev_maybeskey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybes{$rev_maybeskey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $mayberesis_key (sort keys %maybe_resis) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_resis{$mayberesis_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_mayberserovar) = grep { grep { $_ eq $mayberesis_key } @{$mayberesis_serovar{$_}} } keys %mayberesis_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_mayber{$reverse_key}},$rev_mayberserovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_mayberkey (sort keys %rev_mayber) { #for each serovar in resistance gene profile
        my $count_revmayber=@{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "$rev_mayberkey\tTotal = $count_revmayber\n";
        #print OUTFILE "$rev_mayberkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_mayber{$rev_mayberkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    foreach my $maybedecrease_key (sort keys %maybe_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_decrease{$maybedecrease_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_maybedserovar) = grep { grep { $_ eq $maybedecrease_key } @{$maybedecrease_serovar{$_}} } keys %maybedecrease_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_maybed{$reverse_key}},$rev_maybedserovar); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_maybedkey (sort keys %rev_maybed) { #for each serovar in resistance gene profile
        my $count_revmaybed=@{$rev_maybed{$rev_maybedkey}};
        #print OUTFILE "$rev_maybedkey\tTotal = $count_revmaybed\n";
        print OUTFILE "$rev_maybedkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybed{$rev_maybedkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "*****"x30;
}

if ($opt_o){ #didnt want to create new hash with resistance genes: serovar because there is no unique id to use for key in hash -hence could have duplicates
    print OUTFILE "\nThere are $num_pheno_discordant samples that are phenotypically sensitive but genotypically resistant at least once.\n\n";
    foreach my $mispheno_key (sort keys %mismatched_pheno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno{$mispheno_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            my ($rev_mispserovar) = grep { grep { $_ eq $mispheno_key } @{$misp_serovar{$_}} } keys %misp_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_misp{$rev_mispserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_mispkey (sort keys %rev_misp) { #for each serovar in resistance gene profile
        my $count_revmisp=@{$rev_misp{$rev_mispkey}}; #count total
        print OUTFILE "$rev_mispkey\t\tTotal = $count_revmisp\n";
        #print OUTFILE "$rev_mispkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp{$rev_mispkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno2_discordant samples that are phenotypically sensitive but genotypically resistant at least twice.\n\n";
    foreach my $mispheno2_key (sort keys %mismatched_pheno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno2{$mispheno2_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr2
            my ($rev_misp2serovar) = grep { grep { $_ eq $mispheno2_key } @{$misp2_serovar{$_}} } keys %misp2_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_misp2{$rev_misp2serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_misp2key (sort keys %rev_misp2) { #for each serovar in resistance gene profile
        my $count_revmisp2=@{$rev_misp2{$rev_misp2key}}; #count total
        print OUTFILE "$rev_misp2key\t\tTotal = $count_revmisp2\n";
        #print OUTFILE "$rev_misp2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp2{$rev_misp2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno3_discordant samples that are phenotypically sensitive but genotypically resistant at least thrice.\n\n";
    foreach my $mispheno3_key (sort keys %mismatched_pheno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$mismatched_pheno3{$mispheno3_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash i.e mcr1 mcr3
            my ($rev_misp3serovar) = grep { grep { $_ eq $mispheno3_key } @{$misp3_serovar{$_}} } keys %misp3_serovar; # finds the value (sample ID) in misp_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_misp3{$rev_misp3serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_misp3key (sort keys %rev_misp3) { #for each serovar in resistance gene profile
        my $count_revmisp3=@{$rev_misp3{$rev_misp3key}}; #count total
        print OUTFILE "$rev_misp3key\t\tTotal = $count_revmisp3\n";
        #print OUTFILE "$rev_misp3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_misp3{$rev_misp3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_match samples that are phenotypically resistant and concurred genotypically at least once.\n\n";
    foreach my $matchgeno_key (sort keys %matched_geno) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno{$matchgeno_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matgserovar) = grep { grep { $_ eq $matchgeno_key } @{$matg_serovar{$_}} } keys %matg_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matg{$rev_matgserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matgkey (sort keys %rev_matg) { #for each serovar in resistance gene profile
        my $count_revmatg=@{$rev_matg{$rev_matgkey}};
        print OUTFILE "$rev_matgkey\t\tTotal = $count_revmatg\n";
        #print OUTFILE "$rev_matgkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg{$rev_matgkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno2_match samples that are phenotypically resistant and concurred genotypically at least twice.\n\n";
    foreach my $matchgeno2_key (sort keys %matched_geno2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno2{$matchgeno2_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matg2serovar) = grep { grep { $_ eq $matchgeno2_key } @{$matg2_serovar{$_}} } keys %matg2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matg2{$rev_matg2serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matg2key (sort keys %rev_matg2) { #for each serovar in resistance gene profile
        my $count_revmatg2=@{$rev_matg2{$rev_matg2key}};
        print OUTFILE "$rev_matg2key\t\tTotal = $count_revmatg2\n";
        #print OUTFILE "$rev_matg2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg2{$rev_matg2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno3_match samples that are phenotypically resistant and concurred genotypically at least thrice.\n\n";
    foreach my $matchgeno3_key (sort keys %matched_geno3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_geno3{$matchgeno3_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matg3serovar) = grep { grep { $_ eq $matchgeno3_key } @{$matg3_serovar{$_}} } keys %matg3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matg3{$rev_matg3serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matg3key (sort keys %rev_matg3) { #for each serovar in resistance gene profile
        my $count_revmatg3=@{$rev_matg3{$rev_matg3key}};
        print OUTFILE "$rev_matg3key\t\tTotal = $count_revmatg3\n";
        #print OUTFILE "$rev_matg3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matg3{$rev_matg3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_match samples that are with reduced susceptibility and concurred genotypically at least once.\n\n";
    foreach my $matchdecrease_key (sort keys %matched_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease{$matchdecrease_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matdserovar) = grep { grep { $_ eq $matchdecrease_key } @{$matd_serovar{$_}} } keys %matd_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matd{$rev_matdserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matdkey (sort keys %rev_matd) { #for each serovar in resistance gene profile
        my $count_revmatd=@{$rev_matd{$rev_matdkey}};
        print OUTFILE "$rev_matdkey\t\tTotal = $count_revmatd\n";
        #print OUTFILE "$rev_matdkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd{$rev_matdkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease2_match samples that are with reduced susceptibility and concurred genotypically at least twice.\n\n";
    foreach my $matchdecrease2_key (sort keys %matched_decrease2) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease2{$matchdecrease2_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matd2serovar) = grep { grep { $_ eq $matchdecrease2_key } @{$matd2_serovar{$_}} } keys %matd2_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matd2{$rev_matd2serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matd2key (sort keys %rev_matd2) { #for each serovar in resistance gene profile
        my $count_revmatd2=@{$rev_matd2{$rev_matd2key}};
        print OUTFILE "$rev_matd2key\t\tTotal = $count_revmatd2\n";
        #print OUTFILE "$rev_matd2key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd2{$rev_matd2key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease3_match samples that are with reduced susceptibility and concurred genotypically at least thrice.\n\n";
    foreach my $matchdecrease3_key (sort keys %matched_decrease3) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$matched_decrease3{$matchdecrease3_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_matd3serovar) = grep { grep { $_ eq $matchdecrease3_key } @{$matd3_serovar{$_}} } keys %matd3_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_matd3{$rev_matd3serovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_matd3key (sort keys %rev_matd3) { #for each serovar in resistance gene profile
        my $count_revmatd3=@{$rev_matd3{$rev_matd3key}};
        print OUTFILE "$rev_matd3key\t\tTotal = $count_revmatd3\n";
        #print OUTFILE "$rev_matd3key\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_matd3{$rev_matd3key}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_pheno_maybe samples that are sensitive phenotypically and maybe resistant genotypically.\n\n";
    foreach my $maybesens_key (sort keys %maybe_sens) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_sens{$maybesens_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_maybesserovar) = grep { grep { $_ eq $maybesens_key } @{$maybesens_serovar{$_}} } keys %maybesens_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_maybes{$rev_maybesserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_maybeskey (sort keys %rev_maybes) { #for each serovar in resistance gene profile
        my $count_revmaybes=@{$rev_maybes{$rev_maybeskey}};
        print OUTFILE "$rev_maybeskey\t\tTotal = $count_revmaybes\n";
        #print OUTFILE "$rev_maybeskey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybes{$rev_maybeskey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_geno_maybe samples that are resistant phenotypically and maybe resistant genotypically.\n\n";
    foreach my $mayberesis_key (sort keys %maybe_resis) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_resis{$mayberesis_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_mayberserovar) = grep { grep { $_ eq $mayberesis_key } @{$mayberesis_serovar{$_}} } keys %mayberesis_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_mayber{$rev_mayberserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_mayberkey (sort keys %rev_mayber) { #for each serovar in resistance gene profile
        my $count_revmayber=@{$rev_mayber{$rev_mayberkey}};
        print OUTFILE "$rev_mayberkey\t\tTotal = $count_revmayber\n";
        #print OUTFILE "$rev_mayberkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_mayber{$rev_mayberkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "-----"x20;
    print OUTFILE "\nThere are $num_decrease_maybe samples with reduced susceptibility and maybe resistant genotypically.\n\n";
    foreach my $maybedecrease_key (sort keys %maybe_decrease) { #goes through the resistance genes profile for each key sample ID
        my $reverse_key=();
        foreach (@{$maybe_decrease{$maybedecrease_key}}){
            $reverse_key=$column_index{$_}." ";  #string to store the resistance genes profile to put as key in reverse hash
            my ($rev_maybedserovar) = grep { grep { $_ eq $maybedecrease_key } @{$maybedecrease_serovar{$_}} } keys %maybedecrease_serovar; # finds the value (sample ID) in matg_servovar hash to print OUTFILE out the key(serovar)
            push (@{$rev_maybed{$rev_maybedserovar}},$reverse_key); #push each serovar to the corresponding sample ID for that resistance profile
        }
    }
    for my $rev_maybedkey (sort keys %rev_maybed) { #for each serovar in resistance gene profile
        my $count_revmaybed=@{$rev_maybed{$rev_maybedkey}};
        #print OUTFILE "$rev_maybedkey\tTotal = $count_revmaybed\n";
        print OUTFILE "$rev_maybedkey\n";
        my %sums;
        $sums{$_}++ foreach @{$rev_maybed{$rev_maybedkey}}; #then count the occurance of each serovar for that resistance profile
        print OUTFILE "\t$_:\t$sums{$_}\n" foreach keys %sums;
        print OUTFILE "\n";
    }
    print OUTFILE "*****"x30;
}

close;
