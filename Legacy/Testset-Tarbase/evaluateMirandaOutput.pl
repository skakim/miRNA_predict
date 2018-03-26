#!/usr/bin/perl  
#use strict;
#use warnings;

############ PATHS ###########
my $rnafold_binary = "RNAfold";


############ VARIABLES ############
my $usage = "
Usage: 
 perl $0 <miranda output> <class:Positive/Negative>

";

my $fileIN = shift or die $usage;
my $class = shift or die $usage; #defines the real class
my $instance = 0;
my $correct = 0;
my $wrong = 0; 


############# MAIN ############

open(INPUT, $fileIN);

while(<INPUT>)
{
	chomp;
	if (/^Performing Scan:/) {
		$instance++;
	}
    elsif (/^Score for this Scan:/) {
        $_ = <INPUT>;
        chomp;
        #print $_."\n";
        if ($_ eq "No Hits Found above Threshold") { #prediction:false
            if ($class eq "Positive") {
                $wrong++; #if its a true target than prediction is wrong
            }
            else {
                $correct++; #if it is not a target than prediction is correct
            }
            #print $_."\n";
        }
        else {
            #print "Not a target\n";
            if ($class eq "Positive") {
                $correct++; #if its a true target than prediction is correct
            }
            else {
                $wrong++; #if it is not a target than prediction is wrong
            }
        }
	}
}
    print "\n###################\nClass: ".$class,"\n";
    print "Total instances in input file: ".$instance."\n";
    print "Correct predictions: ".$correct."\n"; 
    print "Wrong predictions : ".$wrong."\n###################\n\n"; 
