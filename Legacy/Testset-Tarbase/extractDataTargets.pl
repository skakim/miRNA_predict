#!/usr/bin/perl

my $usage = "Usage: $0 <file name>\n";
my $file = shift or die $usage;
my $TP = 0;
my $FN = 0;

open(FILE, $file);		#open reference file
open(OUT, ">outputTargets.txt");

while (<FILE>) {	
    #print $_."\n";
    #my @info = split(/\t/, $_);
    print OUT $_;
    print OUT "\n\n";
    #print $info[3];
    #   my @info = split(/:/, <FILE>);
    #    if($info[0] eq 'Prediction result') {
    #        chomp ($info[1]);
    #        $mipred_result = $info[1];
    #        $mipred_result =~ s/^\t+//;
    #        print STDERR $mipred_result."\n";
    #        if ($mipred_result eq 'It is a real microRNA precursor'){
    #            $TP++;
    #        }
    #        elsif (($mipred_result eq 'It is a pseudo microRNA precursor') or  ($mipred_result eq 'It is not a pre-miRNA-like hairpin.')){
    #            $FN++;
    #        }
    #    }	
    #    }
}



#print "
##Results:
#\t\t- TP: $TP 
#\t\t- FP: $FP 
#\t\t- TN: $TN 
#\t\t- FN: $FN\n";

