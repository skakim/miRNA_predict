#!/usr/bin/perl -w
use strict;
use warnings;
use Time::HiRes;
use Getopt::Std;
use vars qw($opt_m $opt_t $opt_o);
getopts("m:t:o:");


my $start_time = [Time::HiRes::gettimeofday()];


###########################################################################
# DEFINE PATH VARIABLES
###########################################################################

# set this parameter to the location of the file modelRF.Rdata
#my $Rmodel = "~/Dropbox/UFRGS/Disciplinas/2011-02/CMP569 Algoritmos em Bioinformatica/Projeto/modelRF.RData";
my $Rmodel = "/home/mariana/Dropbox/UFRGS/Disciplinas/2011-02/CMP569 Algoritmos em Bioinformatica/Projeto/modelRFseed.RData";

# set this parameter to point to your R binary
#my $Rbinary = "R64";
my $Rbinary = "R";


###########################################################################
# DEFINE PROGRAM VARIABLES
###########################################################################

my $Rscript = "temp.R";
my $usage = "\nmiRNA target prediction by Random Forests.\nUsage: perl $0 -m file1 -t file2

Where:
	'file1' is a FASTA file with a microRNA query
        'file2' is a FASTA file containing the sequence(s) to be scanned.\n\n";
my $mirnaFile = $opt_m || die "$usage";
my $targetFile = $opt_t || die "$usage";
my $Routput = "Rprediction.csv";
my $mirandaOutput = "tempMiranda.out";
my $idFile = "${mirandaOutput}.ids";
my $attributesFile = "${mirandaOutput}.att";
my $alignFile = "${mirandaOutput}.ali";
my $finalOutput = "TargetsPredicted.csv";
my $finalOutputAli = "TargetsPredicted.ali";

# variables for storing results
my @predictionProbability = ();
my $predictionOutput = "miRNA ID\ttarget ID\tmiRNA first pos\tmiRNA last pos\tTarget first pos\tTarget miRNA last pos\tPrediction Confidence\n";

print STDERR "\n>> miRNA file: $mirnaFile\n>> Candidate targets file: $targetFile\n>> Starting miRNA target prediction...\n";


###########################################################################
# RUN MIRANDA ALIGNER
# EXTRACT FEATURES FOR CLASSIFICATION
###########################################################################
`miranda $mirnaFile $targetFile > $mirandaOutput`;
`perl extractAttributes.pl $mirandaOutput`;


###########################################################################
# MAKE PREDICTION
###########################################################################

#predict class using the trained model (random Forest)
my @Rlines = ( 'library(randomForest)',
           'load("' . $Rmodel . '")',
           'test<-read.table("' . $attributesFile . '",header = TRUE)',
           'pred<-predict(model,test,type="prob")',
           'decision<-pred[,2]',
           'decision<-as.matrix(decision)',
           'write.table(decision,"' . $Routput . '",quote=FALSE,sep=" ",col.names = FALSE,row.names = FALSE)',
           'q("no")');


open SCRIPT, ">$Rscript";
for my $line (@Rlines){
print SCRIPT "$line\n";
}
close SCRIPT;

#run R
my $aa = system("$Rbinary CMD BATCH $Rscript ");


###########################################################################
# PREPARE OUTPUT
# WRITE RESULTS ON FILE
###########################################################################

#read prediction probabilities
open(DATA, $Routput) || die "   ERROR. Can't open $Routput: $!\n";
@predictionProbability = <DATA>;
close DATA;

#read info about miRNA and target IDs' and alignment's positions (first and last)
open (DATA, $idFile) || die "   ERROR. Can't open $idFile: $!\n";
my $index = 0;
my @alignmentIDs = ();
my @line = ();
while (<DATA>) {
    chomp;
    if ($predictionProbability[$index] > 0.5) {
	$predictionOutput = $predictionOutput.$_."\t".($predictionProbability[$index]*100)."\n";
	@line = split(/\t/,$_);
	$alignmentIDs[$index] = $line[0]." vs ".$line[1];
    }
    $index++;
}
open(OUTPUT, ">$finalOutput");
print OUTPUT $predictionOutput;
close(OUTPUT); 

#read info about alignment between miRNA and target
open (DATA, $alignFile) || die "   ERROR. Can't open $alignFile: $!\n";
$index = 0;
my $alignOutput="";
@line = ();
while (<DATA>) {
    chomp;
    if ($predictionProbability[$index] > 0.5) {
	@line = split(/\t/,$_);
	$alignOutput = $alignOutput.$alignmentIDs[$index]."\n".$line[0]."\n".$line[1]."\n".$line[2]."\n\n";
    }
    $index++;
}
open(OUTPUT, ">$finalOutputAli");
print OUTPUT $alignOutput;
close(OUTPUT); 


###########################################################################
# FINAL COMMANDS
###########################################################################

`rm $Rscript ${Rscript}out $Routput`;
`rm $mirandaOutput $idFile $attributesFile $alignFile`;

print STDERR ">> Predictions saved on file $finalOutput\n";

my $diff = Time::HiRes::tv_interval($start_time);
print "\n## Elapsed time: $diff seconds\n\n";

