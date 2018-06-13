#!/usr/bin/perl  
use strict;
#use warnings;

############ PATHS ###########
my $rnaduplex_binary = "RNAduplex";


############ VARIABLES ############
my $usage = "
Usage: 
 perl $0 <miranda output> <class>

";

my $fileIN = shift or die $usage;
my $class_str = shift;
my $instance = 0;
my $read_size = 20;
my $seed_size = 8;
my $class;
my $score;
my $alignLen;
my $query_str;
my $ref_str;
my $energy;
my $totalGC; 
my $totalAT;
my $totalGT;
my $totalGaps;
my $totalMismatches;
my $totalMFE;
my $seedGC; 
my $seedAT;
my $seedGT;
my $seedGaps;
my $seedMismatches;
my $seedMFE;
my $mirnaID;
my $targetID;
my $mirnaFirstPos;
my $mirnaLastPos;
my $targetFirstPos;
my $targetLastPos;
my $linker = "AAAGGGLLLLLLCCCTTT";
my $full_str;
my $length_query;
my $length_ref;
my $RNAduplexInputFile = "InputRNAduplex.txt";
my $RNAduplexOutputFile = "OutputRNAduplex.txt";
my @query;
my @ref;
my @alignRes;
my @alignPos;


#print out all the features extracted from alignment 


############# MAIN ############

open(INPUT, $fileIN);
open(OUTPUT_FEAT, ">$fileIN.att"); #file with attributes
open(OUTPUT_IDS, ">$fileIN.ids");  #file with IDs and alignment positions
open(OUTPUT_ALI, ">$fileIN.ali");  #file with alignment

my $outputFeatures;
my $outputIDs;
my $outputAlignment;

if ($class_str) {
	print OUTPUT_FEAT "Score\tAlignLen\ttotalMFE\tseedMFE\ttotalGC\ttotalAU\ttotalGU\ttotalGaps\ttotalMismatches\tseedGC\tseedAU\tseedGU\tseedGaps\tseedMismatches\tPos_1\tPos_2\tPos_3\tPos_4\tPos_5\tPos_6\tPos_7\tPos_8\tPos_9\tPos_10\tPos_11\tPos_12\tPos_13\tPos_14\tPos_15\tPos_16\tPos_17\tPos_18\tPos_19\tPos_20\tClass\n";
}
else {
	print OUTPUT_FEAT "Score\tAlignLen\ttotalMFE\tseedMFE\ttotalGC\ttotalAU\ttotalGU\ttotalGaps\ttotalMismatches\tseedGC\tseedAU\tseedGU\tseedGaps\tseedMismatches\tPos_1\tPos_2\tPos_3\tPos_4\tPos_5\tPos_6\tPos_7\tPos_8\tPos_9\tPos_10\tPos_11\tPos_12\tPos_13\tPos_14\tPos_15\tPos_16\tPos_17\tPos_18\tPos_19\tPos_20\n";
}

while(<INPUT>)
{
	chomp;
	if (/^Performing Scan:/) {
		#new alignment! increment instance
		$instance++;

		#get sequences IDs
		#$_ = <INPUT>;
		($mirnaID,$targetID) = getIDs($_);
		$outputIDs = "$mirnaID\t$targetID";

		#get score and align length
		$score = getScore($_);
		$alignLen = getAlignLength($_);
        $outputFeatures = $score."\t".$alignLen;
        
        #get alignment position between query and reference
		($mirnaFirstPos,$mirnaLastPos,$targetFirstPos,$targetLastPos) = getAlignPos($_);
	}
	elsif (/^   Forward:/) {
		#new alignment! increment instance
		#$instance++;
		
		#get score and align length
		$score = getScore($_);
		$alignLen = getAlignLength($_);
        $outputFeatures = $score."\t".$alignLen;
        
        #get alignment position between query and reference
		($mirnaFirstPos,$mirnaLastPos,$targetFirstPos,$targetLastPos) = getAlignPos($_);
	}
	elsif (/^   Query:/) {
		#get query sequence
		$query_str = uc(getQuerySequence($_));      #string with query sequence
        $query_str =~ tr/U/T/;                      #change U into U 
        @query = split(//,$query_str);              #array with query sequence
		
		#get alignment
		$_ = <INPUT>;
		$length_query = scalar @query;              #length of query sequence
		@alignRes = getAlignRes($_,$length_query);  #array with alignment result
	}
	elsif (/^   Ref:/) {
		#get refeference sequence
		$ref_str = uc(getRefSequence($_));          #string with reference sequence
		@ref = split(//,$ref_str);                  #array with query sequence
        $length_ref = scalar @ref;
        
		#compute total MFE with RNAduplex
		$full_str = $query_str."_".$ref_str;
		$totalMFE = getMFE($full_str);  
		$outputFeatures = $outputFeatures."\t".$totalMFE;
        
        #compute seed MFE with RNAduplex
        my $seed_query = substr($query_str,($length_query-$seed_size),-1);
        my $seed_ref = substr($ref_str,($length_ref-$seed_size),-1);
        #print STDERR $query_str."\n".$seed_query."\n".$seed_ref."\n";
		$full_str = $seed_query."_".$seed_ref;
		$seedMFE = getMFE($full_str);  
		$outputFeatures = $outputFeatures."\t".$seedMFE;
	}

	elsif(/^Scores for this hit:/){
		#compute TOTAL matchs, wobbles, mismatches and gaps
		($totalGC,$totalAT,$totalGT) = computeMatches(\@query, \@ref);
		$totalGaps = computeGaps(\@query, \@ref);
		$totalMismatches = $read_size - ($totalGaps + $totalGC + $totalAT + $totalGT);
		$outputFeatures = $outputFeatures."\t".$totalGC."\t".$totalAT."\t".$totalGT."\t".$totalGaps."\t".$totalMismatches;
        
        #compute SEED matchs, wobbles, mismatches and gaps (seed: 2 to 8, 5')
		($seedGC,$seedAT,$seedGT) = computeMatchesSeed(\@query, \@ref);
		$seedGaps = computeGapsSeed(\@query, \@ref);
		$seedMismatches = $seed_size - ($seedGaps + $seedGC + $seedAT + $seedGT + 1);
		$outputFeatures = $outputFeatures."\t".$seedGC."\t".$seedAT."\t".$seedGT."\t".$seedGaps."\t".$seedMismatches;

		
		#compute alignment type per position 
		@alignPos = computePositionAlign(\@query, \@ref);
		my @alignPos = grep /\S/, @alignPos; 
		foreach (@alignPos) {
		    $outputFeatures = $outputFeatures."\t".$_;
		} 

		print OUTPUT_IDS $outputIDs."\t".$mirnaFirstPos."\t".$mirnaLastPos."\t".$targetFirstPos."\t".$targetLastPos."\n";
		if ($class_str) {
			print OUTPUT_FEAT $outputFeatures."\t".$class_str."\n";
		}
		else {
			print OUTPUT_FEAT $outputFeatures."\n";
		}
		my $alignString = join('', @alignRes);
		print OUTPUT_ALI $query_str."\t".$alignString."\t".$ref_str."\n";
		
	}
}
close(INPUT);
close(OUTPUT_FEAT);
close(OUTPUT_IDS);
close(OUTPUT_ALI);


############# FUNCTIONS ##############

sub getScore{
	my $string = shift;
    my $value;
    if($string =~ m/Score: (\d+.\d+)/){
        $value = $1;
    }
    return $value;
}
	
sub getAlignLength{
    my $string = shift;
    my $value;
    if($string =~ m/Align Len \((\d+)\)/){
        $value = $1;
    }
    return $value;
}

sub getAlignPos{
    my $string = shift;
    my $mirnaFPos;
    my $mirnaLPos;
    my $targetFPos;
    my $targetLPos;
    if($string =~ m/Q:(\d+) to (\d+)/){
        $mirnaFPos = $1;
	$mirnaLPos = $2;
    }
    if($string =~ m/R:(\d+) to (\d+)/){
        $targetFPos = $1;
	$targetLPos = $2;
    }
    return ($mirnaFPos,$mirnaLPos,$targetFPos,$targetLPos);
}

sub getQuerySequence{
    my $string = shift;
    my $value;
    if($string =~ m/3\'(\s+)(\S+)/){
        $value = $2;
    }
    return $value;
}

sub getAlignRes{
    my $string = shift;
    my $length = shift;
    my @line = ();
    @line = split(//,$string);
    my @align= @line[16..(16+$length-1)];
    return @align;

}
        
sub getRefSequence{
    my $string = shift;
    my $value;
    if($string =~ m/5\'(\s+)(\S+)/){
        $value = $2;
    }
    return $value;
}

sub getMFE{
    my $string = shift;
    my @string_arr = split("_",$string);
    my $value = 0;
    my $line;
    my @line = ();
    open IN,">$RNAduplexInputFile";
    print IN ">mirna\n".$string_arr[0]."\n>target\n".$string_arr[1];
    close IN;
    `$rnaduplex_binary <$RNAduplexInputFile >$RNAduplexOutputFile\n`;
    open IN,"$RNAduplexOutputFile";
    while(<IN>) {
	chomp;
	if(!/^>/) {
		$line = $_;
		@line = split(/  /,$line);
		$value = substr($line[4], 2, -1)
        #$value = $line[4];
	}
    }
    close IN;
    `rm $RNAduplexInputFile $RNAduplexOutputFile`;
    return $value;
}

#sub getMFEseed{
#    my $string = shift;
#    my $value = 0;
#    my $line;
#    my @line = ();
#    open IN,">$RNAduplexInputFile";
#    print IN $string;
#    print IN "\@";# @ is needed by "RNAduplex"
#    close IN;
#    `$rnaduplex_binary -noPS <$RNAduplexInputFile >$RNAduplexOutputFile\n`;
#    open IN,"$RNAduplexOutputFile";
#    while(<IN>) {
#        chomp;
#        if(/^[.(]/) {
#            $line = $_;
#            @line = split(/ /,$line);
#            $value = substr($line[1], 1, -1)
#        }
#    }
#    close IN;
#    `rm $RNAduplexInputFile $RNAduplexOutputFile`;
#    return $value;
#}


sub getEnergy{
    my $string = shift;
    my $value;
    if($string =~ m/Energy:  (\S+)/){
        $value = $1;
    }
    return $value;
}

sub getIDs{
    my $string = shift;
    my @line = ();
    @line = split(/\t/,$string);
    my $mirna = $line[0];
    my $target = $line[1];
    $mirna =~ tr/>//d;
    $target =~ tr/>//d;
    return ($mirna,$target);
}

sub computeMatches{
    my ($query_ref, $ref_ref) = @_;
    my @query = @{$query_ref};
    my @ref = @{$ref_ref};
    my $length = scalar @query;
    my $limit = $length - $read_size;
    my $GC = 0;
    my $AT = 0;
    my $GT = 0;
    
    #from last position, read $read_size nucleotides (value: 20)
    for(my $i = $length-1; $i >= $limit; $i--) {
        if (($query[$i] eq 'G' and $ref[$i] eq 'C') or ($ref[$i] eq 'G' and $query[$i] eq 'C'))
        {
            $GC++;
        }
        elsif (($query[$i] eq 'A' and $ref[$i] eq 'T') or ($ref[$i] eq 'A' and $query[$i] eq 'T'))
        {
            $AT++;
        }
        elsif (($query[$i] eq 'G' and $ref[$i] eq 'T') or ($ref[$i] eq 'G' and $query[$i] eq 'T'))
        {
            $GT++;
        }
    }
    return ($GC,$AT,$GT);    
}

sub computeMatchesSeed{
    my ($query_ref, $ref_ref) = @_;
    my @query = @{$query_ref};
    my @ref = @{$ref_ref};
    my $length = scalar @query;
    my $limit = $length - $seed_size;
    my $GC = 0;
    my $AT = 0;
    my $GT = 0;
    
    #from second position (5'), read $read_size nucleotides (value: 7)
    for(my $i = $length-2; $i >= $limit; $i--) {
        if (($query[$i] eq 'G' and $ref[$i] eq 'C') or ($ref[$i] eq 'G' and $query[$i] eq 'C'))
        {
            $GC++;
        }
        elsif (($query[$i] eq 'A' and $ref[$i] eq 'T') or ($ref[$i] eq 'A' and $query[$i] eq 'T'))
        {
            $AT++;
        }
        elsif (($query[$i] eq 'G' and $ref[$i] eq 'T') or ($ref[$i] eq 'G' and $query[$i] eq 'T'))
        {
            $GT++;
        }
    }
    return ($GC,$AT,$GT);    
}


sub computeGaps{
    my ($query_ref, $ref_ref) = @_;
    my @query = @{$query_ref};
    my @ref = @{$ref_ref};
    my $length = scalar @query;
    my $limit = $length - $read_size;
    my $gap = 0;
    
    #from last position, read $read_size nucleotides (value: 20)
    for(my $i = $length-1; $i >= $limit; $i--) {
        if (($query[$i] eq '-') or ($ref[$i] eq '-'))
        {
            $gap++;
        }
    }
    return ($gap);    
}

sub computeGapsSeed{
    my ($query_ref, $ref_ref) = @_;
    my @query = @{$query_ref};
    my @ref = @{$ref_ref};
    my $length = scalar @query;
    my $limit = $length - $seed_size;
    my $gap = 0;
    
    #from second position (5'), read $read_size nucleotides (value: 7)
    for(my $i = $length-2; $i >= $limit; $i--) {
        if (($query[$i] eq '-') or ($ref[$i] eq '-'))
        {
            $gap++;
        }
    }
    return ($gap);    
}

sub computePositionAlign{
    my ($query_ref, $ref_ref) = @_;
    my @query = @{$query_ref};
    my @ref = @{$ref_ref};
    my $length = scalar @query;
    my $limit = $length - $read_size;
    my @positions = ();
    my $pos = 0;
    #from last position, read $read_size nucleotides (value: 20)
    for(my $i = $length-1; $i >= $limit; $i--) {
        if (($query[$i] eq 'G' and $ref[$i] eq 'C') or ($ref[$i] eq 'G' and $query[$i] eq 'C'))
        {
            $positions[$pos] = 1;  #GC alignment
        }
        elsif (($query[$i] eq 'A' and $ref[$i] eq 'T') or ($ref[$i] eq 'A' and $query[$i] eq 'T'))
        {
            $positions[$pos] = 2;  #AT alignment
        }
        elsif (($query[$i] eq 'G' and $ref[$i] eq 'T') or ($ref[$i] eq 'G' and $query[$i] eq 'T'))
        {
            $positions[$pos] = 3;  #GT alignment
        }
        elsif (($query[$i] eq '-') or ($ref[$i] eq '-'))
        {
            $positions[$pos] = 4;  #gap
        }
        else
        {
            $positions[$pos] = 5;  #mismatch
        }
	$pos++;
    }    
 
    return (@positions);    
}

