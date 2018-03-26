#!/usr/bin/perl

my @targets = `ls neg_*_target.fas`;
my $num_targets = scalar(@targets);


@targets = sort(@targets);

for($i = 0; $i < $num_targets; $i++){
    my %hashMapSeq = ();
    chomp($targets[$i]);
    my $file_targets = $targets[$i];
    $file_targets =~ m/neg_(\d+)_target.fas/;
    my $index = $1;

    `head -n 2 $targets[$i] > target_$index.fas`
}