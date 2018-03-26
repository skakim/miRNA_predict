#!/usr/bin/perl

my @targets = `ls target_*.txt`;
my @mirnas = `ls mir_*.txt`;
my $num_mirnas = scalar(@mirnas);
my $num_targets = scalar(@targets);

@targets = sort(@targets);
@mirnas = sort(@mirnas);

print "miRNAs\tTargets\n";
print $num_mirnas."\t".$num_targets."\n";

if ($num_mirnas != $num_targets) {
	print "ERROR: Different quantity for miRNAs and Targets\n";}
else {
	print "OK!\n";}

for($i = 0; $i < $num_mirnas; $i++){
	my $file_mirnas = $mirnas[$i];
	my $file_targets = $targets[$i];
	#save the file index (different from array index)
	$file_mirnas =~ m/mir_(\d+).txt/;
	my $index = $1;
	chomp($file_mirnas);
	chomp($file_targets);
	print $file_mirnas."\t".$file_targets."\n";
	`miranda $file_mirnas $file_targets > "miranda_$index _results.txt"`;
}
