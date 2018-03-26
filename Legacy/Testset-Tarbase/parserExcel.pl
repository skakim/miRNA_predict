#!/usr/bin/perl -w
#run cpan Spreadsheet::ParseExcel to install packages
use strict;
use Spreadsheet::ParseExcel;


my $usage = "PERL script to parse a worksheet from an excel file.\nUsage: $0 <excel file> <worksheet index>\n\n";
my $file = shift or die $usage;
my $index = shift or die $usage;

my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->Parse($file);
my $worksheet = $workbook->worksheet($index-1);



#for my $worksheet ( $workbook->worksheets() ) {  

my ( $row_min, $row_max ) = $worksheet->row_range();
my ( $col_min, $col_max ) = $worksheet->col_range();
for my $row (1 .. $row_max ) {
#	open(OUT, ">mirna_$row");
#	my $mirna_id = $worksheet->get_cell( $row, 3);
#	my $mirna_seq = $worksheet->get_cell( $row, 4);
#	print OUT ">hsa-".$mirna_id->value()."\n".$mirna_seq->value();
#	close(OUT);
	
 	open(OUT, ">target_$row");
 	my $target_id = $worksheet->get_cell( $row, 7);
 	my $target_seq = $worksheet->get_cell( $row, 8);
 	#my $sequence = $target_seq->value();
 	#$sequence =~ tr{\n}{""};
 	print OUT $target_id->value()." \n".$target_seq->value()."\n";
 	close(OUT);
}
#}
