#!/usr/bin/perl
use strict;
use warnings;
use lib '.';  

use rnapredictor qw(predict_rna_types);
use gene_epigenetic_identifier qw(analyze_gene_and_epigenetics);
use exon_amino_pathway qw(dna_to_rna extract_exons translate_rna map_amino_acids get_pathways 
 validate_dna open_output print_both process_exons_to_aminoacids);
use seq_mutations qw(read_sequence_from_file validate_sequence compare_sequences);

sub main {
print "Welcome to the Gene Pipeline\n";
print "-----------------------------\n";
print "Enter DNA sequence manually or load from file? (m/f):- ";
chomp(my $input_type = <STDIN>);
my $dna_seq = '';
if ($input_type eq 'm') {
print "Enter DNA sequence:- ";
chomp($dna_seq = <STDIN>);
$dna_seq = uc($dna_seq);
unless (validate_dna($dna_seq)) {
die "Invalid DNA sequence. Only A,C,G,T allowed.\n";
}
} elsif ($input_type eq 'f') {
print "Enter filename:- ";
chomp(my $filename = <STDIN>);
$dna_seq = read_sequence_from_file($filename);
unless (defined $dna_seq && validate_sequence($dna_seq)) {
die "Invalid sequence or file.\n";
}
} else {
die "Invalid input choice: '$input_type'. Use 'm' or 'f'.\n";
}
while (1) {
print"------------------------------------------------------------------------\n";
print "Choose an analysis:\n";
print "1. Gene identification, epigenetic modification and gene activation and silencing\n";
print "2. mRNA, tRNA, miRNA, eRNA, rRNA sequences\n";
print "3. RNA sequence, exons, amino acids and their metabolic pathway\n";
print "4. Substitutions, INDELs, percent identity, and visual alignment\n";
print "5. Exit\n";
print"------------------------------------------------------------------------\n";
print "Choice:- ";
chomp(my $choice = <STDIN>);
if ($choice eq '1') {
my $result = analyze_gene_and_epigenetics($dna_seq);
print "\nGene & Epigenetics Report:-\n$result\n";
}
elsif ($choice eq '2') {
my $rna_report = predict_rna_types($dna_seq);
print "\nRNA Type Prediction:-\n$rna_report\n";
}
elsif ($choice eq '3') {
my $fh = open_output("exon_amino_output.txt");
my $rna = dna_to_rna($dna_seq);
print_both("\nRNA Sequence:-\n$rna\n\n", $fh);
my @exons = extract_exons($rna);
print_both("Exons Detected:-\n", $fh);
print_both("$_\n", $fh) for @exons;
print_both("\n", $fh);
foreach my $exon (@exons) {
print_both("-----Exon Information-----\n", $fh);
print_both("Exon Sequence:- $exon\n", $fh);
for my $frame (0..2){
my $aa_seq = translate_rna($exon, $frame);
my @aa_names = map_amino_acids($aa_seq);
my %pathways = get_pathways(@aa_names);
print_both("\nReading Frame $frame:-\n", $fh);
print_both("Amino Acid Sequence:- $aa_seq\n", $fh);
print_both("Amino Acids:- " . join(", ", @aa_names) . "\n", $fh);
print_both("Metabolic Pathways:-\n", $fh);
print_both("  $_ :- $pathways{$_}\n", $fh) for keys %pathways;
}
print_both("\n--------------------------\n\n", $fh);
}
print_both("Results also saved to 'exon_amino_output.txt'\n", $fh);
close $fh;
}
elsif ($choice eq '4') {
print "Enter sample sequence or FASTA file? (s/f): ";
chomp(my $cmp_input = <STDIN>);
my $sample_seq;
if ($cmp_input eq 'f') {
print "Enter sample FASTA file:- ";
chomp(my $file2 = <STDIN>);
$sample_seq = read_sequence_from_file($file2);
unless (defined $sample_seq && validate_sequence($sample_seq)) {
print "Invalid sequence or file.\n";
next;
}
} elsif ($cmp_input eq 's') {
print "Enter sample DNA sequence:- ";
chomp($sample_seq = <STDIN>);
$sample_seq = uc($sample_seq);
unless (validate_sequence($sample_seq)) {
print "Invalid sequence entered.\n";
next;
}
} else {
print "Invalid input. Choose 's' or 'f'.\n";
next;
}
my ($report, $alignment, $identity) = compare_sequences($dna_seq, $sample_seq);
my $outfile = "comparison_output.txt";
open my $fh, '>', $outfile or die "Cannot open $outfile: $!";
print "\nSequence Comparison Report:-\n";
print $report;
print $alignment;
print "Results also saved to '$outfile'\n";
print $fh "Sequence Comparison Report:-\n";
print $fh $report;
print $fh $alignment;
close $fh;
}
elsif ($choice eq '5') {
print "\nExiting pipeline. Thank you!\n";
last;
}
else {
print "Invalid choice. Please try again.\n";
}
}
}
main();
