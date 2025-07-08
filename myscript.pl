#!/usr/bin/perl
use strict;
use warnings;

# Codon table 
my %codon_table = (
    'UUU' => 'F', 'UUC' => 'F', 'UUA' => 'L', 'UUG' => 'L',
    'CUU' => 'L', 'CUC' => 'L', 'CUA' => 'L', 'CUG' => 'L',
    'AUU' => 'I', 'AUC' => 'I', 'AUA' => 'I', 'AUG' => 'M',
    'GUU' => 'V', 'GUC' => 'V', 'GUA' => 'V', 'GUG' => 'V',
    'UCU' => 'S', 'UCC' => 'S', 'UCA' => 'S', 'UCG' => 'S',
    'CCU' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
    'ACU' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
    'GCU' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
    'UAU' => 'Y', 'UAC' => 'Y', 'UAA' => '*', 'UAG' => '*',
    'CAU' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
    'AAU' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K',
    'GAU' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
    'UGU' => 'C', 'UGC' => 'C', 'UGA' => '*', 'UGG' => 'W',
    'CGU' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
    'AGU' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R',
    'GGU' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',
);

# Amino acid full names
my %amino_acid_map = (
  'A' => 'Alanine',
  'C' => 'Cysteine',
  'D' => 'Aspartic Acid',
  'E' => 'Glutamic Acid',
  'F' => 'Phenylalanine',
  'G' => 'Glycine',
  'H' => 'Histidine',
  'I' => 'Isoleucine',
  'K' => 'Lysine',
  'L' => 'Leucine',
  'M' => 'Methionine',
  'N' => 'Asparagine',
  'P' => 'Proline',
  'Q' => 'Glutamine',
  'R' => 'Arginine',
  'S' => 'Serine',
  'T' => 'Threonine',
  'V' => 'Valine',
  'W' => 'Tryptophan',
  'Y' => 'Tyrosine'
);

# Metabolic pathways
my %amino_acid_pathways = (
    'Alanine'      => 'Glucogenic (converted to pyruvate for gluconeogenesis or enters TCA cycle). Enzyme: Alanine aminotransferase (ALT) converts alanine to pyruvate.',
    'Arginine'     => 'Glucogenic (converted to ornithine, enters urea cycle). Enzyme: Arginase converts arginine to ornithine.',
    'Asparagine'   => 'Glucogenic (converted to aspartate, enters TCA cycle). Enzyme: Asparaginase converts asparagine to aspartate.',
    'Aspartic Acid'=> 'Glucogenic (converted to oxaloacetate, enters TCA cycle). Enzyme: Aspartate transaminase (AST) converts aspartate to oxaloacetate.',
    'Cysteine'     => 'Glucogenic (converted to pyruvate). Enzyme: Cystathionine β-synthase converts cysteine to serine, and serine is then converted to pyruvate.',
    'Glutamine'    => 'Glucogenic (converted to glutamate, enters TCA cycle). Enzyme: Glutaminase converts glutamine to glutamate. Glutamate dehydrogenase converts glutamate to α-ketoglutarate.',
    'Glutamic Acid'=> 'Glucogenic (converted to α-ketoglutarate, enters TCA cycle). Enzyme: Glutamate dehydrogenase converts glutamate to α-ketoglutarate.',
    'Glycine'      => 'Glucogenic (converted to serine and then pyruvate). Enzyme: Serine hydroxymethyltransferase converts glycine to serine, and serine is then converted to pyruvate.',
    'Histidine'    => 'Glucogenic (converted to glutamate, enters TCA cycle). Enzyme: Histidase converts histidine to urocanate, which is converted to glutamate.',
    'Isoleucine'   => 'Both Glucogenic and Ketogenic (converted to propionyl-CoA and acetyl-CoA). Enzyme: Branched-chain amino acid aminotransferase and isomerase convert isoleucine to propionyl-CoA and acetyl-CoA.',
    'Leucine'      => 'Ketogenic (converted to acetoacetate and acetyl-CoA). Enzyme: Branched-chain amino acid transferase converts leucine to acetoacetate.',
    'Lysine'       => 'Ketogenic (converted to acetoacetyl-CoA). Enzyme: Lysine ketoglutarate pathway converts lysine to acetoacetyl-CoA.',
    'Methionine'   => 'Glucogenic (converted to homocysteine, enters other pathways). Enzyme: Methionine adenosyltransferase converts methionine to S-adenosylmethionine (SAM), and homocysteine is derived from SAM.',
    'Phenylalanine'=> 'Both Glucogenic and Ketogenic (converted to fumarate and acetoacetate). Enzyme: Phenylalanine hydroxylase converts phenylalanine to tyrosine.',
    'Proline'      => 'Glucogenic (converted to glutamate, enters TCA cycle). Enzyme: Proline dehydrogenase converts proline to glutamate.',
    'Serine'       => 'Glucogenic (converted to pyruvate). Enzyme: Serine dehydratase converts serine to pyruvate.',
    'Threonine'    => 'Both Glucogenic and Ketogenic (converted to propionyl-CoA and acetyl-CoA). Enzyme: Threonine dehydratase converts threonine to α-ketobutyrate, which is converted to propionyl-CoA and acetyl-CoA.',
    'Tryptophan'   => 'Both Glucogenic and Ketogenic (converted to acetyl-CoA and acetoacetate). Enzyme: Tryptophan dioxygenase converts tryptophan to kynurenine, which leads to acetyl-CoA and acetoacetate.',
    'Tyrosine'     => 'Both Glucogenic and Ketogenic (converted to fumarate and acetoacetate). Enzyme: Tyrosine aminotransferase converts tyrosine to p-hydroxyphenylpyruvate, which is then converted to fumarate and acetoacetate.',
    'Valine'       => 'Glucogenic (converted to succinyl-CoA, enters TCA cycle). Enzyme: Branched-chain amino acid transaminase converts valine to succinyl-CoA.'
);

# Convert DNA to RNA
sub dna_to_rna {
    my ($dna) = @_;
    my $rna = $dna;
    $rna =~ tr/T/U/;
    return $rna;
}

# Extract exons by splitting RNA on GU...AG splice sites (introns)
sub extract_exons {
    my ($rna) = @_;
    my @exons;
    my $pos = 0;
    while ($rna =~ /GU(.*?)AG/ig) {
        my $start = $-[0];
        my $end = $+[0];
        my $exon = substr($rna, $pos, $start - $pos);
        push @exons, $exon if length($exon) > 0;
        $pos = $end;
    }
    my $last_exon = substr($rna, $pos);
    push @exons, $last_exon if length($last_exon) > 0;
    return @exons;
}

# Translate RNA to amino acid sequence in a given frame (0,1,2)
sub translate_rna {
    my ($rna, $frame) = @_;
    my $aa_seq = '';
    for (my $i = $frame; $i < length($rna) - 2; $i += 3) {
        my $codon = substr($rna, $i, 3);
        $codon =~ tr/tu/TU/;
        $codon =~ tr/T/U/;
        my $aa = $codon_table{$codon} // 'X';
        last if $aa eq '*';
        $aa_seq .= $aa;
    }
    return $aa_seq;
}

# Map single-letter amino acids to full names
sub map_amino_acids {
    my ($seq) = @_;
    my @names;
    foreach my $aa (split //, $seq) {
        if (exists $amino_acid_map{$aa}) {
            push @names, $amino_acid_map{$aa};
        } else {
            push @names, 'Unknown';
        }
    }
    return @names;
}

# Get metabolic pathways for amino acids (unique)
sub get_pathways {
    my (@aa_names) = @_;
    my %seen;
    my %result;
    foreach my $aa (@aa_names) {
        next if $seen{$aa}++;
        if (exists $amino_acid_pathways{$aa}) {
            $result{$aa} = $amino_acid_pathways{$aa};
        } else {
            $result{$aa} = "Unknown pathway";
        }
    }
    return %result;
}

# Open output filehandle
open(my $fh, '>', 'output.txt') or die "Cannot open output.txt: $!";

# Helper function to print to both terminal and file
sub print_both {
    my ($text) = @_;
    print $text;
    print $fh $text;
}

# Main program

print_both("Enter a DNA sequence: ");
my $dna = <STDIN>;
chomp $dna;
$dna = uc $dna;

unless ($dna =~ /^[ACGT]+$/) {
    die "Invalid DNA sequence. Only A,C,G,T allowed.\n";
}

my $rna = dna_to_rna($dna);
my $length = length($dna);

print_both("\nInput DNA length: $length\n");
print_both("RNA sequence:\n$rna\n\n");

# Extract exons
my @exons = extract_exons($rna);

# Filter exons by length (10 to 300 nt)
@exons = grep { length($_) >= 10 && length($_) <= 300 } @exons;

# If no exons found, use whole RNA as one exon
if (!@exons) {
    print_both("No exons found with length 10-300 nt. Using whole RNA sequence.\n");
    @exons = ($rna);
}

print_both("Exons found (length between 10 and 300): " . scalar(@exons) . "\n");

# For each exon, find ORFs in all 3 frames, translate, map names, and pathways
for my $idx (0..$#exons) {
    my $exon = $exons[$idx];
    print_both("\nExon " . ($idx+1) . " (length " . length($exon) . "):\n");
    for my $frame (0..2) {
        my $aa_seq = translate_rna($exon, $frame);
        next unless length($aa_seq) >= 3;
        print_both("  Frame " . ($frame+1) . " AA sequence: $aa_seq\n");
        my @aa_names = map_amino_acids($aa_seq);
        print_both("  AA full names: " . join(", ", @aa_names) . "\n");
        my %pathways = get_pathways(@aa_names);
        print_both("  Metabolic pathways:\n");
        foreach my $aa (keys %pathways) {
            print_both("    $aa: $pathways{$aa}\n");
        }
    }
}

close($fh);
print_both("\nOutput also saved to output.txt\n");
