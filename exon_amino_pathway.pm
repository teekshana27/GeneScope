package exon_amino_pathway;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(
    dna_to_rna extract_exons translate_rna map_amino_acids get_pathways
    validate_dna open_output print_both process_exons_to_aminoacids
);
# Codon table
my %codon_table = (
    'UUU'=>'F','UUC'=>'F','UUA'=>'L','UUG'=>'L','CUU'=>'L','CUC'=>'L','CUA'=>'L','CUG'=>'L',
    'AUU'=>'I','AUC'=>'I','AUA'=>'I','AUG'=>'M','GUU'=>'V','GUC'=>'V','GUA'=>'V','GUG'=>'V',
    'UCU'=>'S','UCC'=>'S','UCA'=>'S','UCG'=>'S','CCU'=>'P','CCC'=>'P','CCA'=>'P','CCG'=>'P',
    'ACU'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T','GCU'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A',
    'UAU'=>'Y','UAC'=>'Y','UAA'=>'*','UAG'=>'*','CAU'=>'H','CAC'=>'H','CAA'=>'Q','CAG'=>'Q',
    'AAU'=>'N','AAC'=>'N','AAA'=>'K','AAG'=>'K','GAU'=>'D','GAC'=>'D','GAA'=>'E','GAG'=>'E',
    'UGU'=>'C','UGC'=>'C','UGA'=>'*','UGG'=>'W','CGU'=>'R','CGC'=>'R','CGA'=>'R','CGG'=>'R',
    'AGU'=>'S','AGC'=>'S','AGA'=>'R','AGG'=>'R','GGU'=>'G','GGC'=>'G','GGA'=>'G','GGG'=>'G',
);
# Amino acid map
my %amino_acid_map = (
  'A'=>'Alanine','C'=>'Cysteine','D'=>'Aspartic Acid','E'=>'Glutamic Acid','F'=>'Phenylalanine',
  'G'=>'Glycine','H'=>'Histidine','I'=>'Isoleucine','K'=>'Lysine','L'=>'Leucine','M'=>'Methionine',
  'N'=>'Asparagine','P'=>'Proline','Q'=>'Glutamine','R'=>'Arginine','S'=>'Serine','T'=>'Threonine',
  'V'=>'Valine','W'=>'Tryptophan','Y'=>'Tyrosine'
);
# Metabolic pathways
my %amino_acid_pathways = (
    'Alanine'=>'Glucogenic (→ pyruvate). ALT enzyme.',
    'Arginine'=>'Glucogenic (→ ornithine). Arginase enzyme.',
    'Asparagine'=>'Glucogenic (→ aspartate). Asparaginase enzyme.',
    'Aspartic Acid'=>'Glucogenic (→ oxaloacetate). AST enzyme.',
    'Cysteine'=>'Glucogenic (→ pyruvate). CBS enzyme.',
    'Glutamine'=>'Glucogenic (→ glutamate). Glutaminase enzyme.',
    'Glutamic Acid'=>'Glucogenic (→ α-ketoglutarate). GDH enzyme.',
    'Glycine'=>'Glucogenic (→ serine → pyruvate). SHMT enzyme.',
    'Histidine'=>'Glucogenic (→ glutamate). Histidase enzyme.',
    'Isoleucine'=>'Both (→ propionyl-CoA & acetyl-CoA).',
    'Leucine'=>'Ketogenic (→ acetoacetate).',
    'Lysine'=>'Ketogenic (→ acetoacetyl-CoA).',
    'Methionine'=>'Glucogenic (→ homocysteine → SAM).',
    'Phenylalanine'=>'Both (→ fumarate & acetoacetate).',
    'Proline'=>'Glucogenic (→ glutamate).',
    'Serine'=>'Glucogenic (→ pyruvate).',
    'Threonine'=>'Both (→ propionyl-CoA & acetyl-CoA).',
    'Tryptophan'=>'Both (→ acetyl-CoA & acetoacetate).',
    'Tyrosine'=>'Both (→ fumarate & acetoacetate).',
    'Valine'=>'Glucogenic (→ succinyl-CoA).'
);
sub dna_to_rna {
    my ($dna) = @_;
    $dna =~ tr/T/U/;
    return $dna;
}
sub validate_dna {
    my ($dna) = @_;
    return ($dna =~ /^[ACGT]+$/i);
}
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
sub translate_rna {
    my ($rna, $frame) = @_;
    my $aa_seq = '';
    for (my $i = $frame; $i < length($rna) - 2; $i += 3) {
        my $codon = substr($rna, $i, 3);
        $codon = uc $codon;
        $codon =~ tr/T/U/;
        my $aa = $codon_table{$codon} // 'X';
        last if $aa eq '*';
        $aa_seq .= $aa;
    }
    return $aa_seq;
}
sub map_amino_acids {
    my ($seq) = @_;
    my @names;
    foreach my $aa (split //, $seq) {
        push @names, $amino_acid_map{$aa} // 'Unknown';
    }
    return @names;
}
sub get_pathways {
    my (@aa_names) = @_;
    my %seen;
    my %result;
    foreach my $aa (@aa_names) {
        next if $seen{$aa}++;
        $result{$aa} = $amino_acid_pathways{$aa} // "Unknown pathway";
    }
    return %result;
}
sub open_output {
    my $file = shift // 'exon_amino_output.txt';
    open(my $fh, '>', $file) or die "Cannot open $file: $!";
    return $fh;
}
sub print_both {
    my ($text, $fh) = @_;
    print $text;
    print $fh $text if $fh;
}
sub process_exons_to_aminoacids {
    my ($dna) = @_;
    my $rna = dna_to_rna($dna);
    my @exons = extract_exons($rna);
    my @final_output;
    foreach my $exon (@exons) {
for my $frame (0..2) {
    my $aa_seq = translate_rna($exon, $frame);
    my @names = map_amino_acids($aa_seq);
    my %pathways = get_pathways(@names);
    push @final_output, {
        exon => $exon,
        frame => $frame,
        aa_sequence => $aa_seq,
        names => \@names,
        pathways => \%pathways,
    };
}
}
    return @final_output;
}
1;

