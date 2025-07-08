package mRNA;
use strict;
use warnings;
sub score {
    my ($seq) = @_;
    my %scores = (
        length_score          => length_check($seq),
        start_codon           => start_codon_check($seq),
        stop_codon            => stop_codon_check($seq),
        orf_length            => orf_length_check($seq),
        polyA_tail            => polyA_tail_check($seq),
        gc_content            => gc_content_check($seq),
        promoter_motif        => promoter_motif_check($seq),
        splice_site           => splice_site_check($seq),
        codon_bias            => codon_bias_check($seq),
        exon_intron_structure => exon_intron_structure_check($seq),
    );
    my $total = 0; $total += $_ for values %scores;
    $scores{total} = (scalar(keys %scores) ? $total / scalar(keys %scores) : 0);
    return \%scores;
}
sub length_check {
    my $s = shift;
    my $len = length($s);
    return ($len >= 500 && $len <= 10000) ? 1 : ($len > 10000 ? 10000/$len : $len/500);
}
sub start_codon_check {
    my $s = shift;
    return ($s =~ /^ATG/) ? 1 : 0;
}
sub stop_codon_check {
    my $s = shift;
    my $count = () = $s =~ /TAA|TAG|TGA/g;
    return $count > 0 ? 1 : 0;
}
sub orf_length_check {
    my $s = shift;
    my $max_orf_len = 0;
    for (my $i=0; $i<3; $i++) {
        my $pos = index($s, 'ATG', $i);
        while ($pos != -1) {
            for (my $j=$pos+3; $j < length($s); $j+=3) {
                my $codon = substr($s, $j, 3);
                if ($codon =~ /TAA|TAG|TGA/) {
                    my $orf_len = $j + 3 - $pos;
                    $max_orf_len = $orf_len if $orf_len > $max_orf_len;
                    last;
                }
            }
            $pos = index($s, 'ATG', $pos+1);
        }
    }
    return length($s) ? $max_orf_len / length($s) : 0;
}
sub polyA_tail_check {
    my $s = shift;
    return ($s =~ /A{10,}$/) ? 1 : 0;
}
sub gc_content_check {
    my $s = shift;
    my $gc = ($s =~ tr/GC//);
    return length($s) ? $gc / length($s) : 0;
}
sub promoter_motif_check {
    my $s = shift;
    my $subseq = substr($s, 0, 50);
    return ($subseq =~ /TATA[AT]A[AT]/) ? 1 : 0;
}
sub splice_site_check {
    my $s = shift;
    my $count = () = $s =~ /GT.{0,100}?AG/g;
    return $count > 0 ? 1 : 0;
}
sub codon_bias_check {
    my $s = shift;
    my %common_codons = map {$_ => 1} qw(ATG TTT TTC GAA GAG GGT GGC);
    my $count = 0;
    for (my $i=0; $i+3 <= length($s); $i+=3) {
        my $codon = substr($s, $i, 3);
        $count++ if exists $common_codons{$codon};
    }
    return length($s) ? $count / (length($s)/3) : 0;
}
sub exon_intron_structure_check {
    my $s = shift;
    my $count = () = $s =~ /GT.{10,100}?AG/g;
    return $count > 0 ? 1 : 0;
}
1;
