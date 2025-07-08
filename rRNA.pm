package rRNA;
use strict;
use warnings;
sub score {
    my ($seq) = @_;
    my %scores = (
        length_score          => length_check($seq),
        gc_content            => gc_content_check($seq),
        conserved_motif       => conserved_motif_check($seq),
        secondary_structure   => secondary_structure_check($seq),
        high_stop_codons      => high_stop_codons_check($seq),
        no_polyA_tail         => no_polyA_tail_check($seq),
        repeated_regions      => repeated_regions_check($seq),
        pyrimidine_rich       => pyrimidine_rich_check($seq),
        sequence_complexity   => sequence_complexity_check($seq),
        conserved_regions     => conserved_regions_check($seq),
    );
    my $total = 0; $total += $_ for values %scores;
    $scores{total} = (scalar(keys %scores) ? $total / scalar(keys %scores) : 0);
    return \%scores;
}
sub length_check {
    my $s = shift;
    my $len = length($s);
    return ($len >= 1200 && $len <= 5000) ? 1 : ($len > 5000 ? 5000/$len : $len/1200);
}
sub gc_content_check {
    my $s = shift;
    my $gc = ($s =~ tr/GC//);
    my $gc_ratio = length($s) ? $gc / length($s) : 0;
    return 1 - abs(0.62 - $gc_ratio)/0.2 > 0 ? 1 - abs(0.62 - $gc_ratio)/0.2 : 0;
}
sub conserved_motif_check {
    my $s = shift;
    return ($s =~ /CGGAGG/) ? 1 : 0;
}
sub secondary_structure_check {
    my $s = shift;
    my $count = 0;
    for my $i (0..length($s)-8) {
        my $sub = substr($s, $i, 8);
        my $revcomp = reverse $sub;
        $revcomp =~ tr/ATGC/TACG/;
        $count++ if index($s, $revcomp, $i+8) != -1;
    }
    return $count > 2 ? 1 : 0;
}
sub high_stop_codons_check {
    my $s = shift;
    my $stop_count = () = $s =~ /TAA|TAG|TGA/g;
    my $ratio = length($s) ? $stop_count / (length($s)/3) : 0;
    return $ratio > 0.3 ? 1 : $ratio;
}
sub no_polyA_tail_check {
    my $s = shift;
    return ($s =~ /A{10,}$/) ? 0 : 1;
}
sub repeated_regions_check {
    my $s = shift;
    my $count = () = $s =~ /([^A])\1{4,}/g;
    return $count > 0 ? 1 : 0;
}
sub pyrimidine_rich_check {
    my $s = shift;
    my $py = ($s =~ tr/CTct//);
    my $ratio = length($s) ? $py / length($s) : 0;
    return $ratio > 0.4 ? 1 : $ratio/0.4;
}
sub sequence_complexity_check {
    my $s = shift;
    my %freq;
    $freq{$_}++ for split //, $s;
    my $len = length($s);
    my $entropy = 0;
    foreach my $b (keys %freq) {
        my $p = $freq{$b} / $len;
        $entropy -= $p * log($p);
    }
    return $entropy && $len ? $entropy / log(4) : 0;
}
sub conserved_regions_check {
    my $s = shift;
    return ($s =~ /GGG/) ? 1 : 0;
}
1;
