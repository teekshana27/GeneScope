package miRNA;
use strict;
use warnings;
sub score {
    my ($seq) = @_;
    my %scores = (
        length_score       => length_check($seq),
        hairpin_score      => hairpin_check($seq),
        au_rich            => au_rich_check($seq),
        no_long_orf        => no_long_orf_check($seq),
        no_polyA_tail      => no_polyA_tail_check($seq),
        seed_region        => seed_region_check($seq),
        dicer_site         => dicer_site_check($seq),
        mirbase_motif      => mirbase_motif_check($seq),
        mismatch_profile   => mismatch_profile_check($seq),
        precursor_loop     => precursor_loop_check($seq),
    );
    my $total = 0; $total += $_ for values %scores;
    $scores{total} = (scalar(keys %scores) ? $total / scalar(keys %scores) : 0);
    return \%scores;
}
sub length_check {
    my $s = shift;
    my $l = length($s);
    return ($l >= 18 && $l <= 25) ? 1 : ($l > 25 ? 25/$l : $l/18);
}
sub hairpin_check {
    my $s = shift;
    my $count = 0;
    for my $i (0..length($s)-6) {
        my $sub = substr($s, $i, 6);
        my $revcomp = reverse $sub;
        $revcomp =~ tr/ATGC/TACG/;
        $count++ if index($s, $revcomp, $i+6) != -1;
    }
    return $count > 0 ? 1 : 0;
}
sub au_rich_check {
    my $s = shift;
    my $au = ($s =~ tr/AT//);
    return length($s) ? ($au / length($s) > 0.6 ? 1 : 0) : 0;
}
sub no_long_orf_check {
    my $s = shift;
    my $max_orf = 0;
    for my $frame (0..2) {
        for (my $i=$frame; $i+3 <= length($s); $i+=3) {
            my $codon = substr($s, $i, 3);
            if ($codon eq 'ATG') {
                my $j = $i + 3;
                while ($j+3 <= length($s)) {
                    my $stop = substr($s, $j, 3);
                    last if $stop =~ /TAA|TAG|TGA/;
                    $j += 3;
                }
                my $orf_len = $j - $i;
                $max_orf = $orf_len if $orf_len > $max_orf;
            }
        }
    }
    return $max_orf <= 15 ? 1 : 0;
}
sub no_polyA_tail_check {
    my $s = shift;
    return ($s =~ /A{5,}$/) ? 0 : 1;
}
sub seed_region_check {
    my $s = shift;
    my $seed = substr($s,1,7);
    my $au = ($seed =~ tr/AT//);
    return $au / 7 > 0.7 ? 1 : 0;
}
sub dicer_site_check {
    my $s = shift;
    return 0 unless defined $s && length($s) >= 15;
    my $site = substr($s, 9, 6);
    return ($site =~ /GU|AG/) ? 1 : 0;
}
sub mirbase_motif_check {
    my $s = shift;
    return ($s =~ /UGA|UAG|UAA/) ? 1 : 0;
}
sub mismatch_profile_check {
    my $s = shift;
    my $count = () = $s =~ /(.)\1{3,}/g;
    return $count == 0 ? 1 : 0;
}
sub precursor_loop_check {
    my $s = shift;
    my $mid = substr($s, int(length($s)/3), int(length($s)/3));
    return ($mid =~ /TTT|AAA/) ? 1 : 0;
}
1;
