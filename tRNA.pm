package tRNA;
use strict;
use warnings;
sub score {
    my ($seq) = @_;
    my %scores = (
        length_score         => length_check($seq),
        cca_tail             => cca_tail_check($seq),
        d_loop_motif         => d_loop_check($seq),
        tpsic_motif          => tpsic_check($seq),
        gc_content           => gc_content_check($seq),
        anticodon_loop       => anticodon_loop_check($seq),
        discriminator_base   => discriminator_base_check($seq),
        conserved_bases      => conserved_bases_check($seq),
        acceptor_stem        => acceptor_stem_check($seq),
        tertiary_interaction => tertiary_interaction_check($seq),
    );
    my $total = 0; $total += $_ for values %scores;
    $scores{total} = (scalar(keys %scores) ? $total / scalar(keys %scores) : 0);
    return \%scores;
}
sub length_check {
    my $s = shift;
    my $len = length($s);
    return ($len >= 70 && $len <= 95) ? 1 : ($len > 95 ? 95/$len : $len/70);
}
sub cca_tail_check {
    my $s = shift;
    return ($s =~ /CCA$/) ? 1 : 0;
}
sub d_loop_check {
    my $s = shift;
    return ($s =~ /GGU/) ? 1 : 0;
}
sub tpsic_check {
    my $s = shift;
    return ($s =~ /TTC/) ? 1 : 0;
}
sub gc_content_check {
    my $s = shift;
    my $gc = ($s =~ tr/GC//);
    my $gc_ratio = length($s) ? $gc / length($s) : 0;
    return 1 - abs(0.5 - $gc_ratio)/0.2 > 0 ? 1 - abs(0.5 - $gc_ratio)/0.2 : 0;
}
sub anticodon_loop_check {
    my $s = shift;
    return ($s =~ /GAA|GUA/) ? 1 : 0;
}
sub discriminator_base_check {
    my $s = shift;
    my $len = length($s);
    my $base = substr($s, $len-2, 1);
    return ($base =~ /[AG]/) ? 1 : 0;
}
sub conserved_bases_check {
    my $s = shift;
    return ($s =~ /T.{35}G/) ? 1 : 0;
}
sub acceptor_stem_check {
    my $s = shift;
    my $start = substr($s, 0, 7);
    my $end = substr($s, -7);
    my $paired = 0;
    for (my $i=0; $i<7; $i++) {
        my $b1 = substr($start, $i,1);
        my $b2 = substr($end, 6-$i,1);
        $paired++ if (($b1 eq 'A' && $b2 eq 'T') || ($b1 eq 'T' && $b2 eq 'A') || ($b1 eq 'G' && $b2 eq 'C') || ($b1 eq 'C' && $b2 eq 'G'));
    }
    return $paired / 7;
}
sub tertiary_interaction_check {
    my $s = shift;
    return ($s =~ /G.{1}C.{4}G/) ? 1 : 0;
}
1;
