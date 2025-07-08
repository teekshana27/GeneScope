package eRNA;
use strict;
use warnings;
sub score {
    my ($seq) = @_;
    my %scores = (
        length_score                 => length_check($seq),
        enhancer_motif              => enhancer_motif_check($seq),
        transcription_factor_binding => tf_binding_check($seq),
        chromatin_marks             => chromatin_marks_check($seq),
        no_polyA_tail               => no_polyA_tail_check($seq),
        gc_content                  => gc_content_check($seq),
        bidirectional_transcription => bidirectional_trans_check($seq),
        conserved_regions           => conserved_regions_check($seq),
        eRNA_stability              => erna_stability_check($seq),
        secondary_structure         => secondary_structure_check($seq),
    );
    my $total = 0; $total += $_ for values %scores;
    $scores{total} = (scalar(keys %scores) ? $total / scalar(keys %scores) : 0);
    return \%scores;
}
sub length_check {
    my $s = shift;
    my $l = length($s);
    return ($l >= 200 && $l <= 2000) ? 1 : ($l > 2000 ? 2000/$l : $l/200);
}
sub enhancer_motif_check {
    my $s = shift;
    my $count = () = $s =~ /CAAT|GGGCGG/g;
    return $count > 0 ? 1 : 0;
}
sub tf_binding_check {
    my $s = shift;
    my $count = () = $s =~ /TGACTCA|GGGACTTTCC/g;
    return $count > 0 ? 1 : 0;
}
sub chromatin_marks_check {
    my $s = shift;
    my $count = () = $s =~ /CG{2,}/g;
    return $count > 3 ? 1 : 0;
}
sub no_polyA_tail_check {
    my $s = shift;
    return ($s =~ /A{10,}$/) ? 0 : 1;
}
sub gc_content_check {
    my $s = shift;
    my $gc = ($s =~ tr/GC//);
    return length($s) ? $gc / length($s) : 0;
}
sub bidirectional_trans_check {
    my $s = shift;
    my $count = 0;
    for my $i (0..length($s)-6) {
        my $sub = substr($s, $i, 6);
        my $revcomp = reverse $sub;
        $revcomp =~ tr/ATGC/TACG/;
        $count++ if index($s, $revcomp, $i+6) != -1;
    }
    return $count > 1 ? 1 : 0;
}
sub conserved_regions_check {
    my $s = shift;
    return ($s =~ /GGG/) ? 1 : 0;
}
sub erna_stability_check {
    my $s = shift;
    return ($s =~ /A{10,}|T{10,}/) ? 0 : 1;
}
sub secondary_structure_check {
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
1;
