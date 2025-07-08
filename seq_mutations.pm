package seq_mutations;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(
    validate_sequence
    read_sequence_from_file
    find_variants
    show_alignment
    compare_sequences
);
sub validate_sequence {
    my ($seq) = @_;
    return $seq =~ /^[ACGTacgt]+$/;
}
sub read_sequence_from_file {
    my ($filename) = @_;
    open my $fh, '<', $filename or return undef;
    my $seq = '';
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^>/;
        $seq .= $line;
    }
    close $fh;
    return uc($seq);
}
sub find_variants {
    my ($ref, $sample) = @_;
    $ref = uc($ref);
    $sample = uc($sample);
    my $r_len = length($ref);
    my $s_len = length($sample);
    my $min_len = $r_len < $s_len ? $r_len : $s_len;
    my (@subs, @dels, @ins);
    my $match_count = 0;
    for my $i (0 .. $min_len - 1) {
        my $r = substr($ref, $i, 1);
        my $s = substr($sample, $i, 1);
        if ($r eq $s) {
            $match_count++;
        } else {
            push @subs, "Substitution at position " . ($i + 1) . ": $r → $s";
        }
    }
    if ($r_len > $s_len) {
        for my $i ($min_len .. $r_len - 1) {
            my $r = substr($ref, $i, 1);
            push @dels, "Deletion at position " . ($i + 1) . ": $r missing in sample";
        }
    } elsif ($s_len > $r_len) {
        for my $i ($min_len .. $s_len - 1) {
            my $s = substr($sample, $i, 1);
            push @ins, "Insertion at position " . ($i + 1) . ": extra $s in sample";
        }
    }
    my $identity = $min_len > 0 ? ($match_count / $min_len) * 100 : 0;
    return (\@subs, \@dels, \@ins, $identity);
}
sub show_alignment {
    my ($ref, $sample) = @_;
    my $r_len = length($ref);
    my $s_len = length($sample);
    my $len = $r_len > $s_len ? $r_len : $s_len;
    my ($line1, $line2, $line3) = ('', '', '');
    for my $i (0 .. $len - 1) {
        my $r = $i < $r_len ? substr($ref, $i, 1) : '-';
        my $s = $i < $s_len ? substr($sample, $i, 1) : '-';
        $line1 .= $r;
        $line2 .= ($r eq $s) ? '|' : ' ';
        $line3 .= $s;
    }
    return (
        "Alignment:\n" .
        "Ref:    $line1\n" .
        "        $line2\n" .
        "Sample: $line3\n"
    );
}
sub compare_sequences {
    my ($ref, $sample) = @_;
    unless (validate_sequence($ref) && validate_sequence($sample)) {
        return ("Invalid sequence(s).", "", "", 0);
    }
    my ($subs_ref, $dels_ref, $ins_ref, $identity) = find_variants($ref, $sample);
    my @subs = @$subs_ref;
    my @dels = @$dels_ref;
    my @ins  = @$ins_ref;
    my $mutation_report;
    if (@subs or @dels or @ins) {
        $mutation_report = join("\n", @subs, @dels, @ins) . "\n";
        $mutation_report .= "Total mutations: " . (scalar(@subs) + scalar(@dels) + scalar(@ins)) . "\n";
        $mutation_report .= sprintf("Percent identity: %.2f%%\n", $identity);
    } else {
        $mutation_report = "No mutations. Sequences match exactly.\nPercent identity: 100%\n";
    }
    my $alignment_output = show_alignment($ref, $sample);
    return ($mutation_report, $alignment_output, $identity);
}
1; 
