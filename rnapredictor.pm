package rnapredictor;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(predict_rna_types);
use mRNA;
use rRNA;
use tRNA;
use miRNA;
use eRNA;
sub predict_rna_types {
    my ($seq) = @_;
    open my $file_out, '>', 'rna_prediction_output.txt' or die "Cannot open output file: $!";
    $seq =~ s/\s+//g;
    $seq = uc($seq);
    if ($seq =~ /[^ATGC]/) {
        print STDOUT "Invalid characters found! Use only A,T,G,C.\n";
        print $file_out "Invalid characters found! Use only A,T,G,C.\n";
        close $file_out;
        return;
    }
    my %results;
    $results{mRNA}  = mRNA::score($seq);
    $results{rRNA}  = rRNA::score($seq);
    $results{tRNA}  = tRNA::score($seq);
    $results{miRNA} = miRNA::score($seq);
    $results{eRNA}  = eRNA::score($seq);
    _print_dual("[Summary Average Scores]\n", $file_out);
    foreach my $type (sort keys %results) {
        _print_dual(sprintf("%-5s : %.3f\n", $type, $results{$type}->{total}), $file_out);
    }
    my ($best) = sort { $results{$b}->{total} <=> $results{$a}->{total} } keys %results;
    _print_dual("\nLikely RNA type: $best\n", $file_out);

    _print_dual("\nDetailed Scores (1 = high confidence)\n", $file_out);
    foreach my $type (keys %results) {
        _print_dual("\n$type:\n", $file_out);
        foreach my $logic (sort keys %{$results{$type}}) {
            next if $logic eq 'total';
            _print_dual(sprintf("  %-30s : %.3f\n", $logic, $results{$type}{$logic}), $file_out);
        }
    }
    _print_dual("\nResults also saved to rna_prediction_output.txt\n", $file_out);
    close $file_out;
}
sub _print_dual {
    my ($text, $file_out) = @_;
    print STDOUT $text;
    print $file_out $text;
}
1; 
