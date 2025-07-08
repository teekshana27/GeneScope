package gene_epigenetic_identifier;
use strict;
use warnings;
use Exporter 'import';
use File::Temp qw(tempfile);
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use IO::File;
our @EXPORT_OK = qw(analyze_gene_and_epigenetics);
sub apply_epigenetic_modification {
    my @modifications = (
        { type => "DNA Methylation", effect => "silenced", weight => 0.4 },
        { type => "Histone Acetylation", effect => "activated", weight => 0.3 },
        { type => "Histone Methylation (H3K27me3)", effect => "silenced", weight => 0.15 },
        { type => "Histone Methylation (H3K4me3)", effect => "activated", weight => 0.15 },
    );
    my $rand = rand();
    my $sum = 0;
    foreach my $mod (@modifications) {
        $sum += $mod->{weight};
        return $mod if $rand <= $sum;
    }
    return $modifications[0];  
}
sub identify_gene {
    my ($seq) = @_;
    print "\nRunning BLAST for Sequence:\n$seq\n\n";
    my $factory = Bio::Tools::Run::StandAloneBlast->new(
        -prog            => 'blastn',
        -db              => '/home/user/bioproject/local_gene_db',
        -expect          => 10,
        -max_target_seqs => 1,
    );
    my $seq_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna');
    my $blast_report = $factory->blastall($seq_obj);
    unless ($blast_report) {
        warn "BLAST failed for sequence: $seq\n";
        return "Unknown_Gene";
    }
    print "Waiting for BLAST results...\n";
    while (my $result = $blast_report->next_result) {
        while (my $hit = $result->next_hit) {
            print "Hit found: ", $hit->name, "\n";
            return $hit->name;
        }
        print "No hits found in this BLAST result.\n";
    }

    return "Unknown_Gene";
}
sub analyze_gene_and_epigenetics {
    my ($input) = @_;
    my $output_file = "epigenetics_report.csv";  # fixed output file name
    my $seqio;
    my $temp_fh;
    my $temp_filename;
    if (-e $input) {
        $seqio = Bio::SeqIO->new(-file => $input, -format => "fasta");
    }
    else {
        unless ($input =~ /^[ACGT]+$/i) {
            die "Error: Input string is not a valid DNA sequence or existing file: '$input'\n";
        }
        ($temp_fh, $temp_filename) = tempfile(SUFFIX => '.fasta');
        print $temp_fh ">InputSequence\n$input\n";
        close $temp_fh;
        $seqio = Bio::SeqIO->new(-file => $temp_filename, -format => "fasta");
    }
    my $out_fh = IO::File->new(">$output_file") or die "Cannot open output file '$output_file': $!\n";
    print $out_fh "Sequence_ID,Gene_Name,Epigenetic_Modification,Effect\n";
    while (my $seq_obj = $seqio->next_seq) {
        my $seq = $seq_obj->seq;
        my $seq_id = $seq_obj->id;
        print "\nProcessing sequence: $seq_id\n";
        my $gene_name = identify_gene($seq);
        print "Identified gene: $gene_name\n";
        my $mod = apply_epigenetic_modification();
        print "Epigenetic Modification Applied: $mod->{type}\n";
        print "$gene_name is $mod->{effect}!\n";
        print $out_fh "$seq_id,\"$gene_name\",$mod->{type},$mod->{effect}\n";
    }
    $out_fh->close;
    if ($temp_filename) {
        unlink $temp_filename or warn "Could not delete temp file '$temp_filename': $!\n";
    }
return "Result also saved to '$output_file'\n" 
}
1; 
 
