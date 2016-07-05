#!/usr/bin/perl

use Storable qw(dclone);

# Read annotation file, and store gene information in two hashes according to their boundaries.
my $byFivePrime = {};
my $byThreePrime = {};

my $annotationFile = $ARGV[0];
my $regionLimit = $ARGV[1];
open(ANNOTATION, "<", $annotationFile);

while(my $line = <ANNOTATION>) {
    chomp($line);
    
    my @data = split("\t", $line);
    my $chr = $data[0];
    my $start = $data[1];
    my $end = $data[2];
    my $strand = $data[5];
    my %info = (
        "gene" => $data[3],
        "chr" => $chr,
        "start" => $start,
        "end" => $end,
        "strand" => $strand
    );
    
    if(!exists($byFivePrime->{$chr})) {
        $byFivePrime->{$chr} = {};
        $byThreePrime->{$chr} = {};
    }
    
    my $fiveHash = dclone(\%info);
    my $threeHash = dclone(\%info);        
    
    if($strand eq "+") {
        $fiveHash->{"type"} = "promoter";
        $threeHash->{"type"} = "downstream";
    } elsif($strand eq "-") {
        $fiveHash->{"type"} = "downstream";
        $threeHash->{"type"} = "promoter";
    } else {
        $fiveHash->{"type"} = "unknown";
        $threeHash->{"type"} = "unknown";    
    }
    
    $byFivePrime->{$chr}->{$start} = $fiveHash;
    $byThreePrime->{$chr}->{$end} = $threeHash;        
}

close(ANNOTATION);

sub max {
    return $_[0] >= $_[1] ? $_[0] : $_[1];
}

sub min {
    return $_[0] <= $_[1] ? $_[0] : $_[1];
}

# Loop over main file on stdin
while(my $line = <STDIN>) {
    chomp($line);
    my @data = split("\t", $line);
    
    my $chr = $data[0];
    my $start = $data[1];
    my $end = $data[2];
    
    my $midPoint=int($start + (($end - $start)/2));
    my $leftLimitPoint = $start + $regionLimit;
    my $rightLimitPoint = max(0, $end - $regionLimit);
    
    my $fivePrimeType = "unknown";
    my $fivePrimeGene = "unknown";
    my $fivePrimeStrand = "unknown";
    # What's at the upstream limit?
    # Wrapped in an if, since start of chromosome might not have an associated element.
    if(exists($byThreePrime->{$chr}->{$start})) {
        $fivePrimeType = $byThreePrime->{$chr}->{$start}->{"type"};
        $fivePrimeGene = $byThreePrime->{$chr}->{$start}->{"gene"};
        $fivePrimeStrand = $byThreePrime->{$chr}->{$start}->{"strand"};
    }

    # Is the interval long enough to have a middle/"true intergenic" part?
    if($leftLimitPoint >= $midPoint - 1) {
        print $chr . "\t" . $start . "\t" . $midPoint . "\t" . $fivePrimeType . "\t" . $fivePrimeGene . "\t" . $fivePrimeStrand . "\n";
    } else {
        print $chr . "\t" . $start . "\t" . $leftLimitPoint . "\t" . $fivePrimeType . "\t" . $fivePrimeGene . "\t" . $fivePrimeStrand . "\n";
        print $chr . "\t" . ($leftLimitPoint + 1) . "\t" . ($rightLimitPoint - 1) . "\t" . "intergenic" . "\t" . "intergenic" . "\t" . "." . "\n";
    }
    
    my $threePrimeType = "unknown";
    my $threePrimeGene = "unknown";
    my $threePrimeStrand = "unknown";
    # What's at the upstream limit?
    # Wrapped in an if, since start of chromosome might not have an associated element.
    if(exists($byFivePrime->{$chr}->{$end})) {
        $threePrimeType = $byFivePrime->{$chr}->{$end}->{"type"};
        $threePrimeGene = $byFivePrime->{$chr}->{$end}->{"gene"};
        $threePrimeStrand = $byFivePrime->{$chr}->{$end}->{"strand"};
    }
    if($rightLimitPoint <= $midPoint + 1) {
        print $chr . "\t" . ($midPoint + 1) . "\t" . $end . "\t" . $threePrimeType . "\t" . $threePrimeGene . "\t" . $threePrimeStrand . "\n";
    } else {
        print $chr . "\t" . ($leftLimitPoint + 1) .  "\t" . ($rightLimitPoint - 1) . "\t" . "intergenic" . "\t" . "intergenic" . "\t" . "." . "\n";
    }
}