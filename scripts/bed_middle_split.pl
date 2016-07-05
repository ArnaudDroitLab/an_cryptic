#!/usr/bin/perl

# Loop over input file, splitting all bed entries right down the middle.
while(my $line=<STDIN>) {
    chomp($line);
    my @data = split("\t", $line);
    
    my $midPoint=($data[1] + (($data[2] - $data[1])/2))
    print $data[0] . "\t" . $data[1] . "\t" . $midPoint . "\n";
    print $data[0] . "\t" . ($midPoint + 1) . "\t" . $data[2] . "\n";
}