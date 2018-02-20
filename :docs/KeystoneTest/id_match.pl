#!/usr/bin/perl

use strict;
use warnings;


# provide file and open it
my $in_otu = $ARGV[0] || die "Please provide table tab delimited";
my $map = $ARGV[1] || die "Please provide table tab delimited";

open(FILE, "<$in_otu") || die "File $in_otu doesn't exist!!";
open(FILE2, "<$map") || die "File $in_otu doesn't exist!!";

my $line1 = <FILE>;
my $line2 = <FILE>;

my @maps = <FILE2>;

chomp(@maps);

my @otus = split(/\t/, $line1);

my @otus_found = ();

foreach my $otu (@otus){
    
    my @otus_found_line = grep(/$otu\s+/, @maps);
    #push(@otus_found, @otus_found_line);
    chomp(@otus_found_line);
    print "@otus_found_line\n";
     
 
}


