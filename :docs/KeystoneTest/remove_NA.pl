#!/usr/bin/perl

use strict;
use warnings;


# provide file and open it
my $in_otu = $ARGV[0] || die "Please provide table tab delimited";

open(FILE, "<$in_otu") || die "File $in_otu doesn't exist!!";


my @array = <FILE>;

my $element = NA;
#delete @array[$element];

@array = grep { $_ != $element} @array;




print @array;