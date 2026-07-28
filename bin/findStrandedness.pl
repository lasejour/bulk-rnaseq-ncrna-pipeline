#!/usr/bin/perl
use strict;
use warnings;
use feature qw(say);

# genomes
my %FeatureCounts_Strandness = (
    "ISR" => 2,
    "ISF" => 1,
    "IU" => 0,

    "OSR" => 2,
    "OSF" => 1,
    "OU" => 0,

    "MSR" => 2,
    "MSF" => 1,
    "MU" => 0,

    "SR" => 2,
    "SF" => 1,
    "U"  => 0
);

# Check for 1 file as a command line parameter
my $filename = shift or die "Usage: $0 <salmon_quant.log>\n";

# Check to see if the presented file can be opened
open(my $fh, '<:encoding(UTF-8)', $filename)
   or die "Could not open file '$filename' $!\n";

# Loop through the file. Extract the library type
my $value;
while(my $line = <$fh>) {
    chomp $line;

    if($line =~ /(library type as\s)(.+)/) {
        my $lib_type = $2;
        $value = $FeatureCounts_Strandness{$lib_type};
    } else {
        next;
    }

}

if ($value) {
   print $value;
} else {
   print 0;
}


close $fh;
