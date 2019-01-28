#!/usr/bin/perl

use warnings;
use strict;

unless ( @ARGV == 2 )
{
	die "Usage: perl $0 <Classifier Type(SMO, LR, J48)> <Input File Name>\n";
}

my %classifier = (
	"SMO"	=>	"java weka.classifiers.functions.SMO -C 1.0 -L 0.001 -P 1.0E-12 -N 0 -M -V -1 -W 1 -K \"weka.classifiers.functions.supportVector.PolyKernel -E 1.0 -C 250007\"",
	"LR"	=>	"java weka.classifiers.functions.Logistic -R 1.0E-8 -M -1 -num-decimal-places 4",
	"J48"	=>	"java weka.classifiers.trees.J48"
);

my %predictor = (
	"SMO"	=>	"java weka.classifiers.functions.SMO",
	"LR"	=>	"java weka.classifiers.functions.Logistic",
	"J48"	=>	"java weka.classifiers.trees.J48",
);

my $type = $ARGV[0];
my $input = $ARGV[1];
my $commandTR = $classifier{$type} . " -t $input\_TR.arff -d $input\.$type\.model > $input\.$type\.TR\.results"; 
my $commandTS = $predictor{$type} . " -l $input\.$type\.model -T $input\_TS.arff -p 0 > $input\.$type\.results";

print $commandTR, "\n";
`$commandTR`;

print $commandTS, "\n";
`$commandTS`;
