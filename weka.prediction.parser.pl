#!/usr/bin/perl

use warnings;
use strict;

unless ( @ARGV == 4 )
{
	die "Usage: perl $0 <Test Set> <Weka Prediction Result> <Output File Prefix> <Confidence Cutoff>\nperl weka.prediction.parser.pl GSE90116_testset.list GSE90116.LR.results GSE90116 0\n";
}

my $conf = $ARGV[-1];
my $output = $ARGV[-2];

my %test = ();
my %diff = ();

testSetParser();
wekaParser();

#--------------------------------------------------

sub testSetParser
{
	my $_index = 0;
	
	open ( IN, $ARGV[0] ) or die "$!\n";

	while ( <IN> )
	{
		chomp $_;

		my @_token = split ( /\t/, $_ );

		$_index ++;
		$test{$_index} = $_token[0];
		$diff{$_index} = $_token[1];
	}

	close ( IN );
}
sub wekaParser
{
	open ( IN, $ARGV[1] ) or die "$!\n";
	open ( OUT, ">$output\_predicted.result" ) or die "$!\n";

	while ( <IN> )
	{
		chomp $_;
		$_ =~ s/\s$//g;

		next if $_ !~ /\d$/;

		$_ =~ s/^\s+//g;
		$_ =~ s/\s+/\t/g;

		my @_token = split ( /\t/, $_ );

		if ( $_token[-1] >= $conf )
		{
			if ( $_token[2] =~ /1/ )
			{
				print OUT $test{$_token[0]}, "\t$diff{$_token[0]}\tPOS\t$_token[-1]\n";
#			print $_, "\n";
			}
			else
			{
				print OUT $test{$_token[0]}, "\t$diff{$_token[0]}\tNEG\t$_token[-1]\n";
			}
		}
	}

	close ( IN );
	close ( OUT );
}
