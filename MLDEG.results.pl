#!/usr/bin/perl

use strict;
use warnings;

unless ( @ARGV == 3 )
{
	die "Usage: perl $0 <Training Set> <Prediction Result> <Output File Name>\n";
}

my %result = ();

my $outFile = $ARGV[2].".final.result";

trParser();
prParser();
getResult();

##################################################

sub trParser
{
	open ( TR, $ARGV[0] ) or die "$!\n";

	while ( <TR> )
	{
		if ( $_ =~ /POS/ )
		{
			chomp $_;
			my @_token = split ( /\t/, $_ );
			push @{$result{$_token[0]}}, $_token[2], 1;
		}
	}

	close ( TR );
}

sub prParser
{
	open ( PR, $ARGV[1] ) or die "$!\n";

	while ( <PR> )
	{
		if ( $_ =~ /POS/ )
		{
			chomp $_;
			my @_token = split ( /\t/, $_ );
			push @{$result{$_token[0]}}, $_token[1], $_token[3];
		}
	}

	close ( PR );
}

sub getResult
{
	open ( OUT, ">$outFile" ) or die "$!\n";
	foreach my $_gene ( sort { $result{$b}[0] cmp $result{$a}[0] or $result{$b}[1] <=> $result{$a}[1] } keys %result )
	{
		print OUT $_gene, "\t", $result{$_gene}[0], "\t", $result{$_gene}[1], "\n";
	}
	close ( OUT );
}
