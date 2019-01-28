#!/usr/bin/perl

use warnings;
use strict;
use Statistics::Basic qw(:all);

unless ( @ARGV == 6 )
{
	die "Usage: perl $0 <DEG Profile> <PCC stats> <Network Propagation> <log2FC> <p-value> <File Name>\n";
}

my $nTools = 0;
my $zero = 1e-200;
my $log2fc = $ARGV[-3];
my $pval = $ARGV[-2];
my $outFile = $ARGV[-1];
my $cut = 0;

my $FC_CUT = 0.5;
my $PVAL_CUT = 0.05;

my %data = ();
my %net = ();
my %pos = ();
my %neg = ();
my %netProp = ();

my %fisher = ();
my %expFC = ();
my %isSameDirection = ();
my %isDEG = ();
my %diff = ();

my $pos_count = 0;
my $header = "\@relation \'$outFile\'
\@attribute \'log2-fold-change\' real
\@attribute \'combined-p-value\' real
\@attribute \'condition-specific-degree\' real
\@attribute \'condition-specific-degree-ratio\' real
\@attribute \'correlation-coefficient-mean\' real
\@attribute \'correlation-coefficient-std\' real
\@attribute \'correlation-p-value-mean\' real
\@attribute \'correlation-p-value-std\' real
\@attribute \'network-propagation\' real
\@attribute \'class\' {1, 0}
\@data
";

my $training = "";
my $test = "";

setCut();
degParser();
setFisher();
netParser();
netPropParser();
setTraining();

#--------------------------------------------------

sub setCut
{
	my $_pi = 0;
	for ( my $_i = 0; $_i < 4; $_i ++ )
	{
		$_pi += log ( $pval )
	}
	$cut = -2 * $_pi;
}

sub degParser
{
	open ( DEG, $ARGV[0] ) or die "$!\n";

	my @_token = ();
	
	while ( <DEG> )
	{
		@_token = split ( /\t/, $_ );
		for ( my $_i = 1; $_i <= $#_token; $_i ++ )
		{
			push @{$data{$_token[0]}}, $_token[$_i];
		}
	}

	close ( DEG );
	$nTools = $#_token / 2;
}

sub setFisher
{
	foreach my $_gene ( sort keys %data )
	{
		my $_pi = 0;
		my $_fc = 0;
		my $_pval = 0;
		my $_f = 0;
		my $_d = 0;
		my $_score = 0;

		for ( my $_i = 0; $_i < $nTools * 2; $_i += 2 )
		{
			$_pval = $data{$_gene}[$_i];
			$_pi += log ( $_pval );

			if ( $_pval < $PVAL_CUT )
			{
				$_score ++;
			}
		}

		for ( my $_i = 1; $_i < $nTools * 2; $_i += 2 )
		{
			$_f = $data{$_gene}[$_i];
			$_fc += $_f;

			if ( $_f > 0 )
			{
				$_d ++;
			}
			elsif ( $_f < 0 )
			{
				$_d --;
			}
		}

		$fisher{$_gene} = -2 * $_pi;
		$_fc /= $nTools;
		$expFC{$_gene} = $_fc;
		$isDEG{$_gene} = $_score;
		if ( abs ( $_d ) == $nTools )
		{
			$isSameDirection{$_gene} = 1;
		}
		else
		{
			$isSameDirection{$_gene} = 0;
		}
		if ( $expFC{$_gene} > 0 )
		{
			$diff{$_gene} = "U";
		}
		elsif ( $expFC{$_gene} < 0 )
		{
			$diff{$_gene} = "D";
		}
		else
		{
			$diff{$_gene} = "E";
		}
	}
}

sub netParser
{
	open ( IN, $ARGV[1] ) or die "$!\n";

	while ( <IN> )
	{
		chomp $_;

		my @_token = split ( /\t/, $_ );

		for ( my $_i = 1; $_i <= $#_token; $_i ++ )
		{
			push @{$net{$_token[0]}}, $_token[$_i];
		}
	}

	close ( IN );
}

sub netPropParser
{
	open ( NP, $ARGV[2] ) or die "$!\n";

	while ( <NP> )
	{
		chomp $_;

		my @_token = split ( /\t/, $_ );

		$netProp{$_token[0]} = $_token[1];
	}

	close ( NP );
}

sub setTraining
{
	my $_lines = "";
	my $_test_lines = "";

	my $_rank = 0;
	open ( TRSET, ">$outFile\_trainingset.list" ) or die "$!\n";

	foreach my $_gene ( sort { $fisher{$b} <=> $fisher{$a} } keys %fisher )
	{
		next if !defined $net{$_gene};
		if ( $fisher{$_gene} > $cut && abs ( $expFC{$_gene} ) > $log2fc && $isSameDirection{$_gene} )
		{
			$pos_count ++;
			$pos{$_gene} = 1;

			my $_fc		= $expFC{$_gene};
			my $_pval	= $fisher{$_gene};
			my $_deg	= defined ( ${$net{$_gene}}[0] ) ? ${$net{$_gene}}[0] : 0;
			my $_deg_ratio	= defined ( ${$net{$_gene}}[1] ) ? ${$net{$_gene}}[1] : 0;
			my $_corr_mean	= defined ( ${$net{$_gene}}[2] ) ? ${$net{$_gene}}[2] : 0;
			my $_corr_std	= defined ( ${$net{$_gene}}[3] ) ? ${$net{$_gene}}[3] : 0;
			my $_p_mean	= defined ( ${$net{$_gene}}[4] ) ? ${$net{$_gene}}[4] : 1;
			my $_p_std	= defined ( ${$net{$_gene}}[5] ) ? ${$net{$_gene}}[5] : 1;
			my $_np		= defined ( $netProp{$_gene} ) ? $netProp{$_gene} : 0;
			$_lines .= "$_fc\,$_pval\,$_deg\,$_deg_ratio\,$_corr_mean\,$_corr_std\,$_p_mean\,$_p_std\,$_np\,1\n";
			print TRSET $_gene, "\tPOS\t$diff{$_gene}\n";
		}
	}

	open ( TSET, ">$outFile\_testset.list" ) or die "$!\n";

	my $_neg_count = 0;

	foreach my $_gene ( sort { $fisher{$a} <=> $fisher{$b} } keys %fisher )
	{
		next if defined ( $pos{$_gene} );

		my $_fc		= $expFC{$_gene};
		my $_pval	= $fisher{$_gene};
		my $_deg	= defined ( ${$net{$_gene}}[0] ) ? ${$net{$_gene}}[0] : 0;
		my $_deg_ratio	= defined ( ${$net{$_gene}}[1] ) ? ${$net{$_gene}}[1] : 0;
		my $_corr_mean	= defined ( ${$net{$_gene}}[2] ) ? ${$net{$_gene}}[2] : 0;
		my $_corr_std	= defined ( ${$net{$_gene}}[3] ) ? ${$net{$_gene}}[3] : 0;
		my $_p_mean	= defined ( ${$net{$_gene}}[4] ) ? ${$net{$_gene}}[4] : 1;
		my $_p_std	= defined ( ${$net{$_gene}}[5] ) ? ${$net{$_gene}}[5] : 1;
		my $_np		= defined ( $netProp{$_gene} ) ? $netProp{$_gene} : 0;

		if ( $isDEG{$_gene} > 0 )
		{
			$_test_lines .= "$_fc\,$_pval\,$_deg\,$_deg_ratio\,$_corr_mean\,$_corr_std\,$_p_mean\,$_p_std\,$_np\,\?\n";
			print TSET $_gene, "\t$diff{$_gene}\n";
		}
		else
		{
			if ( $_neg_count < $pos_count )
			{
				$_lines .= "$_fc\,$_pval\,$_deg\,$_deg_ratio\,$_corr_mean\,$_corr_std\,$_p_mean\,$_p_std\,$_np\,0\n";
				$_neg_count ++;
				print TRSET $_gene, "\tNEG\t$diff{$_gene}\n";			
			}
		}
	}

	close ( TRSET );
	close ( TSET );


	$training = $header . $_lines;
	$test = $header . $_test_lines;

	open ( TR, ">$outFile\_TR.arff" ) or die "$!\n";
	print TR $training;
	close ( TR );

	open ( TS, ">$outFile\_TS.arff" ) or die "$!\n";
	print TS $test;
	close ( TS );
}

sub log2
{
	return log ( $_[0] ) / log ( 2 );
}
