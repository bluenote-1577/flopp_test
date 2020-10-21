#!/usr/bin/perl -w
use strict;

my $usage =  "usage: perl mk_index.pl whole_pathway_datebase_dir/ output_dir/";
my $index; 
my $index_new;
my $line;
my $path="";
my $gi;
my $number=0;
my @array;
my $j=0;
my @len;
my $lensum=0;
my $percentage;
my $per;
my @pa=();
my @name=();
my $i;
if(@ARGV != 2) 
{
	print $usage."\n"; exit(0);
}

if($ARGV[1]=~/\/$/)
{
	$index=$ARGV[1]."index";
	$percentage=$ARGV[1]."percentage.txt";
}
else
{
    $index=$ARGV[1]."/"."index";
	$percentage=$ARGV[1]."/"."percentage.txt";
}

for($i=1;$i<100;$i++)
{
	$len[$i] =0;
	$name[$i]="";
	$pa[$i]="";
}

if($ARGV[0] =~ /^(.+\/)[^\/]+$/)
{
	$path=$1;
}
if(-d "$path/reftemp")
{
	system("rm -r $path/reftemp");
}
system("mkdir $path/reftemp");
open(IN,"<$ARGV[0]")or die "error reading $ARGV[0] for reading\n";
while(<IN>)
{
	chomp;
	if($_ =~ /^\>(.+)/)
	{		
		if($j>1)
		{
			close OUT;
		}
		$j++;
		open(OUT,">$path/reftemp/$j.fa")or die "error reading $ARGV[0] for reading\n";
		print OUT"$_\n";
		$name[$j]=$1;
		#print"$name[$j]\n";
		$pa[$j]=$path."reftemp".'/'."$j.fa";		
	}
	else
	{
		$len[$j] += length($_);
		$lensum += length($_);
		$_ =uc($_);
		print OUT"$_\n";	
	}

}
close OUT;

open(OUT,">$percentage");
for($i=1;$i<=$j;$i++)
{
	$per=$len[$i]/$lensum;
	print OUT"$name[$i]\t";
	print OUT sprintf"%0.4f\n",$per;
}
close OUT;
open(OUT,">$index");
for($i=1;$i<=$j;$i++)
{
	print OUT"$len[$i]\t$name[$i]\t$pa[$i]\n";
}
close OUT;
