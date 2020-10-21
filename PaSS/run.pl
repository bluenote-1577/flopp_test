#!/usr/bin/perl -w
use strict;
use threads;
my $usage =  "usage: perl run.pl example.fq example.blasr RS/sequel.";
my $line="";
my $temp=1;
my $i;
my @thread=();
my @arg;
if(@ARGV != 3) 
{
	print $usage."\n"; exit(0);
}

for($i=0;$i<13;$i++)
{
	@arg=($ARGV[0],$ARGV[1],$ARGV[2],$i);
	$thread[$i]=threads->new(\&deal,@arg);
}
for($i=0;$i<13;$i++)
{
	$thread[$i]->join();
}
sub deal{
my ($par1,$par2,$par3,$id)=@_;
if($id==0)
{
	system("perl script/eachread_unalign.pl $par2 > unalign.txt ");
}
elsif($id==1)
{
	system("perl script/errsize_config.pl $par2 > errsize.txt ");
}
elsif($id==2)
{
	system("perl script/onebin_config.pl $par2  > onebin.txt ");
}
elsif($id==3)
{
	system("perl script/readlen_qv.pl $par1 $par3 ");
}
}
open(OUT,">sim.config")or die "Can't open the output file: sim.config!\n";
open(IN,"<read_length.txt")or die "error reading read_length.txt for reading\n";
while(<IN>)
{
	
	print OUT "$_";
}
close IN;
open(IN,"<errsize.txt")or die "error reading errsize.txt for reading\n";
while(<IN>)
{
	print OUT "$_";
}
open(IN,"<unalign.txt")or die "error reading unalign.txt for reading\n";
while(<IN>)
{
	print OUT "$_";
}
close IN;

open(IN,"<onebin.txt")or die "error reading head.txt for reading\n";
while(<IN>)
{
	print OUT "$_";
}
close IN;


if($ARGV[2] eq 'RS')
{
	system("rm onehole.fq");
}
system("rm read_length.txt errsize.txt unalign.txt onebin.txt");
close OUT;

