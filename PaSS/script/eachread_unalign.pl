use strict;use warnings;
use POSIX;
my $line;
my @query=();
my $queryline;
my $sbjctline;
my @sbjct=();
my @queryarray;
my @sbjctarray;
my $i;
my @qustart=();
my @quend=();
my @refstart=();
my @refend=();
my $k;
my $temp=0;
my @mark;
my $x=0;
my $linenumber=0;

my %name;
my $y=0;
my $j;
my $hsp=0;
my $pro;
my $quname;
my @strand=();
my @stat=();
my $holequ="";
my $holeref="";
my $holename;
my %hole=();
my $key;
my $sum=0;
my $len;
my @holelen=();
my $p=0;
my @flag=();
my $bound;
my @holepos=();
my $templen;
my @tempqu;
my $fragment;
my $count=0;
my $match=0;
my $del=0;
my $numelse=0;
my @stat1=();
my @stat2=();
my @win1=();
my @win2=();
my $sum1=0;
my $sum2=0;
my $per=0;
my $head;
my $tail;
my @stathead=();
my @stattail=();

for($i=0;$i<300;$i++)
{
	$query[$i]="";
	$sbjct[$i]="";
}
sub convert2
{
	 my ($key)=@_;
     my @temp=split("",$key);
	 my $i;
	 my $j=0;
	 my $x;
	 for($i=0;$i<@temp;$i++)
	 {
	 	if($temp[$i] eq 'A')
	  	{
	    	$x=0;
	 	}
	    elsif($temp[$i] eq 'T')
	    {
	        $x=1;
	    }
	    elsif($temp[$i] eq 'G')
	   	{
	        $x=2;
	    }
	    elsif($temp[$i] eq 'C')
	    {
	        $x=3;
	    }
	    $j += $x*(4**$i);
	 }
	 return ($j);
}

sub convert
{
	my ($number)=@_;
    my $temp;
    my $qu="";
    $temp=$number % 4;
    $number=($number-$temp)/4;

    while($number != 0 || length($qu) <3)
    {
    	if($temp==0)
		{
        	$qu=$qu.'A';
        }
        elsif($temp==1)
        {
           	$qu=$qu.'T';
		}
        elsif($temp==2)
        {
            $qu=$qu.'G';
        }
        elsif($temp==3)
        {
            $qu=$qu.'C';
        }
        $temp=$number%4;
        $number=($number-$temp)/4;
   	}
   	return ($qu);
}

open(IN,"<$ARGV[0]")or die "error reading $ARGV[0] for reading\n";
while(<IN>)
{
    chomp;
    $x++;
    $line = $_;
    if($line =~ /(Query:.+)/)
    {
        $y++;
        $quname=$1;
        $temp=0;
        if(exists $name{$quname})
        {
            $k++;
			$name{$quname}++;      #12.28 change
        }
        else
        {   
			$name{$quname}=1;
			if($line =~ /(Query\:.+\/\d+\/).+/)
			{
				$holename=$1;
				if(exists $hole{$holename})
				{
					$hole{$holename}++;
					$p++;
				}
				else
				{
					$hole{$holename}=1;
					if($y!=1)
            		{

            			$stathead[$count]=$head/$holelen[$p];
					 	$stattail[$count]=$tail/$holelen[$p];
						
						for($i=0;$i<=$p;$i++)
						{
							$holequ=$holequ.$query[$i];
							$holeref=$holeref.$sbjct[$i];
						}
						$holequ=uc($holequ);
						$holeref=uc($holeref);
						@queryarray = split("",$holequ);
						@sbjctarray = split("",$holeref);
						for($i=0;$i<@sbjctarray;$i++)
      		          	{	
							if($sbjctarray[$i] !~ /[ATGC\-]/)
							{
								$numelse++;
							}
							else
							{
                				if($queryarray[$i] eq "-")
							    {   
							       	$del++;
							    }
							    if($queryarray[$i] eq $sbjctarray[$i])
							    {  	 
							    	$match++;
							    }
							}
						}
							
            			$stat1[$count]=(length($holequ)-$numelse-$match)/(length($holequ)-$numelse-$del);
						$sum1 += $stat1[$count];
						#$stat2[$count]=((length($holequ)-$numelse-$match)+($holelen[$p]-(length($holequ)-$numelse-$del))*$ARGV[1])/$holelen[$p];
						#$sum2 += $stat2[$count];
						$count++;

						$holequ="";
						$holeref="";
						@sbjctarray=();
						@queryarray=();
								
			
						@holelen=();
						$p=0;
						@flag=();
						@holepos=();
						$fragment="";
						@sbjct=();
						@query=();
						for($i=0;$i<300;$i++)
						{
							$query[$i]="";
							$sbjct[$i]="";
						}
						$match=0;
						$del=0;
						$numelse=0;
					}
				}
			}
			
			@qustart=();
		    @quend=();
		    @refstart=();
		    @refend=();
		    @mark=();
		    @strand=();
		    $mark[0]=1;
			$k=0;
        }
        next;
    }
	if($line =~ /Target strand: (\d)/)
	{
		$strand[$k]=$1;
	}
    if($line =~ /QueryRange:\s+(\d+).+\s+(\d+)\s+of\s+(\d+)/)
    {
        $qustart[$k]=$1;
        $quend[$k]=$2;
		
    	if($p==0)
		{
			$head=$1;
			$flag[$p]=$1;
			$holelen[$p]=$3;
		}
		else
		{
			$holelen[$p]=$holelen[$p-1]+$3;
			$flag[$p] =$1+$holelen[$p-1];
		}
		$tail =$3-$2;
	}
    if($line =~ /TargetRange:\s+(\d+).+\s+(\d+)\s+of/)
    {
        $refstart[$k]=$1;
        $refend[$k]=$2;
        $linenumber=$x;
        if($k>0)
        {
			if($strand[$k] == $strand[0])
			{
            	for($j=0;$j<(@qustart-1);$j++)
            	{   
                	if($mark[$j]==1)
                	{
                    	if((($qustart[$k] >= $qustart[$j]) && ($qustart[$k]<$quend[$j])) || (($quend[$k]>$qustart[$j])&&($quend[$k]<=$quend[$j])) || (($refstart[$k] >= $refstart[$j]) && ($refstart[$k] < $refend[$j])) || (($refend[$k] > $refstart[$j]) && ($refend[$k] <= $refend[$j])))
                    	{
                        	$temp++;
							if($temp>0)
							{
								last;
							}
                    	}
                	}
            	}
			}
			else
			{
				$temp=1;
			}
            if($temp==0)
            {
                $hsp++;
                $mark[$k]=1; 
            }
            else
            {
                $mark[$k]=0;
            }
        }
    }
    if(($line =~ /\d+\s+([\w\-]+)$/) && (($x-$linenumber)%4==1) && ($temp==0))
    {
        $queryline=$1;
        $query[$p]=$query[$p].$queryline;    
        next;
    }
    if(($line =~ /\d+\s+([\w\-]+)$/) && (($x-$linenumber)%4==3) && ($temp==0))
    {
        $sbjctline=$1;
        $sbjct[$p]=$sbjct[$p].$sbjctline;
        next;
    }
}
close IN;



##########unalign
$stathead[$count]=$head/$holelen[$p];
$stattail[$count]=$tail/$holelen[$p];


for($j=0;$j<=1000;$j++)
{
	$win1[$j]=0;
	$win2[$j]=0;
	
}

for($i=0;$i<$count;$i++)
{
	for($j=0;$j<1000;$j++)
	{
		if($stathead[$i]==0)
		{
			$win1[$j]++;
			last;
		}
		elsif(($stathead[$i]>$j*0.001) && ($stathead[$i]<=($j+1)*0.001))
		{
			$win1[$j+1]++;
			last;
		}
	}
	for($j=0;$j<1000;$j++)
	{
		if($stattail[$i]==0)
		{
			$win2[$j]++;
			last;
		}
		elsif(($stattail[$i]>$j*0.001) && ($stattail[$i]<=($j+1)*0.001))
		{
			$win2[$j+1]++;
			last;
		}
	}
}
print"unalign_head=";
for($j=0;$j<1000;$j++)
{
	$per += $win1[$j]/$count;
	if($j<999)
	{
		printf("%0.4f:",$per);
	}
	else
	{
		printf("1.0000\n");
	}
}
$per=0;
printf("unalign_tail=");
for($j=0;$j<1000;$j++)
{
	$per += $win2[$j]/$count;
	if($j<999)
	{
		printf("%0.4f:",$per);
	}
	else
	{
		printf("1.0000\n");
	}
}

#################

for($i=0;$i<=$p;$i++)
{
	$holequ=$holequ.$query[$i];
	$holeref=$holeref.$sbjct[$i];
}
$holequ=uc($holequ);
$holeref=uc($holeref);
@queryarray = split("",$holequ);
@sbjctarray = split("",$holeref);
for($i=0;$i<@sbjctarray;$i++)
{
	if($sbjctarray[$i] !~ /[ATGC\-]/)
	{
		$numelse++;
	}
	else
	{
		if($queryarray[$i] eq "-")
    	{
			$del++;
		}
		if($queryarray[$i] eq $sbjctarray[$i])
		{
			$match++;
		}
	}
}
$stat1[$count]=(length($holequ)-$numelse-$match)/(length($holequ)-$numelse-$del);
$sum1 += $stat1[$count];
#$stat2[$count]=((length($holequ)-$numelse-$match)+($holelen[$p]-(length($holequ)-$numelse-$del))*$ARGV[1])/$holelen[$p];
#$sum2 += $stat2[$count];
$count++;

for($j=0;$j<=1000;$j++)
{
	$win1[$j]=0;	
}
for($i=0;$i<$count;$i++)
{
	for($j=0;$j<100;$j++)
	{
		if(($stat1[$i]>$j*0.01) && ($stat1[$i]<=($j+1)*0.01))
		{
			$win1[$j]++;
			last;
		}
		#if(($stat2[$i]>$j*0.01) && ($stat2[$i]<=($j+1)*0.01))
		#{
	#		$win2[$j]++;
	#	}
	}
}
$per=0;
print"eachread=";
for($j=0;$j<50;$j++)
{
	$per += ($win1[$j]/$count);
	if($j<49)
	{
		printf("%0.5f:",$per);	
	}
	else
	{
		printf("1.000\n");
	}
}

