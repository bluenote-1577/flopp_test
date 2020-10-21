use strict;use warnings;

my $line;
my $query;
my $queryline;
my $sbjctline;
my $sbjct;
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
my $count1=0;
my $count2=0;
my $count3=0;
my $x=0;
my $linenumber=0;
my $ins;
my $del;
my $sub;
my %name;
my $y=0;

my $hsp=0;
my %length1=();
my %lengthnum1=();
my @position1=();
my @pos1=();
my $leng1;
my $lengthnumsum1=0;

my %length2=();
my %lengthnum2=();
my @position2=();
my @pos2=();
my $leng2;
my $lengthnumsum2=0;

my %length3=();
my %lengthnum3=();
my @position3=();
my @pos3=();
my $leng3;
my $lengthnumsum3=0;
my $key1;
my $key2;
my $key3;
my $j1;
my $j2;
my $j3;
my $pro;
my $quname;
my @errorlength;
my $prosum=0;
my @strand=();
my $j;
####12 types of substitution
my $numAT=0;
my $numAG=0;
my $numAC=0;
my $numTA=0;
my $numTG=0;
my $numTC=0;
my $numGA=0;
my $numGT=0;
my $numGC=0;
my $numCA=0;
my $numCT=0;
my $numCG=0;

my $prosubAT;
my $prosubAG;
my $prosubAC;
my $prosubTA;
my $prosubTG;
my $prosubTC;
my $prosubGA;
my $prosubGT;
my $prosubGC;
my $prosubCA;
my $prosubCT;
my $prosubCG;

my $subtotal=0;
####insdetail_config
my @keys1=();
my $order;
my @stat=();
my $key;
my $per=0;
my $count4=0;

for($i=0;$i<64;$i++)
{
    for($j=0;$j<4;$j++)
    {
        $stat[$i][$j]=0;
    }
}
sub convert
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
            if($y!=1)
            {
                $j1=0;
                $count1=0;
                @position1=();
                @pos1=();
                %length1=();
			 	$j2=0;
		        $count2=0;
		        @position2=();
		        @pos2=();
		        %length2=();

			 	$j3=0;
                $count3=0;
                @position3=();
               	@pos3=();
			   	%length3=();

                $query = uc($query);
                @queryarray = split("",$query);
                $sbjct = uc($sbjct);
                @sbjctarray = split("",$sbjct);
                for($i=0;$i<@queryarray;$i++)
                {
                    if($sbjctarray[$i] =~ /[ATGC\-]/)
                    {
                        
                        if($queryarray[$i] ne $sbjctarray[$i] && $queryarray[$i] ne "-" && $sbjctarray[$i] ne "-")
                        {
                            $subtotal++;  
                            if($queryarray[$i] eq "A")
                            {
                                if($sbjctarray[$i] eq "T")
                                {
                                    $numTA++;
                                }
                                elsif($sbjctarray[$i] eq "G")
                                {
                                   $numGA++;
                                }
                                elsif($sbjctarray[$i] eq "C")
                                {
                                    $numCA++;
                                }
                            }
                            elsif($queryarray[$i] eq "T")
                            {
                                if($sbjctarray[$i] eq "A")
                                {
                                    $numAT++;
                                }
                                elsif($sbjctarray[$i] eq "G")
                                {
                                   $numGT++;
                                }   
                                elsif($sbjctarray[$i] eq "C")
                                {   
                                    $numCT++;
                                }
                            }
                            elsif($queryarray[$i] eq "G")
                            {
                                if($sbjctarray[$i] eq "A")
                                {
                                    $numAG++;
                                }
                                elsif($sbjctarray[$i] eq "T")
                                {
                                    $numTG++;
                                }
                                elsif($sbjctarray[$i] eq "C")
                                {
                                    $numCG++;
                                }
                            }
                            elsif($queryarray[$i] eq "C")
                            {
                                if($sbjctarray[$i] eq "A")
                                {
                                    $numAC++;
                                }
                                elsif($sbjctarray[$i] eq "T")
                                {
                                    $numTC++;
                                }
                                elsif($sbjctarray[$i] eq "G")
                                {   
                                    $numGC++;
                                }
                            }   
                       }
                    }
                }

                for($i=0;$i<@sbjctarray;$i++)
                {
                    if($sbjctarray[$i] eq "-")
                    {
                        $position1[$j1]=$i;
                        $count1++;
                        $j1++;
                    }
                    if($sbjctarray[$i] =~ /[ATGC\-]/)
                    {
                        if($queryarray[$i] eq "-")
                        {
                            $position2[$j2]=$i;
                            $count2++;
                            $j2++;
                        }
                    }
                    if($sbjctarray[$i] =~ /[ATGC\-]/)
                    {
                        if(($queryarray[$i] ne $sbjctarray[$i]) && ($sbjctarray[$i] ne "-") && ($queryarray[$i] ne "-"))
                        {
                            $position3[$j3]=$i;
                            $count3++;
                            $j3++;
                        }
                    }

                }
                if($count1!=0)
                {
                    $pos1[0]=$position1[0];
                    $length1{$pos1[0]}=1;
                    for($i=1;$i<$count1;$i++)
                    {
                        if($position1[$i]==($position1[$i-1]+1))
                        {
                            $pos1[$i]=$pos1[$i-1];
                        }
                        else
                        {
                            $pos1[$i]=$position1[$i];
                        }
                        if(exists $length1{$pos1[$i]})
                        {
                            $length1{$pos1[$i]}++;
                        }
                        else
                        {
                            $length1{$pos1[$i]}=1;
                        }
                    }
                    foreach $key1(keys %length1)
                    {
                        $leng1=$length1{$key1};
                        if(exists $lengthnum1{$leng1})
                        {
                            $lengthnum1{$leng1}++;
                        }
                        else
                        {
                            $lengthnum1{$leng1}=1;
                        }
                        $lengthnumsum1++;
                    }
                    @keys1=keys(%length1);
                    for($j=0;$j<@keys1;$j++)
                    {
                        if(($keys1[$j]+$length1{$keys1[$j]}+2)<@sbjctarray)
                        {
                            $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2];    
                            if($key =~ /[^ATGC\-]+/)                 #7.9 add
                            {
                                next;
                            }
                            elsif($key =~ /^[ATGC]+$/)
                            {
                                ($order)=convert($key);
                                for($i=0;$i<$length1{$keys1[$j]};$i++)
                                {
                                    if($queryarray[$keys1[$j]+$i] eq 'A')
                                    {
                                        $stat[$order][0]++;
                                    }
                                    elsif($queryarray[$keys1[$j]+$i] eq 'T')
                                    {
                                        $stat[$order][1]++;
                                    }
                                    elsif($queryarray[$keys1[$j]+$i] eq 'G')
                                    {
                                        $stat[$order][2]++;
                                    }
                                    elsif($queryarray[$keys1[$j]+$i] eq "C")
                                    {
                                        $stat[$order][3]++;
                                    }
                                }
                            }
                            else
                            {
                                if($sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1] eq '-')
                                {   
                                    if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1})<@sbjctarray)
                                    {
                                        $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}];
                                        if($key =~ /[^ATGC\-]+/)
                                        {
                                            next;
                                        }
                                        elsif($key =~  /^[ATGC]+$/)
                                        {
                                            ($order)=convert($key);
                                            for($i=0;$i<$length1{$keys1[$j]};$i++)
                                            {
                                                if($queryarray[$keys1[$j]+$i] eq 'A')
                                                {
                                                    $stat[$order][0]++;
                                                }
                                                elsif($queryarray[$keys1[$j]+$i] eq 'T')
                                                {
                                                    $stat[$order][1]++;
                                                }
                                                elsif($queryarray[$keys1[$j]+$i] eq 'G')
                                                {
                                                    $stat[$order][2]++;
                                                }
                                                elsif($queryarray[$keys1[$j]+$i] eq "C")
                                                {
                                                    $stat[$order][3]++;
                                                }
                                            }
                                        }
                                        else
                                        {
                                            if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}})<@sbjctarray)
                                            {
                                                $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}}];
                                                if($key =~ /[^ATGC\-]+/)
                                                {
                                                    next;
                                                }
                                                elsif($key =~  /^[ATGC]+$/)
                                                {
                                                    ($order)=convert($key);
                                                    for($i=0;$i<$length1{$keys1[$j]};$i++)
                                                    {
                                                        if($queryarray[$keys1[$j]+$i] eq 'A')
                                                        {
                                                            $stat[$order][0]++;
                                                        }
                                                        elsif($queryarray[$keys1[$j]+$i] eq 'T')
                                                        {
                                                            $stat[$order][1]++;
                                                        }
                                                        elsif($queryarray[$keys1[$j]+$i] eq 'G')
                                                        {
                                                            $stat[$order][2]++;                                  
                                                        }   
                                                        elsif($queryarray[$keys1[$j]+$i] eq "C")
                                                        {
                                                            $stat[$order][3]++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+2})<@sbjctarray)
                                    {
                                        $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+2}];
                                        if($key =~ /[^ATGC\-]+/)
                                        {
                                            next;
                                        }
                                        elsif($key =~  /^[ATGC]+$/)
                                        {
                                            ($order)=convert($key);
                                            for($i=0;$i<$length1{$keys1[$j]};$i++)
                                            {
                                                if($queryarray[$keys1[$j]+$i] eq 'A')
                                                {
                                                    $stat[$order][0]++;
                                                }
                                                elsif($queryarray[$keys1[$j]+$i] eq 'T')
                                                {
                                                    $stat[$order][1]++;
                                                }
                                                elsif($queryarray[$keys1[$j]+$i] eq 'G')
                                                {
                                                    $stat[$order][2]++;
                                                }
                                                elsif($queryarray[$keys1[$j]+$i] eq "C")
                                                {
                                                    $stat[$order][3]++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if($count2!=0)
                {
                    $pos2[0]=$position2[0];
                    $length2{$pos2[0]}=1;
                    for($i=1;$i<$count2;$i++)
                    {
                        if($position2[$i]==($position2[$i-1]+1))
                        {
                            $pos2[$i]=$pos2[$i-1];
                        }
                        else
                        {
                            $pos2[$i]=$position2[$i];
                        }
                        if(exists $length2{$pos2[$i]})
                        {
                            $length2{$pos2[$i]}++;
                        }
                        else
                        {
                            $length2{$pos2[$i]}=1;
                        }
                    }
                    foreach $key2(keys %length2)
                    {
                        $leng2=$length2{$key2};
                        if(exists $lengthnum2{$leng2})
                        {
                            $lengthnum2{$leng2}++;
                        }
                        else
                        {
                            $lengthnum2{$leng2}=1;
                        }
                        $lengthnumsum2++;
                    }
                }
                if($count3!=0)
                {
                    $pos3[0]=$position3[0];
                    $length3{$pos3[0]}=1;
                    for($i=1;$i<$count3;$i++)
                    {
                        if($position3[$i]==($position3[$i-1]+1))
                        {
                            $pos3[$i]=$pos3[$i-1];
                        }
                        else
                        {
                            $pos3[$i]=$position3[$i];
                        }
                        if(exists $length3{$pos3[$i]})
                        {
                            $length3{$pos3[$i]}++;
                        }
                        else
                        {
                            $length3{$pos3[$i]}=1;
                        }
                    }
                    foreach $key3(keys %length3)
                    {
                        $leng3=$length3{$key3};
                        if(exists $lengthnum3{$leng3})
                        {
                            $lengthnum3{$leng3}++;
                        }
                        else
                        {
                            $lengthnum3{$leng3}=1;
                        }
                        $lengthnumsum3++;
                    }
                }
            }
            $query="";
            $sbjct="";
            @queryarray=();
            @sbjctarray=();
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
    if($line =~ /QueryRange:\s+(\d+).+\s+(\d+)\s+of/)
    {
        $qustart[$k]=$1;
        $quend[$k]=$2;
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
        $query=$query.$queryline;    
        next;
    }
    if(($line =~ /\d+\s+([\w\-]+)$/) && (($x-$linenumber)%4==3) && ($temp==0))
    {
        $sbjctline=$1;
        $sbjct=$sbjct.$sbjctline;
        next;
    }
}
close IN;

$j1=0;
$count1=0;
%length1=();
@position1=();
@pos1=();
$j2=0;
$count2=0;
%length2=();
@position2=();
@pos2=();
$j3=0;
$count3=0;
%length3=();
@position3=();
@pos3=();

$query = uc($query);
@queryarray = split("",$query);
$sbjct = uc($sbjct);
@sbjctarray = split("",$sbjct);

for($i=0;$i<@queryarray;$i++)
{
    if($sbjctarray[$i] =~ /[ATGC\-]/)
    {
        if($queryarray[$i] ne $sbjctarray[$i] && $queryarray[$i] ne "-" && $sbjctarray[$i] ne "-")
        {
            $subtotal++;
            if($queryarray[$i] eq "A")
            {
                if($sbjctarray[$i] eq "T")
                {
                    $numTA++;                           #1.15 correct the bug should be numTA rather than numAT
                }
                elsif($sbjctarray[$i] eq "G")
                {
                    $numGA++;
                }
                elsif($sbjctarray[$i] eq "C")
                {
                    $numCA++;
                }
            }
            elsif($queryarray[$i] eq "T")
            {
                if($sbjctarray[$i] eq "A")
                {
                    $numAT++;
                }
                elsif($sbjctarray[$i] eq "G")
                {
                    $numGT++;
                }
                elsif($sbjctarray[$i] eq "C")
                {
                    $numCT++;
                }
            }
            elsif($queryarray[$i] eq "G")
            {
                if($sbjctarray[$i] eq "A")
                {
                    $numAG++;
                }
                elsif($sbjctarray[$i] eq "T")
                {
                   $numTG++;
                }
                elsif($sbjctarray[$i] eq "C")
                {
                    $numCG++;
                }   
            }
            elsif($queryarray[$i] eq "C")
            {
                if($sbjctarray[$i] eq "A")
                {
                    $numAC++;
                }
                elsif($sbjctarray[$i] eq "T")
                {
                    $numTC++;
                }
                 elsif($sbjctarray[$i] eq "G")
                {
                    $numGC++;
                }
            }
        }
    }
}
for($i=0;$i<@sbjctarray;$i++)
{
    if($sbjctarray[$i] eq "-")
    {
        $position1[$j1]=$i;
        $count1++;
        $j1++;
    }
    if($sbjctarray[$i] =~ /[ATGC\-]/)
    {
        if($queryarray[$i] eq "-")
        {
            $position2[$j2]=$i;
            $count2++;
            $j2++;
        }
    }
    if($sbjctarray[$i] =~ /[ATGC\-]/)
    {
        if(($queryarray[$i] ne $sbjctarray[$i]) && ($sbjctarray[$i] ne "-") && ($queryarray[$i] ne "-"))
        {
            $position3[$j3]=$i;
            $count3++;
            $j3++;
        }
    }
}
if($count1!=0)
{
    $pos1[0]=$position1[0];
    $length1{$pos1[0]}=1;
    for($i=1;$i<$count1;$i++)
    {
        if($position1[$i]==($position1[$i-1]+1))
        {
            $pos1[$i]=$pos1[$i-1];
        }
        else
        {
            $pos1[$i]=$position1[$i];
        }
        if(exists $length1{$pos1[$i]})
        {
            $length1{$pos1[$i]}++;
        }
        else
        {
            $length1{$pos1[$i]}=1;
        }
    }
    foreach $key1(keys %length1)
    {
        $leng1=$length1{$key1};
        if(exists $lengthnum1{$leng1})
        {
            $lengthnum1{$leng1}++;
        }
        else
        {
            $lengthnum1{$leng1}=1;
        }
        $lengthnumsum1++;
    }
    @keys1=keys(%length1);
    for($j=0;$j<@keys1;$j++)
    {
        if(($keys1[$j]+$length1{$keys1[$j]}+2)<@sbjctarray)
        {
            $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2];
            if($key =~ /[^ATGC\-]+/)
            {
                next;
            }
            elsif($key =~ /^[ATGC]+$/)
            {   
                ($order)=convert($key);
                for($i=0;$i<$length1{$keys1[$j]};$i++)
                {
                    if($queryarray[$keys1[$j]+$i] eq 'A')
                    {
                        $stat[$order][0]++;
                    }
                    elsif($queryarray[$keys1[$j]+$i] eq 'T')
                    {
                        $stat[$order][1]++;
                    }
                    elsif($queryarray[$keys1[$j]+$i] eq 'G')
                    {
                        $stat[$order][2]++;
                    }
                    elsif($queryarray[$keys1[$j]+$i] eq "C")
                    {
                        $stat[$order][3]++;
                    }
                }
            }
            else
            {
                if($sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1] eq '-')
                {   
                    if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1})<@sbjctarray)
                    {
                        $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}];
                        if($key =~ /[^ATGC\-]+/)
                        {
                            next;
                        }
                        elsif($key =~  /^[ATGC]+$/)
                        {
                            ($order)=convert($key);
                            for($i=0;$i<$length1{$keys1[$j]};$i++)
                            {
                                if($queryarray[$keys1[$j]+$i] eq 'A')
                                {
                                    $stat[$order][0]++;
                                }
                                elsif($queryarray[$keys1[$j]+$i] eq 'T')
                                {
                                    $stat[$order][1]++;
                                }
                                elsif($queryarray[$keys1[$j]+$i] eq 'G')
                                {
                                    $stat[$order][2]++;
                                }
                                elsif($queryarray[$keys1[$j]+$i] eq "C")
                                {
                                    $stat[$order][3]++;
                                }
                            }
                        }
                        else
                        {
                            if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}})<@sbjctarray)
                            {
                                $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}}];
                                if($key =~ /[^ATGC\-]+/)
                                {
                                    next;
                                }
                                elsif($key =~  /^[ATGC]+$/)
                                {
                                    ($order)=convert($key);
                                    for($i=0;$i<$length1{$keys1[$j]};$i++)
                                    {
                                        if($queryarray[$keys1[$j]+$i] eq 'A')
                                        {
                                            $stat[$order][0]++;
                                        }
                                        elsif($queryarray[$keys1[$j]+$i] eq 'T')
                                        {
                                            $stat[$order][1]++;
                                        }
                                        elsif($queryarray[$keys1[$j]+$i] eq 'G')
                                        {
                                            $stat[$order][2]++;                                  
                                        }   
                                        elsif($queryarray[$keys1[$j]+$i] eq "C")
                                        {
                                            $stat[$order][3]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+2})<@sbjctarray)
                    {
                        $key=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+2}];
                        if($key =~ /[^ATGC\-]+/)
                        {
                            next;
                        }
                        elsif($key =~  /^[ATGC]+$/)
                        {
                            ($order)=convert($key);
                            for($i=0;$i<$length1{$keys1[$j]};$i++)
                            {
                                if($queryarray[$keys1[$j]+$i] eq 'A')
                                {
                                    $stat[$order][0]++;
                                }
                                elsif($queryarray[$keys1[$j]+$i] eq 'T')
                                {
                                    $stat[$order][1]++;
                                }
                                elsif($queryarray[$keys1[$j]+$i] eq 'G')
                                {
                                    $stat[$order][2]++;
                                }
                                elsif($queryarray[$keys1[$j]+$i] eq "C")
                                {
                                    $stat[$order][3]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
if($count2!=0)
{
    $pos2[0]=$position2[0];
    $length2{$pos2[0]}=1;
    for($i=1;$i<$count2;$i++)
    {
        if($position2[$i]==($position2[$i-1]+1))
        {
            $pos2[$i]=$pos2[$i-1];
        }
        else
        {
            $pos2[$i]=$position2[$i];
        }
        if(exists $length2{$pos2[$i]})
        {
            $length2{$pos2[$i]}++;
        }
        else
        {
            $length2{$pos2[$i]}=1;
        }
    }
    foreach $key2(keys %length2)
    {
        $leng2=$length2{$key2};
        if(exists $lengthnum2{$leng2})
        {
            $lengthnum2{$leng2}++;
        }
        else
        {
            $lengthnum2{$leng2}=1;
        }
        $lengthnumsum2++;
    }
}
if($count3!=0)
{
    $pos3[0]=$position3[0];
    $length3{$pos3[0]}=1;
    for($i=1;$i<$count3;$i++)
    {
        if($position3[$i]==($position3[$i-1]+1))
        {
            $pos3[$i]=$pos3[$i-1];
        }
        else
        {
            $pos3[$i]=$position3[$i];
        }
        if(exists $length3{$pos3[$i]})
        {
             $length3{$pos3[$i]}++;
        }
        else
        {
            $length3{$pos3[$i]}=1;
        }
    }
    foreach $key3(keys %length3)
    {
        $leng3=$length3{$key3};
        if(exists $lengthnum3{$leng3})
        {
            $lengthnum3{$leng3}++;
        }
        else
        {
            $lengthnum3{$leng3}=1;
        }
        $lengthnumsum3++;
    }
}
###########
$prosubAT=$numAT/$subtotal;
$prosubAG=$numAG/$subtotal;
$prosubAC=$numAC/$subtotal;
$prosubTA=$numTA/$subtotal;
$prosubTG=$numTG/$subtotal;
$prosubTC=$numTC/$subtotal;
$prosubGA=$numGA/$subtotal;
$prosubGT=$numGT/$subtotal;
$prosubGC=$numGC/$subtotal;
$prosubCA=$numCA/$subtotal;
$prosubCT=$numCT/$subtotal;
$prosubCG=$numCG/$subtotal;

print"#AC:AG:AT:CG:CT:CA:GT:GA:GC:TA:TC:TG\n";
print"pacbio_sub_probability=";
printf("%0.5f:%0.5f:%0.5f:%0.5f:%0.5f:%0.5f:%0.5f:%0.5f:%0.5f:%0.5f:%0.5f:%0.5f\n",$prosubAC,$prosubAG,$prosubAT,$prosubCG,$prosubCT,$prosubCA,$prosubGT,$prosubGA,$prosubGC,$prosubTA,$prosubTC,$prosubTG);
##############

print"ins_detail=";
for($i=0;$i<@stat;$i++)
{
    $temp=0;
    $per=0;
    for($j=0;$j<4;$j++)
    {
        $temp += $stat[$i][$j];
    }
    if($temp >0)
    {
        for($j=0;$j<4;$j++)
        {
            $per += ($stat[$i][$j]/$temp);
            if($count4<255)
            {
                printf("%0.3f:",$per);
                $count4++;
            }
            else
            {
                printf("%0.3f\n",$per);
            }

        }
    }
    else
    {
        for($j=0;$j<4;$j++)
        {
            $per += 0.25;
            if($count4<255)
            {
                printf("%0.3f\n",$per);
                $count4++;
            }
            else
            {
                printf("%0.3f\n",$per);
            }
        }
    }
}

#################################

$i=0;
for($i=0;$i<5000;$i++)
{
    $errorlength[$i]=0;
}
foreach $leng1(sort {$a <=> $b} keys %lengthnum1)
{
    $pro=$lengthnum1{$leng1}/$lengthnumsum1;
    $prosum =$prosum+$pro;
    $errorlength[$leng1-1]=$prosum;
}
for($i=0;$i<5000;$i++)
{
    if($errorlength[$i]==0)
    {
        $errorlength[$i]=$errorlength[$i-1];
    }
}
print"inslength=";
for($i=0;$i<5000;$i++)
{
    if($i<4999)
    {
        printf("%0.4f:",$errorlength[$i]);
    }
    else
    {
        print"1.0000\n";
    }
}
for($i=0;$i<5000;$i++)
{
    $errorlength[$i]=0;
}
$prosum=0;
foreach $leng2(sort {$a <=> $b} keys %lengthnum2)
{
    $pro=$lengthnum2{$leng2}/$lengthnumsum2;
    $prosum =$prosum+$pro;
    $errorlength[$leng2-1]=$prosum;
}
for($i=0;$i<5000;$i++)
{
    if($errorlength[$i]==0)
    {
        $errorlength[$i]=$errorlength[$i-1];
    }
}
print"dellength=";
for($i=0;$i<5000;$i++)
{
    if($i<4999)
    {
        printf("%0.4f:",$errorlength[$i]);
    }
    else
    {
        print"1.0000\n";
    }
}
for($i=0;$i<5000;$i++)
{
    $errorlength[$i]=0;
}
$prosum=0;
foreach $leng3(sort {$a <=> $b} keys %lengthnum3)
{
    $pro=$lengthnum3{$leng3}/$lengthnumsum3;
    $prosum =$prosum+$pro;
    $errorlength[$leng3-1]=$prosum;
}
for($i=0;$i<5000;$i++)
{
    if($errorlength[$i]==0)
    {
        $errorlength[$i]=$errorlength[$i-1];
    }
}
print"sublength=";
for($i=0;$i<5000;$i++)
{
    if($i<4999)
    {
        printf("%0.4f:",$errorlength[$i]);
    }
    else
    {
        print"1.0000\n";
    }
}
