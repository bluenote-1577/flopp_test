use strict;use warnings;
my $i=0;
my @count=();
my @hole;
my %pass=();
my $holenum=-1;
my $len;
my $total=0;
my @sum;
my @percentage;
my @length;
my $lengthmin;
my $lengthmax;
my @array;
my @win;
my $window;
my $j=0;
my $percentage=0;
my $start=1;
my $end;
my $name;
my %temp=();

open(FASTQ,"<$ARGV[0]")or die "error reading $ARGV[1] for reading\n";
open(OUT,">read_length.txt") or die "Can't open the output file: read_length.txt!\n";
while(<FASTQ>)
{
    if($_ =~ /^\@(.+\/\d+)\/(\d+)\_(\d+)/)                        #zm:i not np:i
    {
		$name=$1; 
        $len=$3-$2;
        if(exists $temp{$name})
        {
			$temp{$name}++;
            $pass{$holenum}++;
            $count[$holenum]++;
        }
        else
        {
			$temp{$name}=1;
			$holenum++;
            $pass{$holenum}=1;
            $count[$holenum]=0;
        }
        $array[$holenum][$count[$holenum]]=$len;
    }
}

foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}==$start)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passone_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
close FASTQ;

$start=2;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}==$start)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passtwo_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=3;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}==$start)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passthree_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=4;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}==$start)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passfour_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=5;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}==$start)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passfive_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=6;
$end=10;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>=$start && $pass{$holenum}<=$end)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passsix_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}

$start=11;
$end=15;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>=$start && $pass{$holenum}<=$end)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passseven_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=16;
$end=20;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>=$start && $pass{$holenum}<=$end)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passeight_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=21;
$end=25;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>=$start && $pass{$holenum}<=$end)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passnine_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}

$start=26;
$end=30;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>=$start && $pass{$holenum}<=$end)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passten_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=31;
$end=35;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>=$start && $pass{$holenum}<=$end)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passeleven_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}
$start=36;
$end=40;
$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>=$start && $pass{$holenum}<=$end)
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passtwelve_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}

$start=40;

$total=0;
@length=();
@win=();
$percentage=0;
foreach $holenum(keys %pass)
{
    $length[$total]=$array[$holenum][0];
    if($pass{$holenum}>$start )
    {
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$length[$total])
            {
                $length[$total]=$array[$holenum][$j];
            }
        }
        $total++;
    }
}

$lengthmin=$length[0];
$lengthmax=$length[0];
for($i=0;$i<@length;$i++)
{
    if($length[$i]>$lengthmax)
    {
        $lengthmax=$length[$i];
    }
    if($length[$i]<$lengthmin)
    {
        $lengthmin=$length[$i];
    }
}
$window=$lengthmax;
for($j=0;$j<=$window;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<@length;$i++)
{
    for($j=1;$j<=$window;$j++)
    {
        if($length[$i]==$j)
        {
            $win[$j]++;
            last;
        }
    }
}

print OUT "passthirteen_length_rand=";
for($j=0;$j<=$window;$j++)
{
    $percentage += $win[$j]/$total;
    if($j<$window)
    {
        printf OUT ("%0.5f:",$percentage);
    }
    else
    {   
        printf OUT ("%0.5f\n",$percentage);
    }
}




my $lenmax;
my $lenmin;
my @ratio=();
my $per=0;
$total=0;


foreach $holenum(keys %pass)
{
    if($pass{$holenum}==2)
    {
        $lenmax=$array[$holenum][0];
        $lenmin=$array[$holenum][0];
        if($array[$holenum][1]>$lenmax)
        {
            $lenmax=$array[$holenum][1];
        }
        else
        {
            $lenmin=$array[$holenum][1];
        }
        $ratio[$total]=$lenmin/$lenmax;
        $total++;
    }
}
for($j=0;$j<1000;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<$total;$i++)
{   
    for($j=0;$j<1000;$j++)
    {
        if(($ratio[$i]>$j*0.001) && ($ratio[$i]<=($j+1)*0.001))
        {
            $win[$j]++;
            last;
        }
    }
}
print OUT "passtwo_ratio=";
for($j=0;$j<1000;$j++)
{
    $per += $win[$j]/$total;
    if($j<999)
    {
        printf OUT ("%0.5f:",$per);
    }
    else
    {
        printf OUT ("%0.5f\n",$per);
    }
}


@ratio=();
$per=0;
$total=0;
foreach $holenum(keys %pass)
{
    if($pass{$holenum}>2)
    {
        $lenmax=$array[$holenum][0];
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$lenmax)
            {
                $lenmax=$array[$holenum][$j];
            }
        }
        $ratio[$total]=$array[$holenum][0]/$lenmax;               #to statistic the ratio of the first and the last with the lengthmax in the same holenumber
        $total++;
        $ratio[$total]=$array[$holenum][$count[$holenum]]/$lenmax;
        $total++;
    }
}
for($j=0;$j<1000;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<$total;$i++)
{   
    for($j=0;$j<1000;$j++)
    {
        if(($ratio[$i]>$j*0.001) && ($ratio[$i]<=($j+1)*0.001))
        {
            $win[$j]++;
            last;
        }
    }
}
print OUT "passmore_ratio_head=";
for($j=0;$j<1000;$j++)
{
    $per += $win[$j]/$total;
    if($j<999)
    {
        printf OUT ("%0.5f:",$per);
    }
    else
    {
        printf OUT ("%0.5f\n",$per);
    }
}

@ratio=();
$per=0;
$total=0;
foreach $holenum(keys %pass)
{
    if($pass{$holenum}>2)
    {
        $lenmax=$array[$holenum][0];
        for($j=1;$j<=$count[$holenum];$j++)
        {
            if($array[$holenum][$j]>$lenmax)
            {
                $lenmax=$array[$holenum][$j];
            }
        }
       
        for($j=1;$j<$count[$holenum];$j++)
        {

            $ratio[$total]=$array[$holenum][$j]/$lenmax;        #statistics of the ratio of everyreads except the first and last with the lengthmax in the same holenumber 
           
            $total++;
        }
    }
}

for($j=0;$j<1000;$j++)
{
    $win[$j]=0;
}
for($i=0;$i<$total;$i++)
{   
    for($j=0;$j<1000;$j++)
    {
        if(($ratio[$i]>$j*0.001) && ($ratio[$i]<=($j+1)*0.001))
        {
            $win[$j]++;
            last;
        }
    }
}
print OUT "passmore_ratio_mid=";
for($j=0;$j<1000;$j++)
{
    $per += $win[$j]/$total;
    
    if($j<999)
    {
        printf OUT ("%0.5f:",$per);
    }
    else
    {
        printf OUT ("%0.5f\n",$per);
    }
}

my $holenumber=0;
my %coun=();
my $np=0;
my %countnumber=();
my $number=0;
my @percen;
my $prosum=0;
$total=0;

open(FASTQ,"<$ARGV[0]")or die "error reading $ARGV[0] for reading\n";
while(<FASTQ>)
{
    chomp;
    if($_ =~ /^\@(.+\/\d+)\//)
    {
        $holenumber=$1;
        if(exists $coun{$holenumber})
        {
            $coun{$holenumber}++;
        }
        else 
        {
            $coun{$holenumber}=1;
        }
    }
}

foreach $holenumber(keys %coun)
{
    $total++;
    $number=$coun{$holenumber};
    if(exists $countnumber{$number})
    {
       $countnumber{$number}++;
    }
    else
    {
        $countnumber{$number}=1;
    }
}
for($i=0;$i<500;$i++)
{
    $percen[$i]=0;
}
foreach $number(sort {$a <=> $b} keys %countnumber)
{
    $prosum += ($countnumber{$number}/$total);
    $percen[$number]=$prosum;
}
for($i=1;$i<500;$i++)                           #in case that some certain passnumber not exist 
{
    if($percen[$i]==0)
    {
        $percen[$i]=$percen[$i-1];
    }
}
print OUT "passnumber=";
for($i=0;$i<500;$i++)
{
    if($i<499)
    {
        printf OUT ("%0.5f:",$percen[$i]);
    }
    else
    {
        print OUT "1.00000\n";
    }
}

if($ARGV[1] eq "sequel")
{
	exit(0);
}


my $x=0;
my $y=0;
my $query="";
my $qv="";
my $a=0;
my %value=();
my $name1;
my $start1;
my $end1;

open(FASTQ,"<$ARGV[0]")or die "error reading $ARGV[0] for reading\n";
open(OUT1,">onehole.fq") or die "Can't open the output file: onehole.fq!\n";
while(<FASTQ>)
{
    chomp;
    $x++;
    if($_ =~ /^\@(.+\/\d+)\/(\d+)\_(\d+)/)
    {
            $a++;
            $y=$x;
            if(exists $value{$1})
            {
                $value{$1}++;
                $start1=$2;
                for($i=0;$i<($start1-$end1-1);$i++)
                {
                    $query=$query.'?';
                    $qv=$qv.'?';
                }
                $end1=$3;
            }
            else
            {
                $value{$1}=1;
                $end1=$3;
                if($a != 1)
                {
                    print OUT1 "$name1\n";
                    print OUT1 "$query\n";
                    print OUT1 "+\n";
                    print OUT1 "$qv\n";
                }
                $query="";
                $qv="";
                $name1=$_;
            }
        
    }
    if($y>0)
    {
        if($x==($y+1))
        {
            $query=$query.$_;
        }
        if($x==($y+3))
        {
            $qv=$qv.$_;
        }
    }
}
print OUT1 "$name1\n";
print OUT1 "$query\n";
print OUT1 "+\n";
print OUT1 "$qv\n";
close FASTQ;




my $lengthmax1=0;
my $len1;
my @qualityline;
my @quality;
my $linenumber=0;
my @number;

my @rand;

my $k=0;
my $mark1=0;
my $mark2=0;
my $flag=0;
my $tailsize1=2000;
my $m=0;
my @merge;
my $temp=0;
my $bound;
my $interval;
my @per2;
my @count2;
my $t;

my $qualitymin=10;
my $qualitymax=10;

my $tailsize2=1500;                              #bug!!need to be changede because if there is no certian length that agree with the tailsize,there will be some problems
my $s=0;
my @short;
my @long;
my $midmax=1;
my $shortmax=1;
my @tailshort;
my @tailmid;
my @mergeshort;
my @mergelong;
my @count3;
my @countshort;
my @countlong;
my $hh=0;
for($i=0;$i<$tailsize1;$i++)
{
    for($j=0;$j<20;$j++)
    {
        $tailshort[$i][$j]=0;
    }
}

open(IN,"<onehole.fq")or die "error reading $ARGV[0] for reading\n";
while(<IN>)
{
    chomp;
    $linenumber++;
    if ($linenumber%4 == 0)
    {
        @qualityline  = split("",$_);
        $len1 = @qualityline;
        for($i=0;$i<$len1;$i++)
        {
            if($qualityline[$i] ne '?')
            {
                $quality[$i]=ord($qualityline[$i])-33;
   
                if($quality[$i]<$qualitymin)
                {
                    $qualitymin=$quality[$i];
                }
                if($quality[$i]>$qualitymax)
                {
                    $qualitymax=$quality[$i];
                }
            }
        }
        if($len1>$lengthmax1)
        {
            for($i=$lengthmax1;$i<$len1;$i++)
            {
                for($j=0;$j<20;$j++)
                {
                    $long[$i][$j]=0;
                    $short[$i][$j]=0;
                }
            } 
            $lengthmax1=$len1;
        }
        if($len1>=4000)
        {
            $bound=$len1-$tailsize1;
            if($bound>$shortmax)
            {
                $shortmax=$bound;
            }
            for($i=0;$i<$bound;$i++)
            {
                if($qualityline[$i] ne '?')
                {
                    $short[$i][$quality[$i]]++;
                }
            }
            for($i=$bound;$i<$len1;$i++)
            {
                if($qualityline[$i] ne '?')
                {
                    $tailshort[$m][$quality[$i]]++;
                }
                $m++;
            }
            $m=0;
        }
        elsif($len1>2000)
        {
            $m=1000;
            $bound=$len1-1000;
            for($i=0;$i<$bound;$i++)
            {
                if($qualityline[$i] ne '?')
                {
                    $short[$i][$quality[$i]]++;
                }
            }
            for($i=$bound;$i<$len1;$i++)
            {
                if($qualityline[$i] ne '?')
                {
                    $tailshort[$m][$quality[$i]]++;
                }
                $m++;
            }
            $m=0;
        }
        else
        {
            for($i=0;$i<$len1;$i++)
            {
                if($qualityline[$i] ne '?')
                {
                    $short[$i][$quality[$i]]++;
                }
            }
        } 
    }    
}

close IN;

for($j=0;$j<=$qualitymax;$j++)
{
    $mergeshort[0][$j]=0;
    $merge[0][$j]=0;
    $mergelong[0][$j]=0;
}
for($i=0;$i<$shortmax;$i++)
{
    $countshort[$i][0]=$short[$i][0];
    for($j=1;$j<=$qualitymax;$j++)
    {
        $countshort[$i][$j]=$countshort[$i][$j-1]+$short[$i][$j];
    }
    if($countshort[$i][$qualitymax]>=20)
    {
        print OUT "$i-short_rand=";
        for($j=0;$j<20;$j++)
        {
            if($j<=$qualitymax)
            {
                $rand[$j]=$countshort[$i][$j]/$countshort[$i][$qualitymax];
                if($j<19)
                {
                    printf OUT ("%0.5f:",$rand[$j]);
                }
                else
                {
                    printf OUT ("%0.5f\n",$rand[$j]);
                }
            }
            else
            {
                if($j==19)
                {
                    print OUT "1.00000\n";
                }
                else
                {
                    print OUT "1.00000:";
                }
            }
        }
    }
    else
    {
        if($mark1<20)
        {
            for($j=0;$j<=$qualitymax;$j++)
            {
                $mergeshort[$k][$j]+=$countshort[$i][$j];
            }
            $mark1 += $countshort[$i][$qualitymax];
            $mark2++;
            if(($mark1>=20) ||($i==$shortmax-1 ))
            {
                for($t=($i-$mark2+1);$t<=$i;$t++)
                {
                    print OUT "$t-short_rand=";
                    for($j=0;$j<20;$j++)
                    {
                        if($j<=$qualitymax)
                        {
                            $rand[$j]=$mergeshort[$k][$j]/$mergeshort[$k][$qualitymax];
                            if($j<19)
                            {
                                printf OUT ("%0.5f:",$rand[$j]);
                            }
                            else
                            {
                                printf OUT ("%0.5f\n",$rand[$j]);
                            }
                        }
                        else
                        {
                            if($j==19)
                            {
                                print OUT "1.00000\n";
                            }
                            else
                            {
                                print OUT "1.00000:";
                            }
                        }
                    }
                }
                $mark1=0;
                $mark2=0;
                $k++;
                for($j=0;$j<=$qualitymax;$j++)
                {
                    $mergeshort[$k][$j]=0;
                }
            }

        }
    }
}
for($i=0;$i<$tailsize1;$i++)
{
    $count2[$i][0]=$tailshort[$i][0];
    for($j=1;$j<=$qualitymax;$j++)
    {
        $count2[$i][$j] = $count2[$i][$j-1]+$tailshort[$i][$j];
    }
    print OUT "$i-tailshort_rand=";
    for($j=0;$j<20;$j++)
    {
        if($j<=$qualitymax)
        {   
            $rand[$j]=$count2[$i][$j]/$count2[$i][$qualitymax];
            if($j<19)
            {
                printf OUT ("%0.5f:",$rand[$j]);
            }
            else
            {
                printf OUT ("%0.5f\n",$rand[$j]);
            }
        }
        else
        {
            if($j==19)
            {
                print OUT "1.00000\n";
            }
            else
            {
                print OUT "1.00000:";
            }
        }
    }
}


