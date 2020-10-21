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
my %length1=();
my %length2=();
my %length3=();
my @position1=();
my @position2=();
my @position3=();
my @pos1=();
my @pos2=();
my @pos3=();
my $j1;
my $j2;
my $j3;
my @keys1=();
my @keys2=();
my @keys3=();
my $count1=0;
my $count2=0;
my $count3=0;
my $leng;
my $pro;
my $quname;
my @errorlength;
my $prosum=0;
my @strand=();
my @stat=();
my $order;
my %per;
my @keys=();
my %num1=();
my %num2=();
my %num3=();
my $readsum1=0;
my $readsum2=0;
my $readsum3=0;
my $key1;
my $key2;
my $key3;
my $holequ="";
my $holeref="";
my %num4=();
my $key4;
my $holename;
my %hole=();
my $key;
my $sum=0;
my $len;
my $a;
my $b;
my $c;
my $d;
my $e;
my $allevent1=0;
my $allevent2=0;
my $allevent3=0;
my $temp1=0;
my $temp2=0;
my $temp3=0;
my $temp4=0;
my $allevent4=0;
my @holelen=();
my $p=0;
my @flag=();
my @holepos=();
my $templen;
my @tempqu;
my $fragment;
my $hhh=0;
my $tempj1;
my $tempj2;
my $bound1;
my $bound2;
my %countmark=();
my $mapratio;
my @tem=();
my @per=();
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
for($i=0;$i<64;$i++)
{
	($key)=convert($i);
 	$num1{$key}=0;
    $num2{$key}=0;
	$num3{$key}=0;
	$num4{$key}=0;
	$num5{$key}=0;
}
=pod
open(FASTQ,"<$ARGV[1]")or die "error reading $ARGV[0] for reading\n";
 
while(<FASTQ>)
{
	chomp;
    $fastqnum++;
	if($_ =~ /\@.+\/(\d+)\//)
	{
		$order=$1;
		if(exists $countmark{$order})
		{
			$countmark{$order}++;
		}
		else
		{
			$countmark{$order}=1;
			$fastq[$order]="";
		}
	}
	if($fastqnum % 4 ==2)
	{
		$fastq[$order]=$fastq[$order].$_;
	}
}
close FASTQ;
=cut
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
					@tempqu=split("",$query[$p]);
					$holepos[$p][0]=$flag[$p];
					for($i=1;$i<@tempqu;$i++)
					{	
						if($tempqu[$i] ne '-')
						{
							$holepos[$p][$i] = $holepos[$p][$i-1]+1;
						}
						else
						{
							$holepos[$p][$i]=$holepos[$p][$i-1];
						}
					}
					@tempqu=();
					$p++;
				}
				else
				{	
					$hole{$holename}=1;
					if($y!=1)
            		{
					
				        @tempqu=split("",$query[$p]);
				        $holepos[$p][0]=$flag[$p];
				        for($i=1;$i<@tempqu;$i++)
				        {
				            if($tempqu[$i] ne '-')
				            {
				                $holepos[$p][$i]=$holepos[$p][$i-1]+1;
				            }
				            else
				            {
				                $holepos[$p][$i]=$holepos[$p][$i-1];
				            }
				        }
				        @tempqu=();
                		$j1=0;
						$j2=0;
						$j3=0;
						$count1=0;
						$count2=0;
						$count3=0;
						@pos1=();
						@pos2=();
						@pos3=();
						@keys1=();
						@keys2=();
						@keys3=();
						@position1=();
						@position2=();
						@position3=();
						%length1=();
						%length2=();
						%length3=();
						
				         	for($i=0;$i<=$p;$i++)
				          	{
					               	$holequ=$holequ.$query[$i];
					            	$holeref=$holeref.$sbjct[$i];
							}

							$holequ=uc($holequ);
							$holeref=uc($holeref);
							@queryarray = split("",$holequ);
							@sbjctarray = split("",$holeref);

							for($i=0;$i<@queryarray;$i++)
							{
								if($queryarray[$i] ne '-')
								{
									$hhh++;
								}
							}
							#print"11111hhh=$hhh\n";
							$hhh=0;
							for($i=0;$i<@sbjctarray;$i++)
      		          		{	
								if($sbjctarray[$i] =~ /[ATGC\-]/)
								{
		                    	if($sbjctarray[$i] eq "-")      #ins
             					{
                    				$position1[$j1]=$i;
									$count1++;
									$j1++;
								}
								elsif($queryarray[$i] eq "-") #del
								{
									$position2[$j2]=$i;
									$count2++;
									$j2++;
								}
								elsif(($queryarray[$i] ne "-")&&($sbjctarray[$i] ne "-") && ($queryarray[$i] ne $sbjctarray[$i]))  #sub
								{
									$position3[$j3]=$i;
									$count3++;
									$j3++;
								}
								}
                			}
							if($count1 != 0)
							{
								$pos1[0]=$position1[0];
								$length1{$pos1[0]}=1;
								for($i=1;$i<$count1;$i++)
								{
									if(($position1[$i])==($position1[$i-1]+1))
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
								@keys1=keys(%length1);
								for($j=0;$j<@keys1;$j++)
								{
									if(($keys1[$j]+$length1{$keys1[$j]}+2)<@sbjctarray)
									{
										$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2];	
										if($key1 =~ /[^ATGC\-]+/)
										{
											next;
										}
										elsif($key1 =~ /^[ATGC]+$/)
										{
											$readsum1++;
											if(exists $num1{$key1})
											{
												$num1{$key1}++;
											}
											else
											{
												$num1{$key1}=1;
											}
										}
										else
										{
											$temp1++;
											if($sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1] eq '-')
											{
												if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1})<@sbjctarray)
												{
													$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}];
													if($key1 =~ /[^ATGC\-]+/)
													{
														next;
													}
													elsif($key1 =~  /^[ATGC]+$/)
													{
														$allevent1++;
														$readsum1++;
														if(exists $num1{$key1})
														{
															$num1{$key1}++;
														}	
														else
														{
															$num1{$key1}=1;
														}
													}
													else
													{
														if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}})<@sbjctarray)
														{
															$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}}];
															if($key1 =~ /[^ATGC\-]+/)
															{
																next;
															}
															elsif($key1 =~ /^[ATGC]+$/)
															{
																$allevent1++;
																$readsum1++;
																if(exists $num1{$key1})
																{
																	$num1{$key1}++;
																}
																else
																{
																	$num1{$key1}=1;
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
													$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+2}];
													if($key1 =~ /[^ATGC\-]+/)
													{
														next;
													}
													elsif($key1 =~ /^[ATGC]+$/)
													{
														$allevent1++;
														$readsum1++;
														if(exists $num1{$key1})
														{
															$num1{$key1}++;
														}
														else
														{
															$num1{$key1}=1;
														}
													}
												}
											}
										}			
									}
								}
							}
							for($i=0;$i<@sbjctarray;$i++)
							{
								if($queryarray[$i] eq $sbjctarray[$i])
						        {
						        	if(($i+2)<@sbjctarray)
						           	{
						            	$key4=$sbjctarray[$i].$sbjctarray[$i+1].$sbjctarray[$i+2];
						            	if($key4 =~ /[^ATGC\-]+/)
										{
											next;
										}
										elsif($key4 =~ /^[ATGC]+$/)
						              	{
											$allevent4++;
					    	          		if(exists $num4{$key4})
					        	      		{
					            	    		$num4{$key4}++;
					               			}
						               		else
						               		{
						                  		$num4{$key4}=1;
						             		}
						            	}
										else
										{
											$temp4++;
											if($sbjctarray[$i+1] eq '-')
											{
												if(($i+2+$length1{$i+1})<@sbjctarray)
												{
													$key4=$sbjctarray[$i].$sbjctarray[$i+1+$length1{$i+1}].$sbjctarray[$i+2+$length1{$i+1}];
													if($key4 =~ /[^ATGC\-]+/)
													{
														next;
													}
													elsif($key4 =~ /^[ATGC]+$/)
													{
														if(exists $num4{$key4})
														{
															$num4{$key4}++;
														}
														else
														{
															$num4{$key4}=1;
														}
													}
													else
													{
														if(($i+2+$length1{$i+1}+$length1{$i+2+$length1{$i+1}})<@sbjctarray)
														{
															$key4=$sbjctarray[$i].$sbjctarray[$i+1+$length1{$i+1}].$sbjctarray[$i+2+$length1{$i+1}+$length1{$i+2+$length1{$i+1}}];
															if($key4 =~ /[^ATGC\-]+/)
															{
																next;
															}
															elsif($key4 =~  /^[ATGC]+$/)
															{
																if(exists $num4{$key4})
																{
																	$num4{$key4}++;
																}
																else
																{
																	$num4{$key4}=1;
																}
															}
														}
													}
												}
											}
											else
											{
												if(($i+2+$length1{$i+2})<@sbjctarray)
												{
													$key4=$sbjctarray[$i].$sbjctarray[$i+1].$sbjctarray[$i+2+$length1{$i+2}];
													if($key4 =~ /[^ATGC\-]+/)
													{
														next;
													}
													elsif($key4 =~  /^[ATGC]+$/)
													{
														if(exists $num4{$key4})
														{
															$num4{$key4}++;
														}
														else
														{
															$num4{$key4}=1;
														}
													}
												}
											}
										}
									}
						         }
						   	 }
							if($count2 != 0)
							{
								$pos2[0]=$position2[0];
				    		$length2{$pos2[0]}=1;

							for($i=1;$i<$count2;$i++)
			        		{
			         			if(($position2[$i])==($position2[$i-1]+1))
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
							@keys2=keys(%length2);
	              			for($j=0;$j<@keys2;$j++)
	              			{
	                			if(($keys2[$j]+2)<@sbjctarray)
	                 			{
	                    			$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1].$sbjctarray[$keys2[$j]+2];
	                      			if($key2 =~ /[^ATGC\-]+/)
									{
										next;
									}
									elsif($key2 =~ /^[ATGC]+$/)
	                  				{
										$readsum2++;
	                   					if(exists $num2{$key2})
	                     				{
	                      					$num2{$key2}++;
	                     				}
	                  					else
	                   					{
	                    					$num2{$key2}=1;
	                   					}
	               					}
	               					else
									{
										$temp2++;
										if($sbjctarray[$keys2[$j]+1] eq '-')
										{
											if(($keys2[$j]+2+$length1{$keys2[$j]+1})<@sbjctarray)
											{
												$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1+$length1{$keys2[$j]+1}].$sbjctarray[$keys2[$j]+2+$length1{$keys2[$j]+1}];
												if($key2 =~ /[^ATGC\-]+/)
												{
													next;
												}
												elsif($key2 =~  /^[ATGC]+$/)
												{
													$allevent2++;
													$readsum2++;
													if(exists $num2{$key2})
													{
														$num2{$key2}++;
													}	
													else
													{
														$num2{$key2}=1;
													}
												}
												else
												{
													if(($keys2[$j]+2+$length1{$keys2[$j]+1}+$length1{$keys2[$j]+2+$length1{$keys2[$j]+1}})<@sbjctarray)
													{
														$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1+$length1{$keys2[$j]+1}].$sbjctarray[$keys2[$j]+2+$length1{$keys2[$j]+1}+$length1{$keys2[$j]+2+$length1{$keys2[$j]+1}}];
														if($key2 =~ /[^ATGC\-]+/)
														{
															next;
														}
														elsif($key2 =~ /^[ATGC]+$/)
														{
															$allevent2++;
															$readsum2++;
															if(exists $num2{$key2})
															{
																$num2{$key2}++;
															}
															else
															{
																$num2{$key2}=1;
															}
														}
													}
												}
											}
										}
										else
										{
											if(($keys2[$j]+2+$length1{$keys2[$j]+2})<@sbjctarray)
											{
												$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1].$sbjctarray[$keys2[$j]+2+$length1{$keys2[$j]+2}];
												if($key2 =~ /[^ATGC\-]+/)
												{
													next;
												}
												elsif($key2 =~ /^[ATGC]+$/)
												{
													$allevent2++;
													$readsum2++;
													if(exists $num2{$key2})
													{
														$num2{$key2}++;
													}
													else
													{
														$num2{$key2}=1;
													}
												}
											}
										}
									}
	            				}
	           				}
						}
						if($count3 != 0)
			    		{
			    			$pos3[0]=$position3[0];
							$length3{$pos3[0]}=1;
							for($i=1;$i<$count3;$i++)
							{
			        			if(($position3[$i])==($position3[$i-1]+1))
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
			      			@keys3=keys(%length3);
			  				for($j=0;$j<@keys3;$j++)
			      			{
			        			if(($keys3[$j]+2)<@sbjctarray)
			     				{
			            			$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1].$sbjctarray[$keys3[$j]+2];
			           				if($key3 =~ /[^ATGC\-]+/)
									{
										next;
									}
									elsif($key3 =~ /^[ATGC]+$/)
			          				{
										$readsum3++;
			               				if(exists $num3{$key3})
			                			{
			                   				$num3{$key3}++;
			                  			}
			               				else
			              				{
			                 				$num3{$key3}=1;
			               				}
			          				}
			          				else
									{
										$temp3++;
							 			if($sbjctarray[$keys3[$j]+1] eq '-')
										{
											if(($keys3[$j]+2+$length1{$keys3[$j]+1})<@sbjctarray)
											{
												$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1+$length1{$keys3[$j]+1}].$sbjctarray[$keys3[$j]+2+$length1{$keys3[$j]+1}];
												if($key3 =~ /[^ATGC\-]+/)
												{
													next;
												}
												elsif($key3 =~  /^[ATGC]+$/)
												{
													$allevent3++;
													$readsum3++;
													if(exists $num3{$key3})
													{
														$num3{$key3}++;
													}	
													else
													{
														$num3{$key3}=1;
													}
												}
												else
												{
													if(($keys3[$j]+2+$length1{$keys3[$j]+1}+$length1{$keys3[$j]+2+$length1{$keys3[$j]+1}})<@sbjctarray)
													{
														$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1+$length1{$keys3[$j]+1}].$sbjctarray[$keys3[$j]+2+$length1{$keys3[$j]+1}+$length1{$keys3[$j]+2+$length1{$keys3[$j]+1}}];
														if($key3 =~ /[^ATGC\-]+/)
														{
															next;
														}
														elsif($key3 =~ /^[ATGC]+$/)
														{
															$allevent3++;
															$readsum3++;
															if(exists $num3{$key3})
															{
																$num3{$key3}++;
															}
															else
															{
																$num3{$key3}=1;
															}
														}
													}
												}
											}
										}
										else
										{
											if(($keys3[$j]+2+$length1{$keys3[$j]+2})<@sbjctarray)
											{
												$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1].$sbjctarray[$keys3[$j]+2+$length1{$keys3[$j]+2}];
												if($key3 =~ /[^ATGC\-]+/)
												{
													next;
												}
												elsif($key3 =~ /^[ATGC]+$/)
												{
													$allevent3++;
													$readsum3++;
													if(exists $num3{$key3})
													{
														$num3{$key3}++;
													}
													else
													{
														$num3{$key3}=1;
													}
												}
											}
										}
									}
								}
							}
						}
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
							for($i=0;$i<100;$i++)
							{
								$query[$i]="";
								$sbjct[$i]="";
							}
						
					}
					if($holename =~ /.+\/(\d+)\//)
				    {
				        $holenum=$1;
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
			$flag[$p]=$1;
			$holelen[$p]=$3;
		}
		else
		{
			$holelen[$p]=$holelen[$p-1]+$3;
			$flag[$p] =$1+$holelen[$p-1];
		}
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

$j1=0;
$j2=0;
$j3=0;
$count1=0;
$count2=0;
$count3=0;
%length1=();
%length2=();
%length3=();
@position1=();
@position2=();
@position3=();
@pos1=();
@pos2=();
@pos3=();
#print"111111p=$p\n";
@tempqu=split("",$query[$p]);
$holepos[$p][0]=$flag[$p];
for($i=1;$i<@tempqu;$i++)
{
	if($tempqu[$i] ne '-')
    {
        $holepos[$p][$i]=$holepos[$p][$i-1]+1;
    }
    else
    {
        $holepos[$p][$i]=$holepos[$p][$i-1];
    }
}
	for($i=0;$i<=$p;$i++)
	{
		  $holequ=$holequ.$query[$i];
			$holeref=$holeref.$sbjct[$i];
 	}

	#print"holequ=$holequ\n";
	$holequ=uc($holequ);
	$holeref=uc($holeref);
	@queryarray = split("",$holequ);
	@sbjctarray = split("",$holeref);
	for($i=0;$i<@queryarray;$i++)
	{
		if($queryarray[$i] ne '-')
		{
			$hhh++;
		}
	}
	#print"22222hhhh=$hhh\n";
for($i=0;$i<@sbjctarray;$i++)
{
	if($sbjctarray[$i] =~ /[ATGC\-]/)
	{
	if($sbjctarray[$i] eq "-")
    {
    	$position1[$j1]=$i;
		$count1++;
		$j1++;
	}
	elsif($queryarray[$i] eq "-")
	{
		$position2[$j2]=$i;
		$count2++;
		$j2++;
	}
	elsif(($queryarray[$i] ne "-")&&($sbjctarray[$i] ne "-") && ($queryarray[$i] ne $sbjctarray[$i]))
	{
		$position3[$j3]=$i;
		$count3++;
		$j3++;
	}
	}
}
if($count1 != 0)
{
	$pos1[0]=$position1[0];
	$length1{$pos1[0]}=1;
	for($i=1;$i<$count1;$i++)
	{
		if(($position1[$i])==($position1[$i-1]+1))
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
	@keys1=keys(%length1);
	for($j=0;$j<@keys1;$j++)
	{
		if(($keys1[$j]+$length1{$keys1[$j]}+2)<@sbjctarray)
		{
			$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2];
			if($key1 =~ /[^ATGC\-]+/)
			{
				next;
			}
			elsif($key1 =~ /^[ATGC]+$/)
			{
				$readsum1++;
				if(exists $num1{$key1})
				{
					$num1{$key1}++;
				}
				else
				{
					$num1{$key1}=1;
				}
			}
			else
			{
				$temp1++;
				if($sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1] eq '-')
				{
					if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1})<@sbjctarray)
					{
						$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}];
						if($key1 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key1 =~  /^[ATGC]+$/)
						{
							$allevent1++;
							$readsum1++;
							if(exists $num1{$key1})
							{
								$num1{$key1}++;
							}	
							else
							{
								$num1{$key1}=1;
							}
						}
						else
						{
							if(($keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}})<@sbjctarray)
							{
								$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}+$length1{$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+1}}];
								if($key1 =~ /[^ATGC\-]+/)
								{
									next;
								}
								elsif($key1 =~ /^[ATGC]+$/)
								{
									$allevent1++;
									$readsum1++;
									if(exists $num1{$key1})
									{
										$num1{$key1}++;
									}
									else
									{
										$num1{$key1}=1;
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
						$key1=$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+1].$sbjctarray[$keys1[$j]+$length1{$keys1[$j]}+2+$length1{$keys1[$j]+$length1{$keys1[$j]}+2}];
						if($key1 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key1 =~ /^[ATGC]+$/)
						{
							$allevent1++;
							$readsum1++;
							if(exists $num1{$key1})
							{
								$num1{$key1}++;
							}
							else
							{
								$num1{$key1}=1;
							}
						}
					}
				}
			}
		}
	}
}
for($i=0;$i<@sbjctarray;$i++)
{
	if($queryarray[$i] eq $sbjctarray[$i])
   	{
    	if(($i+2)<@sbjctarray)
     	{
       		$key4=$sbjctarray[$i].$sbjctarray[$i+1].$sbjctarray[$i+2];
         	if($key4 =~ /[^ATGC\-]+/)
			{
				next;
			}
			elsif($key4 =~ /^[ATGC]+$/)
         	{
				$allevent4++;
           		if(exists $num4{$key4})
            	{
                	$num4{$key4}++;
              	}
           		else
          		{
              		$num4{$key4}=1;
           		}
         	}
         	else
       		{
				$temp4++;
           		if($sbjctarray[$i+1] eq '-')
          		{
                	if(($i+2+$length1{$i+1})<@sbjctarray)
               		{
                   		$key4=$sbjctarray[$i].$sbjctarray[$i+1+$length1{$i+1}].$sbjctarray[$i+2+$length1{$i+1}];
                    	if($key4 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key4 =~ /^[ATGC]+$/)
                   		{
                      		if(exists $num4{$key4})
                    		{
                       			$num4{$key4}++;
                     		}
                      		else
                     		{
                        		$num4{$key4}=1;
                			}
                 		}
                   		else
                   		{
                      		if(($i+2+$length1{$i+1}+$length1{$i+2+$length1{$i+1}})<@sbjctarray)
                       		{
								$key4=$sbjctarray[$i].$sbjctarray[$i+1+$length1{$i+1}].$sbjctarray[$i+2+$length1{$i+1}+$length1{$i+2+$length1{$i+1}}];
								if($key4 =~ /[^ATGC\-]+/)
								{
									next;
								}
								elsif($key4 =~ /^[ATGC]+$/)
								{
									if(exists $num4{$key4})
									{
										$num4{$key4}++;
									}
									else
									{
										$num4{$key4}=1;
									}
								}
							}
						}
					}
				}
				else
				{
					if(($i+2+$length1{$i+2})<@sbjctarray)
					{
						$key4=$sbjctarray[$i].$sbjctarray[$i+1].$sbjctarray[$i+2+$length1{$i+2}];
						if($key4 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key4 =~ /^[ATGC]+$/)
						{
							if(exists $num4{$key4})
							{
								$num4{$key4}++;
							}
							else
							{
								$num4{$key4}=1;
							}
						}
					}
				}
			}
		}
	}
}


if($count2 != 0)
{
	$pos2[0]=$position2[0];
    $length2{$pos2[0]}=1;

	for($i=1;$i<$count2;$i++)
   	{
    	if(($position2[$i])==($position2[$i-1]+1))
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
  	@keys2=keys(%length2);
	for($j=0;$j<@keys2;$j++)
    {
    	if(($keys2[$j]+2)<@sbjctarray)
        {
        	$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1].$sbjctarray[$keys2[$j]+2];
            if($key2 =~ /[^ATGC\-]+/)
			{
				next;
			}
			elsif($key2 =~ /^[ATGC]+$/)
          	{
				$readsum2++;
             	if(exists $num2{$key2})
                {
                	$num2{$key2}++;
               	}
               	else
           		{
             		$num2{$key2}=1;
           		}
			}
			else
			{
				$temp2++;
				if($sbjctarray[$keys2[$j]+1] eq '-')
				{
					if(($keys2[$j]+2+$length1{$keys2[$j]+1})<@sbjctarray)
					{
						$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1+$length1{$keys2[$j]+1}].$sbjctarray[$keys2[$j]+2+$length1{$keys2[$j]+1}];
						if($key2 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key2 =~  /^[ATGC]+$/)
						{
							$allevent2++;
							$readsum2++;
							if(exists $num2{$key2})
							{
								$num2{$key2}++;
							}	
							else
							{
								$num2{$key2}=1;
							}
						}
						else
						{
							if(($keys2[$j]+2+$length1{$keys2[$j]+1}+$length1{$keys2[$j]+2+$length1{$keys2[$j]+1}})<@sbjctarray)
							{
								$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1+$length1{$keys2[$j]+1}].$sbjctarray[$keys2[$j]+2+$length1{$keys2[$j]+1}+$length1{$keys2[$j]+2+$length1{$keys2[$j]+1}}];
								if($key2 =~ /[^ATGC\-]+/)
								{
									next;
								}
								elsif($key2 =~ /^[ATGC]+$/)
								{
									$allevent2++;
									$readsum2++;
									if(exists $num2{$key2})
									{
										$num2{$key2}++;
									}
									else
									{
										$num2{$key2}=1;
									}
								}
							}
						}
					}
				}
				else
				{
					if(($keys2[$j]+2+$length1{$keys2[$j]+2})<@sbjctarray)
					{
						$key2=$sbjctarray[$keys2[$j]].$sbjctarray[$keys2[$j]+1].$sbjctarray[$keys2[$j]+2+$length1{$keys2[$j]+2}];
						if($key2 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key2 =~ /^[ATGC]+$/)
						{
							$allevent2++;
							$readsum2++;
							if(exists $num2{$key2})
							{
								$num2{$key2}++;
							}
							else
							{
								$num2{$key2}=1;
							}
						}
					}
				}
			}
		}
	}
}
if($count3 != 0)
{
	$pos3[0]=$position3[0];
    $length3{$pos3[0]}=1;
    for($i=1;$i<$count3;$i++)
    {
    	if(($position3[$i])==($position3[$i-1]+1))
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
    @keys3=keys(%length3);
   	for($j=0;$j<@keys3;$j++)
   	{
      	if(($keys3[$j]+2)<@sbjctarray)
      	{
       		$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1].$sbjctarray[$keys3[$j]+2];
        	if($key3 =~ /[^ATGC\-]+/)
			{
				next;
			}
			elsif($key3 =~ /^[ATGC]+$/)
       		{
				$readsum3++;
           		if(exists $num3{$key3})
        		{
          			$num3{$key3}++;
       			}
       			else
       			{
            		$num3{$key3}=1;
      			}
			}
			else
			{
				$temp3++;
				if($sbjctarray[$keys3[$j]+1] eq '-')
				{
					if(($keys3[$j]+2+$length1{$keys3[$j]+1})<@sbjctarray)
					{
						$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1+$length1{$keys3[$j]+1}].$sbjctarray[$keys3[$j]+2+$length1{$keys3[$j]+1}];
						if($key3 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key3 =~  /^[ATGC]+$/)
						{
							$allevent3++;
							$readsum3++;
							if(exists $num3{$key3})
							{
								$num3{$key3}++;
							}	
							else
							{
								$num3{$key3}=1;
							}
						}
						else
						{
							if(($keys3[$j]+2+$length1{$keys3[$j]+1}+$length1{$keys3[$j]+2+$length1{$keys3[$j]+1}})<@sbjctarray)
							{
								$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1+$length1{$keys3[$j]+1}].$sbjctarray[$keys3[$j]+2+$length1{$keys3[$j]+1}+$length1{$keys3[$j]+2+$length1{$keys3[$j]+1}}];
								if($key3 =~ /[^ATGC\-]+/)
								{
									next;
								}
								elsif($key3 =~ /^[ATGC]+$/)
								{
									$allevent3++;
									$readsum3++;
									if(exists $num3{$key3})
									{
										$num3{$key3}++;
									}
									else
									{
										$num3{$key3}=1;
									}
								}
							}
						}
					}
				}
				else
				{
					if(($keys3[$j]+2+$length1{$keys3[$j]+2})<@sbjctarray)
					{
						$key3=$sbjctarray[$keys3[$j]].$sbjctarray[$keys3[$j]+1].$sbjctarray[$keys3[$j]+2+$length1{$keys3[$j]+2}];
						if($key3 =~ /[^ATGC\-]+/)
						{
							next;
						}
						elsif($key3 =~ /^[ATGC]+$/)
						{
							$allevent3++;
							$readsum3++;
							if(exists $num3{$key3})
							{
								$num3{$key3}++;
							}
							else
							{
								$num3{$key3}=1;
							}
						}
					}
				}
			}					
		}
	}
}

#$mapratio=(1-$ARGV[4])*0.8935/0.8734;
for($i=0;$i<64;$i++)
{
	($key)=convert($i);
	$sum=$num1{$key}+$num2{$key}+$num3{$key}+$num4{$key};
	if($sum>0)
	{
		$stat[$i][0]=$num4{$key}/$sum*0.981;   #match
		$stat[$i][1]=($num4{$key}+$num1{$key})/$sum*0.991;   #ins
		$stat[$i][2]=($num4{$key}+$num1{$key}+$num2{$key})/$sum*1.003;      #del
		$stat[$i][3]=1;    #sub
	}
	else
	{
		$stat[$i][0]=0.85;
		$stat[$i][1]=0.09;
		$stat[$i][2]=0.045;
		$stat[$i][3]=0.015;
	}
	$sum=0;
}
for($i=0;$i<@stat;$i++)
{
	for($j=1;$j<4;$j++)
  	{
    	$tem[$i][$j]=$stat[$i][$j]-$stat[$i][$j-1];
		$sum += $tem[$i][$j];
	}
	$per[$i][1]=$tem[$i][1]/$sum;
    $per[$i][2]=($tem[$i][1]+$tem[$i][2])/$sum;
    $per[$i][3]=1;
    $sum=0;
}

print"fourper_head=";
for($i=0;$i<64;$i++)
{
	printf("%0.5f:",$stat[$i][0]);
}
for($i=0;$i<64;$i++)
{
	for($j=1;$j<4;$j++)
	{
		if($i==63 && $j==3)
		{
			printf("%0.5f\n",$per[$i][$j]);
		}
		else
		{
			printf("%0.5f:",$per[$i][$j]);
		}
	}
}


