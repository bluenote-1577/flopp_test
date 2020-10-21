#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#define MAX_LINE 100000
#define NN 200000

//======================================================================================
//function:usage 
//======================================================================================
//int bp();
void read_print();
void input_right();
void minus_strand();
int species_number();
int get_qv_type();
int first_line();
int remain_lines();
int quality();
int length_distribution();
int start_position();
char convert();
int error_length();
int pass_number();
float length_ratio();
void *multi_pass_thread();
//int ffun_per();
//int decide_event();
//float get_eachread();
int bp(int,float *,float *,float *,float *,int *,int,float [],int,int,int,float [],char *,char *,char *,float [],char *,int,int,float *,int,float [],int,int,int,float *,float *,float *,float *,float *,float [],float [],float *);     
int ffun_per(float *,char [],int *,float [],float [],float [],int,int);
int decide_event(float *,int,int,int *,float []);
float get_eachread(int *,float []);
float get_unalign(int *,float []);



void usage()
{
	printf("This is a sequencing simulator for PacBio sequencing: PaSS.\n");
	//printf("By combining completed genomes currently available and some initial metagenome sequencing information, it can create a simulated metagenome closer-than-ever to the reality.\n");
	//printf("With this system, an optimized metagenomics analysis strategy for a specific project can be created according to the sequencing platforms available to the user, the complexity of targeted metagenome, and the goal of metagenome sequencing. Currently, NeSSM supports 454,sange and Illumina sequencing platforms.\n"); 
	printf("PaSS can be helpful to evaluate or develop tools for PacBio sequencing.\n\n");
	printf("Usage: ./PaSS  [options]\n");
	printf("[options]:\n");
	printf("-list <input_file>             percentage.txt\n");
	printf("-index <index_file>            index\n");
	printf("-m <sequencing_method>         'pacbio_RS' or 'pacbio_sequel'.\n");
 	printf("-c <error_model_file>          error model file. e.g. 'sim.config'.\n");
	printf("-r <reads_number>              number of reads to generate.\n");
	printf("-t <threads_number>            number of threads to use.default is 1.\n");
	printf("-o <output_file>               output file.\n");
	printf("-d                             If '-d' is set, the ground truth of simulation will output concurrently.\n");
	//printf("-l <sequencing length>          length of the reads, default for illumina sequencing 36(36bp).\n");
	//printf("-exact <exact simulation>       0 means the length of read is decided by the parameter -l, 1 means the length of read is decided by the distribution of length according to a real data, default is 0.\n");
}

void runNeSSM(int argc,char **argv);   //declare a function:runNeSSM
//======================================================================================
//function:rand_length--generate a read's length according to average and standard deviation,  
//		   the result obeys normal distribution. the length doesn't exceed the max.
//=====================================================================================
int rand_length(float average,float sd,int max,int min,int mark)
{
	float out,random,test,tmp;
	int flag=0;
    
	if(sd==0)
	{
		return (int)average;
	}

	while(flag==0)
	{
		random=(float)((float)rand()/RAND_MAX-0.5)*50*sd+average;
		test=(float)rand()/RAND_MAX*2.0/sqrt(2.0*3.1415926)/sd;
		out=random;
		random=1.0/sqrt(2.0*3.1415926)/sd*(float)exp((-1)*(random-average)*(random-average)/2.0/sd/sd);

		if((mark==1)&&(test<random)&&(out>min)&&(out<=max))    //round-off
		{
			flag=1;
		}
        else if(((mark==5)||(mark==4))&&(test<random)&&(exp(out)>min)&&(exp(out)<=max))
        {
           	flag=1;
            tmp=out;
            out=exp(tmp);
        }	
    }
	return (int)(out);
}
//=====================================================================================
//
//======================================================================================
int minus_bp(int truth,float *rate,float *insbias,float *delbias,float *subbias,int *seed,float sublength[],int holenum,int start,int end,int sum_len,int last_len,int x,int flag_pn,int every_length,int point,char *buff_sequence,char *sequence_name,int buff_point,int len_sequence_name,int lable_method,float type[],char *line2,char *line4,char *minus_buffle,float subratio[],float *err_rand,float qual[],int len,int length_max,int input_max,float *short_rand,float *mid_rand,float *long_rand,float *tailshort,float *tailmid,float inslength[],float dellength[],float *ins_less,int buff,int flag_buff,int species_reads)
{
    memset(line2,'\0',(input_max*2)*sizeof(char));         //attention memeset line2 and line4 before the simulation of every read 
    memset(line4,'\0',(input_max*2)*sizeof(char));
    
    buff_point=first_line(holenum,start,end,point,buff_sequence,sequence_name,buff_point,len_sequence_name,flag_pn);
	buff_point=bp(truth,rate,insbias,delbias,subbias,seed,holenum,sublength,sum_len,last_len,lable_method,type,line2,line4,minus_buffle,subratio,buff_sequence,buff_point,every_length,err_rand,(len-point-1),qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less);
    buff_point=remain_lines(buff_sequence,buff_point,line2,line4,every_length);
    flag_buff++;
    return buff_point;
}
//=======================================================================================
//
//=======================================================================================
int plus_bp(int truth,float *rate,float *insbias,float *delbias,float *subbias,int *seed,float sublength[],int holenum,int start,int end,int sum_len,int last_len,int x,int flag_pn,int every_length,int point,char *buff_sequence,char *sequence_name,int buff_point,int len_sequence_name,int lable_method,float type[],char *line2,char *line4,char *buffle,float subratio[],float *err_rand,float qual[],int len,int length_max,int input_max,float *short_rand,float *mid_rand,float *long_rand,float *tailshort,float *tailmid,float inslength[],float dellength[],float *ins_less,int buff,int flag_buff,int species_reads)
{
    memset(line2,'\0',(input_max*2)*sizeof(char));  //attention:memset line2 and line4 before the simulation of every read
    memset(line4,'\0',(input_max*2)*sizeof(char));

    buff_point=first_line(holenum,start,end,point,buff_sequence,sequence_name,buff_point,len_sequence_name,flag_pn);
    buff_point=bp(truth,rate,insbias,delbias,subbias,seed,holenum,sublength,sum_len,last_len,lable_method,type,line2,line4,buffle,subratio,buff_sequence,buff_point,every_length,err_rand,point,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less);
    buff_point=remain_lines(buff_sequence,buff_point,line2,line4,every_length);
    flag_buff++;
    return buff_point;
}
//======================================================================================
//function:split--split a string according with a small string and return how many segments
//          after spliting.
//======================================================================================
int split(char **arr,char *str,char *del)
{
	int i=0;
	char *s=NULL;
	s=strtok(str,del);
	while(s!=NULL)
	{
		*arr++=s;
		s=strtok(NULL,del);
		i++;
	}
	return i;
}
//======================================================================================
//function:sub_single--substitute one base
//======================================================================================
char sub_single(int *seed,char b1,char b2,char b3,float p1,float p2,float p3)
{
	float per1,per2,per;
	per=(float)(rand_r(seed))/(float)(RAND_MAX);
	per1=p1/(p1+p2+p3);
	per2=p2/(p1+p2+p3);
    if(per<=per1)
		return b1;
	else if(per>per1 && per<=(per1+per2))
		return b2;
	else 
		return b3;	
} 
//======================================================================================
//function:sub_all--substitute the inputed base
//======================================================================================
char sub_all(int *seed,char seq,float ac,float ag,float at,float cg,float ct,float ca,float gt,float ga,float gc,float ta,float tc,float tg)
{ 
    if(seq=='A')
		return sub_single(seed,'C','G','T',ac,ag,at);
	else if(seq=='C')
		return sub_single(seed,'G','T','A',cg,ct,ca);
	else if(seq=='G')
		return sub_single(seed,'T','A','C',gt,ga,gc);
	else if(seq=='T')
		return sub_single(seed,'A','C','G',ta,tc,tg);
}
//======================================================================================
//function:ins_all--generate a base for insertion
//======================================================================================
char ins_all(int *seed)
{
	float random=(float)(rand_r(seed))/(float)(RAND_MAX);
    if(random<=0.25)
    	return 'A';
	else if(random<=0.5)
       	return 'T';
    else if(random<=0.75)
        return 'G';
   	else
      	return 'C';
}

main(int argc,char **argv)
{
	if(argc>=10)
	{
    	int t,lable;
		lable=0;
		for(t=1;t<=argc-1;t+=2)
		{
			if(strcmp(argv[t],"-index")==0){lable+=1;}
			else if(strcmp(argv[t],"-list")==0){lable+=1;}
			else if(strcmp(argv[t],"-m")==0){lable+=1;}
			else if(strcmp(argv[t],"-c")==0){lable+=1;}
			else if(strcmp(argv[t],"-o")==0){lable+=1;}

		}
		if(lable==5)
		{
			time_t start,end;
			start=time(NULL);
			printf("This is a PacBio sequencing simulator PaSS\n");
			printf("Start...\n");
			runNeSSM(argc,argv);
			printf("The end!\n");
			end=time(NULL);
			printf("print %f s\n",difftime(end,start));
			return 0;
		}
	    else
		{
			printf("You must set the parameter '-list','-index','-m','-c','-o'.\n\n");
			usage();
			return -1;
		}
	}
	else
	{
		usage();
		return -1;
	}
}

struct fun_para
{
    int ranseed;
    int id;
	float *lenone;
	float *lentwo;
	float *lenthree;
	float *lenfour;
	float *lenfive;
	float *lensix;
	float *lenseven;
	float *leneight;
	float *lennine;
    float *lenten;
    float *leneleven;
    float *lentwelve;
    float *lenthirteen;		
    //FILE *filefs;
    int current_reads;
    int *passnum;
    float *sublen;
  	int hole_start;
    int hole_num;
 	char *buffseq;
	char *seqname;
  	//int buffpoint;
  	int lenseqname;
  	int lablemethod;
  	float *errtype;
  	char *linetwo;
 	char *linefour;
	char *minus;
  	char *plus;
  	float *subrate;
  	float *errrand;
  	float *quality;
 	int leng;
  	int lengthmax;
  	int inputmax;
  	float *shortrand;
	float *midrand;
  	float *longrand;
  	float *tail_short;
  	float *tail_mid;
  	float *ins_length;
  	float *del_length;
 	float *insless;
  	float *tworatio;
	float *headratio;
	float *midratio;
	float *insarray;
	float *delarray;
	float *subarray;
	int buffnum;
  	int flagbuff;
  	int speciesreads;			
	float *eacharray;
	float *unheadarray;
	float *untailarray;
	int iftruth;
	float *averarray;
};
//======================================================================================
//function:runNeSSM--the main function to simulate
//======================================================================================

void runNeSSM(int argc,char **argv)
{
	//some variants
	int reads,gap,length,circle,lable,local,flag_buff,buff,gap_lable,lable_method,len,len_sequence_name,everynumber[2000],min,exact,flag_bias;
	int species_reads,len_buf,i,x,y,flag,total,point1,point2,point,flag_pn,buff_point,buff_point1,buff_point2,length_max,every_length;
	buff=1; 
	gap_lable=0;      
	total=0; 
	reads=1000;
	gap=200;
	length=50;
	lable=0;  
	exact=1;
	char index[MAX_LINE],list[MAX_LINE],name[MAX_LINE],*index_name[MAX_LINE],*path,*ref_name,*number[2],*tmp_buffle;
	char fastq[MAX_LINE];
    char *buffle,*sequence_name,name1[2000][200],*line2_1,*line2_2,*line4_1,*line4_2,*line2,*line4,*minus_buffle,bias[MAX_LINE];
	char config[MAX_LINE]="sim.config";
	char method[20]="pacbio_RS";
	FILE *fp,*findex,*fpath,*fs1,*fs2,*fc,*fbias;
	float random_point,reads_number,every_percent,sd_ratio[1],type[3],sd,subratio[12],*err_rand,length_rand[50000],inslength[5000],dellength[5000],sublength[5000],ins_less[500],*lengthone_rand,*lengthtwo_rand,*lengththree_rand,*lengthfour_rand,*lengthfive_rand,*lengthsix_rand,*lengthseven_rand,lengtheight_rand[50000],lengthnine_rand[50000];
	float lengthten_rand[50000],lengtheleven_rand[50000],lengthtwelve_rand[50000],lengththirteen_rand[50000];
    float *short_rand,*mid_rand,*long_rand,*tailshort,*tailmid;
    int gi,max_t,get_gi,flag_gi;
	char *all,*part_t[2],*value_t[NN];
	float value_f[NN],qual[50];
    int count=0;
    int test=0;    
	flag_bias=0;   //don't use the coverage bias file
    int input_max=-1;                              //10.13 change
    float passnumber[100];
	int cir_minus,cir_plus;   //the number of circulation on minus strand or minus strand
    float two_ratio[1000],head_ratio[1000],mid_ratio[1000];  //the ratio of the read with the lengthmax in the same holenumber,two_ratio is for pass=2,head_ratio is for first and last for pass>2,the tail_ratio is for mid that pass>2
    int tempmax;   //record the lengthmax when pass>2 
    float per; 
    int temp_count=0;
	int minus_point,plus_point;
	int last_len=0;     //1.19 change 
	int sum_len=0;
	int read_count=0;
	int len_arr[200];
	int j=0;
	int order_flag;
	time_t now;
	struct tm *tm_now;
	int holenum=0;
	int end;
	char temp_tm[80];
	void *rst;
    int ret;
    int temp_hole=0;
    int THREADS_NUM=1;
    int holesta=0;
    int readcir[1000],cirnum;             //not readcir[THREADS_NUM]!!!!!!
    float insbias[256],delbias[2048],subbias[256],eachread[50],unhead[1000],untail[1000];
	int truth=-1;
	float avererr[1];
	avererr[0]=0.14;


    for(i=0;i<50;i++)
		qual[i]=pow(10,(i/(-10.0)));
    
	for(circle=1;circle<=argc-1;circle+=2)					//read the input commend
	{     
		if(strcmp(argv[circle],"-r")==0){reads=atoi(argv[circle+1]);}
		//else if(strcmp(argv[circle],"-l")==0){length=atoi(argv[circle+1]);}
		else if(strcmp(argv[circle],"-index")==0){strcpy(index,argv[circle+1]);}
		//else if(strcmp(argv[circle],"-e")==0){lable=atoi(argv[circle+1]);}
		else if(strcmp(argv[circle],"-list")==0){strcpy(list,argv[circle+1]);}
		else if(strcmp(argv[circle],"-m")==0){strcpy(method,argv[circle+1]);}
		else if(strcmp(argv[circle],"-o")==0){local=circle;}  
		//else if(strcmp(argv[circle],"-w")==0){gap=atoi(argv[circle+1]);gap_lable=1;}
		else if(strcmp(argv[circle],"-c")==0){strcpy(config,argv[circle+1]);}
		//else if(strcmp(argv[circle],"-buff")==0){buff=atoi(argv[circle+1]);}
		//else if(strcmp(argv[circle],"-exact")==0){exact=atoi(argv[circle+1]);}
		//else if(strcmp(argv[circle],"-b")==0){strcpy(bias,argv[circle+1]);flag_bias=1;}
        else if(strcmp(argv[circle],"-input_max")==0){input_max=atoi(argv[circle+1]);}           //10.13 change:add a new parameter:input_max
		else if(strcmp(argv[circle],"-d")==0){truth=1;}
		else if(strcmp(argv[circle],"-t")==0){THREADS_NUM=atoi(argv[circle+1]);}
		else {printf("Waring:there is no this parameter '%s',please see the Usage again!\n",argv[circle]);exit(1);}  
	    
    }
	char *line02[THREADS_NUM],*line04[THREADS_NUM],*buff_sequence0[THREADS_NUM];
	int create_num=0;
    int *pass[THREADS_NUM];
    struct fun_para para[THREADS_NUM];
    FILE *tempfs[THREADS_NUM],*fs;
    int readnum[THREADS_NUM],tempseed[THREADS_NUM];
	pthread_t threads[THREADS_NUM];
	for(i=0;i<THREADS_NUM;i++)
	{
		readnum[i]=0;
        tempseed[i]=0;
        pass[i]=(int *)malloc(50000*sizeof(int));
        memset(pass[i],'\0',50000*sizeof(int));
    }

	if(lable==0)								//open the file to reserve the simulation
	{	
		memset(fastq,'\0',MAX_LINE*sizeof(char));
        strcpy(fastq,argv[local+1]);
        strcat(fastq,".fq"); 
        fs=fopen(fastq,"w");   
        if(fs==NULL){printf("Warning:can not make the fastq file in this directory.\n");exit(1);}
		printf("the fastq file:%s.\n",fastq);
	}

    fp=fopen(list,"r");							//open the species' list
	if(fp==NULL){printf("Warning:can not find the percentage.txt in this directory.\n");exit(1);}

    findex=fopen(index,"r");					//open the index file
	if(findex==NULL){printf("Warning:can not find the index file in this directory.\n");exit(1);}

    fc=fopen(config,"r");						//open the config file
	if(fc==NULL){printf("Warning:can not find the configure file in this directory.\n");exit(1);}

	sd_ratio[0]=2;	
    
  	err_rand=(float *)malloc(5000*50*sizeof(float));
   	memset(err_rand,'\0',5000*50*sizeof(float));
    short_rand=(float *)malloc(200000*20*sizeof(float));              //11.3 change
  	memset(short_rand,'\0',200000*20*sizeof(float));
    mid_rand=(float *)malloc(30000*20*sizeof(float));
    memset(mid_rand,'\0',30000*20*sizeof(float));
    long_rand=(float *)malloc(50000*20*sizeof(float));
    memset(long_rand,'\0',50000*20*sizeof(float));
    tailshort=(float *)malloc(3000*20*sizeof(float));
    memset(tailshort,'\0',3000*20*sizeof(float));
    tailmid=(float *)malloc(1500*20*sizeof(float));
    memset(tailmid,'\0',1500*20*sizeof(float));
    lengthone_rand=(float *)malloc(200000*sizeof(float));
	memset(lengthone_rand,'\0',200000*sizeof(float));
	lengthtwo_rand=(float *)malloc(200000*sizeof(float));
	memset(lengthtwo_rand,'\0',200000*sizeof(float));
	lengththree_rand=(float *)malloc(200000*sizeof(float));
	memset(lengththree_rand,'\0',200000*sizeof(float));
	lengthfour_rand=(float *)malloc(200000*sizeof(float));
	memset(lengthfour_rand,'\0',200000*sizeof(float));
	lengthfive_rand=(float *)malloc(200000*sizeof(float));
	memset(lengthfive_rand,'\0',200000*sizeof(float));
	lengthsix_rand=(float *)malloc(200000*sizeof(float));
	memset(lengthsix_rand,'\0',200000*sizeof(float));
	lengthseven_rand=(float *)malloc(200000*sizeof(float));
	memset(lengthseven_rand,'\0',200000*sizeof(float));
	min=1;
        
    if(strstr(method,"pacbio_RS")!=NULL)
	{
		lable_method=4;
    }
    else if(strstr(method,"pacbio_sequel")!=NULL)
   	{   
        lable_method=5;
    }  
	else
	{
		printf("This system only can simulate PacBio datas!\n");
		exit(1);
	}
    
	length_rand[0]=2; //as a flag
		length_max=get_qv_type(avererr,unhead,untail,eachread,insbias,delbias,subbias,sublength,lengththirteen_rand,lengthtwelve_rand,lengthten_rand,lengtheleven_rand,lengthfive_rand,lengthsix_rand,lengthseven_rand,lengtheight_rand,lengthnine_rand,two_ratio,head_ratio,mid_ratio,fc,type,subratio,sd_ratio,err_rand,lable_method,length_rand,exact,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,lengthone_rand,lengthtwo_rand,lengththree_rand,lengthfour_rand,passnumber);//function and get the simulation type
        if(input_max==-1)
        {
            input_max=length_max;
        }
        else
        {
            input_max=input_max;
        }
                                       
    //input_right(length,reads,index,lable,gap_lable,gap,length_max,method,exact,length_rand);
    
	if(exact==0)
		if(sd_ratio[0]!=2)    //get the sd_ratio from config file
			sd=sd_ratio[0]*length;
		else if(lable_method==1)   //454,initiation
			sd=0.02*length;
		else if(lable_method==4)
            sd=1.0454;
        else if(lable_method==5)
            sd=0.5966;
        else
			sd=0;

	srand((unsigned)time(0));				//seed the rand			
	//srand(0);
    total=species_number(fp,everynumber,name1,reads);				//function and get the whole number of species
   	//printf("total=%d\n",total);
    
	line2=(char *)malloc((input_max+1)*sizeof(char));
	line4=(char *)malloc((input_max+1)*sizeof(char));
	for(i=0;i<THREADS_NUM;i++)
  	{
    	line02[i]=(char *)malloc((input_max*2)*sizeof(char));
       	line04[i]=(char *)malloc((input_max*2)*sizeof(char));

 	}

	path=(char *)malloc(2000*sizeof(char));
    
	for(circle=0;circle<total;circle++)				 //in a specie,assign reads number
	{  
        species_reads=everynumber[circle];
		//printf("species_reads=%d\n",species_reads);
		if(species_reads == 0)
		{
			continue;                     //7.11 add
		}
		flag=0;                         //flag:index has this species or not
		
			 while(fgets(name,MAX_LINE,findex)!=NULL)                
			 {   
			 	name[strlen(name)-1]='\0'; 
			   	split(index_name,name,"\t");
				//printf("index_name[0]=%s\tindex_name[1]=%s\t index_name[2]=%s\n",index_name[0],index_name[1],index_name[2]);
				if(memcmp(name1[circle],index_name[1],strlen(name1[circle]))==0)  //add one condition
				{  
					flag=1;
				   	len=atoi(index_name[0]);
					memset(path,'\0',2000*sizeof(char));
					strcat(path,index_name[2]);
					//break;
				}
			}
		

		rewind(findex);

		if(flag==0)
		{
			printf("Warning:%s can't find in the index!!!",name1[circle]);
			continue;
		}

		//simulation
		buffle=(char *)malloc(len*2*sizeof(char));
		memset(buffle,'\0',len*2*sizeof(char));
		fpath=fopen(path,"r");
		if(fpath==NULL)
		{
			printf("Warning:can not find the index file pathway in this directory.\n");
			exit(1);
		}

		tmp_buffle=(char *)malloc(len*2*sizeof(char));
		ref_name=(char *)malloc(10000*sizeof(char));       //9.26 400 change to 10000 
		memset(tmp_buffle,'\0',len*2*sizeof(char));
		memset(ref_name,'\0',10000*sizeof(char));          //9.26
		fread(tmp_buffle,len*2,1,fpath);
		fclose(fpath);	  
		flag=0;
		char *p1,*p2,*p3;
		p1=tmp_buffle;p2=ref_name;p3=buffle;

		while(*p1 != '\0')			//get the specie's name and its sequence
		{
			if(*p1 == '\n') {p1++;flag+=1;continue;}
			else
			{
				if(flag==0){*p2++=*p1++;}
				else{*p3++=*p1++;}
			}
		}
		free(tmp_buffle);
		len_buf=strlen(ref_name);
		sequence_name=(char *)malloc(100*sizeof(char));         //len_buf+2 change to 100
		memset(sequence_name,'\0',100*sizeof(char));
		sequence_name[0]='@';
		time(&now);                  //1.25 add normal name
		tm_now=localtime(&now);
		//strftime(temp_tm,80,"m%y%m%d_%H%M%S_00000_cbarcode_s1_p0/",tm_now);
		strftime(temp_tm,80,"m%y%m%d_%H%M%S/",tm_now);
		strcat(sequence_name,temp_tm);

		//strcat(sequence_name,ref_name);
		len_sequence_name=strlen(sequence_name);
		minus_buffle=(char *)malloc(len*2*sizeof(char));
        memset(minus_buffle,'\0',len*2*sizeof(char));
		minus_strand(buffle,minus_buffle,len);

		flag_buff=0;           //buff's flag
		buff_point=0;		   //record the positon in the sequence
		buff_point1=0;
		buff_point2=0;	
///////
		if(flag_bias==1)
		{
			flag_gi=0;
			fbias=fopen(bias,"r");
			if(fbias==NULL){printf("Warning:can not find the coverage bias file.\n");exit(1);}
			all=(char *)malloc(NN*11*sizeof(char));
			memset(all,'\0',NN*11*sizeof(char));
			while(fgets(all,NN*11*sizeof(char),fbias)!=NULL)
			{
				all[strlen(all)-1]='\0';
				split(part_t,all,"=");
		
				get_gi=atoi(part_t[0]);
				if(gi==get_gi)
				{
					max_t=split(value_t,part_t[1],":");
					flag_gi=1;
					break;
				}
			}
			fclose(fbias);
			if(flag_gi==0)
			{
				printf("Warning:can not find the coverage bias information about %s\n",name1[circle]);
				exit(1);
			}
			for(x=0;x<max_t;x++)
				value_f[x]=atof(value_t[x]);
			free(all);
		}
///////////	
       if(lable==0)											//simulation   single
		{
            for(i=0;i<THREADS_NUM;i++)
            {
                para[i].id=i;
                para[i].lenone=lengthone_rand;
                para[i].lentwo=lengthtwo_rand;
                para[i].lenthree=lengththree_rand;
                para[i].lenfour=lengthfour_rand;
                para[i].lenfive=lengthfive_rand;
                para[i].lensix=lengthsix_rand;
                para[i].lenseven=lengthseven_rand;
                para[i].leneight=lengtheight_rand;
                para[i].lennine=lengthnine_rand;
                para[i].lenten=lengthten_rand;
                para[i].leneleven=lengtheleven_rand;
                para[i].lentwelve=lengthtwelve_rand;
                para[i].lenthirteen=lengththirteen_rand;
                //para[i].filefs=tempfs[i];
                para[i].seqname=sequence_name;
                para[i].sublen=sublength;
                //para[i].buffpoint=buff_point;
                para[i].lenseqname=len_sequence_name;
                para[i].lablemethod=lable_method;
                para[i].errtype=type;
                para[i].plus=buffle;
                para[i].minus=minus_buffle;
                para[i].subrate=subratio;
                para[i].errrand=err_rand;
                para[i].quality=qual;
                para[i].leng=len;
                para[i].lengthmax=length_max;
                para[i].inputmax=input_max;
                para[i].shortrand=short_rand;
                para[i].midrand=mid_rand;
                para[i].longrand=long_rand;
                para[i].tail_short=tailshort;
                para[i].tail_mid=tailmid;
                para[i].ins_length=inslength;
                para[i].del_length=dellength;
                para[i].insless=ins_less;
                para[i].buffnum=buff;
                para[i].flagbuff=flag_buff;
                para[i].speciesreads=species_reads;
                para[i].tworatio=two_ratio;
                para[i].headratio=head_ratio;
                para[i].midratio=mid_ratio;
            	para[i].insarray=insbias;
				para[i].delarray=delbias;
				para[i].subarray=subbias;
				para[i].eacharray=eachread;
				para[i].unheadarray=unhead;
				para[i].untailarray=untail;
				para[i].iftruth=truth;
				para[i].averarray=avererr;
			}
            if(species_reads % 5000==0)
            {
                cirnum=(species_reads/5000);
            }
            else
            {
                cirnum=species_reads/5000+1;
            }
            
            //printf("cirnum=%d\n",cirnum);
            for(x=0;x<cirnum;x++)
			{
                readcir[x]=0;
                if(x<cirnum-1)
                {
                    readcir[x]=5000;
                }
                else
                {
                    readcir[x]=species_reads-5000*(cirnum-1);
                }      
                if(exact==1)
                {
                    for(i=0;i<THREADS_NUM-1;i++)
                    {
                        readnum[i]=readcir[x]/THREADS_NUM;
                    }
                    readnum[THREADS_NUM-1]=readcir[x]-readnum[0]*(THREADS_NUM-1);
                    for(i=0;i<THREADS_NUM;i++)
                    {
                        memset(line02[i],'\0',(input_max*2)*sizeof(char));
                        memset(line04[i],'\0',(input_max*2)*sizeof(char));
                  
                        if(readnum[0]>readnum[THREADS_NUM-1])
						{
							if(truth==1)
							{
								buff_sequence0[i]=(char *)malloc((readnum[0]*5*input_max)*sizeof(char));
								memset(buff_sequence0[i],'\0',(readnum[0]*5*input_max)*sizeof(char));
							}
							else
							{
								buff_sequence0[i]=(char *)malloc((readnum[0]*3*input_max)*sizeof(char));
								buff_sequence0[i]=(char *)malloc((readnum[0]*3*input_max)*sizeof(char));
							}
						}
						else
						{
							if(truth==1)
							{
								buff_sequence0[i]=(char *)malloc((readnum[THREADS_NUM-1]*5*input_max)*sizeof(char));
								memset(buff_sequence0[i],'\0',(readnum[THREADS_NUM-1]*5*input_max)*sizeof(char));
							}
							else
							{
								buff_sequence0[i]=(char *)malloc((readnum[THREADS_NUM-1]*3*input_max)*sizeof(char));
								memset(buff_sequence0[i],'\0',(readnum[THREADS_NUM-1]*3*input_max)*sizeof(char));
							}
						}
                    }   
                    time_t start,end;
                    int realtemp=0;
                    for(i=0;i<THREADS_NUM;i++)
                	{
                        start=time(NULL);
                        tempseed[i]=rand();
                        para[i].ranseed=tempseed[i];
                        
                        para[i].hole_start=holesta;
                        temp_hole=0;
                        for(j=0;j<readnum[i];j++)
                        {
                    	    pass[i][temp_hole]=pass_number(passnumber);
                            realtemp += pass[i][temp_hole];
                            j=j+pass[i][temp_hole]-1;
                            if(j>=readnum[i])
                            {
                                j=j-pass[i][temp_hole];    //j=j-pass[i][temp_hole] not j=j-pass[i][temp_hole]+1 !!
                                realtemp -= pass[i][temp_hole];
                                continue;
                            }
                            temp_hole++;
                            holesta++;
                    	}
                      
                        realtemp=0;

                        para[i].hole_num=temp_hole;
                        para[i].current_reads=readnum[i];
						para[i].passnum=pass[i];
                        para[i].linetwo=line02[i];
                        para[i].linefour=line04[i];
                        para[i].buffseq=buff_sequence0[i];
							
                    	ret=pthread_create(&threads[i],NULL,multi_pass_thread,&para[i]);
                        if(ret != 0)
                        {
                            printf("error!\n");
                        }
                        end=time(NULL);
                    
                    }
                    for(i=0;i<THREADS_NUM;i++)
                    {
                        start=time(NULL);
						ret=pthread_join(threads[i],(void*)&rst);
					
						if(ret != 0)
                        {
                            printf("error2\n\n");
                        }
                        end=time(NULL);
                        start=time(NULL);
						fprintf(fs,"%s",buff_sequence0[i]);
                        end=time(NULL);     
                        free(buff_sequence0[i]);
                    }
                }
            }
        }
		free(buffle);     
		free(ref_name); 
		free(sequence_name);
		free(minus_buffle);
	}
	fclose(fp);
	fclose(fs);
	free(line2);
	free(line4);

	fclose(findex);
	free(path);
	if(lable_method!=3)
	{
        free(err_rand);
        free(short_rand);
        free(mid_rand);
        free(long_rand);
        free(tailshort);
        free(tailmid);
    }
}

//======================================================================================
//
//======================================================================================
void *multi_pass_thread(void *arg)
{ 
	int i;
	int len_arr[200];
	
	struct fun_para *para;
	para=(struct fun_para *)arg;
    int seeda=para->ranseed;
    int *seed= &seeda;
	int threadid=para->id;
   
    float *lengthone_rand=para->lenone;
	float *lengthtwo_rand=para->lentwo;;
	float *lengththree_rand=para->lenthree;
	float *lengthfour_rand=para->lenfour;
	float *lengthfive_rand=para->lenfive;
	float *lengthsix_rand=para->lensix;
	float *lengthseven_rand=para->lenseven;
	float *lengtheight_rand=para->leneight;
	float *lengthnine_rand=para->lennine;
	float *lengthten_rand=para->lenten;
	float *lengtheleven_rand=para->leneleven;
	float *lengthtwelve_rand=para->lentwelve;
	float *lengththirteen_rand=para->lenthirteen;
	//FILE *fs=para->filefs;
	char *buff_sequence=para->buffseq;
	char *sequence_name=para->seqname;
   	float *sublength= para->sublen;

	int len_sequence_name=para->lenseqname; 
    int lable_method=para->lablemethod;
	float *type=para->errtype;
    char *line2=para->linetwo;
    char *line4=para->linefour;
    char *buffle=para->plus;
    char *minus_buffle=para->minus;
    float *subratio=para->subrate;
	float *err_rand=para->errrand;
 	float *qual=para->quality;
   	int len=para->leng;
   	int length_max=para->lengthmax;
    int input_max=para->inputmax;
 	float *short_rand=para->shortrand;
  	float *mid_rand=para->midrand;
   	float *long_rand=para->longrand;
    float *tailshort=para->tail_short;
	float *tailmid= para->tail_mid;
	float *inslength=para->ins_length;
	float *dellength= para->del_length;
	float *ins_less= para->insless;
	float *two_ratio=para->tworatio;
	float *head_ratio=para->headratio;
	float *mid_ratio=para->midratio;
	float *insbias=para->insarray;
	float *delbias=para->delarray;
	float *subbias=para->subarray;
	int buff=para->buffnum;
    int flag_buff=para->flagbuff;
    int species_reads=para->speciesreads;
	int x=para->current_reads;
	int *passarr=para->passnum;
	int holesum=para->hole_num;
    int holenum=para->hole_start;
	float *eachread=para->eacharray;	
	float *unhead=para->unheadarray;
	float *untail=para->untailarray;
    float *avererr=para->averarray;
	int truth=para->iftruth;
	int j=0;
	int end,cir_minus,cir_plus,minus_point,plus_point,point;
	float per;
	int tempmax;
	int every_length;
	int flag_pn;
	int sum_len=0;
	int read_count=0;
	int temp_count=0;        
	int last_len=0;        //1.19 change another holenumber initialise the last_len and temp_count
    int order_flag=0;	
    flag_buff=0;
    int buff_point=0;
    int a=0;
    int pass;
    char name[50];
    int temphh=0;
    time_t timestart,timeend;
    float temprate,rate[3];
	memset(rate,0,3*sizeof(float));
	/*memset(name,'\0',50*sizeof(char));
    sprintf(name,"%d.txt",threadid);
    fs=fopen(name,"w");*/
    //printf("%dholesum=%d\n",threadid,holesum);
    for(i=0;i<200;i++)
    {
        len_arr[i]=0;
    }
    for(a=0;a<holesum;a++)
    {
        temphh += passarr[a];
    }

    //printf("first threadid%d\ttemphh=%d\n",threadid,temphh);
    temphh=0;
    timestart=time(NULL);
	//printf("first!!!threadid=%d\t,timestart=%ld\n",threadid,timestart);
    for(a=0;a<holesum;a++)
    {
		//timestart=time(NULL);
		//printf("threadid=%d\t,timestart=%ld\n",threadid,timestart);
        //printf("a=%d\n",a);
        holenum++; 
		//printf("holenum=%d\n",holenum);
        sum_len=0;
        read_count=0;
        temp_count=0;
        last_len=0;
        order_flag=0;
		temprate=get_eachread(seed,eachread);
		//printf("temprate=%f\n",temprate);
		
		rate[0]=(float)((1-temprate)/(1-avererr[0]));
		//printf("avererr[0]=%f\n",avererr[0]);
		//rate[0]=(float)((1-temprate)/(1-0.1428));
		rate[1]=get_unalign(seed,unhead);
		rate[2]=get_unalign(seed,untail);
		memset(line2,'\0',(input_max*2)*sizeof(char));
        memset(line4,'\0',(input_max*2)*sizeof(char));
        pass=passarr[a];
        temphh += pass;
        //printf("%dpass=%d\n",threadid,pass);
        
        if(pass==0)
        {
            printf("hahahahaerror!!!!!!!!!!!!!\tholenum=%d\n",holenum);
        }

        if(pass==1)
   	    {
    	    every_length=length_distribution(seed,lengthone_rand,lable_method,length_max,input_max);    
		    tempmax=every_length;
        }
        else if(pass==2)
        {
            every_length=length_distribution(seed,lengthtwo_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if(pass==3)
        {   
            every_length=length_distribution(seed,lengththree_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if(pass==4)
        {
            every_length=length_distribution(seed,lengthfour_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if(pass==5)
        {
            every_length=length_distribution(seed,lengthfive_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if((pass>=6) && (pass<=10))
        {
            every_length=length_distribution(seed,lengthsix_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if((pass>=11) && (pass<=15))
        {
            every_length=length_distribution(seed,lengthseven_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if((pass>=16) && (pass<=20))
        {
            every_length=length_distribution(seed,lengtheight_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if((pass>=21) && (pass<=25))
        {
            every_length=length_distribution(seed,lengthnine_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if((pass>=26) && (pass<=30))
        {
            every_length=length_distribution(seed,lengthten_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if((pass>=31) && (pass<=35))
        {
            every_length=length_distribution(seed,lengtheleven_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else if((pass>=36) && (pass<=40))
        {   
            every_length=length_distribution(seed,lengthtwelve_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
        else
        {   
            every_length=length_distribution(seed,lengththirteen_rand,lable_method,length_max,input_max);
            tempmax=every_length;
        }
		if(tempmax>len)   //7.15 add,maybe the length of chromosome is so short
		{
			a--;
			holenum--;
			temphh -= pass;
			continue;
		}
        if(pass==1)
	    {
	        sum_len = tempmax;
					
        }
	    else if(pass==2)
	    {
            sum_len=0;
		    sum_len += (tempmax+44);	
		    per=length_ratio(seed,two_ratio);
		    len_arr[1]=(int)(tempmax*per);
		    if(len_arr[1]<0)
			{
				printf("len_arr[1]=%d\tper=%f\n",len_arr[1],per);
			}
			if(len_arr[1]==0)
		    {
		        len_arr[1]=1;
		    }
		    //printf("holenum=%d\n",holenum);
		    //printf("length=%d\n",len_arr[1]);
		    //if(len_arr[1]==0)
		    //{
		    	//printf("per=%f\n",per);
		    //}
		    sum_len += len_arr[1];
		    if(sum_len > (length_max+2000))
		    {
                //printf("error3\tthreadid=%d\tholenum=%d\n",threadid,holenum);    
                a--;
                holenum--;
                temphh -= pass;
                continue;
            }
		}
	    if(pass<=2)
        {
	        every_length=tempmax;
            end=last_len+every_length-1;
            if ((float)(rand_r(seed))/(float)(RAND_MAX)<=0.5)           //minus strand
            {
                flag_pn=0;
                point=(int)(rand_r(seed)/(float)(RAND_MAX)*(len-every_length)+every_length);
                buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
            }
            else                             //plus strand
            {
                flag_pn=1;
                    point=(int)(rand_r(seed)/(float)(RAND_MAX)*(len-every_length));
                buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
            }
        
            if(pass==2)
            {
                last_len += (tempmax+44);
		        every_length=len_arr[1];
		        end=last_len+every_length-1;
                if (flag_pn==1)       //the first one is plus 
                {
			        flag_pn=0;                                             //1.17 correct the bug!! should change the flag_pn
                    point=point+tempmax-1;          //!!not +every_length because every_length is changed !!
                    buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
		        }
                else               //the first one is minus 
                {
			        flag_pn=1;
                    point=point-tempmax+1;
                    buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);                  
                }
            }
        }
        else if(pass>2)
        {
            sum_len=0;
            if ((float)(rand_r(seed))/(float)(RAND_MAX)<=0.4286)
		    {
		        order_flag=0;	
			    len_arr[0]=tempmax;
		        //printf("length=%d\n",len_arr[0]);
                sum_len += (len_arr[0]+44);
	            for(j=1;j<=(pass-2);j++)
		        {
		            per=length_ratio(seed,mid_ratio);
			        len_arr[j]=(int)(per*tempmax);
					if(len_arr[j]<0)
					{
						printf("len_arr[%d]=%d\tper=%f\n",j,len_arr[j],per);
					}
			        if(len_arr[j]==0)
			        {
			            len_arr[j]=1;
			        }
					sum_len += (len_arr[j]+44);
		        }
		        per=length_ratio(seed,head_ratio);
		        len_arr[pass-1]=(int)(per*tempmax);
				if(len_arr[pass-1]<0)
				{
					printf("len_arr[%d]=%d\tper=%f\n",pass-1,len_arr[pass-1],per);
				}
		        if(len_arr[pass-1]==0)
		        {
			        len_arr[pass-1]=1;
		        }   
		        sum_len += len_arr[pass-1];
	        }
	        else
	        {
                order_flag=1;
    	        per=length_ratio(seed,head_ratio);
		        len_arr[0]=(int)(per*tempmax);
				if(len_arr[0]<0)
				{
					printf("len_arr[0]=%d\tper=%f\n",len_arr[0],per);
				}
		        if(len_arr[0]==0)
		        {
			        len_arr[0]=1;
		        }
		        sum_len += (len_arr[0]+44);
		        len_arr[1]=tempmax;
		        sum_len += (len_arr[1]+44);
		        for(j=2;j<=(pass-2);j++)
		        {
			        per=length_ratio(seed,mid_ratio);
			        len_arr[j]=(int)(per*tempmax);
					if(len_arr[j]<0)
					{
						printf("len_arr[%d]=%d\tper=%f\n",j,len_arr[j],per);
					}
			        if(len_arr[j]==0)
			        {
			    	    len_arr[j]=1;
			        }
			        sum_len += (len_arr[j]+44);
		        }
		        per=length_ratio(seed,head_ratio);
		        len_arr[pass-1]=(int)(per*tempmax);
				if(len_arr[pass-1]<0)
				{
					printf("len_arr[%d]=%d\tper=%f\n",pass-1,len_arr[pass-1],per);
				}
		        if(len_arr[pass-1]==0)
		        {
			        len_arr[pass-1]=1;
		        }
		        sum_len += len_arr[pass-1];
	        }
	        if(sum_len > (length_max+2000))
	        {
                a--;
                holenum--;
                temphh -= pass;
                continue;
            }
            if(sum_len<0)
            {
                printf("sum_len<0 error\tholenum=%d\tx=%d\n",holenum);
            }
		    //printf("pass=%d\n",pass);
           
			if(order_flag==0)
            {                          //1.24 correct the bug:cannot use (float)(rand())/(float)(RAND_MAX)<=0.4286 again!!!!
							
			    every_length=len_arr[0];
			    end=last_len+every_length-1;
			    if ((float)(rand_r(seed))/(float)(RAND_MAX)<=0.5)           //minus strand
			    {
				    flag_pn=0;
	    		        point=(int)(rand_r(seed)/(float)(RAND_MAX)*(len-every_length)+every_length); 
				    buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                }
                else                             //plus strand
			    {
			        flag_pn=1;
			            point=(int)(rand_r(seed)/(float)(RAND_MAX)*(len-every_length));                           
				    buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
			    }
			    read_count++;
			    last_len += (every_length+44);
                if(pass % 2 == 1)        //passnumber is odd
                {
                    cir_minus=(pass-1)/2;
                    cir_plus=(pass-1)/2;
                    if(flag_pn==0)      //the first one is minus
                    {
                        for(i=0;i<cir_plus;i++)
                        {
				    	    every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    point=point-tempmax+1;           //start position of plus strand 
						    flag_pn=1;
                            buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                            last_len += (every_length+44);
						    read_count++;
                            every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    point=point+tempmax-1;
						    flag_pn=0;
						    temp_count++;
						    buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                            last_len += (every_length+44);
						    read_count++;
					    }
				    }
                    else          //the first one is plus 
				    {
                        for(i=0;i<cir_minus;i++)
                        {
				    	    every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    point=point+tempmax-1;          //start of the minus strand
						    flag_pn=0;
                            buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                            read_count++;
						    last_len += (every_length+44);   
									
						    every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    flag_pn=1;
						    point=point-tempmax+1;          //start position of the plus strand
						    buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
						    read_count++;
						    last_len += (every_length+44);
					    }
			        }		
			    }
                else
			    {    
                    if(flag_pn==0)     //the first one is minus 
                    {
                        cir_plus=pass/2;
                        cir_minus=pass/2-1;
                        for(i=0;i<(cir_plus-1);i++)
                        {
						    every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    flag_pn=1;
						    point=point-tempmax+1;   //transfer the position
                            buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads); 
                            read_count++;
						    last_len += (every_length+44);
						    flag_pn=0;
						    every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    point=point+tempmax-1;
						    buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
						    read_count++;
						    last_len += (every_length+44);
					    }
					    every_length=len_arr[read_count];
					    end=last_len+every_length-1;
					    flag_pn=1;
					    point=point-tempmax+1;
                        buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
				    }
                    else                //the first one is plus 
                    {
                        cir_minus=pass/2;
                        cir_plus=pass/2-1;
                        for(i=0;i<(cir_minus-1);i++)
                        {
					        every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    flag_pn=0;
						    point=point+tempmax-1;
                            buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                            read_count++;
						    last_len += (every_length+44);
						    flag_pn=1;
						    every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    point=point-tempmax+1;
						    buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
						    read_count++;
						    last_len += (every_length+44);
					    }
					    flag_pn=0;
					    every_length=len_arr[read_count];
					    end=last_len+every_length-1;
					    point=point+tempmax-1;
                        buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
				    }
			    }
		    }
            else          //the lengthest one is a mid one 
            {
			    if ((float)(rand_r(seed))/(float)(RAND_MAX)<=0.5)           //minus strand       
			    {
				    flag_pn=0;
					    point=(int)(rand_r(seed)/(float)(RAND_MAX)*(len-every_length)+every_length); 
					    minus_point=point;
					    plus_point=point-tempmax+1;
			    }
			    else                             //plus strand
			    {
		            flag_pn=1;
			            point=(int)(rand_r(seed)/(float)(RAND_MAX)*(len-every_length)); 
					    plus_point=point;
					    minus_point=point+tempmax-1;
                }
			    if(pass % 2==1)
			    {
			        cir_minus=(pass-1)/2;
				    cir_plus=(pass-1)/2;
                    if(flag_pn==0)      //the first one is minus
                    {
					    flag_pn=0;
					    every_length=len_arr[0];
					    end=last_len+every_length-1;
					    point=minus_point;
			    	    buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                        read_count++;
					    last_len += (every_length+44);
								
			    	    flag_pn=1;
					    point=plus_point;
					    every_length=tempmax;         //1.19 the lengthest one is the second one defalut 
					    end=last_len+every_length-1;
					    buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
					    read_count++;
					    last_len += (every_length+44);
					    for(i=0;i<(cir_plus-1);i++)
                        {
						    point=minus_point;
						    flag_pn=0;
						    every_length=len_arr[read_count];
						    end=last_len+every_length-1;
						    buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
						    read_count++;
						    last_len += (every_length+44);	
										
    						point=plus_point;
	    					flag_pn=1;
                            every_length=len_arr[read_count];
                            end=last_len+every_length-1;
                             buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                            read_count++;
						    last_len += (every_length+44); 
	    				}
		    			flag_pn=0;
			    		every_length=len_arr[read_count];
				    	end=last_len+every_length-1;
					    point=minus_point;
                        buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
			    	}
    				else
	    			{
			    		flag_pn=1;
				    	every_length=len_arr[0];
    					end=last_len+every_length-1;
	    				point=plus_point;
				    	buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
				    	read_count++;
	    				last_len += (every_length+44);
			    				
				    	flag_pn=0;
    					point=minus_point;
	    				every_length=tempmax;        //the lengthest one 
	    				end=last_len+every_length-1;
	    				buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
		    			read_count++;
			    		last_len += (every_length+44);
					    for(i=0;i<(cir_minus-1);i++)
    					{
	    					point=plus_point;
		    	       	    flag_pn=1;
			    		    every_length=len_arr[read_count];
				    		end=last_len+every_length-1;
                            buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
		    				read_count++;
    						last_len += (every_length+44);
	    						
		    				point=minus_point;
			    	        flag_pn=0;
				            every_length=len_arr[read_count];
					        end=last_len+every_length-1;
		    		        buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
			    			read_count++;
    						last_len += (every_length+44);
	    					//printf("last_len=%d\n",last_len);
		    			}
    
	    				flag_pn=1;
		    			point=plus_point;
			    		every_length=len_arr[read_count];
    					end=last_len+every_length-1;
	    				buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
	    			}
		    	}
                else        
                {
    	    		if(flag_pn==0)
			        {
    					cir_minus=(pass-2)/2;
	    				flag_pn=0;
		    			every_length=len_arr[0];
			    		end=last_len+every_length-1;
    					point=minus_point;
			    		buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                        read_count++;
					    last_len += (every_length+44);
	    					
		    			flag_pn=1;
			    		point=plus_point;
    			        every_length=tempmax;         //1.19 the lengthest one is the second one defalut
	    		        end=last_len+every_length-1;
			            buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
				    	read_count++;
    					last_len += (every_length+44);
		    			for(i=0;i<cir_minus;i++)
                        {
				    		point=minus_point;
				            flag_pn=0;
    						every_length=len_arr[read_count];
	    			        end=last_len+every_length-1;
		    		        buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
			    		    read_count++;
				    		last_len += (every_length+44);
	    					temp_count++;
		    		    	if(temp_count==(cir_minus*2-1))
			    			{
				    			continue;
					    	}
    						flag_pn=1;
	    					point=plus_point;
                           	every_length=len_arr[read_count];
                            end=last_len+every_length-1;
                            buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                            read_count++;
    						last_len += (every_length+44);
		    				temp_count++;
			    		}
				    			
    					every_length=len_arr[read_count];
	    				end=last_len+every_length-1;
		    			flag_pn=1;
		    			point=plus_point;
                       buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
				    }
    				else
	    	        {
				    	cir_plus=(pass-2)/2;
    					flag_pn=1;
	    			    every_length=len_arr[read_count];
	    				end=last_len+every_length-1;
		    			point=plus_point;
    				    buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
	    				read_count++;
		    			last_len += (every_length+44);
    
	    				flag_pn=0;
		    	        point=minus_point;
			            every_length=tempmax;        //the lengthest one
			            end=last_len+every_length-1;
	    				buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
		    			read_count++;
			    		last_len += (every_length+44);	
				    	for(i=0;i<cir_plus;i++)
                        {
    						flag_pn=1;
	    					every_length=len_arr[read_count];
		    			    end=last_len+every_length-1;
			    			point=plus_point;
	    				    buff_point=plus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
		    				read_count++;
			    			last_len += (every_length+44);
    						temp_count++;
	    					if(temp_count==(cir_plus*2-1))
		    				{
			    				continue;
				    		}
    						flag_pn=0;
	    					point=minus_point;
                            every_length=len_arr[read_count];
			    	    	end=last_len+every_length-1;
                            buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
                            read_count++;
				    		last_len += (every_length+44);
    						temp_count++;
	    				}
    
	    				flag_pn=0;
		    			every_length=len_arr[read_count];
			    		end=last_len+every_length-1;
				    	point=minus_point;
                        buff_point=minus_bp(truth,rate,insbias,delbias,subbias,seed,sublength,holenum,last_len,end,sum_len,last_len,x,flag_pn,every_length,point,buff_sequence,sequence_name,buff_point,len_sequence_name,lable_method,type,line2,line4,minus_buffle,subratio,err_rand,qual,len,length_max,input_max,short_rand,mid_rand,long_rand,tailshort,tailmid,inslength,dellength,ins_less,buff,flag_buff,species_reads);
			        }
			    }
		    }
        }
	}
    timeend=time(NULL);
}
//======================================================================================
//function:buff_point_bp--add error position and type to result
//======================================================================================
int buff_point_bp(int position_err,int buff_point,char a1,char a2,char *buff_sequence)
{
	char s[6];
	sprintf(s,"%d",position_err);
	buff_sequence[buff_point]='|';
	buff_point++;
	strcat(&buff_sequence[buff_point],s);
	buff_point+=strlen(s);
	//printf("position_err=%d\n",position_err);
	//printf("strlen(position_err)=%d\n",strlen(s));
		
	buff_sequence[buff_point]=':';
	buff_point++;
	buff_sequence[buff_point]=a1;
	buff_point++;
	buff_sequence[buff_point]=':';
	buff_point++;
	buff_sequence[buff_point]=a2;
	buff_point++;
	return buff_point;
}
//======================================================================================
//function:bp--simulate a read
//======================================================================================
int bp(int truth,float *rate,float *insbias,float *delbias,float *subbias,int *seed,int holenum,float sublength[],int sum_len,int last_len,int lable_method,float type[],char *line2,char *line4,char *buffle,float subratio[],char *buff_sequence,int buff_point,int length,float *err_rand,int position,float qual[],int len_seq,int length_max,int input_max,float *short_rand,float *mid_rand,float *long_rand,float *tailshort,float *tailmid,float inslength[],float dellength[],float *ins_less)
{
    //printf("uuuuuuuu sum_len=%d\tlast_len=%d\tlength=%d\n",sum_len,last_len,length);
	int a,turn,add_mark=0;
	float random,err_random,p;
	int quality();
    int error_length();
    float ccsqual[100];
    float subqual[20];
    int i=0;
    int del_len,ins_len,sub_len;
    char ins_error[5];        // 11.19 the detail insertion when ins_len < 4
    int temp;
	char convert();
	char key[5];
	memset(key,'\0',5*sizeof(char));
    int flag=0;
    int cir_num=0;
    int lentemp;
	float p0,p1,p2;

    for(i=0;i<20;i++)
    {
        subqual[i]=pow(10,(i/(-10.0)));
    }
	for(a=0;(a+add_mark<length)&&((a+position)<len_seq);a++)          //length is every_length of every read,len_seq is the length of reference
	{
		if(lable_method==4)
		{
            turn=quality(seed,sum_len,(a+add_mark+last_len),err_rand,lable_method,input_max,length_max,length,short_rand,mid_rand,long_rand,tailshort,tailmid);  
      	}
		else
		{
			turn=28;
		}
        if((a+position+2)<len_seq)
		{
					key[0]=buffle[a+position];
            		key[1]=buffle[a+position+1];
					key[2]=buffle[a+position+2];
					if(((key[0] =='A')||(key[0]=='T')||(key[0]=='G')||(key[0]=='C')) && ((key[1] =='A')||(key[1]=='T')||(key[1]=='G')||(key[1]=='C')) && ((key[2] =='A')||(key[2]=='T')||(key[2]=='G')||(key[2]=='C')))
					{
						flag=ffun_per(rate,key,seed,insbias,delbias,subbias,a+add_mark,sum_len);
					}
					else
					{
						flag=4;
					}
				}
				else
				{
					if((buffle[a+position]=='A') || (buffle[a+position]=='T') || (buffle[a+position]=='G') || (buffle[a+position]=='C'))
					{
						
						random=(float)(rand_r(seed))/(float)(RAND_MAX);
						if(random < 0.15)
						{
							if(random<=0.6)
							{
								flag=1;
							}
							else if(random<=0.9)
							{
								flag=2;
							}
							else
							{
								flag=3;
							}
						}
						else
						{
							flag=0;
						}
					}
					else
					{
						flag=4;
					}
		}
				
			if(flag==3)                            //substitution
			{
            	sub_len=error_length(seed,sublength);
				if(((a+position+sub_len-1)<len_seq) && ((a+add_mark+sub_len-1)<length))
				{
					for(i=0;i<sub_len;i++)
					{
						line2[a+add_mark]=sub_all(seed,buffle[a+position],subratio[0],subratio[1],subratio[2],subratio[3],subratio[4],subratio[5],subratio[6],subratio[7],subratio[8],subratio[9],subratio[10],subratio[11]);
						if(truth==1)
						{
							buff_point=buff_point_bp(a+add_mark,buff_point,buffle[a+position],line2[a+add_mark],buff_sequence);
						}
						line4[a+add_mark]=(char)(turn+33);
						a++;
					}
					a--;
				}
				else
				{
					if((length-a-add_mark) <= (len_seq-a-position))
					{
						temp=length-a-add_mark;
					}
					else
					{
						temp=len_seq-a-position;
					}
					for(i=0;i<temp;i++)
					{
						line2[a+add_mark]=sub_all(seed,buffle[a+position],subratio[0],subratio[1],subratio[2],subratio[3],subratio[4],subratio[5],subratio[6],subratio[7],subratio[8],subratio[9],subratio[10],subratio[11]);
						if(truth==1)
						{
							buff_point=buff_point_bp(a+add_mark,buff_point,buffle[a+position],line2[a+add_mark],buff_sequence);
						}
						line4[a+add_mark]=(char)(turn+33);
						a++;
					}
					a--;
				}
			}
			else if(flag==1)                    //insertion
			{
                ins_len=error_length(seed,inslength);
        
				if((a+add_mark+ins_len-1)<length)
                {
                	for(i=0;i<ins_len;i++)
                    {
                        if((a+position+2)<len_seq)
                        {
						    line2[a+add_mark]=convert(seed,key,ins_less);
                        }
                        else
                        {
                        	line2[a+add_mark]=ins_all(seed);
                        }
						if(truth==1)
						{
                        	buff_point=buff_point_bp(a+add_mark,buff_point,'-',line2[a+add_mark],buff_sequence);
                        }
						line4[a+add_mark]=(char)(turn+33);
                        
                        add_mark++;
                    }
        		}
                else
                {
                  	for(i=0;i<(length-a-add_mark);i++)
                    {
                        if((a+position+2)<len_seq)
						{
							line2[a+add_mark]=convert(seed,key,ins_less);
						}
						else
						{
							line2[a+add_mark]=ins_all(seed);
                        }
						if(truth==1)
						{
							buff_point=buff_point_bp(a+add_mark,buff_point,'-',line2[a+add_mark],buff_sequence);
                        }
						line4[a+add_mark]=(char)(turn+33);
                        add_mark++;
                    }
            	}
				if((a+add_mark)<length)
				{
                    line2[a+add_mark]=buffle[a+position];
					if(lable_method==4)
					{
						turn=quality(seed,sum_len,(a+add_mark+last_len),err_rand,lable_method,input_max,length_max,length,short_rand,mid_rand,long_rand,tailshort,tailmid);
					}
					else
					{
						turn=28;
					}
					line4[a+add_mark]=(char)(turn+33);
				}
			}
			else if(flag==2)                                    //deletion
			{
				del_len=error_length(seed,dellength); 
                if((a+position+del_len-1)<len_seq)
              	{
                  	for(i=0;i<del_len;i++)
                   	{
						if(truth==1)
						{
                        	buff_point=buff_point_bp(a+add_mark,buff_point,buffle[a+position],'-',buff_sequence);
					  	}
						add_mark--;
                       	a++;
                	}
                    a--;
				}
              	else
                {
                   	for(i=0;i<(len_seq-a-position);i++)
                    {
						if(truth==1)
						{
                        	buff_point=buff_point_bp(a+add_mark,buff_point,buffle[a+position],'-',buff_sequence);
                        }
						add_mark--;
                        a++;
               		}
                    a--;
            	}
          	}
			else if(flag==0)
			{   
            	line2[a+add_mark]=buffle[a+position];
				line4[a+add_mark]=(char)(turn+33); 
			}
			else if(flag==4)
			{
				 random=(float)(rand_r(seed))/(float)(RAND_MAX);
				 if(random<=0.25)
				 {
				 	line2[a+add_mark]='A';
				 }
				 else if(random<=0.5)
				 {
				 	line2[a+add_mark]='T';
				 }
				 else if(random<=0.75)
				 {
				 	line2[a+add_mark]='G';
				 }
				 else
				 {
				 	line2[a+add_mark]='C';
				 }
				 line4[a+add_mark]=(char)(turn+33);
			}
		}
	
	return buff_point;
}
//======================================================================================
//function:input_right--judge the inputed parameters right or wrong
//======================================================================================
void input_right(int length,int reads,char index[],int lable,int gap_lable,int gap,int length_max,char method[],int exact,float length_rand[])
{
	if((length>length_max)&&(exact==0)) 
	{
		printf("Warning:The max read's length is:%d bps under the config file, so your reads' length can't exceed the max length!\n",length_max);
		exit(1);
	}

	printf("the reads number: %d.\n",reads);
	printf("the index file: %s.\n",index);
	printf("the method: %s.\n",method);
	if(exact==1)
	{
		if(length_rand[0]==2)
		{
			printf("this simulation-config-file can't allow exact simulation!\n");
			exit(1);
		}
		printf("this simulation is exact.\n");
	}
	else
		printf("the read length: %dbp.\n",length);

	if(lable==0) 
	{
		if(gap_lable==1) 
		{
			printf("Warning:You want to get the pair-end file,but this is the single model.Please use the parameter '-e 1'.\n");
			exit(1);
		}
		printf("this is the single model.\n");
	}
	if(lable==1) 
	{ 
		printf("this is a pair-end model.\n");
		printf("the width of pairs %d.\n",gap);
	}
}
//======================================================================================
//function:minus_strand--get the "-" strand according to the "+" strand
//======================================================================================
void minus_strand(char *buffle,char *minus_buffle,int len)
{
	int circle_f;
	for(circle_f=0;circle_f<len;circle_f++)
	{
		if (buffle[circle_f]=='A')
			minus_buffle[len-1-circle_f]='T';
		else if (buffle[circle_f]=='T')
			minus_buffle[len-1-circle_f]='A';
		else if (buffle[circle_f]=='C')
			minus_buffle[len-1-circle_f]='G';
		else if(buffle[circle_f]=='G')
			minus_buffle[len-1-circle_f]='C';
		else 
			minus_buffle[len-1-circle_f]=buffle[circle_f];   //7.12 add 
		
	}
}
//======================================================================================
//function:get_qv_type--get error config and return the simulation type
//======================================================================================
int get_qv_type(float *avererr,float *unhead,float *untail,float *eachread,float *insbias,float *delbias,float *subbias,float sublength[],float lengththirteen_rand[],float lengthtwelve_rand[],float lengthten_rand[],float lengtheleven_rand[],float lengthfive_rand[],float lengthsix_rand[],float lengthseven_rand[],float lengtheight_rand[],float lengthnine_rand[],float two_ratio[],float head_ratio[],float mid_ratio[],FILE *fc,float type[],float subratio[],float sd_ratio[],float *err_rand,int lable_method,float length_rand[],int exact,float *short_rand,float *mid_rand,float *long_rand,float *tailshort,float *tailmid,float inslength[],float dellength[],float ins_less[],float lengthone_rand[],float lengthtwo_rand[],float lengththree_rand[],float lengthfour_rand[],float passnumber[])
{
    
	int circle_f,turn;
	char *qv_config;
	qv_config=(char *)malloc(2000000*sizeof(char));
	memset(qv_config,'\0',2000000*sizeof(char));

	float value[20];              //11.3 change to 20 from 50 
    float ccsvalue[100];
	int length_max=0;
	int get_value();
    turn =0;
    int short_max=0;
    int mid_max=0;
    int long_max=0;
    int tailshort_max=0;
    int tailmid_max=0;
    int test1=0;
    int test2=0;
    int test3=0;
	while(fgets(qv_config,1000000,fc)!=NULL)
	{
		qv_config[strlen(qv_config)-1]='\0';
        if((strstr(qv_config,"-short_rand")!=NULL))
        { 
            get_value(qv_config,value);
            for(circle_f=0;circle_f<20;circle_f++)
            {
                short_rand[length_max*20+circle_f]=value[circle_f];
            }
            length_max++;
            continue;
        }
        else if((strstr(qv_config,"tailshort_rand")!= NULL))
        {
            get_value(qv_config,value);
            for(circle_f=0;circle_f<20;circle_f++)
            {
                tailshort[tailshort_max*20+circle_f]=value[circle_f];
            }
            tailshort_max++;
            continue;
        }
        else if((strstr(qv_config,"tailmid_rand")!= NULL))
        {
            get_value(qv_config,value);
            for(circle_f=0;circle_f<20;circle_f++)
            {
                tailmid[tailmid_max*20+circle_f]=value[circle_f];
            }
            tailmid_max++;
            continue;
        }
		if(strstr(qv_config,"pacbio_type")!=NULL)
		{ 
			get_value(qv_config,type);
			get_value(qv_config,subratio);
			continue;
        }
        if((strstr(qv_config,"inslength")!=NULL))
        {
            test1=get_value(qv_config,inslength);
            continue;
        }
        if((strstr(qv_config,"dellength")!= NULL))
        {
            test1=get_value(qv_config,dellength);
            continue;
        }
		if((strstr(qv_config,"sublength")!= NULL ))
		{
			test2=get_value(qv_config,sublength);
			continue;
		}
        if((strstr(qv_config,"ins_detail")!= NULL))
        {
            test3=get_value(qv_config,ins_less);
            continue;
        }
        if((strstr(qv_config,"passtwo_ratio")!= NULL ))
        {
            test1=get_value(qv_config,two_ratio);
			continue;
        }
        else if((strstr(qv_config,"passmore_ratio_head")!= NULL ))
        {
            test1=get_value(qv_config,head_ratio);
			continue;
        }
        else if((strstr(qv_config,"passmore_ratio_mid")!= NULL ))
        {
            test1=get_value(qv_config,mid_ratio);
			continue;
        }
		//if((exact==0)&&(strstr(qv_config,"sd_ratio")!=NULL))
        //{
		//	get_value(qv_config,sd_ratio);
		//}
        //else if((exact==1)&&(strstr(qv_config,"length_rand")!=NULL))
		//{
		//	length_max=get_value(qv_config,length_rand);
        //    if(lable_method==1)
		//		length_max=length_max+39;
            //else if(lable_method==4)               //10.13change : to agree with in the length_distribution( ) length=1 and length=772
                //length_max=length_max;
        //    else if(lable_method==5)
        //        length_max=length_max+771;
		//	else if(lable_method==2)
		//		length_max=length_max+29;
        //    continue;
		//}
        else if((exact==1)&&(strstr(qv_config,"passone_length_rand")!=NULL))
        {
		
            get_value(qv_config,lengthone_rand);
            continue;
        }
        else if((exact==1)&&(strstr(qv_config,"passtwo_length_rand")!=NULL))
        {
            get_value(qv_config,lengthtwo_rand);
           
			continue;
		}
        else if((exact==1)&&(strstr(qv_config,"passthree_length_rand")!=NULL))
        {
           	get_value(qv_config,lengththree_rand);
			continue;
		}
        else if((exact==1)&&(strstr(qv_config,"passfour_length_rand")!=NULL))
        {   
            get_value(qv_config,lengthfour_rand);
        	continue;
		}
        else if((exact==1)&&(strstr(qv_config,"passfive_length_rand")!=NULL))
        {
            get_value(qv_config,lengthfive_rand);
        	continue;
		}
        else if((exact==1)&&(strstr(qv_config,"passsix_length_rand")!=NULL))
        {
            get_value(qv_config,lengthsix_rand);
       		continue;
	   	}
        else if((exact==1)&&(strstr(qv_config,"passseven_length_rand")!=NULL))
        {
            get_value(qv_config,lengthseven_rand);
            continue;
        }
        else if((exact==1)&&(strstr(qv_config,"passeight_length_rand")!=NULL))
        {
            get_value(qv_config,lengtheight_rand);
            continue;
        }
        else if((exact==1)&&(strstr(qv_config,"passnine_length_rand")!=NULL))
        {
            get_value(qv_config,lengthnine_rand);
            continue;
        }
        else if((exact==1)&&(strstr(qv_config,"passten_length_rand")!=NULL))
        {
            get_value(qv_config,lengthten_rand);
            continue;
        }
        else if((exact==1)&&(strstr(qv_config,"passeleven_length_rand")!=NULL))
        {
            get_value(qv_config,lengtheleven_rand);
            continue;
        }
        else if((exact==1)&&(strstr(qv_config,"passtwelve_length_rand")!=NULL))
        {
            get_value(qv_config,lengthtwelve_rand);
            continue;
        }
        else if((exact==1)&&(strstr(qv_config,"passthirteen_length_rand")!=NULL))
        {
            get_value(qv_config,lengththirteen_rand);
            continue;
        }
        else if((strstr(qv_config,"passnumber")!= NULL ))
	    {
            get_value(qv_config,passnumber);
            continue;
		}
		else if((strstr(qv_config,"fourper_head")!= NULL ))
		{
			get_value(qv_config,insbias);
			continue;
		}
		else if((strstr(qv_config,"fourper_mid")!= NULL ))
		{
			get_value(qv_config,delbias);
			continue;
		}
		else if((strstr(qv_config,"fourper_tail")!= NULL ))
		{
			get_value(qv_config,subbias);
			continue;
		}
		else if((strstr(qv_config,"eachread")!= NULL ))
		{
			get_value(qv_config,eachread);
			continue;
		}
		else if((strstr(qv_config,"unalign_head")!= NULL ))
		{
			get_value(qv_config,unhead);
			continue;
		}
		else if((strstr(qv_config,"unalign_tail")!= NULL))
		{
			get_value(qv_config,untail);
			continue;
		}
		else if((strstr(qv_config,"aver_err")!= NULL))
		{
			get_value(qv_config,avererr);
			continue;
		}
    }
	
	fclose(fc);
	if(lable_method==5)
	{
		length_max=79542;
	}
	return length_max;
}
///======================================================================================
//function:get_value
//======================================================================================
int get_value(char qv_config[],float need[])
{
	int circle_f,length;
	char *part[2],*value[100000];
	split(part,qv_config,"=");
	length=split(value,part[1],":");
    //printf("length=%d\n",length);
	for(circle_f=0;circle_f<length;circle_f++)
	{
		need[circle_f]=atof(value[circle_f]);
        
    }
	return length;
}
//======================================================================================
//function:species_number--how many reads in every species and return how many species
//======================================================================================
int species_number(FILE *fp,int everynumber[],char name1[2000][200],int reads)
{
	int i;
    int count=0;
	int total=0;
	float percentage1[2000][1];
	float phase[2000],random;
	int circle1,circle2;
	char species_percent[MAX_LINE],*species_name[MAX_LINE];
	while(fgets(species_percent,MAX_LINE,fp)!=NULL)
	{
		species_percent[strlen(species_percent)-1]='\0';
		count=split(species_name,species_percent,"\t");
        //printf("count=%d",count);
		strcpy(name1[total],species_name[0]);
		percentage1[total][0]=atof(species_name[1]);
		everynumber[total]=0;
		total++;
	}
    //printf("total=%d",total);    
	for(circle1=0;circle1<total-1;circle1++)
	{
        
		phase[circle1]=0;
		for(circle2=0;circle2<=circle1;circle2++)
		{
			phase[circle1]=phase[circle1]+percentage1[circle2][0];
		}
    }
	for(circle1=0;circle1<reads;circle1++)
	{
		random=(float)rand()/(float)(RAND_MAX);
		
		if(random>phase[total-2])
		{
			everynumber[total-1]++;
		}
		else if(random<=phase[0])
		{
			everynumber[0]++;
		}
		else
		{
			for(circle2=0;circle2<total-2;circle2++)
			{
				if(random>phase[circle2]&&random<=phase[circle2+1])
				{
					everynumber[circle2+1]++;
					continue;
				}
			}
		}
	}
	return total;
}
//======================================================================================
//function:first_line--add the reference name,start position and "+" "-"strain to the sequence
//======================================================================================
int first_line(int holenum,int start,int end,int point,char *buff_sequence,char *sequence_name,int buff_point,int len_sequence_name,int flag_pn)//flag=1:the first in pair-end and 2
{
	char temp1[10];
	char temp2[10];
	char temp3[10];
	char prints[20];
    strcat(&buff_sequence[buff_point],sequence_name);
	//printf("sequence_name=%s\n",sequence_name);
	//printf("buff_sequence[buff_point]=%s\n",&buff_sequence[buff_point]);
	//printf("buff_point=%d\n",buff_point);
	buff_point=buff_point+len_sequence_name;
	sprintf(temp1,"%d",holenum);
	strcat(&buff_sequence[buff_point],temp1);
	buff_point += strlen(temp1);
	buff_sequence[buff_point]='/';
	buff_point++;
	sprintf(temp2,"%d",start);
	strcat(&buff_sequence[buff_point],temp2);
	buff_point += strlen(temp2);
	buff_sequence[buff_point]='_';
	buff_point++;
	sprintf(temp3,"%d",end);
	strcat(&buff_sequence[buff_point],temp3);
	buff_point += strlen(temp3);

	buff_sequence[buff_point]='|';
	buff_point++;
	sprintf(prints,"%d",point);
	strcat(&buff_sequence[buff_point],prints);
	buff_point=buff_point+strlen(prints);
	
	//printf("point=%d\n",point);
	//printf("strlen(point)=%d\n",strlen(prints));
	
	
	buff_sequence[buff_point]='|';
	if (flag_pn==1)
		buff_sequence[buff_point+1]='+';
	else
		buff_sequence[buff_point+1]='-';
	buff_point+=2;
	//printf("%s\n",buff_sequence);
	return buff_point;
}
//======================================================================================
//function:remain_lines--add simulated sequence and it QV to the result:buff_sequence
//======================================================================================
int remain_lines(char *buff_sequence,int buff_point,char *line2,char *line4,int every_length)
{
	buff_sequence[buff_point]='\n';
	buff_point++;
	strcat(&buff_sequence[buff_point],line2);
	buff_point+=every_length;
	buff_sequence[buff_point]='\n';
	buff_sequence[buff_point+1]='+';
	buff_sequence[buff_point+2]='\n';
	buff_point+=3;
	strcat(&buff_sequence[buff_point],line4);
	buff_point+=every_length;
	buff_sequence[buff_point]='\n';
	buff_point++;
	return buff_point;
}
//======================================================================================
//function:quality--return an error value
//======================================================================================
int quality(int *seed,int sum_len,int pos_base,float *err_rand,int lable_method,int input_max,int length_max,int every_length,float *short_rand,float *mid_rand,float *long_rand,float *tailshort,float *tailmid)  //pos_base:from zero 
{
    //10.30 add parameter every_length and *tail_rand
	float random;
	int circle,step,turn;
    int temp;                         //10.21 change temp is the position in the length_max
    int boundary;
    int templength;
	step=20;
	random=(float)(rand_r(seed))/(float)(RAND_MAX);
	boundary=sum_len-2000;
	
	if(pos_base<boundary)
    {
    	temp=(int)(1.0*pos_base/input_max*length_max);
		if(random<=short_rand[step*temp+10])
        {
        	for(circle=(step*temp+10);circle>step*temp;circle--)
            {
            	if((random>short_rand[circle-1])&&(random<=short_rand[circle]))
                {
                	turn=circle-step*temp;
                    return turn;
                }
            }
           	return 0;
         }
         else
         {
         	for(circle=(step*temp+11);circle<step*(temp+1);circle++)
            {
            	if((random>short_rand[circle-1])&&(random<=short_rand[circle]))
                {
                	turn=circle-step*temp;
                	return turn;
                }
            }
        }
	}  
    else
    {
        temp=pos_base-boundary;         //200 is the size of the tail_rand
        if(random<=tailshort[step*temp+10])
	    {
	 		for(circle=(step*temp+10);circle>step*temp;circle--)
		    {
			  	if((random>tailshort[circle-1])&&(random<=tailshort[circle]))
			    {   
					turn=circle-step*temp;
                    return turn;
			    }
		    }
		    return 0;
		}
	    else
       	{
			for(circle=(step*temp+11);circle<step*(temp+1);circle++)
		    {
				if((random>tailshort[circle-1])&&(random<=tailshort[circle]))
			    {
					turn=circle-step*temp;
                    return turn;
                }
            }
		}
	}
}
//======================================================================================
//function:length_distribution--return an random length
//======================================================================================
int length_distribution(int *seed,float length_rand[],int lable_method,int length_max,int input_max)
{
    float random;
	int i,length;
	random=(float)(rand_r(seed))/(float)(RAND_MAX);

        for(i=0; ;i++)
        {
            if((random>length_rand[i])&&(random<=length_rand[i+1]))
            {    
                if(input_max==length_max)
                {
                    length=i+1;
                    return length;
                }
                else
                {
                    length=(int)((i+(random-length_rand[i])/(length_rand[i+1]-length_rand[i]))*(input_max-1)/(length_max-1)+1);       //10.13change
                    return length;
                }
            }
        }
}
//======================================================================================
//function:start_position--return an random start position
//======================================================================================
int start_position(int max_t,int flag_pn,int len,int every_length,float value_f[])
{
	int i,pos,flag,start;
	float random;

	flag=0;
	pos=-1;
	while(flag==0)
	{
		random=(float)(rand())/(float)(RAND_MAX);
		//random=0.2;
		start=(int)(max_t*random);
		if(random>value_f[start])
		{
			for(i=start;i<max_t;i++)
			{
				if(random<value_f[i])
				{
					pos=i;
					break;
				}
			}
		}
		else
		{
			for(i=start;i>=0;i--)
			{
				if(random>value_f[i])
				{
					pos=i+1;
					break;
				}
			}
			if(pos==-1)
			{
				pos=0;
				break;
			}
		}
	
		pos=pos*100+rand()%100;
		if(flag_pn==1)
		{
			if((len-every_length)>=pos)
				flag=1;
			else
				flag=0;
		}
		else
		{
			if((pos>=every_length)&&(pos<len))
				flag=1;
			else
				flag=0;
        }
    }
}
//===========================================================================
//
//===========================================================================
int error_length(int *seed,float *error_length)
{
    double random;
    int len;
    int i;
    double temp;
    random=(double)(rand_r(seed))/(double)(RAND_MAX);
	if(random<=error_length[0])
    {
        return 1;
    }
    else if(random<=error_length[2500])
    {
        for(i=1;i<=2500;i++)
        {
            if((random>error_length[i-1]) && (random<=error_length[i]))
            {
                len=i+1;
                return len;
            }
        }
    }
    else if(random<=error_length[4999])
    {
        for(i=2501;i<=4999;i++)
        {
            if((random>error_length[i-1]) && (random<=error_length[i]))
			{
				len=i+1;
				return len;
			}
        }
    }
}
//=======================================================================
//
//=======================================================================
int pass_number(float *passnumber)
{
	int i;
	float random;
    random=(float)(rand())/(float)(RAND_MAX);

	if(random<=passnumber[50])
	{	
		for(i=1;i<=50;i++)
		{
			if((random>passnumber[i-1]) && (random<=passnumber[i]))
			{
				return i;
			}
		}
	}
	else
	{
		for(i=51;i<100;i++)
		{
			if((random>passnumber[i-1]) && (random<=passnumber[i]))
			{
				return i;
			}
		}
	}
}
//=======================================================================
//
//=======================================================================
float length_ratio(int *seed,float *ratio)
{
    float random;
    float temp;
    int i;
    float ran;
    random=(float)(rand_r(seed))/(float)(RAND_MAX);
    // random=0.2;
	if(random<=ratio[0])
    {
        ran=(float)(rand_r(seed))/(float)(RAND_MAX);
        temp=(float)(0.001*ran);
        return temp;
    }
    else
    {
        if(random<=ratio[499])
        {
            for(i=1;i<500;i++)
            {
                if((random > ratio[i-1]) && (random <= ratio[i]))
                {
                    ran=(float)(rand_r(seed))/(float)(RAND_MAX);
                    temp=(float)(0.001*(i+ran));
                    return temp;
                }
            }
        }
        else
        {
            for(i=500;i<1000;i++)
            {
                if((random>ratio[i-1]) && (random<=ratio[i]))
                {
                    ran=(float)(rand_r(seed))/(float)(RAND_MAX);
                    temp=(float)(0.001*(i+ran));
                    return temp;
                }
            }
        }
    }
}
//============================================================================
//
//=============================================================================
char convert(int *seed,char key[],float ins_less[])
{
	int i;
	int number;
	int temp=0;
	float random;
	for(i=0;i<3;i++)
	{
    	if(key[i]=='A')
    	{
       	 	number=0;
    	}
    	else if(key[i] == 'T')
    	{
      		number=1;
    	}
    	else if(key[i] == 'G')
    	{
        	number=2;
    	}
    	else if(key[i]=='C')
    	{
        	number=3;
    	}
        temp += number*pow(4,i);
	}
	random=(float)(rand_r(seed))/(float)(RAND_MAX);
	if(random <= ins_less[temp*4])
	{
		return 'A';
	}
	else if(random <= ins_less[temp*4+1])
	{
		return 'T';
	}
	else if(random <= ins_less[temp*4+2])
	{
		return 'G';
	}
	else
	{
		return 'C';
	}
}
//=====================================================================
//
//=====================================================================
int ffun_per(float *rate,char key[],int *seed,float headper[],float midper[],float tailper[],int pos,int sum_len)
{
	int i,number;
	int temp=0;
	float random;
	int event;
	int start;
	float unalign=0.6;
	if((pos<(rate[1]*sum_len)) || (pos> ((1-rate[2])*sum_len)))
	{
		 random=(float)(rand_r(seed))/(float)(RAND_MAX);
		 if(random<= unalign*0.8935/0.8734*rate[0])
		 {
		 	event=0;
		 }
		 else
		 {
		 	random=(float)(rand_r(seed))/(float)(RAND_MAX);
		 	if(random<=0.6)
		 	{	
		 		event=1;
		 	}
		 	else if(random<=0.9)
		 	{
		 		event=2;
		 	}
		 	else
		 	{
		 		event=3;
		 	}
		}
	}
	else
	{
		for(i=0;i<3;i++)
		{
			if(key[i]=='A')
			{	
				number=0;
			}
			else if(key[i]=='T')
			{
				number=1;
			}
			else if(key[i]=='G')
			{
				number=2;
			}
			else if(key[i]=='C')
			{
				number=3;
			}
			temp += number*pow(4,i);
		}
		start=0;
		event=decide_event(rate,start,temp,seed,headper);
	}
	return event;
}
//========================================================

//========================================================
int decide_event(float *rate,int start,int temp,int *seed,float array[])
{
	float random;
	float hhh;
	random=(float)(rand_r(seed))/(float)(RAND_MAX);
	hhh=(float)(array[start+temp]*rate[0]);
	if(random<=hhh)
	{
		return 0;              //match
	}
	else
	{
		random=(float)(rand_r(seed))/(float)(RAND_MAX);
		if(random<=array[start+64+temp*3])
		{
			return 1;
		}
		else if(random <= array[start+64+temp*3+1])
		{
			return 2;
		}
		else
		{
			return 3;
		}
	}
}
//=============================================================
//
//=============================================================
float get_eachread(int *seed,float eachread[])
{
	int i;
	float value;
	float random;
	random=(float)(rand())/(float)(RAND_MAX);
	if(random<=eachread[0])
	{
		return 0.005;
	}
	else if(random <= eachread[25])
	{
		for(i=1;i<=25;i++)
		{
			if((random>eachread[i-1]) && (random<=eachread[i]))
			{
				value=(float)(0.01*i+0.005);
				return value;
			}
		}
	}
	else
	{
		for(i=26;i<50;i++)
		{
			if((random>eachread[i-1]) && (random<=eachread[i]))
			{
				value=(float)(0.01*i+0.005);
				return value;
			}
		}
	}
}
float get_unalign(int *seed,float *unalign)
{
	int i;
	float random,value;
	random=(float)(rand())/(float)(RAND_MAX);
	if(random<=unalign[0])
	{
		return 0.0;
	}
	else if(random <= unalign[500])
	{
		for(i=1;i<=500;i++)
		{
			if((random>unalign[i-1]) && (random<=unalign[i]))
			{
				value=0.001*i;
				return value;
			}
		}
	}
	else
	{
		for(i=501;i<1000;i++)
		{
			if((random>unalign[i-1]) && (random<=unalign[i]))
			{
				value=0.001*i;
				return value;
			}
		}
	}
}
