/* program symulatora  sieci  neuronowej         */
/* na wejsciu wczytuje tablice neuronow i synaps */
/* przygotowane wczesniej                        */
/* Piotr Franaszczuk & Pawel Kudela              */
/* JHU Baltimore, MD                             */
/******** wersja 2.5 Feb 2002 *********************/

#define RISC
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <float.h>
#include <signal.h>
#include "lnet.h"
#include "clust_cn.h"

#define randoms(num) (int)floor((double)num*drand48())

FILE *config;              /* zbior tekstowy zawierajacy konfiguracje */

/*double nn=0,nv=0,nw=0,nx=0,nc=0,nb=0;*/
struct SYNAPSE * synapses;
struct NEURON  * neurons;
struct PARAMETERS   * params;
struct SYNAPSESD * synapsesd;
struct SYNAPSESD * ssz, * stsz;
//struct NEUROND * neuonsd;
//double * synapsesd;
//long int jd;
int *szum_file,szum_i,stim_i;

SPIKE_INTERVAL *szum_interval;
SPIKE_INTERVAL zero_interval = 0;

int            neuron_nr;
int            no_neurons;  /* liczba neuronow w sieci  */
int            no_kinds;    /* liczba rodzajow neuronow */
int            no_synapses; /* liczba synaps w sieci    */
VAR	       timet;
int            lambda_int;  /* integer lambda for noise */
int            out_i;
int            his_nr;      /* his_nr - liczba neuronow dla ktorych obliczany jest histogram */
int            out_nr;      /* out_nr - liczba neuronow dla ktorych zapisywany jest output   */

char  *p=NULL,*path;                /* working path */

int histo_f,out_f,vout_f,noise_f,iext_f,link_f,stim_f,stim_flag;  /* flagi sterujace*/
int firststimflag=0;
int stimpulsecounter=0;
int licz=0;
double lambda;
double lambda1=0.00020068;

/* struct SYNAPSE ** synaps_j; */

long stim_pause,stim_len,stim_counter=1;
long stimulus_period, stim_period=0;
VAR iext;               /* To provide first time through information to stimulus subroutine */
// int stimflag,pauseflag; /* stimulation on/off flags for subroutine stimulate */
FILE* stim_file, *Vfile=NULL;

struct LINK **link_in;
struct LINKD **linkw;
//struct LINKD **linkstimw;
//double **linkw;
NEURON_INDEX ** link_out;
int link_in_nproc,link_out_nproc,*link_in_no,*link_out_no;

extern char **proc_name;
extern int *conn_sock;                  /* array of  sockets descriptors */ 
extern struct CONN cn_inp;                    /* input conn.sockets  fo host */
extern struct CONN cn_out;                    /* out conn.sockets  for host */
extern int nproc;

extern int nfd;                          /* largest file descriptor to check in select */
extern fd_set fdr0;                      /* descriptor set for select */ 
extern int *buf_in;     /* Max buffer there are always two leading shorts */  //was short
extern int this_proc;                     /* this proc no */
extern char *this_name;
unsigned short  PORT;
int his_licz = 0, his_out=0; /* his_out - liczba krokow po ktorych histogram jest zapisywany na dysk */
int Vave_licz=0, Vave_out=0; /* same variables for Vave sampling */
int count = 0;
float   dt = DELTA_T; 
int     COUNT_NO;
char    par_name[80], syn_name[80], neur_name[80], cfg_name[80],*argv1;
FILE* his_file=NULL;
FILE* Ca_file=NULL;
FILE* Noise_file=NULL;
FILE* Vave_file=NULL;
FILE* VPSP_file=NULL;
FILE **out_file=NULL;
char    buf1[400],buf2[80];
double his_fout;
double testpulsecount;
int Noise_sum=0;

int sigflag=-2,sz=0,stz=0;

int apcount=0; /* Keeping track of stimulation efficacy */

#ifdef DEB
FILE *lWfile, *Wfile,*Cfile,*Sfile;
#endif
//VAR I_NMDA, G_NMDA, R_NMDA, t1_NMDA, t2_NMDA, Co;
//VAR V_NMDA;

void sigproc(int sig)
{
  sigflag=sig;
 
  sprintf(buf2,"received signal %i;his_out-his_licz=%i",sig,his_out-his_licz);
  if(this_proc==0)printf("%s\n",buf2);
  print_err(INFO,"sigproc",buf2,NULL,NULL);
  return;
}
  
int main(int argc, char *argv[] )
{

  VAR     Amp;
   
  int new_sim;
  char * pname="main";
  char erbuf[200]; 

  struct NEURON  * neuron_i, * neuron_j;
  struct SYNAPSE * synapse_j;
  struct PARAMETERS *par_i;
  SYNAPSE_NUMBER syn_j;
  SPIKE_INTERVAL time_i, delay_j;

  SPIKE_INTERVAL *delay_i;

  int index, end_index;

  NEURON_INDEX his_i=0;
  HIS *histo=NULL;
  double nsec;  /* dlugosc symulacji */
  VAR decr=1,indecr; //,decrR=1;
   
  char *fopen_flag;

  double ca_sum,ca_ave;
  double V_sum, V_ave;
  double VPSP_sum, VPSP;
  double testpulsecount;
  double I_SYNE, I_SYNI;
  int kk;

  COUNT_NO = 0.01/dt; 

  if(argc < 4)
    {
      printf("\nUsage: net <name> <nsec> <decr> <port> [flag]\n");
      printf("\nwhere: <name> cfg file name  (no ext),");
      printf("\n       <nsec> time of simulation in sec ,");
      printf("\n       <port> port# for communiccations  ,");
      printf("\n       [flag] if present continue previous simulation.\n");
      return 1;
    }
  path=strdup(argv[0]);
  argv1=strdup(argv[1]);

  p = strrchr(path,'/');
  if(p==NULL)strcpy(path,"./");
  else *(p+1)='\0';
  freopen("netclust.log","w",stderr);
  if(signal(SIGUSR1,sigproc)==SIG_ERR )
    print_err(PERROR,"main","signal USR1",NULL,NULL);
  if(signal(SIGPWR,sigproc)==SIG_ERR )
    print_err(PERROR,"main","signal PWR",NULL,NULL);

  nsec = atof(argv[2]);
  indecr=atof(argv[3]);
  PORT=atoi(argv[4]);
  if(argc>5 && atoi(argv[5]))
    new_sim=0;
  else
    new_sim=1;

  unlink("done");

  fopen_flag=new_sim?"w":"a"; /* append or create new */   
    
  if(nsec<=0)nsec=1;

  init_path(argv[0]);   /* initialization  */
  
  read_conn(argv1);


  if(new_sim)
    file_name(cfg_name,argv1,".cfg");
  else
    file_name(cfg_name,argv1,".cfg.bak");

  config = fopen(cfg_name,"r");
  if(config==NULL)
    print_err(FERROR,"main","cfg fopen",cfg_name,config);

  set_flags();

     if(link_f)
	proc_conn();
  
  fscanf(config,"%s",par_name);
  fscanf(config,"%i",&no_kinds);
  read_params(par_name);

  fscanf(config,"%s",syn_name);
  fscanf(config,"%i",&no_synapses);
  fscanf(config,"%s",neur_name);
  fscanf(config,"%i",&no_neurons);
  fscanf(config,"%i,%lf",&his_nr,&his_fout);
  fscanf(config,"%i",&out_nr);
  fscanf(config,"%lf",&lambda);
  fclose(config);

  if(no_synapses < 1 )
    print_err(ERROR,"main","no_synapses < 1",NULL,NULL);

  read_synapses(syn_name);
  
  if(no_neurons < 2 )
    print_err(ERROR,"main","no_neurons < 2",NULL,NULL); 
    
  read_neurons(neur_name,new_sim);  

  if(stim_f)
    {
       if(new_sim)
	 file_name(buf1,argv1,".stim");
       else
	 file_name(buf1,argv1,".stim.bak");

      stim_file=fopen(buf1,"r"); /* same as cfg name */
      if(stim_file==NULL)
	print_err(FERROR,"main","stim file fopen",buf1,stim_file);
      if(!new_sim)
	fscanf(stim_file,"%li,%li\n",&stim_counter,&stim_len);
      else
	stim_counter=1;

      /* stimflag=0;
      pauseflag=0;
      iext=0; */

       if(stz)
        stszalloc(stz);
    }
     
  if(histo_f)
    {    
      sprintf(buf2,"his_nr=%i,his_out=%lf",his_nr,his_fout);
      print_err(INFO,"",buf2,NULL,NULL);
      his_out=rint(his_fout/(1000.*DELTA_T));
      //Vave_out=rint(DELTA_T*10/DELTA_T);
      histo = (HIS*)calloc(his_nr, sizeof(HIS));
      if(histo==NULL)
	print_err(ERROR,"main","calloc(histo)",NULL,NULL);
    }

  if(out_f)
    { 
      sprintf(buf2,"out_nr=%i",out_nr);
      print_err(INFO,"",buf2,NULL,NULL);
      out_file = (FILE**)malloc(out_nr*sizeof(FILE*));
      if(out_file==NULL)
	print_err(ERROR,"main","malloc(out_file)",NULL,NULL);
    }

  if(link_f)
    read_links(argv1,new_sim);  /* use same name as cfg_name */


  if(vout_f)
    { 
      strcpy(buf2,"V_");strcat(buf2,argv1);
      file_name(buf1,buf2,".out");

      Vfile = fopen(buf1,fopen_flag);

      if(Vfile==NULL)
	print_err(FERROR,"main","Vfile fopen",buf1,Vfile);
    }

  
#ifdef DEB
  //Wfile = fopen(file_name(buf1,"Wfile",""),"w");
  //lWfile= fopen(file_name(buf1,"lWfile",""),"w");
  //Cfile = fopen(file_name(buf1,"Cfile",""),"wb");
  //Sfile = fopen(file_name(buf1,"Sfile",""),"w");
#endif
  if(out_f)
    {

      for( out_i = 0; out_i < out_nr;out_i++)
	{
	  sprintf(buf2,"NEURON_%s.%d.%d",argv1,out_i,this_proc);
	  file_name(buf1,buf2,NULL);

	  out_file[out_i] = fopen(buf1,fopen_flag);
 
	  if(out_file[out_i]==NULL)
	    print_err(FERROR,"main","out_file fopen",buf1,out_file[out_i]);
	}
    }

  if(histo_f)
    { 

      strcpy(buf2,"histo_");strcat(buf2,argv1);
      file_name(buf1,buf2,"");
  
      his_file = fopen(buf1,fopen_flag);
      if(his_file==NULL)
	print_err(FERROR,"main","his_file fopen",buf1,his_file);

      /* [Ca] output file setup */
      strcpy(buf2,"Ca_");strcat(buf2,argv1);
      file_name(buf1,buf2,"");
  
      Ca_file = fopen(buf1,fopen_flag);
      if(Ca_file==NULL)
	print_err(FERROR,"main","Ca_file fopen",buf1,Ca_file);

      /* Noise output file setup */
      strcpy(buf2,"Noise_");strcat(buf2,argv1);
      file_name(buf1,buf2,"");
  
      Noise_file = fopen(buf1,fopen_flag);
      if(Noise_file==NULL)
	print_err(FERROR,"main","Noise_file fopen",buf1,Noise_file);

      /* Vave output file setup */
      strcpy(buf2,"Vave_");strcat(buf2,argv1);
      file_name(buf1,buf2,"");

      Vave_file = fopen(buf1,fopen_flag);
      if(Vave_file==NULL)
	print_err(FERROR,"main","Vave_file fopen",buf1,Vave_file);

      /* VPSP output file setup */
      strcpy(buf2,"VPSP_");strcat(buf2,argv1);
      file_name(buf1,buf2,"");

      VPSP_file = fopen(buf1,fopen_flag);
      if(VPSP_file==NULL)
	print_err(FERROR,"main","VPSP_file fopen",buf1,VPSP_file);

    }

  if(noise_f)
    {
      sprintf(buf1,"noise lambda %lf",lambda); 
      print_err(INFO,"",buf1,NULL,NULL);
      srand48(127);
      if(sz)
       sszalloc(sz);
    }


  if(new_sim)init();       /* initial values of NEURON variables */     

  sigflag=0;

  /* glowna petla symulacji */
if(this_proc==0){
printf("\n time    Aep[0]   Aip[0]");
}
 
  for(;;)
    {
      /* petla calkowania i generacji spike'ow */
      if(histo_f)
	{
	  if(++his_licz >= his_out )
	    {
	      if(fwrite(histo,sizeof(HIS),his_nr,his_file)!=his_nr)
		print_err(FERROR,"main","fwrite(histo)",buf1,his_file);
	      if(sigflag)break;
		
	      memset(histo,0,his_nr*sizeof(HIS));
	      his_licz=0;

	    }

	  if(++Vave_licz >= Vave_out )
	  {

               /* Cellular calcium averaging and printing to file */
	       ca_sum=0;
               V_sum=0;
	       VPSP_sum=0;
	       for( neuron_i = neurons,neuron_nr=0; neuron_nr < no_neurons; neuron_i++,neuron_nr++ ){
	            ca_sum+=neuron_i->C[6];
                    
		    par_i = params + neuron_i->param;

		    I_SYNE = 0;
		    for (kk=0;kk<NEXCCLASS;kk++){

			I_SYNE += (neuron_i->Ve_o[kk] + neuron_i->Ve_d[kk])*(par_i->Esyn_e - neuron_i->V);

		    }

		    I_SYNI=0;
		    for (kk=0;kk<NINHCLASS;kk++){

			I_SYNI += (neuron_i->Vi_o[kk] + neuron_i->Vi_d[kk])*(neuron_i->V - par_i->Esyn_i);

		    }

		    VPSP_sum+=I_SYNE+I_SYNI;

		    V_sum+=neuron_i->V;
		    
               }

               ca_ave=ca_sum/no_neurons;
	       V_ave=V_sum/no_neurons;
	       VPSP=VPSP_sum/no_neurons;

	       /* sprintf(erbuf,"Ave outer shell [Ca] = %lf\n",ca_ave);		
		  print_err(INFO,pname,erbuf,NULL,NULL); */
               
	       if(fwrite(&ca_ave,sizeof(double),1,Ca_file)!=1)
		 print_err(FERROR,"main","fwrite(ca_ave)",buf1,Ca_file);

	       if(fwrite(&Noise_sum,sizeof(int),1,Noise_file)!=1)
		 print_err(FERROR,"main","fwrite(Noise_sum)",buf1,Noise_file);

	       if(fwrite(&V_ave,sizeof(double),1,Vave_file)!=1)
		 print_err(FERROR,"main","fwrite(V_ave)",buf1,Vave_file);

	       //sprintf(erbuf,"Made it Vave write, Vave=%lf\n",V_ave);		
	       //      print_err(INFO,pname,erbuf,NULL,NULL);

	       if(fwrite(&VPSP,sizeof(double),1,VPSP_file)!=1)
		 print_err(FERROR,"main","fwrite(VPSP)",buf1,VPSP_file);

	       Vave_licz=0;

	       Noise_sum=0;

	  }


	  else 
	    if(sigflag<0)break;
	  his_i=0;
	}
      else
	if(sigflag) break;

      if(licz-- == 0)
	{
	 timet = count*COUNT_NO;
	  if(this_proc==0){
	    //if(timet!=0)
               printf("\n%5.2f sec %lf %lf",timet*dt,(params+0)->Aep[0],(params+0)->Aip[0]);
            //else
	    //  printf("\nstart");
	  }
	  if((timet*dt >= 1 )&&((params+0)->Aip[0] > 0.000011)) //0.5
	     {
              decr=indecr;
              //decrR=0.9999997;
              //lambda=lambda1;              
             }
          else
           {
                decr=1.0;
           }


	  /* Patch for limiting noise time */
	  //if((timet*dt) >= 0.4)
	  //  noise_f=0; 

	  fflush(stdout);
	  count++;
            
	  licz=COUNT_NO-1;
	  if(count*COUNT_NO*dt > nsec)  
	    break;
	}

      //  (params+0)->Aip*=decr; 
    /* if(timet*dt>=1.0){
       if( (params+0)->Aep < 0.00036 ){
       (params+0)->Aep*=1.0000004999;
       (params+1)->Aep*=1.0000004999;
       }
       (params+0)->R*=0.999999;}*/ /* was 0.9999997 */
      //  R_NMDA *=decrR;
      //if( (lambda <0.04)&&(timet*dt >= 10 ))
      // lambda*=1.000001999;

      firststimflag=0;
      if(stim_f)stimulus();
      if(iext_f) --stim_period;
      out_i=0;
      szum_i=0;
      stim_i=0;
      
      
      for( neuron_i = neurons,neuron_nr=0; neuron_nr < no_neurons; neuron_i++,neuron_nr++ )
	{
	  neuron_i->interval++;
	  if( calkuj( neuron_i ) )
	    {
	      /* nowy spike */
	      if( histo_f && (neuron_i->flags & HISTO)
		  && (++histo[his_i] == 0 ) )
		histo[his_i] = -1;
	      if(neuron_i->syn_num)
                {
		  switch( neuron_i->first )
		    {
		    case EMPTY_LIST:
		      neuron_i->first = ONE_SPIKE;
		      break;
		    case ONE_SPIKE:
		      neuron_i->first = neuron_i->last = 0;
		      neuron_i->spikes[0] = neuron_i->sum_interval = neuron_i->interval;
		      break;
		    default:
		      neuron_i->last = (neuron_i->last+1)%SPIKE_LIST_LENGTH;
		      if( neuron_i->last == neuron_i->first ) 
			  print_err(ERROR,"main","spike list full",NULL,NULL);
		      neuron_i->spikes[ neuron_i->last ] = neuron_i->interval;
		      neuron_i->sum_interval += neuron_i->interval;
		    }
                } 
	      else
                {      
		  switch( neuron_i->first )
		    {  
		    case EMPTY_LIST:
		      break; 
		    case ONE_SPIKE:
		      neuron_i->sum_interval = neuron_i->interval;
		      break;
		    default: 
		      neuron_i->sum_interval += neuron_i->interval;
		    }

                }  
	      if( out_f && neuron_i->flags & OUT )

		if(fwrite(&(neuron_i->interval), sizeof(SPIKE_INTERVAL),1,out_file[out_i])!=1)
		  {
		    sprintf(buf2,"NEURON_%s.%d.%d",argv1,out_i,this_proc);
		    print_err(FERROR,"main","fwrite(neuron interval)",buf2,out_file[out_i]);
		  }
	      neuron_i->interval = 0;

	    }
	  else
	    /* nie ma spike'u */
	    /* sprawdzamy czy dlugo nie ma */
	    if( neuron_i->interval == MAX_INTERVAL )
	      {
		neuron_i->interval = 0;

		if( out_f && neuron_i->flags & OUT )
		  if(fwrite(&zero_interval, sizeof(SPIKE_INTERVAL),1,out_file[out_i])!=1)
		    {
		      sprintf(buf2,"NEURON_%s.%d.%d",argv1,out_i,this_proc);
		      print_err(FERROR,"main","fwrite(zero_interval)",buf2,out_file[out_i]);
		    }
	      }
	  if( out_f && neuron_i->flags & OUT ) out_i++;
	  if( histo_f && (neuron_i->flags & HISTO)) his_i++;
	  
	
	  /* tutaj dodawanie szumu */

           
	  if( (noise_f) && (neuron_i->flags & NOISE) )
	    {
                 //(ssz+szum_i)->w += (om((ssz+szum_i)->C) - (ssz+szum_i)->w)/taul((ssz+szum_i)->C);
               //I_NMDA = G_NMDA*((ssz+szum_i)->V_d+(ssz+szum_i)->V_o)*B_NMDA(neuron_i->V)*(V_NMDA((ssz+szum_i)->C)-neuron_i->V);
               //(ssz+szum_i)->C += I_NMDA-R_NMDA*((ssz+szum_i)->C-Co);
               //(ssz+szum_i)->V_d *=t1_NMDA;
               //(ssz+szum_i)->V_o *=t2_NMDA;
               //if(!(licz%10) && (neuron_nr == 20))   
               //fprintf(Sfile,"%d,%lf,%d,%lf,%lf,%lf\n",(int)neuron_nr,(ssz+szum_i)->w,(int)neuron_nr,(ssz+szum_i)->C,I_NMDA,6*neuron_i->C[0]);
               /* sprintf(erbuf,"Time stamp in noise section licz=%i\n",licz);
		  print_err(INFO,pname,erbuf,NULL,NULL); */
               szum( neuron_i);
               szum_i++;
	    }


          if( neuron_i->flags & IEXT  )
            {
                 //(stsz+stim_i)->w += (om((stsz+stim_i)->C) -(stsz+stim_i)->w)/taul((stsz+stim_i)->C);
              //I_NMDA =G_NMDA*((stsz+stim_i)->V_d+(stsz+stim_i)->V_o)*B_NMDA(neuron_i->V)*(V_NMDA((stsz+stim_i)->C)-neuron_i->V);
              //(stsz+stim_i)->C += I_NMDA-R_NMDA*((stsz+stim_i)->C-Co);
              //(stsz+stim_i)->V_d *=t1_NMDA;
              //(stsz+stim_i)->V_o *=t2_NMDA;
              //if(!(licz%10) && (neuron_nr == 20))   
              //fprintf(Sfile,"%d,%lf,%d,%lf,%lf,%lf\n",(int)neuron_nr,(stsz+stim_i)->w,(int)neuron_nr,(stsz+stim_i)->C,I_NMDA,6*neuron_i->C[0]);
		if( (iext_f==1) && (firststimflag==1) ){
            
	           testpulsecount=stimpulsecounter;

	           if ( (remainder(testpulsecount,2) || 0) && (neuron_i->polarity & POSPOL) ){

                        stim_freq(neuron_i);              
		 
                   }    

		   if ( (remainder(testpulsecount,2) == 0) && (neuron_i->polarity & NEGPOL) ){

                        stim_freq(neuron_i);   

                        /* sprintf(erbuf,"Made it to stim, stimpulsecounter=%i\n",stimpulsecounter);		
			   print_err(INFO,pname,erbuf,NULL,NULL); */ 
		 
                   }    

	      }

             }

          

	} /* nastepny neuron */
      
      if(iext_f && stim_period==0)  
         stim_period=stimulus_period;

      if(link_f)link_update();

	
      /* petla update'u synaps */

      synapse_j = synapses;
      //jd = 0;
      for( neuron_i = neurons,neuron_nr=0; neuron_nr < no_neurons; neuron_i++,neuron_nr++ )
	if( neuron_i->first != EMPTY_LIST )

	  {	   
	    delay_i = (params+neuron_i->param)->del;

	    /* petla po synapsach wyjsciowych i-tego neuronu */
	    /* synapsy w danym neuronie musza byc uporzadkowane wg. */
	    /* rosnacych opoznien */
	    for( syn_j = 0; syn_j < neuron_i->syn_num; syn_j++)
	      {
		//long int jd=synapse_j - synapses;
		neuron_j = neurons + synapse_j->neuron;

                if(synapse_j->weight > 0)  /* synaptic facilitation */
                {
                //(synapsesd+jd)->w += (om((synapsesd+jd)->C) - (synapsesd+jd)->w)/taul((synapsesd+jd)->C);
                //I_NMDA = G_NMDA*((synapsesd+jd)->V_d+(synapsesd+jd)->V_o)*B_NMDA(neuron_j->V)*(V_NMDA((synapsesd+jd)->C)-neuron_j->V);
                //(synapsesd+jd)->C += I_NMDA-R_NMDA*((synapsesd+jd)->C-Co);
                //(synapsesd+jd)->V_d *=t1_NMDA;
                //(synapsesd+jd)->V_o *=t2_NMDA;
#ifdef DEB
                //if (!(licz%10) && (synapse_j->neuron == 20))
                //fprintf(Wfile,"%d,%lf,%d,%lf,%lf,%lf\n",(int)neuron_nr,(synapsesd+jd)->w,(int)synapse_j->neuron,(synapsesd+jd)->C,I_NMDA,6*neuron_j->C[0]);
#endif
		} 
		delay_j = *(delay_i+ neuron_j->param)
		  + synapse_j->delay;
		time_i = neuron_i->interval + neuron_i->sum_interval;
		if( (index = neuron_i->first) != ONE_SPIKE )
		  {
		    /* sprawdzam liste spike'ow */
		    end_index = (neuron_i->last+1)%SPIKE_LIST_LENGTH;
		    while( time_i > delay_j && index != end_index )
		      {
			time_i -= neuron_i->spikes[ index ];
			index = (index+1)%SPIKE_LIST_LENGTH;
		      }
		  }
		if( time_i == delay_j )
		  {

		 
		    /* obliczanie amplitudy do dodania w PSP */

		    if( synapse_j->weight > 0 )
		      {
			//Amp = 2.5*(synapsesd+jd)->w * synapse_j->weight * (params+neuron_j->param)->Aep;
                        Amp =  synapse_j->weight * (params+neuron_j->param)->Aep[synapse_j->preclass_num];
			//Amp = Amp*0.625;

                        neuron_j->Ve_d[synapse_j->preclass_num] += Amp;
			neuron_j->Ve_o[synapse_j->preclass_num] -= Amp;
                        //(synapsesd+jd)->V_d +=Amp;
                        //(synapsesd+jd)->V_o -=Amp;

    
		      }
		    else
		      {
			//printf("weight= %i\n",synapse_j->weight);
			Amp =  synapse_j->weight *(params+neuron_j->param)->Aip[synapse_j->preclass_num-NEXCCLASS];

			neuron_j->Vi_d[synapse_j->preclass_num-NEXCCLASS] += Amp;
			neuron_j->Vi_o[synapse_j->preclass_num-NEXCCLASS] -= Amp;

		      }


	            if( syn_j == neuron_i->syn_num-1  )
		      {

			/* usuwam z listy najstarszy spike */
			if( neuron_i->first == ONE_SPIKE )
			  neuron_i->first = EMPTY_LIST;
			else
			  {
			    index = neuron_i->first;
			    neuron_i->sum_interval -=  neuron_i->spikes[ index ];
			    if( index == neuron_i->last )
			      neuron_i->first = ONE_SPIKE;
			    else
			      neuron_i->first = (index+1)%SPIKE_LIST_LENGTH;
			  }
		      }
		  } /* koniec if( time_i == delay_j ) */
		synapse_j++;
                //jd++;
	      }/* nastepna synapsa */
	  }
	else       /* EMPTY_LIST*/
          {
            for( syn_j = 0; syn_j < neuron_i->syn_num; syn_j++)
	          {
		   /* long int jd=synapse_j - synapses; */
                   neuron_j = neurons + synapse_j->neuron;

                   if(synapse_j->weight > 0)  /* synaptic facilitation */
                   { 
                   //(synapsesd+jd)->w += (om((synapsesd+jd)->C) - (synapsesd+jd)->w)/taul((synapsesd+jd)->C);
                   //I_NMDA = G_NMDA*((synapsesd+jd)->V_d+(synapsesd+jd)->V_o)*B_NMDA(neuron_j->V)*(V_NMDA((synapsesd+jd)->C)-neuron_j->V);
                   //(synapsesd+jd)->C += I_NMDA-R_NMDA*((synapsesd+jd)->C-Co);
                   //(synapsesd+jd)->V_d *=t1_NMDA;
                   //(synapsesd+jd)->V_o *=t2_NMDA;
#ifdef DEB
                   //if (!(licz%10) && (synapse_j->neuron == 20) )
                   //fprintf(Wfile,"%d,%lf,%d,%lf,%lf,%lf\n",(int)neuron_nr,(synapsesd+jd)->w,(int)synapse_j->neuron,(synapsesd+jd)->C,I_NMDA,6*neuron_j->C[0]);
#endif
                   }
                synapse_j++;
                //jd++;
                  }
          } 
      //	  synapse_j +=  /* liczba synaps i-tego neuronu */   
      //            neuron_i->syn_num;

      /* nastepny neuron */
    }/* koniec glownej petli */
  fflush(NULL);
  close_all(sigflag);

  sprintf(erbuf,"apcount = %i\n",apcount);		
  print_err(INFO,pname,erbuf,NULL,NULL); 


#ifdef DEB
  //fflush(Wfile);
  //fclose(Wfile);
  //fflush(lWfile);
  //fclose(lWfile);
  //fflush(Cfile);
  //fclose(Cfile);
  //fclose(Sfile);
#endif

  return 0;
}
void close_all(int sig)
{ int out_i;
  
 
 if(histo_f && his_file){
     fclose(his_file);
     fclose(Ca_file);
 }
 if(out_f)
   {
     for( out_i = 0; out_i < out_nr;out_i++)
       if(out_file[out_i])fclose(out_file[out_i]);
   }
 if(vout_f && Vfile) {fclose(Vfile);}
 if(sigflag>-2)
   {
     write_neurons(file_name(buf1,neur_name,".neu.bak"));
     if(link_f)write_links(file_name(buf1,neur_name,".lin.bak"));
     save_stimulus();
   }
 if(this_proc==0)printf("\nTime: %5.2f sec\n",(count*COUNT_NO-licz)*dt);
 if(histo_f)
   {
     sprintf(buf1,"his_licz=%i,his_out=%i",his_licz,his_out);
     print_err(INFO,"close_all",buf1,NULL,NULL);
   }
 sprintf(buf1,"exit on signal %d",sig);
 print_err(INFO,"",buf1,NULL,NULL);
 close_socks();
 fopen("done","w");
 return;
   
}

void  put_flags(char *buf)
{
  *buf=0;
  if(stim_f)iext_f=1;
  if(noise_f)strcat(buf,"N,");
  if(histo_f)strcat(buf,"H,");
  if(vout_f)strcat(buf,"V,");
  if(out_f)strcat(buf,"O,");
  if(iext_f)strcat(buf,"I,");
  if(stim_f)strcat(buf,"S,");
  strcat(buf,"$");
}

int save_stimulus()
{
  double pause,len;
  char buf[20];
     
  VAR iext=params->i_ext;

  if(stim_f)
    {
      FILE* stim_bak;
      //double pause,len;
      //VAR iext;

      stim_bak=fopen(file_name(buf1,argv1,".stim.bak"),"w");
      if(stim_bak==NULL)
	print_err(FERROR,"save_stimulus","stim.bak  fopen",buf1,stim_bak);
      fprintf(stim_bak,"%li,%li\n",stim_counter,stim_len);
      while(fscanf(stim_file,"%lf %lf %lf\n",&pause,&len,&iext)==3)
	fprintf(stim_bak,"%lf %lf %lf\n",pause,len,iext);
      fclose(stim_bak);  
    }


  config = fopen(file_name(cfg_name,argv1,".cfg.bak"),"w");
  if(config==NULL)
    print_err(FERROR,"save_stimulus","cfg.bak fopen",cfg_name,config);     
  put_flags(buf);
  fprintf(config,"%s\n",buf);
  fprintf(config,"%s",par_name);
  fprintf(config,"%i",no_kinds);
  fprintf(config,"%s",syn_name);
  fprintf(config,"%i",no_synapses);
  fprintf(config,"%s",neur_name);
  fprintf(config,"%i",no_neurons);
   
  fprintf(config,"%i,%lf",his_nr,his_fout); 
  fprintf(config,"%i",out_nr);
  fprintf(config,"%lf",lambda);
  fclose(config);
  return 0;
}

int stimulus()   
{
  double pause,len;
  VAR iext;
  int i;
  // int stimflag,pauseflag;
  
  if(--stim_counter) return 0; // on start stim_count=1 and iext_f  true

  /* if (stimflag){
       stimflag=0;
       pauseflag=1;}
  else {
       stimflag=1;
       pauseflag=0;} */

  if(iext_f)
    {
	// if (iext>=0){

             if(fscanf(stim_file,"%lf %lf %lf\n",&pause,&len,&iext)!=3)
	          {stim_f=0;iext_f=0;return 2;}
             else
	     {
	          if(pause > LONG_MAX*DELTA_T|| len > LONG_MAX*DELTA_T) 
	               print_err(ERROR,"stimulus"," pause or len to large\n",NULL,NULL);  
	          stim_pause=pause/DELTA_T;
	          stim_len=len/DELTA_T;
                  /* if (stim_pause){
                       stimflag=0;
                       pauseflag=1;}
                  else {
                       stimflag=1;
                       pauseflag=0;} */
             }

//	}

      if(iext>0)
      {     
       for(i=0;i<no_kinds;i++)
	 (params+i)->i_ext=iext;
      }
      else
      {
       stimulus_period=-1./iext/DELTA_T;
       stim_period=1;
       for(i=0;i<no_kinds;i++)
         (params+i)->i_ext=0;    
      }
	   
      if(stim_pause)
	{
	  iext_f=0;
	  stim_counter=stim_pause;
          if(this_proc==0)printf("\n- stim_pause=%ld",stim_pause);fflush(stdout);
	}
      else
	{ // this is for initial pause == 0
	  stim_counter=stim_len;
	  iext_f=1;
          firststimflag=1;
	  stimpulsecounter=stimpulsecounter+1;
          if(this_proc==0)printf("\n+ stim_len=%ld",stim_len);fflush(stdout);
	}
    }
  else
    {
      stim_counter=stim_len; // this should never be zero
      iext_f=1;
      firststimflag=1;
      stimpulsecounter=stimpulsecounter+1;
      if(this_proc==0)printf("\n+ stim_len=%ld",stim_len);fflush(stdout);
    }
  return 2;
}
	
   
void add_to_list(struct LINK *l)
{
  char erbuf[100];

    
  if(l->first != EMPTY_LIST)
    {
      l->last=(l->last+1)%LINK_LIST_LENGTH;
      if(l->last==l->first)        
	{     
	  sprintf(erbuf,"\nLink list overflow licz = %i,neuron :%hu ",licz,l->neuron_post);
	  print_err(ERROR,"add_to_list",erbuf,NULL,NULL);
	}
      l->interval[l->last]=l->delay;
    }
  else
    {
      l->first=l->last=0;
      l->interval[0]=l->delay;

      /* if( iext_f==1 && firststimflag==1 ){
           sprintf(erbuf,"\nIn add_to_list l->first=%i,l->interval[0]=%i",l->first,l->interval[0]);
           print_err(INFO,"add_to_list",erbuf,NULL,NULL);
	 } */

    }

    
}
    
    
void link_update()
{
  struct LINK *l_in;
  struct LINKD  *lw;
  //struct LINKD  *lstimw;
  struct timeval tim ={100,0};
  //double *lw;
  NEURON_INDEX * l_out;
  int  *j,index,end_index,i;
  int *out_buf; // was short
  char * pname="link_update";
  int m,n,k,buflen=0,ii;
  fd_set fdr1,fdr2;
  struct NEURON *neuron_j; 
  char erbuf[200];
  int jd;
  //ushort postneuron;
  //int Deelay,Weeight;
  //byte Feirst,Laast;

  int mysend(int,char *, int, int); //was char * for 2nd arg
  int myrecv(int,char *, int, int);   //was char * for 2nd arg

   
  /* update local links */
  if(link_in_nproc>cn_inp.nconn)
     
    { l_in=link_in[cn_inp.nconn]; /* last is for this */
    l_out=link_out[cn_out.nconn];
    for(n=0;n<link_in_no[cn_inp.nconn];n++,l_out++)
      if((neurons+*(l_out))->flags & SPIKE_DET)
	add_to_list(l_in+n);	   
    }

  /* fill  out_buffers and send */ 
  j=cn_out.lconn;
   
  for(m=0;m<cn_out.nconn;m++,j++)
    {

      k=1;
      l_out=link_out[m];
      out_buf=cn_out.buf[m];
      out_buf[0]=(int)licz;   /* time stamp common for all */
      for(n=0;n<link_out_no[m];n++,l_out++)
	{
	  if(k==MAX_SPIKE_OUT+1)
	       print_err(ERROR,pname,"Too many spikes for one packet",NULL,NULL); 
	  if((neurons+*(l_out))->flags&SPIKE_DET)
	      out_buf[++k]=n;           /* !!!! limit links to less then 32768 but uses only 2 bytes*/ //was (short)n
	}
      k--;       
      out_buf[1]=k; //was (short)k
	
      buflen=(k+2)*sizeof(int); //was sizeof(short)
	
        
      if((mysend(conn_sock[*j],(char*)out_buf,buflen,0))!=buflen) //was (char*)out_buf
	print_err(PERROR,pname,"send ",NULL,NULL);

     
    }
  
  /* read input from other proc */

  k=cn_inp.nconn;

  if(k>0){
    j=cn_inp.lconn;
    memcpy(&fdr1,&fdr0,sizeof(fd_set));
  }
  /* check if there are packets read in last iteration */
  for(i=0;i<cn_inp.nconn && k>0;i++,j++)
    {

      if((ii=cn_inp.buf[i][0])>0) /* ii index to begining of next packet */
      { int *b; //was short
	    
	b=cn_inp.buf[i]+ii;   
	    
	if(b[0]!=(int)licz) /* synchro error */
	  { char erbuf[100];
	     
	  sprintf(erbuf,"next packet sync licz=%i,b0=%i,cn[%i]=%i",licz,b[0],i,cn_inp.buf[i][0]);		
	  print_err(ERROR,pname,erbuf,NULL,NULL);
	  }	
	buflen=b[1]+2;

        /* if (buflen>1000){
	     sprintf(erbuf,"processed next licz=%i,buflen=%i,b[950]=%i from %s\n",b[0],buflen,b[950],proc_name[*j]);
	     print_err(DEBUG,pname,erbuf,NULL,NULL);
	     } */
	   
	for(ii=2;ii<buflen;ii++)
	  add_to_list(link_in[i]+b[ii]);
	FD_CLR(conn_sock[*j],&fdr1);
	if(b[ii]>-1)	 
	  cn_inp.buf[i][0]+=ii; /* index to next packet */	 
	else 
	  cn_inp.buf[i][0]=0;        /* no more packets left */
	k--;
	}
    }
 
  while(k)
    {
       
      memcpy(&fdr2,&fdr1,sizeof(fd_set));
  
      if((n=select(nfd,&fdr2,NULL,NULL,&tim))<0)
	print_err(PERROR,pname,"select",NULL,NULL);

      /* sprintf(erbuf,"Select result ,n=%i",n);
	 print_err(INFO,pname,erbuf,NULL,NULL); */

      if(n==0){ /* close gracefully */
        sprintf(erbuf,"Select time-out info,k=%i,licz=%i\n",k,licz);
	print_err(INFO,pname,erbuf,NULL,NULL);
	sigflag=-1;
	break;
      }

      j=cn_inp.lconn;

      for(i=0;n>0 && k>0 && i<cn_inp.nconn;i++,j++)
	if(FD_ISSET(conn_sock[*j],&fdr2))
	  {  int jj;
	  if((ii=myrecv(conn_sock[*j],(char*)buf_in,PACKETBUFFER,0))<0)	//was (char*)
	    print_err(PERROR,pname,"recv",NULL,NULL); 
	

          ii/=sizeof(int); //was (short)

          /* if (buf_in[1]>800){
	     sprintf(erbuf,"processed next licz=%i,buflen=%i,buf_in[790]=%i, ii=%i from %s\n",buf_in[0],buf_in[1]+2,buf_in[790],ii,proc_name[*j]);
	     print_err(INFO,pname,erbuf,NULL,NULL);
	     } */ 


	  if(ii==0){n--;continue;}
        
	  if(ii<2 || ii<(buflen=(buf_in[1]+2)))
	    {
	      sprintf(erbuf,"recv too short ,ii=%i b0=%i b1=%i from %s",ii,buf_in[0],buf_in[1],proc_name[*j]);
	      print_err(ERROR,pname,erbuf,NULL,NULL);
	    }
	  
	    sprintf(erbuf,"recv ii=%i,licz=%i,buflen=%i from %s\n",ii,buf_in[0],buflen,proc_name[*j]);
	    print_err(DEBUG,pname,erbuf,NULL,NULL); 

	  if(buf_in[0]!=licz) /* synchro error */
	    {
	      print_err(ERROR,pname,"sync",NULL,NULL);
	    }	

	  for(jj=2;jj<buflen;jj++)
	    {
	      add_to_list(link_in[i]+buf_in[jj]);
	    }
	   
	  if(ii>buflen) 
	    { 
	      /* recv more packets */


              /* sprintf(erbuf,"Condition ii>buflen ,ii=%i buflen=%i",ii,buflen);
		 print_err(INFO,pname,erbuf,NULL,NULL); */

	      /* check if all packets complete */
	      m=buflen;
	      while(m<ii)m+=buf_in[m+1]+2;
	      if(m!=ii)
		{
		  sprintf(erbuf,"recv next wrong length licz=%i,ii=%i,m=%i from %s\n",buf_in[buflen],ii,m,proc_name[*j]);
		  print_err(ERROR,pname,erbuf,NULL,NULL);
		}
	                                   
	      cn_inp.buf[i][0]=1; /* always one  because it doesn't recv until all previous processed */		    
	      buf_in[ii]=-1;   /* marker of end of buf */
	      //               fprintf(stderr,"%p %p %p\n",path,buf_in+ii,cn_inp.buf[i]+(ii-buflen+1)*2);
	      memcpy(cn_inp.buf[i]+1,buf_in+buflen,(ii-buflen+1)*sizeof(int)); //was sizeof(short)
		 
	    }
	  else
	    /* recv 1 packet */
	    cn_inp.buf[i][0]=0;           


	  FD_CLR(conn_sock[*j],&fdr1);
	  k--;n--;
	  }/* for i */
    }/* while(k) */  

     
      
  /* update synapses */	

  for(m=0;m<link_in_nproc;m++)
    {
      l_in=link_in[m];
      lw=linkw[m];
      for(n=0;n<link_in_no[m];n++,l_in++)
	{
	  neuron_j=neurons+l_in->neuron_post;
	  jd=l_in-link_in[m];
          if( l_in->weight > 0 )
            {
             //(lw+jd)->w += (om((lw+jd)->C) - (lw+jd)->w)/taul((lw+jd)->C);
             //I_NMDA = G_NMDA*((lw+jd)->V_d+(lw+jd)->V_o)*B_NMDA(neuron_j->V)*(V_NMDA((lw+jd)->C)-neuron_j->V);
             //(lw+jd)->C += I_NMDA-R_NMDA*((lw+jd)->C-Co);
             //(lw+jd)->V_d *=t1_NMDA;
             //(lw+jd)->V_o *=t2_NMDA;
#ifdef DEB
             //if( !(licz%10) ) //&& (l_in->neuron_post == 31))
             //  fprintf(lWfile,"%d,%lf,%d,%lf,%lf,%lf\n",(int)l_in->neuron_post,(lw+jd)->w,(int)l_in->neuron_post,(lw+jd)->C,I_NMDA,6*neuron_j->C[0]);
#endif
            }
            
	  if( (index = l_in->first) != EMPTY_LIST )
	    {
	      if(!l_in->interval[index])
		{ double Amp;

		if( l_in->weight > 0 )
		  {
		    Amp = l_in->weight *(params+neuron_j->param)->Aep[l_in->preclass_num];

                    //Amp = 4*lw[jd] * l_in->weight * (params+neuron_j->param)->Aep;
                    //Amp = 2.5*(lw+jd)->w * l_in->weight * (params+neuron_j->param)->Aep;
		    neuron_j->Ve_d[l_in->preclass_num] += Amp;
		    neuron_j->Ve_o[l_in->preclass_num] -= Amp;
                    //(lw+jd)->V_d += Amp;
                    //(lw+jd)->V_o -= Amp;
	
		  }
		else
		  {
		    Amp = l_in->weight *(params+neuron_j->param)->Aip[l_in->preclass_num-NEXCCLASS];

		    neuron_j->Vi_d[l_in->preclass_num-NEXCCLASS] += Amp;
		    neuron_j->Vi_o[l_in->preclass_num-NEXCCLASS] -= Amp;

		  }
			      	
		if(index==l_in->last)
		  l_in->first=EMPTY_LIST;
		else
		  l_in->first=index=(index+1)%LINK_LIST_LENGTH;
				 
		}
	      end_index = (l_in->last+1)%LINK_LIST_LENGTH;
		
	      while( index != end_index )
		{
		  --(l_in->interval[ index ]);
		  index = (index+1)%LINK_LIST_LENGTH;
		}
	    }
	} /* end for n */
    }/* end for m */

}



void read_links(char * name, int new_sim)
{
  FILE* fil;
  NEURON_INDEX ** l_out;
  int i,j,n,k;
  //int tallystim;
  struct LINK **l_in;
  struct LINKD **lw;
  //struct LINKD **lstimw;
  //double **lw;
  char buf[400];
  char *proc="read_links";
  //char erbuf[200];
  //char *pname="read_links";

  /* read output links */
  
  fil=fopen(file_name(buf,name,".lout"),"r");
  if(!fil) print_err(FERROR,proc,"fopen",buf,fil); 
    
  if(fread(&n,sizeof(int),1,fil)!=1 ||
     (n!=cn_out.nconn && n!=cn_out.nconn+1))
    print_err(FERROR,proc,"fread link_out_nproc",buf,fil);
  link_out_nproc=n;
  if(n>0)
    {
      link_out=l_out=(NEURON_INDEX**)malloc(sizeof(NEURON_INDEX*)*n);
      if(!l_out)print_err(ERROR,proc,"malloc(l_out)",NULL,NULL);

      link_out_no=(int*)malloc(sizeof(int)*n);
      if(!link_out_no)print_err(ERROR,proc,"malloc(link_out_no)",NULL,NULL);

      for(i=0,k=0;i<link_out_nproc;i++,l_out++)
        {
	  if(fread(&j,sizeof(int),1,fil)!=1)
	    print_err(FERROR,proc,"fread proc no",buf,fil);
	 
	  if(j!=this_proc && (cn_out.lconn==NULL || j!= cn_out.lconn[k++]))
	    print_err(FERROR,proc,"wrong proc no",buf,fil);
	  if(fread(&j,sizeof(int),1,fil)!=1)
	    print_err(FERROR,proc,"fread link_out_no",buf,fil);
	  link_out_no[i]=j;
	  *l_out=(NEURON_INDEX*)malloc(sizeof(NEURON_INDEX)*j);
	  if(!l_out)print_err(ERROR,proc,"malloc(*l_out)",NULL,NULL);
	  if(fread(*l_out,sizeof(NEURON_INDEX),j,fil)!=j)
	    print_err(FERROR,proc,"fread l_out",buf,fil);
	  /*         
		     for(n=0;n<j;n++)
		     { 
		     fprintf(stderr,"lout %i \n",link_out[i][n]);
		     fflush(stderr);
		     }
	  */
	}

           
    }

  fclose(fil);

  /* read input links */    
  if(new_sim)
    file_name(buf,name,".lin");
  else
    file_name(buf,name,".lin.bak");

  fil=fopen(buf,"r");
  if(!fil) print_err(FERROR,proc,"fopen read no of proc",buf,fil); 

  if(fread(&n,sizeof(int),1,fil)!=1 ||
     (n!=cn_inp.nconn && n!=cn_inp.nconn+1))
    print_err(FERROR,proc,"fread link_in_nproc",buf,fil);
  link_in_nproc=n;

  if(n>0){
    link_in=l_in=(struct LINK**)malloc(sizeof(struct LINK*)*n);
    if(!l_in)print_err(ERROR,proc,"malloc(l_in)",NULL,NULL);
    
    //linkw=lw=(double**)malloc(sizeof(double*)*n);
    linkw=lw= (struct LINKD **)malloc(sizeof(struct LINKD*)*n);

    link_in_no=(int*)malloc(sizeof(int)*n);
    if(!link_in_no)print_err(ERROR,proc,"malloc(link_in_no)",NULL,NULL);

    for(i=0,k=0;i<link_in_nproc;i++,l_in++,lw++)
      {
	if(fread(&j,sizeof(int),1,fil)!=1)
	  print_err(FERROR,proc,"fread proc no",buf,fil);
	 
	if(j!=this_proc && (cn_inp.lconn==NULL || j!= cn_inp.lconn[k++])) 
	  print_err(FERROR,proc,"wrong proc no",buf,fil);

	/*          fprintf(stderr,"link_in proc %i=%i : ",i,j); */
	if(fread(&j,sizeof(int),1,fil)!=1)
	  print_err(FERROR,proc,"fread link_in_no",buf,fil);
	link_in_no[i]=j;
	/*	     fprintf(stderr,"%i\n",j); */
	*l_in=(struct LINK*)malloc(sizeof(struct LINK)*j);
	if(!l_in)print_err(ERROR,proc,"malloc(*l_in)",NULL,NULL);
    
        //*lw=(double *)malloc(sizeof(double)*j);
        *lw=(struct LINKD *)malloc(sizeof(struct LINKD)*j);
    
	for(n=0;n<j;n++)
	  {
	    if(fread((*l_in)+n,sizeof(struct LINK),1,fil)!=1)
	      print_err(FERROR,proc,"fread l_in",buf,fil);	   
	    /*	         fprintf(stderr,"link %p l->neuron %i\n",link_in[i]+n,(link_in[i]+n)->neuron_post);
			 fflush(stderr); */
            //*((*lw)+n)= 0.25;    //((*l_in)+n)->weight;  
	    //((*lw)+n)->w   = 0.25;
	    //((*lw)+n)->V_o = 0;
	    //((*lw)+n)->V_d = 0;
	    //((*lw)+n)->C   = 0.0001;
           
	  }

      }
  }
  fclose(fil);	   
}


void read_params(char *name)
{
  FILE *fil;
  struct PARAMETERS *p;
  int  i,j,kk;
  char buf[400],*proc="read_params";
  params = p = malloc( sizeof( struct PARAMETERS )*no_kinds);
  if(p==NULL)print_err(ERROR,proc,"malloc(params)",NULL,NULL);
  fil = fopen(file_name(buf,name,".par"),"r");

  if(fil==NULL) print_err(FERROR,proc,"fopen",buf,fil);

  for(i = 0; i < no_kinds; i++,p++ )
    {
    
      if(fread(p,sizeof(struct PARAMETERS),1,fil)!=1)
	print_err(FERROR,proc,"fopen",buf,fil);


      if((p->del = malloc( sizeof( SPIKE_INTERVAL) * no_kinds ))==NULL)
	print_err(ERROR,proc,"malloc(p->del)",NULL,NULL);
      for(j = 0; j < no_kinds; j++)
	{
	  if(fread(p->del+j,sizeof( SPIKE_INTERVAL ),1,fil)!=1)
	    print_err(FERROR,proc,"fread(p->del)",buf,fil);
	}

      for (kk=0;kk<NEXCCLASS;kk++){

           p->de_o[kk]=exp(-DELTA_T*1000/p->de_o[kk]);
           p->de_d[kk]=exp(-DELTA_T*1000/p->de_d[kk]);

      }
       

      for (kk=0;kk<NINHCLASS;kk++){

           p->di_o[kk]=exp(-DELTA_T*1000/p->di_o[kk]); 
           p->di_d[kk]=exp(-DELTA_T*1000/p->di_d[kk]);

      }

    }

  fclose(fil);

}
void read_synapses(char *name)
{
  FILE* fil;
  struct SYNAPSE *s;
  struct SYNAPSESD *sd;
  //double *sd;
  int i;
  char buf[400],*proc="read_synapses";
  synapses = s = malloc( no_synapses * sizeof( struct SYNAPSE ) );
  if(!s)  print_err(ERROR,proc,"malloc(synapses)",NULL,NULL);
  synapsesd = sd = (struct SYNAPSESD *)malloc(no_synapses * sizeof( struct SYNAPSESD ) );
//synapsesd = sd = (double *)malloc(no_synapses*sizeof(double));
  if(!sd)  print_err(ERROR,proc,"malloc(synapsesd)",NULL,NULL);
  
  fil = fopen(file_name(buf,name,".syn"),"r");
  if(fil==NULL) print_err(FERROR,proc,"fopen",buf,fil) ;

  for(i = 0; i < no_synapses; i++,s++,sd++)
    {
      if(fread(s,sizeof( struct SYNAPSE ),1,fil)!=1)
	print_err(FERROR,proc,"fread(synapse)",buf,fil);
      sd->w   = 0.25;  //(double)s->weight;
      sd->V_d = 0;
      sd->V_o = 0;
      sd->C   = 0.0001;
      
    }

  fclose(fil);
  
}
void read_neurons(char *name,int new_sim)
{
  FILE *fil;
  int i,h=0,o=0;
  struct NEURON *n;
//  struct NEUROND *nd;
  char buf[400],*proc="read_neurons";
  neurons = n = malloc( no_neurons * sizeof( struct NEURON ) );
  if( neurons==NULL )  print_err(ERROR,proc,"malloc(neurons)",NULL,NULL);
//  neuronsd = nd = (struct NEUROND *) malloc( no_neurons * sizeof( struct NEUROND ) );
//  if( neuronsd==NULL )  print_err(ERROR,proc,"malloc(neuronsd)",NULL,NULL);

  if(new_sim)
    file_name(buf,name,".neu");
  else
    file_name(buf,name,".neu.bak");
  fil = fopen(buf,"r");            ;
  if(fil==NULL) print_err(FERROR,proc,"fopen",buf,fil) ;

  for(i = 0; i < no_neurons; i++,n++ ) /* ,nd++ ) */
    {

      if(fread(n,sizeof( struct NEURON ),1,fil)!=1)
	print_err(FERROR,proc,"fread(neurons)",buf,fil);
//       nd->C0=0.;
//       nd->C1=0.;
//       nd->C2=0.;
//       nd->C3=0.;
//       nd->C4=0.;
//      nd->C5=0.;

      if(histo_f && n->flags & HISTO) h++;
      if(out_f && n->flags & OUT ) o++;
      if(noise_f && n->flags & NOISE) sz++;
      if(iext_f && n->flags & IEXT) stz++;
    }
  if(histo_f)
    {
      if( h != his_nr )
	{
	  char errbuf[40];
	  sprintf(errbuf,proc,"his_nr:%i != no of hist flags: %i",his_nr,h);
	  print_err(ERROR,proc,errbuf,NULL,NULL);	
	}
    }
  if(out_f)
    {
      if( o != out_nr )
	{      char errbuf[40];
	sprintf(errbuf,proc,"out_nr:%i != no of out flags: %i",out_nr,o);
	print_err(ERROR,proc,errbuf,NULL,NULL);
	}
    }
  fclose(fil);
  
}

void write_neurons(char *name)
{
  FILE* fil;
  int i;


  struct NEURON *n = neurons;
  fil = fopen(name,"w");
  if(fil==NULL) 
    {  fprintf(stderr,"%s\tERROR:write_neurons: fopen(%s)\n",this_name,name);
    return;
    }
          
  for(i = 0; i < no_neurons; i++,n++ )
    {
      if(fwrite(n,sizeof( struct NEURON ),1,fil)!=1)
	{
	  fprintf(stderr,"%s\tERROR:write_neurons: fwrite(%s),n=%i\n",this_name,name,i);
	  return;
	}
    }

  fclose(fil);
}


void write_links(char *name)
{
  FILE* fil;
  int i,nproc;
  
  fil = fopen(name,"w");
  if(fil==NULL)
    {  fprintf(stderr,"%s\tERROR:write_links: fopen(%s)\n",this_name,name);
    return;
    }
  for(i = 0; i < link_in_nproc; i++ )
    {
      if(i<cn_inp.nconn)
	nproc=cn_inp.lconn[i];
      else
	nproc=this_proc;
      if(fwrite(&nproc,sizeof(int),1,fil)!=1)goto error;
       
      if(fwrite(link_in_no+i,sizeof(int),1,fil)!=1)goto error;
       
      if(fwrite(link_in+i,sizeof(struct LINK),1,fil)!=1)goto error;
       
    }

  fclose(fil);
  return;
 error:
  {
    fprintf(stderr,"%s\tERROR:write_links: fwrite(%s),n=%i\n",this_name,name,i);
    return;
  }
}


void  set_flags()
{
  /*signed */char z;
  /*unsigned*/ char buf[100],*l;
  l=buf;
  histo_f=out_f=vout_f=noise_f=iext_f=link_f=stim_f=0;  /* flagi sterujace*/

  fscanf(config,"%s",buf);
  do
    {
      sscanf(l,"%c,",&z);
      l+=2;

      switch(z)
	{
	case 'N': noise_f = 1;
	  break;
	case 'H': histo_f = 1;
	  break;
	case 'V': vout_f = 1;
	  break;
	case 'O': out_f = 1;
	  break;
	case 'I': iext_f = 1;
	  break;
	case 'L': link_f = 1;
	  break;
	case 'S': stim_f = 1;
	  iext_f = 1;
	case '$': break;
	default :
	  sprintf(buf,"Unknown option %c in cfg",z);
	  print_err(ERROR,"set_flags",buf,NULL,NULL);
	}
    }while (z != '$');

	  
  strcpy(buf,"Options:");
  if (noise_f) strcat(buf,"NOISE,");
  if (iext_f) strcat(buf,"IEXT,");
  if (histo_f) strcat(buf,"HISTO,");
  if (out_f) strcat(buf,"OUT,");
  if (vout_f) strcat(buf,"VOUT,");
  if (link_f) strcat(buf,"LINK,");
  if (stim_f) strcat(buf,"STIM,");
  print_err(INFO," ",buf,NULL,NULL);
	  
}

int calkuj( struct NEURON *neuron)
{
  struct PARAMETERS *par = params + neuron->param;
  /*	ushort V_int;   indeks do tablic funkcji(V) */
  VAR W_wp,V,W,Vd;//Cd,Ud;
  VAR X,B,I_Ca,I_K,I_Na,I_SYNE, I_SYNI;//,C,U;
  VAR dr;//,C0,C1,C2,C3,C4,C5,U0,U1,U2,U3,U4,U5;
//  VAR U0d,U1d,U2d,U3d,U4d,U5d,C0d,C1d,C2d,C3d,C4d,C5d;
  VAR C[7],U[7],Cd[7],Ud[7];
  /*      	VAR Vd,Wd,Cd,Xd,Bd; */
  VAR G_Na_m(VAR ,struct PARAMETERS *);
  /* VAR tau(VAR); */

  short V_out,i;
//  ushort C_out[7];
  int vcount =(int) (0.0001/DELTA_T);
  int kk;
  //char *pname="calkuj";
  //char erbuf[200];

  V = neuron->V;
  W = neuron->W;
  X = neuron->X;
  B = neuron->B;
//  C = neuron->C;
  for(i=0;i<7;i++)
  {     
   C[i] = neuron->C[i];
   U[i] = neuron->U[i];
  }
//  C0 = neuron->C0;
//  C1 = neuron->C1;
//  C2 = neuron->C2;
//  C3 = neuron->C3;
//  C4 = neuron->C4;
//  C5 = neuron->C5;
//  U0 = neuron->U0; 
//  U1 = neuron->U1; 
//  U2 = neuron->U2; 
//  U3 = neuron->U3; 
//  U4 = neuron->U4; 
//  U5 = neuron->U5;
  dr = par->dr;
       
       
  /*
    Vd=V;
    Wd=W;
    Cd=C;
    Xd=X;
    Bd=B;*/

  /*V_int = (short)rint(V*10.)+V_INT_OFFSET;*/
  if(!(licz%vcount)&& vout_f && neuron->flags & VOUT) 
    {
	if(fabs(V*50.)<SHRT_MAX) /* 32,767 */
	V_out = (short)rint(V*50.);
      else
	/*	if((V_out + V_INT_OFFSET) >= TAB_SIZE) */
	{
	  V_out=(V<0)?-SHRT_MAX:SHRT_MAX;
	  printf("\nV out of range %lf\n",V);
          //sprintf(erbuf,"time= %i, par->D= %11.10lf, par->dr=%11.10lf\n",licz,par->D,par->dr);
          //  print_err(INFO,pname,erbuf,NULL,NULL);
	}
      if(fwrite(&V_out,sizeof(V_out),1,Vfile)!=1)
	print_err(FERROR,"calkuj","fwrite","Vfile",Vfile);
      //for(i=0;i<7;i++)
      //C_out[i]=(ushort)rint(C[i]*100000);
      //fwrite(C_out,sizeof(ushort),7,Cfile);      
      
    }
	       
  W_wp = W*W;
  W_wp *= W_wp; /*  dla wp = 4 */

  I_K=par->G_K_s * W_wp
		
    +((par->G_K_Ca * C[6]) / (par->Kd + C[6]))
		    
    +(par->G_a) * A_inf(V) * B ;

//  I_Ca = (par->G_Ca)*X*X*(par->Kc/(par->Kc+C))*(V - par->V_Ca);
  I_Ca = (par->PCa)*X*X*(par->Kc/(par->Kc+C[6]))*15.05155*V*(C[6] - par->Cao * exp(-0.078*V))/(1-exp(-0.078*V));
  I_Na = G_Na_m(V,par) * ( (VAR)1. - W ) * (V - par->V_Na); 

  I_SYNE = 0;
  for (kk=0;kk<NEXCCLASS;kk++){

       I_SYNE += (neuron->Ve_o[kk] + neuron->Ve_d[kk])*(par->Esyn_e - V);

  }

  I_SYNI=0;
  for (kk=0;kk<NINHCLASS;kk++){

       I_SYNI += (neuron->Vi_o[kk] + neuron->Vi_d[kk])*(V - par->Esyn_i);

  }

  /* if ( (I_SYNE>.05) ||(I_SYNI<-0.05) ){
       sprintf(erbuf,"I_SYNE= %f, I_SYNI= %f time= %i\n",I_SYNE,I_SYNI,licz);
            print_err(INFO,pname,erbuf,NULL,NULL);
	    } */	


//    C += - par->Kp * I_Ca - par->R * C;
  Ud[6] = par->forw * C[6] * (1 - U[6]) - par->back * U[6];
  Cd[6] = - par->Kp * I_Ca + par->back * par->Utot * U[6] - par->forw* C[6] * par->Utot * (1-U[6]) - par->gpump * (C[6]/(C[6]+par->Kpump))+par->D*(par->Rsh-dr)*(C[5]-C[6])/(par->Rsh*dr*dr);
 
  Ud[5] = par->forw * C[5] * (1 -U[5]) - par->back * U[5];
  Cd[5] = par->back*par->Utot5*U[5]-par->forw*C[5]*par->Utot5*(1-U[5])+par->D*((par->Rsh5+dr)*C[6]-2*par->Rsh5*C[5]+(par->Rsh5-dr)*C[4] )/(par->Rsh5*dr*dr);
 
  Ud[4] = par->forw*C[4]*(1-U[4])-par->back*U[4];
  Cd[4] = par->back*par->Utot4*U[4]-par->forw*C[4]*par->Utot4*(1-U[4])+par->D*((par->Rsh4+dr)*C[5]-2*par->Rsh4*C[4]+(par->Rsh4-dr)*C[3] )/(par->Rsh4*dr*dr);
  
  Ud[3] = par->forw*C[3]*(1-U[3])-par->back*U[3];                                                                                                                 
  Cd[3] = par->back*par->Utot3*U[3]-par->forw*C[3]*par->Utot3*(1-U[3])+par->D*((par->Rsh3+dr)*C[4]-2*par->Rsh3*C[3]+(par->Rsh3-dr)*C[2] )/(par->Rsh3*dr*dr);
  
  Ud[2] = par->forw*C[2]*(1-U[2])-par->back*U[2];                                                                                                                
  Cd[2] = par->back*par->Utot2*U[2]-par->forw*C[2]*par->Utot2*(1-U[2])+par->D*((par->Rsh2+dr)*C[3]-2*par->Rsh2*C[2]+(par->Rsh2-dr)*C[1] )/(par->Rsh2*dr*dr);
  
  Ud[1] = par->forw*C[1]*(1-U[1])-par->back*U[1];                                                                                                                   
  Cd[1] = par->back*par->Utot1*U[1]-par->forw*C[1]*par->Utot1*(1-U[1])+par->D*((par->Rsh1+dr)*C[2]-2*par->Rsh1*C[1]+(par->Rsh1-dr)*C[0] )/(par->Rsh1*dr*dr);
  
  Ud[0] = par->forw*C[0]*(1-U[0])-par->back*U[0];                                                                                                                 
  Cd[0] = par->back*par->Utot0*U[0]-par->forw*C[0]*par->Utot0*(1-U[0])+3*par->D*(C[1]-C[0])/(par->Rsh0*par->Rsh0);
 

  
  X += (X_inf(V) - X)/par->tau_x;
	
	

  Vd= -I_Na
    -par->G_L * (V - par->V_L)
    -I_K*(V-par->V_K)
    -I_Ca
       
    + I_SYNI
    + I_SYNE
    ;
                           
  /* if (Vd>60){

       sprintf(erbuf,"Vd= %f,I_SYNE= %f, I_SYNI= %f time= %i\n",Vd,I_SYNE,I_SYNI,licz);
       print_err(INFO,pname,erbuf,NULL,NULL); 
       for (kk=0;kk<NEXCCLASS;kk++){
            sprintf(erbuf,"Ve_o[kk]= %f,Ve_d[kk]= %f, kk= %i\n",neuron->Ve_o[kk],neuron->Ve_d[kk],kk);
            print_err(INFO,pname,erbuf,NULL,NULL); 
       }
           
       } */
   
    
  // if( iext_f && neuron->flags & IEXT) Vd += par->i_ext;

  //Vd += (0.05); /* injecting 5 microamp/cm^2 inward current in all cells */

  W += (W_inf(V) - W)/tau(V);

              
  B += (B_inf(V) - B)/(par->tau_b);
 
  /*if(iext_f && neuron->flags & IEXT)
    fprintf(fout,"\nV=%f,dV=%f,I_Na=%f,I_K=%f,I_KCA=%f,I_A=%f,I_Ca=%f,I_EXT=%f,I_L=%f,I_SYNI=%f,I_SYNE=%f",
    neuron->V,Vd,I_Na,par->G_K_s * W_wp*(V-par->V_K), ((par->G_K_Ca * C) / (par->Kd + C))*(V-par->V_K), (par->G_a) * A_inf(V) * B *(V-par->V_K), I_Ca, par->i_ext, par->G_L * (V - par->V_L), I_SYNI, I_SYNE );*/ 
  neuron->V +=Vd;
  neuron->W = W;
  neuron->X = X;
  neuron->B = B;
//  neuron->C = C;
  for(i=0;i<7;i++)
  {     
     neuron->C[i] += Cd[i];
     neuron->U[i] += Ud[i];
  }
//  neuron->U0 += U0d;
//  neuron->U1 += U1d;
//  neuron->U2 += U2d;
//  neuron->U3 += U3d;
//  neuron->U4 += U4d;
//  neuron->U5 += U5d;
//  neuron->C1 += C1d;
//  neuron->C2 += C2d;
//  neuron->C3 += C3d;
//  neuron->C4 += C4d;
//  neuron->C5 += C5d;
  
  /* zmiana potencjalow Vi i Ve */


       for (kk=0;kk<NEXCCLASS;kk++){

            neuron->Ve_o[kk] *= par->de_o[kk];
            neuron->Ve_d[kk] *= par->de_d[kk];

       }

       for (kk=0;kk<NINHCLASS;kk++){

            neuron->Vi_o[kk] *= par->di_o[kk];
            neuron->Vi_d[kk] *= par->di_d[kk];

       }

  /*       

  if(fabs((V-Vd)/V)<EPS)
  nv++;
  if(fabs((W-Wd)/W)<EPS)
  nw++;
  if(fabs((C-Cd)/C)<EPS)
  nc++;
  if(fabs((X-Xd)/X)<EPS)
  nx++;
  if(fabs((B-Bd)/C)<EPS)
  nb++;  
  nn++; 
  */

  /* jesli jest spike return True*/
  if(neuron->flags & BYL_SPIKE )
    {
      if( V < AP_tr)
	neuron->flags &= ~BYL_SPIKE;
      if(neuron->flags& SPIKE_DET)
	neuron->flags &= ~SPIKE_DET;  
      return 0;
    }
  else
    if( V > AP_tr )
      {
	 
	/* Stim efficacy check routine */
	if( (iext_f==1) && (firststimflag==1) ){
	    
	     if ( (remainder(testpulsecount,2) || 0) && (neuron->polarity & POSPOL) ){
                  apcount=apcount+1;              		 
             }    

        }

	neuron->flags |= BYL_SPIKE;
	neuron->flags |= SPIKE_DET;
	return 1;
      }
    else
      return 0;

}


void szum( struct NEURON *neuron)
{
  double Amp;
  int p;


  if(drand48() < lambda)
    {
      //Amp = 2.5*(ssz+szum_i)->w*(params+neuron->param)->An;
      //Amp = 0.625*(params+neuron->param)->An;
      Amp = (params+neuron->param)->An;
      neuron->Ve_d[NOISECLASS] += Amp;
      neuron->Ve_o[NOISECLASS] -= Amp;
      //(ssz+szum_i)->V_d +=Amp;
      //(ssz+szum_i)->V_o -=Amp;     
      p=neuron - neurons;
      Noise_sum+=1;     
    }

}



void init()
{
    NEURON_INDEX i,kk,j;
  struct NEURON *n=neurons;
  /*   noise_f=0; */
  for(i = 0; i < no_neurons; i++,n++ )
    {
      n->V=V_0;
      n->W=W_0;
      n->X=X_0;
      n->B=B_0;
//      n->C=C_0;

      for (kk=0;kk<NEXCCLASS;kk++){

            n->Ve_o[kk] = 0;
            n->Ve_d[kk] = 0;

       }

       for (kk=0;kk<NINHCLASS;kk++){

            n->Vi_o[kk] = 0;
            n->Vi_d[kk] = 0;

       }

      for(j=0;j<7;j++)
      {     
       n->U[j]=U_0;
       n->C[j]=C_0;
      }
//      n->C1=C1_0;
//      n->C2=C2_0;
//      n->C3=C3_0;
//      n->C4=C4_0;
//      n->C5=C5_0;  
//      n->U0=U0_0;
//      n->U1=U1_0;
//      n->U2=U2_0;
//      n->U3=U3_0;
//      n->U4=U4_0;
//      n->U5=U5_0;
    }
}

int sszalloc(int si)
{
     int i;
     struct SYNAPSESD *szd;
     ssz = szd = (struct SYNAPSESD *)malloc(si * sizeof( struct SYNAPSESD ) );
     for(i=0;i<si;i++)
     {
          (szd+i)->w=0.25;
          (szd+i)->V_d=0;
          (szd+i)->V_o=0;
          (szd+i)->C=0.0001;
     }
          return i;
}


int stszalloc(int si)
{
     int i;
     struct SYNAPSESD *szd;
     stsz = szd = (struct SYNAPSESD *)malloc(si * sizeof( struct SYNAPSESD ) );
     for(i=0;i<si;i++)
     {
          (szd+i)->w=0.25;
          (szd+i)->V_d=0;
          (szd+i)->V_o=0;
          (szd+i)->C=0.0001;
     }
          return i;
}


int stim_freq(struct NEURON *neuron)
{
 double Amp;
 //char * pname="stim_freq";
 //char erbuf[200]; 

 /* sprintf(erbuf,"Made it to stim_freq, stim_period=%li\n",stim_period);		
    print_err(INFO,pname,erbuf,NULL,NULL); */ 

 if(stim_period) return 0;
    else
   {
     // Amp = 2.5*(stsz+stim_i)->w*(params+neuron->param)->An;
     //Amp = (params+neuron->param)->An*10.0;
     //Amp = 0;
     Amp   = .06;
     neuron->Ve_d[NOISECLASS] += Amp;
     neuron->Ve_o[NOISECLASS] -= Amp;
     // (stsz+stim_i)->V_d +=Amp;
     // (stsz+stim_i)->V_o -=Amp; 
     //  stim_period=stimulus_period;
     return 1; 
   }

}

int mysend(int sockcur, char *outbuf,int buflen, int flag) //was char for second arg
{

     int len,ret=0,retcheck,send_count=0;
     char * pname="mysend";
     char erbuf[200];

//     while (buflen>0 && send_count-->0){
     while(buflen>0){

	 send_count+=1;
	 len = (buflen>MAX_PACKET_SIZE)?MAX_PACKET_SIZE:buflen;
         
         if((retcheck=send(sockcur,outbuf,len,flag))!=len)
	 { if(errno==EAGAIN)
	     continue;
            print_err(PERROR,pname,"send ",NULL,NULL);
	 }
        
         ret+=retcheck;
         outbuf+=len;
         buflen-=len;

     }
     if(send_count>1){

          sprintf(erbuf,"send count limit: send_count=%i,licz=%i\n", send_count,licz);		
	      print_err(INFO,pname,erbuf,NULL,NULL);

     }

     return ret;

}     

int myrecv(int sockcur, char *inbuf,int bufmaxlen, int flag) //was char for second arg
{

     int k,len,l,n;
     char * pname="myrecv";
     char erbuf[200];
     struct timeval tim ={100,0};
     fd_set fdr;
  

     /* Read 1st packet */
     if((k=recv(sockcur,inbuf,bufmaxlen,flag))<0)	     
	    print_err(PERROR,pname,"recv",NULL,NULL); 

     if (k==0) return k;
 
     if(((int *)inbuf)[0]!=(int)licz) /* synchro error */ // was ((short *)inbuf)[0]
	    {
              sprintf(erbuf,"sync error in myrecv licz=%i,inbuf[0]=%i,inbuf[1]=%i,k=%i\n", licz,((short *)inbuf)[0],((short *)inbuf)[1],k);		
	      print_err(ERROR,pname,erbuf,NULL,NULL);
	    }	

     len=(((int *)inbuf)[1] + 2)*sizeof(int); //was ((short *)inbuf)[1] and sizeof(short)
     

     while (k<len){

	 FD_ZERO(&fdr);
	 FD_SET(sockcur,&fdr);
	 

         if((n=select(nfd,&fdr,NULL,NULL,&tim))<1)
	      print_err(PERROR,pname,"select",NULL,NULL);

         sprintf(erbuf,"Loop over packets in myrecv, k=%i, n=%i,licz=%i",k,n,licz);		
	      print_err(INFO,pname,erbuf,NULL,NULL);

	 if((l=recv(sockcur,inbuf+k,bufmaxlen,flag))<0)
	     print_err(PERROR,pname,"recv ",NULL,NULL);
	 
	 k+=l;
/*
	 if(k<0)
	 {
	      sprintf(erbuf,"k<0, k=%i, n=%i,licz=%i",k,n,licz);		
	      print_err(ERROR,pname,erbuf,NULL,NULL);
	 }
*/
}

     return k;
}
  
