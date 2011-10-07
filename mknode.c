/* Program do generacji polaczen synaptycznych 
   na wyjsciu ma produkowac plik z polaczeniami z zapisana para neuronow waga i
   delay. Nastenamepny program (lub druga czesc tego samego, moze wykorzystac to
   do generacji struktury neuronow i synaps we wlasciwym formacie.

   05/29/2003 P. Kudela, now multi-class supported  
   07/27/2005 S. Anderson 
*/
#define RISC
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "lnet.h"

#define randomm(num) (int)((double)num*(double)random()/((double)(RAND_MAX+1.0)))
#define randomize()     srandom((unsigned)time(NULL))

int lrodz,histoint;               /* number of neuron classes */
int *NN;                 /* NN[lrodz] total number of neurons within each class */
int *NOF;                /* NOF[lrodz] table of index offsets */ 
int *NB; 

int **N;             /* N[lrodz][3] size of subnetwork for each class   */
int **k;             /* k[lrodz][3] scaling factor to align class networks  */

int ***o;            /* o[lrodz][lrodz][3] neighborhood array  */
int **s;             /* s[lrodz][lrodz] synapsy na we/wy z otoczenie */
int **w,**sdw;    /* w[lrodz][lrodz] sdw[][] synaptic weights and dispersions */
int **d,**sdd;    /* d[lrodz][lrodz] sdd[][] synaptic delays and dispersions */

int *tab;         /* temporary neighborhood array tab[maxot]*/ 

int get_kind(int nr)
{
  int kind=0;
  while(NOF[++kind]<=nr);
  return kind-1;
}

struct SYN_TMP
{
NEURON_INDEX    neuron_pre;
NEURON_INDEX    neuron_post;
NEURON_INDEX    preclass_num;
WEIGHT          weight;
DELAY           delay;
};

int sort_delay(const void *, const void *);
int sort_neu_num(const void *, const void *);
void set_neu(struct NEURON *);
unsigned seed_save();

FILE *cfg;
int histo_f,out_f,vout_f,noise_f,iext_f,link_f,stim_f;  /* flags for simulator */
char flags[100];

void  set_flags()
{
  /*signed*/ char z;
  /*unsigned*/ char buf[100],*l;
  l=flags;
  histo_f=out_f=vout_f=noise_f=iext_f=link_f=stim_f=0;  /* flags for simulator */

  fscanf(cfg,"%s",flags);
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
	  fprintf(stderr,"Unknown option %c in cfg",z);
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
  puts(buf);
	  
}

int main(int argc, char  *argv[])
{
 int i,j,l,m,n,jp;
 int nr,ir,ix,iy,iz;  /* index for neurons */
 int maxot, tmpot;
 int xyz[3]={0,0,0};
 int no,nl;
 ulong ns,is;
 int *neu_pre;
 int Nneu;
 int torus_flag;
 int syn_n;
 char label[]="xyz",name[100],fname[100],buff[10],strdummy[10];
 FILE *inputfile,*outputfile,*par,*synfile,*neufile,*opar,*nodefile;
 struct SYN_TMP *syn_tmp;
 struct SYNAPSE syn;
 struct NEURON *neu;
 struct PARAMETERS p;
 int nproc,nsub,iproc,isub,vsub,nproc2,nsub2;
 SPIKE_INTERVAL delays[MAX_KIND][MAX_KIND];
 int otoczenie(int,int,int*,int);
 int otoczenie2(int,int,int*,int);
 int losuj(int*,int);
 int dummy,gumba;
 double testrandnoise;
 extern int zapisz(int,int,int,int,FILE *);
 extern int sort_neu_num(const void *, const void *);

 int **flag;
 int **flagv; 
 if(argc<5)
   {fprintf(stderr,"USAGE:\n\tmknodes <inp_name> <out_fname> <sqrt(subnets/node)> <sqrt(nodes)> [<seed>]\n");
   return 1;
   }

  if(argc>5)
   srandom((unsigned)atoi(argv[5]));
  else
   printf("seed %u\n",seed_save()); /* seed for srandom */

 strcpy(name,argv[1]);
  inputfile = fopen(strcat(name,".sub"),"r");
  if(inputfile==NULL){fprintf(stderr, "Subnetwork definition file not present %s\n",name);return -1;}

 strcpy(name,argv[1]); /*input file name*/
 strcpy(fname,argv[2]); /*output file name*/
 nsub2=atoi(argv[3]);  /*sqrt subnets per node*/
 nproc2=atoi(argv[4]); /*sqrt nodes*/
 nsub=nsub2*nsub2;    /* number of subnets per node */
 nproc=nproc2*nproc2; /* number of nodes */
 
 fscanf(inputfile,"%d",&torus_flag);
 fscanf(inputfile,"%d",&lrodz);

 cfg=fopen(strcat(name,".cfg"),"r");
 if(cfg==NULL)
    {fprintf(stderr,"Can't open configuration file %s\n",name);return 2;}

  set_flags();
  fscanf(cfg,"%d",&dummy);
  fscanf(cfg,"%d",&dummy);
  fscanf(cfg,"%s",&strdummy);
  fscanf(cfg,"%d",&dummy);
  fscanf(cfg,"%s",&strdummy);
  fscanf(cfg,"%d",&dummy);
  fscanf(cfg,"%d,%d",&dummy,&histoint);
  /* printf("dummy=%d,histoint=%d\n",dummy,histoint); */
 
  fclose(cfg);
 
  flag=(int**)calloc(nsub,sizeof(int*));
  for(i=0;i<nsub;i++)
    flag[i]=(int*)calloc(nproc,sizeof(int));  /* Subnet# X Node# flag array */ 

  flagv=(int**)calloc(nsub,sizeof(int*));     /* Subnet# X Node# flag array */
  for(i=0;i<nsub;i++)
    flagv[i]=(int*)calloc(nproc,sizeof(int));
 
                                       /* table for rodzajow neurons  */
 N=(int**)malloc(lrodz*sizeof(int*));        /* table allocation for size of networks for each class  */
  for(i=0;i<lrodz;i++)
    N[i]=(int*)malloc(3*sizeof(int));
 k=(int**)malloc(lrodz*3*sizeof(int));        /* table allocation for network scaling factor for each class  */
  for(i=0;i<lrodz;i++)
    k[i]=(int*)malloc(3*sizeof(int)); 

 o=(int ***)malloc(lrodz*sizeof(int **));  /* allocating neighborhood tables  */
   for(i=0;i<lrodz;i++)
     o[i]=(int **)malloc(lrodz*sizeof(int *));
    for(i=0;i<lrodz;i++)  
     for(j=0;j<lrodz;j++)
        o[i][j]=(int *)malloc(3*sizeof(int));

 s=(int **)malloc(lrodz*sizeof(int*));           /* synaptic table allocation */
   for(i=0;i<lrodz;i++)
    s[i]=(int*)malloc(lrodz*sizeof(int));

   w=(int **)malloc(lrodz*sizeof(int*));           /* allocation for synaptic weight storage */
   for(i=0;i<lrodz;i++)
    w[i]=(int*)malloc(lrodz*sizeof(int));

 sdw=(int **)malloc(lrodz*sizeof(int*));           /* allocation for synaptic weight dispersion  */
   for(i=0;i<lrodz;i++)
    sdw[i]=(int*)malloc(lrodz*sizeof(int));

 d=(int **)malloc(lrodz*sizeof(int*));           /* allocation for synaptic delay */
   for(i=0;i<lrodz;i++)
    d[i]=(int*)malloc(lrodz*sizeof(int));

 sdd=(int **)malloc(lrodz*sizeof(int*));           /* allocation for synaptic delay dispersion */
   for(i=0;i<lrodz;i++)
    sdd[i]=(int*)malloc(lrodz*sizeof(int));


/* input parameters from files to fill arrays */

strcpy(buff,fname); /* copying output file name to 'buff'*/
nodefile = fopen(strcat(buff,".node"),"wt"); /* file used by jhist program for linking histograms */
fprintf(nodefile,"%d\n",lrodz);

for(i=0;i<lrodz;i++)
  {
   fscanf(inputfile,"%d,",&nr);
   fprintf(nodefile,"%d\n",nr);
     N[i][0]=N[i][1]=nr; /* square mapping */
     N[i][2]=1;
     k[i][0]=k[i][1]=k[i][2]=1;
  }
 fprintf(nodefile,"%d\n",nsub2);
 fprintf(nodefile,"%d\n",nproc2);
 fclose(nodefile);
 fscanf(inputfile,"%d",&nr); /* reading summation check on total number of neurons per subnetwork */

/* setting neighborhood and other input parameters  */

 for(i=0;i<lrodz;i++)
  for(j=0;j<lrodz;j++)
   { 
     fscanf(inputfile,"%d,%d,%d",&o[j][i][0],&o[j][i][1],&o[j][i][2]); /* connection range for given neuron */
     fscanf(inputfile,"%d",&s[j][i]);  /* number of links from given neuron */
     fscanf(inputfile,"%d",&w[j][i]);  /* synaptic weights for given neuron */
     fscanf(inputfile,"%d",&sdw[j][i]);/* weight dispersion */
     fscanf(inputfile,"%d",&d[j][i]);  /* synaptic delay */
     fscanf(inputfile,"%d",&sdd[j][i]);/* delay dispersion */
   }
 
 fclose(inputfile);

/* Setting up noise flags in array */

 iproc=(nproc2/2)*(nproc2+1); /* works for even case */


if(nproc%2)
   isub=0;
 else
   isub=(nsub2/2)*(nsub2+1);


/* Old Central Noise Source Routine */
//if(isub=0) /* central noise source for oddXodd node array including case of one processor */
//   {
//     isub=(nsub2/2)*(nsub2+1);
//     iproc = (nproc-1)/2;
//     flag[isub-nsub2-1][iproc]=flag[isub-nsub2][iproc]=flag[isub-1][iproc]=flag[isub][iproc]=4;  // 4 neurons marked with noise flag
//     printf("4 neurons with NOISE in subnets %d,%d,%d,%d in node %d\n",isub-nsub2-1,isub-nsub2,isub-1,isub,iproc);
//   }
// else  /*sets central four neurons with noise flags in evenXeven case*/
//   {
//   flag[nsub-1][iproc-nproc2-1]=flag[nsub-nsub2][iproc-nproc2]=flag[nsub2-1][iproc-1]=flag[0][iproc]=4; 
//    printf("4 neurons with NOISE in subnet/node %d/%d,%d/%d,%d/%d,%d/%d\n",
//   	   nsub-1,iproc-nproc2-1,nsub-nsub2,iproc-nproc2,nsub2-1,iproc-1,0,iproc);
   //for(i=0;i<nsub;i++)
   //  for(j=0;j<nproc;j++)
   //    flag[i][j]=4;
//   }


/* Setting up vout flags in array for use with SVE program */
   if(nsub2%2)
     vsub=(nsub2*(nsub2-1))/2+(nsub2-1)/2; /* for odd, provides index of middle subnetwork */
   else
   vsub=nsub2*nsub2/2-nsub2/2;
 
   flagv[vsub][0]=flagv[vsub][nproc2-1]=flagv[vsub][nproc2*(nproc2-1)]=flagv[vsub][nproc-1]=nr; /* nr is the total number of neurons/subnet */  

/*

OFFSET CALCULATIONS TO FILL TABLES

*/

 NN=(int*)malloc(lrodz*sizeof(long int));   /* table with numbers of neurons within each class */
 for(i=0;i<lrodz;i++)
  NN[i]=(N[i][0]?N[i][0]:1)*(N[i][1]?N[i][1]:1)*(N[i][2]?N[i][2]:1);
 
 NOF=(int*)malloc((lrodz+1)*sizeof(int));  /* table with offsets */
  for(i=1,NOF[0]=0;i<=lrodz;i++)
    NOF[i]=NOF[i-1]+NN[i-1];

 NB=(int*)calloc(lrodz,sizeof(int));     /* base no of neurons in subnets */
  if (NB==NULL){fprintf(stderr,"Insufficient memory for  NB\n");exit(1);}

 

/* maximal neighborhood calculation */

 maxot=0;
 if(torus_flag!=2)
 for(i=0;i<lrodz;i++)
  for(j=0;j<lrodz;j++)
   {
     tmpot=1;
     for(l=0;l<3;l++)
     tmpot*=2*o[i][j][l]+1;
      if(tmpot > maxot)
       maxot=tmpot;
   }
  else
   for(i=0;i<lrodz;i++)
    {
    if(maxot<NN[i])
       maxot=NN[i];
    } 
 tab=(int*)malloc(maxot*sizeof(int));

 Nneu=0;
 for(i=0;i<lrodz;i++)
     Nneu+=NN[i];   /* Calculate the total number of neurons in the subnetwork */

 if(Nneu!=nr) {fprintf(stderr,"\nerror: check lines 2 & 3 in input file\n");
   return 2;
   }else
     Nneu*=nsub; /* Total number of neurons of all classes per node */

 if(Nneu>MAX_NEURON)
   {fprintf(stderr,"No of Neurons per node : %d > %d\n",Nneu,MAX_NEURON);
   return 2;
   }
 //printf("Nneu=%i\n",Nneu);

 neu=(struct NEURON *)calloc(Nneu,sizeof(struct NEURON));
   if (neu==NULL){fprintf(stderr,"Insufficient memory for neu");exit(1);}

   for(iproc=0;iproc<nproc;iproc++) /* nproc is number of nodes, big loop over nodes */
   {
     ulong nsum=0;

     sprintf(name,"%s.syncomp.%d",fname,iproc);
     if(!(synfile=fopen(name,"w")))
	  {fprintf(stderr,"Can't open %s\n",fname);return 2;}
     sprintf(name,"%s.neutmp.%d",fname,iproc);
     if(!(neufile=fopen(name,"w")))
          {fprintf(stderr,"Can't open %s\n",fname);return 2;}
     for(i=0;i<lrodz;i++)
	 NB[i]=NOF[i]*(nsub-1);
     ns=0;                           /* number of synapses, talleyed below when writing synapse file  */

     sprintf(name,"tmpsyn.%d",iproc);
     if(!(outputfile=fopen(name,"wb")))
          {fprintf(stderr,"Can't open %s\n",name);return 2;}
     //outputfile=fopen("syn.tmp","wb");
     //if(outputfile==NULL){fprintf(stderr,"\nError in opening output file");exit(1);}

     for(isub=0;isub<nsub;isub++)
        {
	 
          //for(i=0;i<lrodz;i++)
          //printf("\nNN[%d]=%d, NOF[%d]=%d",i,NN[i],i,NOF[i]); /* checks */
          //printf("\nmaxot=%d,Nneu=%d\n",maxot,Nneu);
     
          /*

          3.Neighborhood calculations

          */
 

          /*randomize();*/
          nr=0;                           /* neighborhood counting */

          for(ir=0;ir<lrodz;ir++)  /* loop over neuron class */
             for (ix=0,xyz[0]=0;ix<N[ir][0];ix++,xyz[0]+=k[ir][0])     /* x-loop */
                for(iy=0,xyz[1]=0;iy<N[ir][1];iy++,xyz[1]+=k[ir][1])    /* y-loop */
                   for(iz=0,xyz[2]=0;iz<N[ir][2];iz++,xyz[2]+=k[ir][2],nr++)  /* z-loop  */
                      for(jp=0;jp<lrodz;jp++)     /* loop on classes of neurons */
                      {
                              /*     wygeneruj tablice z numerami neuronow rodzaju jp 
                                     w otoczeniu neuronu nr (czyli neuronu rodzaju ir,
                                     w polozeniu ix,iy,iz) np.*/
			  switch(torus_flag) /* tabulating number of neurons available for connection */
                              {
                                  case 0:
				     no=otoczenie(nr,ir,xyz,jp);  /* no torus boundary condition */
                                     break;
                                  case 1:
                                     no=otoczenie2(nr,ir,xyz,jp); /* torus boundary condition */
                                     break;
                                  case 2:                         /* case 2 choosing from whole array  */
                                     for(i=NOF[jp],no=0;i<NN[jp]+NOF[jp];i++)
                                     if(i!=nr)
                                     tab[no++]=i;
                                     break;
                                  default:
                                     fprintf(stderr,"wrong option for torus flag %d\n",torus_flag);
                                     exit(0);
                              }


			  if(!no)continue;
			  l=(no<s[ir][jp])?no:s[ir][jp];
			  /* random connection formation */
			  for(i=0;i<l;i++) 
			  {
			      /* random connection formation */
			      nl=losuj(tab,no); /*printf("%d,%d,%d,%d\n",ir,jp, nr,nl);*/
			      /* filling of tables for delay/synaptic weights */
			      ns+=zapisz(nl,nr,ir,jp,outputfile); /*writing output to file 0.neu, 0.syn  */    

			      /* printf("\npostindex=%d, postclass=%d, preindex=%d, preclass=%d",nr,ir+1,nl,jp+1); */ /* print check variables */
                 
			  } /* end of drawing */
                      } /* jp */
	              for(i=0;i<lrodz;i++)
	              NB[i]+=NN[i];

        } // end of isub

  fclose(outputfile);
 

  sprintf(name,"tmpsyn.%d",iproc);
  if(!(outputfile=fopen(name,"r")))
       {fprintf(stderr,"Can't open %s\n",name);return 2;}
  //outputfile=fopen("syn.tmp","r");

  syn_tmp=(struct SYN_TMP *)calloc(ns,sizeof(struct SYN_TMP)); /* tablica wszystkch syn_tmp */
   if (syn_tmp==NULL){fprintf(stderr,"Can't allocate temporary synapse memory space\n");exit(1);}
  fread(syn_tmp,sizeof(struct SYN_TMP),(size_t)ns,outputfile); /* czytaj syn_tmp z syn.tmp */   
  fclose(outputfile);
  //printf("\nns=%d iproc=%d\n",ns,iproc);/* print check variables */
  //unlink(name);
  //unlink("syn.tmp");  /* syn.tmp nie jest dalej potrzebny */

  qsort( (void *)syn_tmp,(size_t)ns,sizeof(struct SYN_TMP),sort_neu_num ); /* sorting synapses by presynaptic neuron index for simulator use
									    sort_neu_num defined below*/ 
  
  neu_pre=calloc(Nneu,sizeof(int));                 /* table of neurons with number of synapses (presynaptic neurons)*/
  is=0;syn_n=0;
  for(i=0;i<Nneu;i++)
   {
    neu_pre[i]=0;
    while(i==(syn_tmp+syn_n)->neuron_pre){
    syn_n++;
    neu_pre[i]++; }
    qsort((void *)(syn_tmp+is),(size_t)neu_pre[i],sizeof(struct SYN_TMP), sort_delay );   /* sort by delay */
    is=syn_n;
   } 



  for(is=0;is<ns;is++)
   {
     //     int kind=get_kind((syn_tmp+is)->neuron_post);
    syn.neuron=(syn_tmp+is)->neuron_post;
    syn.weight=(syn_tmp+is)->weight;
    syn.delay=(syn_tmp+is)->delay;
    //printf("=%d,iproc=%d,i=%d,j=%d\n",isub,iproc,i,j);
    syn.preclass_num=(syn_tmp+is)->preclass_num;
    fwrite(&syn,sizeof(struct SYNAPSE),1,synfile);
   }

  nsum+=ns;
  free(syn_tmp);
 
  /* teraz neurony */
 
  for(i=0;i<Nneu;i++){
     (neu+i)->syn_num=neu_pre[i];   /* preparing to write neuron files, structure NEURONS */
     set_neu(neu+i);}        /* setting initial voltage and other values */

  

  for(i=0;i<lrodz;i++)
     for(j=0;j<NN[i]*nsub;j++)
        (neu+NOF[i]*nsub+j)->param=i;      /* assigning type of neuron */

  /* noise flags */
  //for(isub=0;isub<nsub;isub++)
  //   { l=isub*NN[0];
  //      if(flag[isub][iproc]){      
  //         for(i=0;i<flag[isub][iproc];i++)
  //         {
 	      //j=randomm(NN[0])+l;
	      //j=l+N[0][0]/2*(N[0][0]+1)+i;
  //	      j=l+i; 
              //j=isub*N[0][0]*N[0][0]+i;
  //            printf("Noise flags ** isub=%d,iproc=%d,i=%d,j=%d\n",isub,iproc,i,j);
  //	      (neu+j)->flags|=NOISE;  /* only 1st kind with NOISE */
  //         } 
  //      }    /* end if */
  //   }


  /* noise out flags */
  gumba=0;
  for(isub=0;isub<nsub;isub++)
  { 

           l=0;
           for(n=0;n<lrodz;n++)
           {
              for(i=0;i<NN[n];i++)  
              {
                 j=l+isub*NN[n]+i;
                 //printf("Noise flags ** isub=%d,iproc=%d,i=%d,j=%d\n",isub,iproc,i,j);
		 testrandnoise=(double)random()/((double)RAND_MAX+1);
		 // printf("testrandnoise=%f\n",testrandnoise);
		 if( (testrandnoise<0.01) && ( (n==0) || (n==2) ) ){
                      (neu+j)->flags|=NOISE;
		 }

              } 
              l+=nsub*NN[n];
            }
         
  }

  /* vout flags */
  for(isub=0;isub<nsub;isub++)
     { 

        //l=isub*NN[0];
        //m=nsub*NN[0]+isub*NN[1];  

        if(flagv[isub][iproc])
        {
           l=0;
           for(n=0;n<lrodz;n++)
           {
              for(i=0;i<NN[n];i++)  
              {
                 j=l+isub*NN[n]+i;
                 printf("Voltage out flags ** vsub=%d,iproc=%d,i=%d,j=%d,class=%d\n",isub,iproc,i,j,n);
                 (neu+j)->flags|=VOUT;  /* only 1st kind VOUT */
                 //(neu+j)->flags|=IEXT; 
              } 
              l+=nsub*NN[n];
            }
         }

         /*
         for(i=0;i<flagv[isub][iproc];i++)
         {
            j=l+i;
            printf("vsub=%d,iproc=%d,i=%d,j=%d\n",isub,iproc,i,j);
            (neu+j)->flags|=VOUT;  // only 1st kind VOUT 
         }

         for(i=0;i<sqrt(flagv[isub][iproc]);i++) // inhibitorty 
         {
            j=m+i;
            printf("vsub=%d,iproc=%d,i=%d,j=%d\n",isub,iproc,i,j);
            (neu+j)->flags|=VOUT;  // 2nd kind VOUT 
         }
         */
     }
 
  fwrite(neu,sizeof(struct NEURON),Nneu,neufile);

 
  fclose(synfile);
  fclose(neufile);
  
  sprintf(name,"%s.cfg.%d",fname,iproc);
  cfg=fopen(name,"w");   /* zbior *.cfg */
  //fprintf(outputfile,"%s\n","ONSET");         /* creating .cfg file */
  fprintf(cfg,"%s\n",flags);                  /* */
  fprintf(cfg,"%s\n",fname);                  /* */
  fprintf(cfg,"%d\n",lrodz);
  fprintf(cfg,"%s\n",fname);                  /* */
  fprintf(cfg,"%lu\n",nsum);             /* */
  fprintf(cfg,"%s\n",fname);                  /* */   
  fprintf(cfg,"%d\n",Nneu);            /* */
  fprintf(cfg,"%d,%d\n",Nneu,histoint);      /* histogram parameters-length of bin  */ 
  fprintf(cfg,"0\n");                /* out_no*/
  fprintf(cfg,"0.00003068");            /* the lambda actually used by simulator, was 3.068X10^-4 */
  fclose(cfg); 
 } // end iproc


    /* copy par file */

  strcpy(name,argv[1]);
  par=fopen(strcat(name,".par"),"r");
  if(par==NULL)
    {fprintf(stderr,"Can't open %s\n",name);return 2;}
 
  for(i=0;i<nproc;i++)
    {
      int ni=0;
      sprintf(name,"%s.par.%d",argv[2],i);
      opar=fopen(name,"w");
      if(opar==NULL)
	{fprintf(stderr,"Can't open %s\n",name);return 2;}
      rewind(par);
      while(ni<lrodz && fread(&p,sizeof(struct PARAMETERS),1,par)==1 && 
	    fread(&delays[ni][0],sizeof( SPIKE_INTERVAL ),lrodz,par)==lrodz)
	{
	     
	  fwrite(&p,sizeof(struct PARAMETERS),1,opar);
	  fwrite(&delays[ni][0],sizeof( SPIKE_INTERVAL ),lrodz,opar);
	  ni++;
	}
      fclose(opar);
    }
  fclose(par);

  /* copy stim file */
  if(stim_f)
     {
        strcpy(name,argv[1]);
        par=fopen(strcat(name,".stim"),"r");
        if(par==NULL)
	   {fprintf(stderr,"Can't open %s\n",name);return 2;}
    
        for(i=0;i<nproc;i++)
	{
	   rewind(par);
	   sprintf(name,"%s.stim.%d",argv[2],i);
	   opar=fopen(name,"w");
	   if(opar==NULL)
	      {fprintf(stderr,"Can't open %s\n",name);return 2;}
	   while(!feof(opar))
	      {  char buf[100];
	         if(fscanf(par,"%[^\n]\n",buf)!=1)break;
	         fprintf(opar,"%s\n",buf);
	      }
	    fclose(opar);
	      
	 }
         fclose(par);
      }
  

printf("\n");
return 0;
} /*koniec main()*/


/*
sort_delay  - needed in qsort
*/

int sort_delay(const void * a, const void * b )
{
char delay1[1],delay2[1];
struct SYN_TMP *syn1=NULL, *syn2=NULL;
syn1=(struct SYN_TMP *)a;
syn2=(struct SYN_TMP *)b;
*delay1=syn1->delay;
*delay2=syn2->delay;

return( strcmp(delay1,delay2) );
}


/*
sort_neu_num - needed in qsort() 
*/

int sort_neu_num(const void * a, const void * b)
{
NEURON_INDEX neu1,neu2;
struct SYN_TMP *syn_tmp1=NULL,*syn_tmp2=NULL;
syn_tmp1=(struct SYN_TMP *)a;
syn_tmp2=(struct SYN_TMP *)b;
neu1=syn_tmp1->neuron_pre;
neu2=syn_tmp2->neuron_pre;
return(neu1-neu2);
}

/*
  
4. Generation of neighborhood table:

*/

int rozrzut(int sd)
{
 int los;
 los=randomm(100);
  if(los < 50 )
   return( -1*randomm(sd) );
  else
   return(    randomm(sd) );
}


int losuj(int *tablica,int index) /* random connections within identified regions of allowed connections */
{
int i,j;
do
{
i=randomm(index);
}while(tablica[i]<0);
j=tablica[i];
tablica[i]=-1;
return j;
}

int zapisz(int inr, int inl, int i, int j, FILE * outputfile)
{
struct SYN_TMP syn_tmp; 
extern  int rozrzut(int);

syn_tmp.neuron_post=inl+NB[i];
syn_tmp.neuron_pre=inr+NB[j]; 
syn_tmp.preclass_num=j;      
syn_tmp.weight=w[i][j]+rozrzut(sdw[i][j]);
   /* obliczamy najpierw wspolrzedne i[3] neuronu na siatce jp, ktory jest najblizej x,y,z */
syn_tmp.delay =d[i][j]+rozrzut(sdd[i][j]);
//printf("pre= %d,post= %d,preclass= %d, weight= %d,delay= %d\n",syn_tmp.neuron_pre,syn_tmp.neuron_post,syn_tmp.preclass_num,syn_tmp.weight,syn_tmp.delay);
return(fwrite(&syn_tmp,sizeof(struct SYN_TMP),1,outputfile));
}


int otoczenie(int nr ,int ir,int xyz[3],int jp)
{
   int nx,ny,no,i[3],ix,iy,iz,mini[3],maxi[3], *ot;
   int j;

   /* tab must be declared earlier with sufficient size */

   ot = o[ir][jp];
   
   /*   The next line calculates the coordinates of neurons from one class to the others */  

  no=0;
  for(j=0;j<3;j++)
         i[j]=xyz[j]/k[jp][j];  

  /* No torus boundary condition */

  for(j=0;j<3;j++)
  {
     mini[j]=i[j]-ot[j];
     maxi[j]=i[j]+ot[j]+1;
     if(mini[j]<0)mini[j]=0;
     if(maxi[j]>N[jp][j])maxi[j]=N[jp][j];
  }
   
   for(ix=mini[0];ix<maxi[0];ix++)
      { nx=ix*N[jp][1]*N[jp][2]+NOF[jp];
      for(iy=mini[1];iy<maxi[1];iy++)
         { ny=nx+iy*N[jp][2];
         for(iz=mini[2];iz<maxi[2];iz++)
         {
            if((j=ny+iz)!=nr)
            tab[no++]=j;   
         }
      }
   }

return no;
}


int otoczenie2(int nr, int ir,int xyz[3],int jp)
{
   int nx,ny,no,i[3],ix,iy,iz,mini[3],maxi[3], *ot;
   int j,ii,NX,NY,NZ;
   /* Torus boundary condition */

   ot = o[ir][jp];
   /* calculating location of neurons used to create connections  x,y,z */
  
  
  no=0;
  for(j=0;j<3;j++)
         i[j]=xyz[j]/k[jp][j];  

  /* teraz otoczenie
     wersja z zawijaniem na torusie */

      
 for(j=0;j<3;j++)
 {
   mini[j]=i[j]-ot[j];
   maxi[j]=i[j]+ot[j]+1;
 }
   NX=N[jp][0];
   NY=N[jp][1];
   NZ=N[jp][2];
 
   for(ix=mini[0];ix<maxi[0];ix++)
   {  ii=(ix+NX)%NX;
      nx=ii*NZ*NY+NOF[jp];
      for(iy=mini[1];iy<maxi[1];iy++)
      {  ii=(iy+NY)%NY;
         ny=nx+ii*NZ;
         for(iz=mini[2];iz<maxi[2];iz++)
         {
             if((j=ny+(iz+NZ)%NZ)!=nr)
             tab[no++]=j;
         }
      }
    }
return no;
} 


void set_neu(struct NEURON * neu)
{
   int j,kk;
   switch(neu->param)
   {
      default:
      neu->flags = HISTO;
      neu->polarity = 0x0;

      for (kk=0;kk<NEXCCLASS;kk++){	

           neu->Ve_o[kk]=0;
           neu->Ve_d[kk]=0;

      }


      for (kk=0;kk<NINHCLASS;kk++){
	
           neu->Vi_o[kk]=0;
           neu->Vi_d[kk]=0;

      }


      neu->interval=0;
      neu->sum_interval=0;
      for(j=0;j<SPIKE_LIST_LENGTH;j++)
         neu->spikes[j]=0;
      neu->first=EMPTY_LIST;
      neu->last=0x1;
      neu->V=V_0;
      neu->W=W_0;
#ifdef CALCIUM
      for(j=1;j<7;j++)
      {     
         neu->C[j]=0.;
         neu->U[j]=0.;
      } 
//      neu->C=C_0;
      neu->X=X_0;
      neu->B=B_0;
#endif
   }
 
}

unsigned seed_save()
{
 unsigned seed;
 seed = (unsigned)time(NULL);
 srandom(seed);
return seed;
}


