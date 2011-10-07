/* Program do generacji polaczen synaptycznych miedzy subnets
   Na wyjsciu ma produkowac plik z polaczeniami z zapisana para neuronow waga i
   delay. Nastepny program (lub druga czesc tego samego) moze wykorzystac to
   do generacji struktury neuronow i synaps we wlasciwym formacie.

   01/29/03 P. Kudela, losuj2 prawdopodobienstwo polaczenia zanika exp.  
                         a w zapisz2 delay wzrasta liniowo  
   05/29/03 P. Kudela, now multi-class supported 
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

int lrodz;               /* number of neuron classes */
int *NN;                 /* NN[lrodz] total number of neurons within each class */
int *NOF;                /* NOF[lrodz] table of index offsets */ 

int **N;             /* N[lrodz][3] size of subnetwork for each class */
int **k;             /* k[lrodz][3] scaling factor to align class networks */

int ***o;            /* o[lrodz][lrodz][3] neighborhood array */
int **s;             /* s[lrodz][lrodz] synapsy na we/wy z otoczenie */
int **w,**sdw;    /* w[lrodz][lrodz] sdw[][] synaptic weights and dispersions  */
int **d,**sdd;    /* d[lrodz][lrodz] sdd[][] synaptic delays and dispersions */
double **Ap,**lambda;/* Ap[lrodz][lrodz] lambda[][] amplitude and spatial decay rate of connection probability  */ 

int *tab;         /* temporary neighborhood table tab[maxot]*/ 
double *dist, *costheta, *sintheta;      /* table of distances between candidate neurons and costheta, sintheta, updated 08/23/05 */
int *postx, *posty; /* Array of postsynaptic neuron coordinates for accepted target */
double *zclass;  /* Array of class z-values within cortical model */
int *horaxonflagpos,*veraxonflagpos, *axonflag; /* Array of flags for stimulated axon portions */
int *horaxonflagneg,*veraxonflagneg;

/* Neuron Class            Type               Firing Pattern
        1           Layer II/III Pyramid             RS
        2           Layer IV Spiny                   RS
        3           Layer V Pyramid                  IB
        4           Layer VI Pyramid                 IB
        5           Basket Layers II-VI              FS
        6           Double Bouquet II-VI             FS
        7           Chandelier                       FS     */

#define P23 1
#define P4  2
#define P5  3
#define P6  4
#define BASKET 5
#define BOUQUET 6
#define CHANDELIER 7
#define BASKETANG 0.0 //1.5708  /* for offset costheta distr for basket cell axons */
#define PIIIANG 0.0   /* for offset costheta distr for PIII cell axons */
#define PVANG 0.0 /*for offset costheta distr for PV cell axons */
#define PVIANG 0.0    /* for offset costheta distr for PVI cell axons */

#define SPATIALSCALE 25 /* center to center column spacing in microns */

#define BASKETLAMBDA 50 /* 1/e decay for basket cell axon connections */
#define BASKETNEAR 50 /* distance in microns within which basket conn'ions are isotropic */

#define PIIINEARRAD 300   /* Near connection radius for PIII neurons */
#define STRIPWIDTH  500   /* 500 micron strip width for domain simulations */
#define PIIIFARAVE 2000   /* Ave dist to long range connection spot */
#define PIIISPOTRAD 200   /* Long range spot radius */

#define BOUNDARYLAYER 300 /* Boundary layer distance, to reduce connections */

#define ISLENGTH    30    /* One half initial segment length (microns) */

/* Definitions for E-field effects */

#define SEG 500  /* spacing btwn nodes of Ranvier for 5.7 micron myelinated axon */
//#define VOLT .2 /* stimulation voltage applied to metal electrode in voltage controlled mode */
#define I0 .010 /* Stimulation current in current controlled mode in Amps */
#define SIGMA 0.3 /* Tissue conductivity in A/Vm */
#define THRESHD2V 2 /* Threshold for activation fct (mV), second difference of voltage along fiber element */ 
#define A 1000 /* Electrode disk radius (microns) */
#define STIMPOS 2340 /* Electrode Z-value (microns), taken with z=0 at the top of the white matter */
#define DISKX 32 /* Center of electrode disk x-value (lattice constants) */
#define DISKY 32 /* Center of electrode disk y-value (lattice constants) */ 
#define PI 3.14159

/* Cortical layer z-values */
#define ZVIB 162 /* VIbeta */
#define ZVIA 488 /* VIalpha */
#define ZV   1130 /* V */
#define ZIVCBETA 1230 /* IVCbeta */
#define ZIVCALPHA 1230 /* IVCalpha */
#define ZIVB 1230 /* IVB */
#define ZIVA 1230 /* IVA */
#define ZIII 1540 /* III */
#define ZII  1940 /* II */
#define ZI   2185 /* I */

#define Gi   7.6e-07 /* Intracompartment conductivity in inverse Ohms */ 

struct SYN_TMP
{
NEURON_INDEX    neuron_pre;
NEURON_INDEX    neuron_post;
NEURON_INDEX    preclass_num;
WEIGHT          weight;
DELAY           delay;
VAR		axdelay;
};

int sort_delay(const void *, const void *);
int sort_neu_num(const void *, const void *);
void set_neu(struct NEURON *);
unsigned seed_save();

int main(int argc, char  *argv[])
{
 int i,j,l,jp;
 int nr,ir,ix,iy,iz;  /* index for neurons*/
 int maxot, tmpot;
 int xyz[3]={0,0,0};
 int no,nl;
 ulong ns,is;
 int *neu_pre;
 int Nneu;
 int torus_flag;
 int syn_n;
 char label[]="xyz",name[10];
 FILE *inputfile,*outputfile,*graphfile,*axflagfile;
 struct SYN_TMP *syn_tmp;
 struct SYNAPSE syn;
 struct NEURON *neu;
 int otoczenie(int,int,int*,int);
 int otoczenie2(int,int,int*,int);
 int losuj(int*,int);
 int losuj2(int*,double*,double*,double*,int,double,double,int,int,int*);
 extern int zapisz(int,int,int,int,FILE *);
 extern int zapisz2(int,int,int,int,double,FILE *,FILE *,FILE *,int);
 extern int sort_neu_num(const void *, const void *);
 int nsub;
 int pre[3];
 double horzmid,horxmid,horymid,horrhomid,horVmid,horrhopre,horVpre,horrhopost,horVpost;
 double horseglength,hordelsquaredV;
 double verzmid,verxmid,verymid,verrhomid,verVmid,verrhopre,verVpre,verrhopost,verVpost;
 double verseglength,verdelsquaredV,aradius,coeff;
 int synapsesum;
 double currentinjsum, avecurrentinj;
 int nps,s_j,s_i,critdex;

/*

1.Reading and checking input file information

*/

/* printf("Made it to file checking\n"); */

 if(argc<3)
   {
     fprintf(stderr,"USAGE:\n\tmk_link <file_name> <sqrt(no of subnets in all)> [<seed>]\n");
     return 1;
   }
 if(argc>3)
   srandom((unsigned)atoi(argv[3]));
  else
   printf("seed %u\n",seed_save()); /* seed for srand */
 srand48(127L);
 //inputfile = fopen(strcpy(name,argv[1]),"r");
 inputfile = fopen(strcat(strcpy(name,argv[1]),".net"),"r");
 if(inputfile==NULL){fprintf(stderr, "Network definition file not present %s\n",name);return -1;}

 fscanf(inputfile,"%d",&torus_flag);
 fscanf(inputfile,"%d",&lrodz);

 nsub=atoi(argv[2]);
                                       /* Memory allocations  */
 N=(int**)malloc(lrodz*sizeof(int*));        /* table allocation for size of networks for each class */
  for(i=0;i<lrodz;i++)
    N[i]=(int*)malloc(3*sizeof(int));
 k=(int**)malloc(lrodz*3*sizeof(int));        /* table allocation for network scaling factor for each class */
  for(i=0;i<lrodz;i++)
    k[i]=(int*)malloc(3*sizeof(int)); 

 o=(int ***)malloc(lrodz*sizeof(int **));  /* allocating neighborhood tables  */
   for(i=0;i<lrodz;i++)
     o[i]=(int **)malloc(lrodz*sizeof(int *));
    for(i=0;i<lrodz;i++)  
     for(j=0;j<lrodz;j++)
        o[i][j]=(int *)malloc(3*sizeof(int));

 s=(int **)malloc(lrodz*sizeof(int*));      /* synaptic table allocation */
   for(i=0;i<lrodz;i++)
    s[i]=(int*)malloc(lrodz*sizeof(int));

 w=(int **)malloc(lrodz*sizeof(int*));           /* allocation for synaptic weights */
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

  Ap=(double **)malloc(lrodz*sizeof(double));      /* used to compute probability of connection */
   for(i=0;i<lrodz;i++)
    Ap[i]=(double*)malloc(lrodz*sizeof(double));

  lambda=(double **)malloc(lrodz*sizeof(double));  /* used to commpute probability of connection */
   for(i=0;i<lrodz;i++)
    lambda[i]=(double*)malloc(lrodz*sizeof(double));

  zclass=(double*)malloc(lrodz*sizeof(double));

/* Load cortical layer z-information for E-field calc, array index by class number */
   zclass[0]= ZII;
   zclass[1]= ZIVB;   
   zclass[2]= ZV;  
   zclass[3]= ZVIA;
   zclass[4]= ZII;  
   zclass[5]= ZII;  
   zclass[6]= ZIVB;
 
/* initialize mapping and scaling arrays */

for(i=0;i<lrodz;i++)
{
   N[i][0]=N[i][1]=nsub; //square mapping
   N[i][2]=1;
   k[i][0]=k[i][1]=k[i][2]=1;
}


/* input parameters from file to fill arrays */


 for(i=0;i<lrodz;i++)
  for(j=0;j<lrodz;j++)
   { 
     fscanf(inputfile,"%d,%d,%d",&o[j][i][0],&o[j][i][1],&o[j][i][2]);
     fscanf(inputfile,"%d",&s[j][i]);
     fscanf(inputfile,"%d",&w[j][i]);
     fscanf(inputfile,"%d",&sdw[j][i]);
     fscanf(inputfile,"%d",&d[j][i]);
     fscanf(inputfile,"%d",&sdd[j][i]);
     fscanf(inputfile,"%lf",&Ap[j][i]);    
     fscanf(inputfile,"%lf",&lambda[j][i]); 
   }
 fclose(inputfile);
 

/*

2. Offset calculations and table filling

*/

 /* printf("Made it to offset calculations\n"); */

 NN=(int*)malloc(lrodz*sizeof(long int));   /* table with number of subnetworks containing each class */
 for(i=0;i<lrodz;i++)
  NN[i]=(N[i][0]?N[i][0]:1)*(N[i][1]?N[i][1]:1)*(N[i][2]?N[i][2]:1);
 
 NOF=(int*)malloc(lrodz*sizeof(int)); /* table with offsets */
  for(i=1,NOF[0]=0;i<lrodz;i++)
    NOF[i]=NOF[i-1]+NN[i-1];


/* maximal neigborhood calculation */

 maxot=0;
 if(torus_flag!=2) /* i.e. not choosing from whole array */
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
 dist=(double*)malloc(maxot*sizeof(double));
 costheta=(double*)malloc(maxot*sizeof(double));
 sintheta=(double*)malloc(maxot*sizeof(double));
 postx=(int*)malloc(maxot*sizeof(int));
 posty=(int*)malloc(maxot*sizeof(int));
 horaxonflagpos=(int*)malloc(maxot*sizeof(int));
 veraxonflagpos=(int*)malloc(maxot*sizeof(int));
 horaxonflagneg=(int*)malloc(maxot*sizeof(int));
 veraxonflagneg=(int*)malloc(maxot*sizeof(int));
 axonflag=(int*)malloc(maxot*sizeof(int));

 Nneu=0;
 for(i=0;i<lrodz;i++)
  Nneu+=NN[i];   /* calculate subnetworks X # of classes */

 for(i=0;i<lrodz;i++)
  printf("\nNN[%d]=%d, NOF[%d]=%d",i,NN[i],i,NOF[i]); /* co by sprawdzic */
  printf("\nmaxot=%d,Nneu=%d\n",maxot,Nneu);
 
/*

3.Setting up for output of links

*/
 
 outputfile=fopen("syn.tmp","wb");
 graphfile=fopen("syn.graph","wb");
 axflagfile=fopen("syn.axflag","wb");

 if(outputfile==NULL){fprintf(stderr,"\nError in opening synapse output file");exit(1);}
 if(graphfile==NULL){fprintf(stderr,"\nError in opening synapse output file");exit(1);}
 if(axflagfile==NULL){fprintf(stderr,"\nError in opening axflagfile file");exit(1);}


    /*randomize();*/
     nr=0;                           /* subnetwork "neuron" number */
     ns=0;                           /* synapse index */
     currentinjsum=0;                /* in microamps */
     synapsesum=0;
     for(ir=0;ir<lrodz;ir++)         /* loop on classes of pre-neurons */
     {
	 //printf("Working on Preclass %d\n",ir+1);
        for (ix=0,xyz[0]=0;ix<N[ir][0];ix++,xyz[0]+=k[ir][0])     /* x-loop */
	{
	    //printf("ix= %d\n",ix+1);
           for(iy=0,xyz[1]=0;iy<N[ir][1];iy++,xyz[1]+=k[ir][1])    /* y-loop */
              for(iz=0,xyz[2]=0;iz<N[ir][2];iz++,xyz[2]+=k[ir][2],nr++)  /* z-loop  */
                  for(jp=0;jp<lrodz;jp++)     /* loop on classes of post-neurons */
                  {
                     /*     generate tables with numbers of neurons of each class jp,in the neighborhood of neuron nr (neuron from class ir
			    connects to jp (used to be the other way around) */
                     switch(torus_flag)
                     {
                        case 0:                           /* no torus boundary condition */
                           no=otoczenie(nr,ir,xyz,jp);
                           break;
                        case 1:                           /* torus boundary */
                           no=otoczenie2(nr,ir,xyz,jp);  
                           break;
                        case 2:                           /* case 2 choosing from whole array */
                           for(i=NOF[jp],no=0;i<NN[jp]+NOF[jp];i++)
                              if(i!=nr)
                                 tab[no++]=i;
                           break;
                        default:
                           fprintf(stderr,"wrong option for torus flag %d\n",torus_flag);
                           exit(0);
                      }

                     if(!no)continue;

                     l=(no<s[jp][ir])?no:s[jp][ir];
                     /* if no is greater than the number of allowed synapses, then set l=# of allowed synapses  */
		     //printf("Allowed/required conn comp no=%d,preclass=%d,postclass=%d,s[jp][ir]=%d\n",no,ir+1,jp+1,s[jp][ir]);/* network linkage checks */


		     /* Boundary effect correction, reducing allowed connections for neurons in specific column number ranges, and row number ranges */

		     nps=nsub; //sqrt(subnets/node)*sqrt(nodes)
                     s_j=(nr)%nps; /* column index of presynaptic subnetwork */
                     s_i=(nr)/nps-ir*nps; /* row index of presynaptic subnetwork */
		     critdex=BOUNDARYLAYER/SPATIALSCALE;

		     //printf("Subnet indices: column=%d row=%d nps=%d critdex=%d nr=%d\n",s_j,s_i,nps,critdex,nr);

		     //printf("Pre-l: l=%d\n",l);

		     if ( (s_j<critdex) && (s_i>=critdex) && (s_i<=((nps-1)-critdex)) ){
		          l=(l*5)/8;
		     }	 

		     if ( (s_j>((nps-1)-critdex)) && (s_i>=critdex) && (s_i<=((nps-1)-critdex)) ){
		          l=(l*5)/8;
		     }	 

		     if ( (s_i<critdex) && (s_j>=critdex) && (s_j<=((nps-1)-critdex)) ){
		          l=(l*5)/8;
		     }	 

		     if ( (s_i>((nps-1)-critdex)) && (s_j>=critdex) && (s_j<=((nps-1)-critdex)) ){
		          l=(l*5)/8;
		     }	 

		     if ( (s_j<critdex) && (s_i<critdex) ){
		          l=(l*3)/8;
		     }	 

		     if ( (s_j<critdex) && (s_i>((nps-1)-critdex)) ){
		          l=(l*3)/8;
		     }	 

		     if ( (s_j>((nps-1)-critdex)) && (s_i<critdex) ){
		          l=(l*3)/8;
		     }	 

		     if ( (s_j>((nps-1)-critdex)) && (s_i>((nps-1)-critdex)) ){
		          l=(l*3)/8;
		     } 

		     //printf("Post-l: l=%d\n",l);

		     /* Form connections and determine stimulation locations */

                     for(i=0;i<l;i++)    /* losuje z calego otoczenia - nowa wersja */
                     {                    /* z losuj2  01/29/03 */ 
                        
                        for(j=0;j<3;j++)
			     pre[j]=xyz[j]/k[jp][j];  

			// printf("no=%d,preclass=%d,postclass=%d,l=%d,i=%d\n",no,ir+1,jp+1,l,i);  /* network linkage checks */
			nl=-1;
			while (nl==-1){
                             nl=losuj2(tab,dist,costheta,sintheta,no,Ap[jp][ir],lambda[jp][ir],jp,ir,pre);
			}
			//printf("i= %d,l= %d,nl= %d\n",i,l,nl);
			if(nl!=-1)
                        {

	                /* Set synapse flag if gradient requirement met */
                        /* Check E field gradient requirement */
			/* Charged balance, bipolar pulse */
   
                        /* Horizontal axon segment */
			aradius=A;
			aradius=aradius/(1000*1000);
			coeff=I0/(4*PI*SIGMA*aradius);
			//printf("\n aradius=%f, coeff=%f",aradius,coeff);


			axonflag[nl]=0;
		
                        horaxonflagpos[nl]=0;
			horaxonflagneg[nl]=0;
			horzmid=zclass[jp];
                        horxmid=(pre[0]+postx[nl])/2;
			horymid=(pre[1]+posty[nl])/2;
			horrhomid=sqrt(pow((DISKX-horxmid),2)+pow((DISKY-horymid),2))*25;
			horVmid=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-horzmid),2)+pow((horrhomid-A),2))+sqrt(pow((STIMPOS-horzmid),2)+pow((horrhomid+A),2))));
			horrhopre=sqrt(pow((DISKX-pre[0]),2)+pow((DISKY-pre[1]),2))*25;
			horVpre=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-zclass[jp]),2)+pow((horrhopre-A),2))+sqrt(pow((STIMPOS-zclass[jp]),2)+pow((horrhopre+A),2))));
			horrhopost=sqrt(pow((DISKX-postx[nl]),2)+pow((DISKY-posty[nl]),2))*25;
			horVpost=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-zclass[jp]),2)+pow((horrhopost-A),2))+sqrt(pow((STIMPOS-zclass[jp]),2)+pow((horrhopost+A),2))));
			horseglength=sqrt(pow((postx[nl]-pre[0]),2)+pow((posty[nl]-pre[1]),2))*25;
			/* hordelsquaredV in mV */ 
			hordelsquaredV = (horVpre+horVpost-2*horVmid)*1000;
                        if (hordelsquaredV > THRESHD2V)
			    horaxonflagpos[nl]=1;
			if (-hordelsquaredV > THRESHD2V)
			    horaxonflagneg[nl]=1;

                        // Vertical axon segment  (Old Version)
                        //veraxonflagpos[nl]=0;
			//veraxonflagneg[nl]=0;
			//verzmid=(zclass[ir]+zclass[jp])/2;
                        //verxmid=pre[0];
			//verymid=pre[1];
			//verrhomid=sqrt(pow((DISKX-verxmid),2)+pow((DISKY-verymid),2))*25;
			//verVmid=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-verzmid),2)+pow((verrhomid-A),2))+sqrt(pow((STIMPOS-verzmid),2)+pow((verrhomid+A),2))));
			//verrhopre=sqrt(pow((DISKX-pre[0]),2)+pow((DISKY-pre[1]),2))*25;
			//verVpre=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-zclass[ir]),2)+pow((verrhopre-A),2))+sqrt(pow((STIMPOS-zclass[ir]),2)+pow((verrhopre+A),2))));
			//verrhopost=sqrt(pow((DISKX-pre[0]),2)+pow((DISKY-pre[1]),2))*25;
			//verVpost=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-zclass[jp]),2)+pow((verrhopost-A),2))+sqrt(pow((STIMPOS-zclass[jp]),2)+pow((verrhopost+A),2))));
			//verseglength=sqrt(pow((zclass[jp]-zclass[ir]),2));
			///* verdelsquaredV in mV */
			//verdelsquaredV = (verVpre+verVpost-2*verVmid)*1000;
                        //if (verdelsquaredV > THRESHD2V)
			//    veraxonflagpos[nl]=1;
			//if (-verdelsquaredV > THRESHD2V)
			//    veraxonflagneg[nl]=1;

                        // Vertical axon segment  (IS Version)
                        veraxonflagpos[nl]=0;
			veraxonflagneg[nl]=0;
			verzmid=(2*zclass[ir]-2*ISLENGTH)/2;
                        verxmid=pre[0];
			verymid=pre[1];
			verrhomid=sqrt(pow((DISKX-verxmid),2)+pow((DISKY-verymid),2))*25;
			verVmid=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-verzmid),2)+pow((verrhomid-A),2))+sqrt(pow((STIMPOS-verzmid),2)+pow((verrhomid+A),2))));
			verrhopre=sqrt(pow((DISKX-pre[0]),2)+pow((DISKY-pre[1]),2))*25;
			verVpre=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-zclass[ir]),2)+pow((verrhopre-A),2))+sqrt(pow((STIMPOS-zclass[ir]),2)+pow((verrhopre+A),2))));
			verrhopost=sqrt(pow((DISKX-pre[0]),2)+pow((DISKY-pre[1]),2))*25;
			verVpost=(coeff)*asin(2*A/(sqrt(pow((STIMPOS-zclass[ir]+2*ISLENGTH),2)+pow((verrhopost-A),2))+sqrt(pow((STIMPOS-zclass[ir]+2*ISLENGTH),2)+pow((verrhopost+A),2))));
			verseglength=2*ISLENGTH; //sqrt(pow((zclass[jp]-zclass[ir]),2));
			/* verdelsquaredV in mV */
			verdelsquaredV = (verVpre+verVpost-2*verVmid)*1000;
                        if (verdelsquaredV > THRESHD2V)
			    veraxonflagpos[nl]=1;
			if (-verdelsquaredV > THRESHD2V)
			    veraxonflagneg[nl]=1;


			/* if (-verdelsquaredV > THRESHD2V){
			     currentinjsum+=(-verdelsquaredV)*Gi;
			     synapsesum=synapsesum+1;
			} */

			//if ( (horaxonflagpos[nl]==1) || (veraxonflagpos[nl]==1) )
			if ( veraxonflagpos[nl]==1 )
			    axonflag[nl]=1;

			//if ( (horaxonflagneg[nl]==1) || (veraxonflagneg[nl]==1) )
			if ( veraxonflagneg[nl]==1 )
			    axonflag[nl]=2;

                        ns+=zapisz2(nr,tab[nl],jp,ir,dist[nl],outputfile,graphfile,axflagfile,axonflag[nl]);
			/* printf("\npostindex=%d, postclass=%d, preindex=%d, preclass=%d",nr,ir+1,tab[nl],jp+1);*/  /* network linkage checks */
			// tab[nl]=-1;

			}
			
                      } /* koniec petli losowania */
                  } /* jp */

	}
     }
  fclose(outputfile);
  fclose(graphfile);
  fclose(axflagfile);

  // avecurrentinj=currentinjsum/synapsesum;
  //printf("currentinjsum=%f,avecurrentinj=%f,synapsesum=%i\n",currentinjsum,avecurrentinj,synapsesum);  /* network linkage checks */


printf("\n");
return 0;
} // koniec main()


//sort_delay  - uzywana w qsort() przy porzadkowaniu delay'ow

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
sort_neu_num - uzywana w qsort() 
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
  
4. Szczegoly generacji tablicy otoczenia:

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


int losuj(int *tablica,int index)
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

int losuj2(int *tablica, double *distance, double *costht,double *sintht,int index, double Ap, double lambda, int classpost, int classpre,int *precoord)    
/* zwraca index w tab, a nie numer neuronu jak losuj() */
{
  int i,j,flagspatial;
  //double Ap=1.0;
  //double lambda=6;
  double los,losang,los1,los2,cutoff;
  double p3nearrad,p3farave,p3spotrad;
  double stripsize,preyy,postycoord,prefloor,postfloor;
  double precheck,postcheck;

  // i=randomm(index);

    do
    {
     i=randomm(index);
     // printf("i=%d,tablica[i]=%d,index=%d\n",i,tablica[i],index);  /* network linkage checks */
    }while(tablica[i]<0); 
  
  switch(classpre)
                     {

                        case 0:                           /* Pre - Layer II/III Pyramid */

			   /* Patch Handling */

			   p3nearrad=PIIINEARRAD/SPATIALSCALE;
                           stripsize=STRIPWIDTH/SPATIALSCALE;
			   //p3farave=PIIIFARAVE/SPATIALSCALE;
			   //p3spotrad=PIIISPOTRAD/SPATIALSCALE;
                           los=drand48();
                           // los2= Ap*exp(-1.0*distance[i]/lambda);
			   //losang=drand48();
                           // los2=1;
			   flagspatial=0;
			   if (distance[i]<=p3nearrad)
			       flagspatial=1;
                           else {
			       /* Stripflag A=even, B=odd */ 
                               preyy=precoord[1];
			       postycoord=posty[i];
                               prefloor=floor(preyy/stripsize);
                               postfloor=floor(postycoord/stripsize);
			       precheck=prefloor-2*floor(prefloor/2.0);
			       postcheck=postfloor-2*floor(postfloor/2.0);
                               //printf("precheck=%f,postcheck=%f,preyy=%f\n",precheck,postcheck,preyy); 
			       if ( precheck == postcheck)
			   	   flagspatial=0;
			       else
                                   flagspatial=0;

			       //printf("precheck=%f,postcheck=%f,flagspatial=%d\n",precheck,postcheck,flagspatial); 

                           }
			 
			   /* Angular isotropy routine */
			   //los1=(costht[i]*cos(PIIIANG)+sintht[i]*sin(PIIIANG)); 
			   //los1*=los1; 
			   //los1*=los1;
			   //los1*=los1; 
                           /*printf("p3nearrad=%lf,p3farave=%lf,p3spotrad=%lf,flagspatial=%i,distance=%lf\n",p3nearrad,p3farave,p3spotrad,flagspatial,distance[i]);*/  

                           if( ( flagspatial ) /*&&  (losang < los1)*/ ) 
                           {
                              return i;
                           }
                           else
                           {
			       // tablica[i]=-1; 
                              return -1;
                           }
                           break;

			   //los=drand48();
                           //// los2= Ap*exp(-1.0*distance[i]/lambda);
                           //los2=1;
                           ////printf("index=%i, neuron=%i, distance=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,distance[i],los,los2);  
                           //if( los < los2 )
                           //{
                           //   return i;
                           //}
                           //else
                           //{
			   //    // tablica[i]=-1; 
                           //   return -1;
                           //}  
                           //break;


                        case 1:                           /* Pre - Layer IV Spiny */
                           los=drand48();
                           // los2= Ap*exp(-1.0*distance[i]/lambda);
                           los2=1;
                           //printf("index=%i, neuron=%i, distance=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,distance[i],los,los2);  
                           if( los < los2 )
                           {
                              return i;
                           }
                           else
                           {
			       // tablica[i]=-1; 
                              return -1;
                           }  
                           break;
                        case 2:                           /* Pre - Layer V Pyramid */
                           los=drand48();
                           // los2= Ap*exp(-1.0*distance[i]/lambda);
			   //losang=drand48();
                           los2=1;
			   //los1=(costht[i]*cos(PVANG)+sintht[i]*sin(PVANG));
			   //los1*=los1;
			   //los1*=los1;
			   //los1*=los1;
                           //printf("index=%i, neuron=%i, distance=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,distance[i],los,los2);  
			   if( ( los < los2 ) /*  && (losang < los1) */ )
                           {
                              return i;
                           }
                           else
                           {
			       // tablica[i]=-1; 
                              return -1;
                           }                           
                           break;
                        case 3:                           /* Pre - Layer VI Pyramid */
                           los=drand48();
                           // los2= Ap*exp(-1.0*distance[i]/lambda);
			   //losang=drand48();
                           los2=1;
			   //los1=(costht[i]*cos(PVIANG)+sintht[i]*sin(PVIANG));
			   //los1*=los1;
			   //los1*=los1;
			   //los1*=los1;
                           //printf("index=%i, neuron=%i, distance=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,distance[i],los,los2);  
			   if( ( los < los2 ) /*  && (losang < los1) */ )
                           {
                              return i;
                           }
                           else
                           {
			       // tablica[i]=-1; 
                              return -1;
                           }                          
                           break;
                        case 4:                           /* Pre - Basket Cell */
			   cutoff = BASKETNEAR/SPATIALSCALE;
			   if ( distance[i] <= cutoff ){       /* Isotropic connections to 50 microns */
			       return i;
                           }
			   else { 
                              los=drand48();
			      //losang=drand48();
			      lambda = BASKETLAMBDA/SPATIALSCALE;
                              los2= Ap*exp(-1.0*distance[i]/lambda); 
			      //los1=(costht[i]*cos(BASKETANG)+sintht[i]*sin(BASKETANG));
			      //los1*=los1;
			      //los1*=los1;
			      //los1*=los1;
			      // los2=1;
                              /* printf("index=%i, neuron=%i, distance=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,distance[i],los,los2); */ 
                              if( (los < los2) /* && (losang < los1) */ )
                              {
				  /* printf("index=%i, neuron=%i, costht=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,costht[i],losang,los1);*/
                                 return i;
                              }
                              else
                              {
				  // tablica[i]=-1; 
                                 return -1;
                              }
                           }                           
                           break; 
                        case 5:                           /* Pre - Double Bouquet */
                           los=drand48();
                           // los2= Ap*exp(-1.0*distance[i]/lambda);
                           los2=1;
                           //printf("index=%i, neuron=%i, distance=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,distance[i],los,los2);  
                           if( los < los2 )
                           {
                              return i;
                           }
                           else
                           {
			       // tablica[i]=-1; 
                              return -1;
                           }                           
                           break;  
                        case 6:                           /* Pre - Chandelier Cell  */
                           los=drand48();
                           // los2= Ap*exp(-1.0*distance[i]/lambda);
                           los2=1;
                           //printf("index=%i, neuron=%i, distance=%lf, drand48=%3.2lf, prob= %3.2lf\n",i,j,distance[i],los,los2);  
                           if( los < los2 )
                           {
                              return i;
                           }
                           else
                           {
			       // tablica[i]=-1; 
                              return -1;
                           }
                           break;
                        default:
                           fprintf(stderr,"Problem with neuron class in losuj2 %d\n",classpre);
                           exit(0);
                      }

}

int zapisz(int inr, int inl, int i, int j, FILE * outputfile)
{

struct SYN_TMP syn_tmp; 
extern  int rozrzut(int);
syn_tmp.neuron_post=inl;
syn_tmp.neuron_pre=inr;
syn_tmp.weight=w[i][j]+rozrzut(sdw[i][j]);
syn_tmp.delay =d[i][j]+rozrzut(sdd[i][j]);
//printf("%i %i %i %i\n",syn_tmp.neuron_pre, syn_tmp.neuron_post , syn_tmp.weight,syn_tmp.delay );
return(fwrite(&syn_tmp,sizeof(struct SYN_TMP),1,outputfile));

}

int zapisz2(int inr, int inl, int i, int j, double distance,  FILE * outputfile,  FILE * graphfile, FILE * axflagfile, int axflag)
{
struct SYN_TMP syn_tmp;
extern  int rozrzut(int);
int axonal_del=5; /* was 5 */
syn_tmp.neuron_post=inl;
syn_tmp.neuron_pre=inr;
syn_tmp.preclass_num=j;
syn_tmp.weight=w[i][j]+rozrzut(sdw[i][j]);
syn_tmp.delay =d[i][j]+rozrzut(sdd[i][j]);
//printf("distance=%3.2lf, axonal del=%3.2lf, delay=%i, delay + axonal del=%i\n", distance, distance*axonal_del, syn_tmp.delay, syn_tmp.delay+(int)rint(distance*axonal_del));
syn_tmp.axdelay=(int)rint(distance*axonal_del)+3*rozrzut(sdd[i][j]);
// printf("%i %i %i %i %i\n",syn_tmp.neuron_pre, syn_tmp.neuron_post , syn_tmp.weight, syn_tmp.delay , syn_tmp.axdelay);
if ((j)==4){
 fwrite(&inr,sizeof(ushort),1,graphfile);
 fwrite(&j,sizeof(ushort),1,graphfile);
 fwrite(&inl,sizeof(ushort),1,graphfile);
 fwrite(&i,sizeof(ushort),1,graphfile);
}
/* printf("axflag= %i\n",axflag); */
 fwrite(&axflag,sizeof(ushort),1,axflagfile);

return(fwrite(&syn_tmp,sizeof(struct SYN_TMP),1,outputfile));
}

int otoczenie(int nr ,int ir,int xyz[3],int jp)
 {
   int nx,ny,no,i[3],ix,iy,iz,mini[3],maxi[3], *ot;
   int j;
   /* tab must be declared earlier with sufficient size */

   ot = o[jp][ir];
   /* obliczamy najpierw wspolrzedne i[3] neuronu na siatce jp, ktory jest najblizej x,y,z */ /*???*/
  
  
  no=0;
  for(j=0;j<3;j++)
         i[j]=xyz[j]/k[jp][j];  

  /* teraz otoczenie
     wersja bez zawijania (brzegowe maja mniejsze otoczenie) */ /*???*/

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
            if((j=ny+iz)!=nr){
              dist[no]=sqrt((ix-i[0])*(ix-i[0])+(iy-i[1])*(iy-i[1])+(iz-i[2])*(iz-i[2]));
              costheta[no]=(ix-i[0])/dist[no];
	      sintheta[no]=(iy-i[1])/dist[no];
	      /* printf("ix=%i, iy=%i, i[0]=%i,i[1]=%i,costheta=%lf,sintheta=%lf,dist=%lf\n",ix,iy,i[0],i[1],costheta[no],sintheta[no],dist[no]);*/
              // printf("%6.1lf,",dist[no]); /* tablica z odleglosciami - nowa wersja 01/29/03 */
	      postx[no]=ix;
	      posty[no]=iy;

              tab[no++]=j;}
          }
        }
     }
//printf("\n\n");
return no;
}


int otoczenie2(int nr, int ir,int xyz[3],int jp)
{
   int nx,ny,no,i[3],ix,iy,iz,mini[3],maxi[3], *ot;
   int j,ii,NX,NY,NZ;
   /* tab must be declared earlier with sufficient size */

   ot = o[jp][ir];
   /* obliczamy najpierw wspolrzedne i[3] neuronu na siatce jp, ktory jest najblizej x,y,z */ /* ??? */
  
  
  no=0;
  for(j=0;j<3;j++)
         i[j]=xyz[j]/k[jp][j];  

  /* Torus case */

      
 for(j=0;j<3;j++)
 {
   mini[j]=i[j]-ot[j];
   maxi[j]=i[j]+ot[j]+1;
 }
   NX=N[jp][0];
   NY=N[jp][1];
   NZ=N[jp][2];
 
   for(ix=mini[0];ix<maxi[0];ix++)
    { ii=(ix+NX)%NX;
      nx=ii*NZ*NY+NOF[jp];
      for(iy=mini[1];iy<maxi[1];iy++)
       { ii=(iy+NY)%NY;
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
   neu->flags|=HISTO;
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
//   neu->C=C_0;
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


