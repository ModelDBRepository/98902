/* Procedure generates text link file for subnetwork connections
   input from syn.tmp (pre- and post- numbers mean subnetworks)
   output in <name>.link file for crl 
   P.J. Franaszcuk JHU ver1. Jan 2002

   05/29/03 P. Kudela, now multi-class supported
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "lnet.h"

#define randomm(num) (int)((double)num*(double)random()/((double)(RAND_MAX+1.0)))

struct SYN_TMP
{
NEURON_INDEX    neuron_pre;
NEURON_INDEX    neuron_post;
NEURON_INDEX    preclass_num;
WEIGHT          weight;
DELAY           delay;
VAR		axdelay;
};

struct SYN_TMPMKNODE
{
NEURON_INDEX     neuron_pre;
NEURON_INDEX     neuron_post;
NEURON_INDEX     preclass_num;
WEIGHT           weight;
DELAY            delay;
};

int sort_neu_num(const void *, const void *);
int sort_delay(const void *, const void *);

int main(int argc, char  *argv[])
{
 #define MAXTIME 500
 char name[100],fname[100];
 int ns,np; /* sqrt of :(no_of subnets in node, no_of nodes)*/
 int m; /* number of connections between subnets */
 int i,s_i,s_j,d_i,d_j,nneu,nneui,nps,nns,nis,d_sub,s_sub,is,id,lrodz;
 int predex, postdex,syncounter,ifile,weirdcount;
 int pospolstimcounter, negpolstimcounter;
 int preclassnum;
 int syntime[MAXTIME],linktime[MAXTIME];
 int possynpospol, negsynpospol, possynnegpol, negsynnegpol;
 int posspikecount, negspikecount;
 int **syntotal,j,synave;
 FILE *inp,*out,*nodefile,*inpax,**neufile,**neuwritefile,*cfg;
 FILE **synfile,**tmpsyn;
 FILE *synt, *linkt;
 struct SYN_TMP syn;
 struct SYN_TMPMKNODE synintranode,syndummy,*syn_tmp;
 struct SYNAPSE synwrite;
 struct NEURON neu;
 int src_neu,dst_neu,src_proc,dst_proc,w,nr,classtest;
 unsigned short axstimflag;
 double del,delay,delintrasyn;
 int  idel;
 unsigned int seed, bytetest;
 int *N,*NOF;                /* NOF[lrodz] table with neuron index offsets */
 int *nsyn;
 int *neu_pre;
 int Nneu=0,syn_n,oldsyn_num;
 int histoint;
 extern int sort_neu_num(const void *, const void *);
 
 if(argc<5)
   {
     fprintf(stderr,"\nUSAGE:\n\tsublink <name> <m> <del> <seed>\n"
	            "\tname - output name (no. ext)\n"
	            "\tm  - no of connections between subnets\n"
	            "\tdel - axonal delay factor (float)\n"
	            "\tseed - seed for random generator\n");
    	             return 3;
   }

 if((nodefile=fopen(strcat(strcpy(name,argv[1]),".node"),"r"))==NULL)
  {fprintf(stderr, "Can't open %s file\n",name);return 1;}
 fscanf(nodefile,"%d",&lrodz);
 N=(int*)malloc(lrodz*sizeof(int)); 
 for(i=0;i<lrodz;i++)
   fscanf(nodefile,"%d",&N[i]);
 fscanf(nodefile,"%d",&ns);
 fscanf(nodefile,"%d",&np); 
 fclose(nodefile);

 syntotal=(int**)malloc(lrodz*sizeof(int*));        /* syntotal[pre][post] total of synapses on post class coming from pre class */
  for(i=0;i<lrodz;i++)
    syntotal[i]=(int*)malloc(lrodz*sizeof(int));

  for(i=0;i<lrodz;i++)
    Nneu+=N[i]*N[i];
  Nneu=Nneu*ns*ns;
  histoint=10;
  //printf("Nneu= %i\n",Nneu);

 for(i=0;i<lrodz;i++)
     for(j=0;j<lrodz;j++){
          syntotal[i][j]=0;;
     }
 /*Initialize timing histograms */
 for(i=0;i<MAXTIME;i++){
      linktime[i]=0;
      syntime[i]=0;
 }

 neufile=(FILE **)malloc(np*np*sizeof(FILE*));
 neuwritefile=(FILE **)malloc(np*np*sizeof(FILE*));
 tmpsyn=(FILE **)malloc(np*np*sizeof(FILE*));
 synfile=(FILE **)malloc(np*np*sizeof(FILE*));

 nsyn=(int*)malloc(np*np*sizeof(int));  /* table of neuron index offsets */

 NOF=(int*)malloc(lrodz*sizeof(int));  /* table of neuron index offsets */
  for(i=1,NOF[0]=0;i<lrodz;i++)
    NOF[i]=NOF[i-1]+N[i-1]*ns*N[i-1]*ns;

 inp=fopen("syn.tmp","r");
 if(!inp){fprintf(stderr, "Can't open syn.tmp file\n");return 1;}

 inpax=fopen("syn.axflag","r");
 if(!inpax){fprintf(stderr, "Can't open syn.axflag file\n");return 1;}

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      sprintf(name,"tmp.neutmp.%d",ifile);
      neufile[ifile]=fopen(name,"r+");
      if(!neufile[ifile]){fprintf(stderr, "Can't open tmp.neutmp.%d file\n",ifile);return 1;}
 }

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      sprintf(name,"tmpsyn.%d",ifile);
      tmpsyn[ifile]=fopen(name,"a+");
      if(!tmpsyn[ifile]){fprintf(stderr, "Can't open tmpsyn.%d file\n",ifile);return 1;}
 }

 m=atoi(argv[2]);   //no of conn  per sub
 delay=atof(argv[3]); //del
 seed=atoi(argv[4]); //seed
 nps=ns*np;  // sqrt(subnets/node)*sqrt(nodes)
 
 out=fopen(strcat(strcpy(name,argv[1]),".link"),"w");
 if(!out){fprintf(stderr, "Can't open link file %s\n",name);return 1;}

 syncounter=0; 
 weirdcount=0;
 pospolstimcounter=0;
 negpolstimcounter=0;
 possynpospol=0;
 negsynpospol=0;
 possynnegpol=0;
 negsynnegpol=0;
 posspikecount=0;
 negspikecount=0;

 while(1)
   {
     if(fread(&syn,sizeof(struct SYN_TMP),1,inp)!=1)
        if(feof(inp)) break;
     else
        { fprintf(stderr,"fread error in syn.tmp\n");return 2;}

     syncounter=syncounter+1;
     /* printf("syncounter= %3i\n",syncounter); */

     if(fread(&axstimflag,sizeof(ushort),1,inpax)!=1){
	 
	 if (feof(inpax))
              fprintf(stderr,"reached eof in syn.axflag\n");
         else
	      { fprintf(stderr,"fread error in syn.axflag\n");return 2;}
     }    
     
     predex=syn.neuron_pre;
     s_j=syn.neuron_pre%nps; /* column index of subnet within syn structure neuron_pre element */
     s_i=syn.neuron_pre/nps; /* row index of subnet within syn structure neuron_pre element */
        
     //printf("Subnet indices: column=%d row=%d nps=%d predex=%d\n",s_j,s_i,nps,predex);

     preclassnum=syn.preclass_num; /* Preclass number taken directly from stucture */

     is=s_i/nps; /* Calculated preclass number */

     /* printf("Stuct preclass= %i Calculated preclass= %i\n",is,preclassnum); */

     src_proc=((s_i%nps)/ns)*np+s_j/ns; /* node no */

     s_sub= NOF[is]+(((s_i%nps)%ns)*N[is]*ns+(s_j%ns)*N[is])*N[is];

     postdex=syn.neuron_post;
     d_j=syn.neuron_post%nps;
     d_i=syn.neuron_post/nps;
 
     id=d_i/nps;  /* klasa neuronu */
     dst_proc=((d_i%nps)/ns)*np+d_j/ns;
     //printf("Destination proc # %3i out of %3i\n",dst_proc,(np*np-1));
     d_sub= NOF[id]+(((d_i%nps)%ns)*N[id]*ns+(d_j%ns)*N[id])*N[id];

//     if(! ((s_j == d_j) && ((s_i%nps) == (d_i%nps))) )  /* co by nie tworzyc linkow in sub */
//     {
     del=syn.delay*1000.*DELTA_T+delay*syn.axdelay;
     delintrasyn=rint(syn.delay+delay*syn.axdelay/(1000*DELTA_T));
     //delintrasyn=syn.delay+delay*syn.axdelay;
     w=syn.weight;

     syntotal[is][id]=syntotal[is][id]+1;

     for(i=0;i<m;i++)
       {
         nr=randomm(N[is]*N[is]);
         src_neu=s_sub+nr;
         nr=randomm(N[id]*N[id]);
         dst_neu=d_sub+nr;

         if( (s_j == d_j) && ((s_i%nps) == (d_i%nps)) )
	    weirdcount=weirdcount+1;
	 /* printf("%i:%i %i:%i %lf %i\n",src_proc,src_neu,dst_proc,dst_neu,del,w); */ 
             
         else {
         /* if( (s_j != d_j) && ((s_i%nps) != (d_i%nps)) ){ */

             if (axstimflag==1){
		  pospolstimcounter=pospolstimcounter+1;

		  if(w>0){
		       possynpospol=possynpospol+1;}
		  else{
		       negsynpospol=negsynpospol+1;}

		  /* printf("axstimflag= %i\n",axstimflag); */     

	          fseek(neufile[src_proc],(src_neu)*sizeof(struct NEURON),SEEK_SET);	
                  if(fread(&neu,sizeof(struct NEURON),1,neufile[src_proc])!=1)
                       if(feof(inp)) break;
                  else
                       { fprintf(stderr,"fread error in tmp.neutmp.%d\n",src_proc);return 2; } 
		  /* printf("%3i:%3i %3i:%3i\n",is,predex,id,postdex);*/
                  bytetest=neu.flags;
		  /* printf("old flagbyte= %3i\n",bytetest); */


		  // if( !(neu.flags & IEXT) && !(neu.polarity & POSPOL)  ){
		  if( !(neu.polarity & POSPOL)  ){

		       posspikecount=posspikecount+1;

                       neu.flags=neu.flags|=IEXT;
		       neu.polarity=neu.polarity|=POSPOL;
                       bytetest=neu.flags;
		       /* printf("new flagbyte= %3i\n",bytetest); */ 
	               fseek(neufile[src_proc],(src_neu)*sizeof(struct NEURON),SEEK_SET);	
                       fwrite(&neu,sizeof(struct NEURON),1,neufile[src_proc]);

		       /* printf("file position= %d\n",ftell(neufile[dst_proc])); */

		       fflush(neufile[src_proc]);

		       /* fseek(neufile[dst_proc],(dst_neu)*sizeof(struct NEURON),SEEK_SET);	
                       if(fread(&neu,sizeof(struct NEURON),1,neufile[dst_proc])!=1)
                            if(feof(inp)) break;
		            else {
			         fprintf(stderr, "sizeof struct neu= %d\n",sizeof(struct NEURON));
                                 fprintf(stderr,"bad read file position= %d\n",ftell(neufile[dst_proc]));
			         fprintf(stderr,"%i:%i %i:dst_neu=%i %lf %i d_sub=%i\n",src_proc,src_neu,dst_proc,dst_neu,del,w,d_sub);
                                 fprintf(stderr,"fread error in tmp.neu.%d %3i:%3i %3i:%3i %lf %3i\n",dst_proc,is,predex,id,postdex,del,w);return 2;
                            } 
                       bytetest=neu.flags;
		       printf("newnew flagbyte= %3i\n",bytetest); */

                  }


             } 

	     if (axstimflag==2){
		  negpolstimcounter=negpolstimcounter+1;

		  if(w>0){
		       possynnegpol=possynnegpol+1;}
		  else{
		       negsynnegpol=negsynnegpol+1;}

		  fseek(neufile[src_proc],(src_neu)*sizeof(struct NEURON),SEEK_SET);	
                  if(fread(&neu,sizeof(struct NEURON),1,neufile[src_proc])!=1)
                       if(feof(inp)) break;
                  else
                       { fprintf(stderr,"fread error in tmp.neutmp.%d\n",src_proc);return 2; } 

		  // if( !(neu.flags & IEXT) && !(neu.polarity & NEGPOL)  ){
		  if( !(neu.polarity & NEGPOL)  ){

		       negspikecount=negspikecount+1;

                       neu.flags=neu.flags|=IEXT;
		       neu.polarity=neu.polarity|=NEGPOL;
                
	               fseek(neufile[src_proc],(src_neu)*sizeof(struct NEURON),SEEK_SET);	
                       fwrite(&neu,sizeof(struct NEURON),1,neufile[src_proc]);

		       fflush(neufile[src_proc]);

                  }

             } 
             
	     if ( (src_proc != dst_proc)){
	     //if ( (src_proc != dst_proc) || (src_proc == dst_proc) ){
                  fprintf(out,"%i:%i %i:%i %i %lf %i\n",src_proc,src_neu,dst_proc,dst_neu,preclassnum,del,w);
		  idel=del;
		  //printf("delay=%lf idel=%i\n",del,idel);
		  linktime[idel]+=1;
	     }
	     else {	
	     // Append (src_proc==dst_proc) Synaptic info to temporary synapse file generated in mknode.c

	          synintranode.neuron_post=dst_neu;
                  synintranode.neuron_pre=src_neu; 
                  synintranode.preclass_num=preclassnum;      
                  synintranode.weight=w;
                  synintranode.delay =delintrasyn;
                  //idel=delintrasyn*1000*DELTA_T;
	          //printf("delintrasyn=%i idel=%i\n",delintrasyn,idel);
                  fwrite(&synintranode,sizeof(struct SYN_TMPMKNODE),1,tmpsyn[src_proc]);
                  //printf("pre= %d,post= %d,preclass= %d, weight= %d,delay= %d, del= %f\n",synintranode.neuron_pre,synintranode.neuron_post,synintranode.preclass_num,synintranode.weight,synintranode.delay, del);

             }

              /* printf("%3i:%3i %3i:%3i %lf %3i\n",is,predex,id,postdex,del,w); */
         }
       }
    } 


 // printf("weirdcount= %3i\n",weirdcount);
 fprintf(stdout,"syncounter=%i \n",syncounter); 
 fprintf(stdout,"pospolstimcounter=%i possynpospol=%i negsynpospol=%i\n",pospolstimcounter,possynpospol,negsynpospol);
 fprintf(stdout,"negpolstimcounter=%i possynnegpol=%i negsynnegpol=%i\n",negpolstimcounter,possynnegpol,negsynnegpol);
 fprintf(stdout,"posspikecount=%i negspikecount=%i\n",posspikecount,negspikecount); 

 for(i=0;i<lrodz;i++)
     for(j=0;j<lrodz;j++){
	  synave = syntotal[i][j]/(N[j]*N[j]*nps*nps);
          fprintf(stdout,"Pre Class=%i Post Class=%i Total Synapses=%i Synapses/PreCell=%i\n",i,j,syntotal[i][j],synave);
     }

 fclose(inp);
 fclose(out);
 fclose(inpax);

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      fclose(tmpsyn[ifile]);
 }

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      fclose(neufile[ifile]);
 }


//Final InterNodal Synapse Sorting, and Synapse File Writing

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      sprintf(name,"tmpsyn.%d",ifile);
      tmpsyn[ifile]=fopen(name,"r+");
      if(!tmpsyn[ifile]){fprintf(stderr, "Can't open tmpsyn.%d file\n",ifile);return 1;}
 }

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      sprintf(name,"%s.syn.%d",argv[1],ifile);
      synfile[ifile]=fopen(name,"w");
      if(!synfile[ifile]){fprintf(stderr, "Can't open %s.synfile.%d file\n",argv[1],ifile);return 1;}
 }

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      sprintf(name,"tmp.neutmp.%d",ifile);
      neufile[ifile]=fopen(name,"r+");
      if(!neufile[ifile]){fprintf(stderr, "Can't open tmp.neutmp.%d file\n",ifile);return 1;}
 }

 for (ifile=0;ifile<(np*np);ifile++){ /* Loop over processors */
      sprintf(name,"tmp.neu.%d",ifile);
      neuwritefile[ifile]=fopen(name,"w+");
      if(!neuwritefile[ifile]){fprintf(stderr, "Can't open tmp.neu.%d file\n",ifile);return 1;}
 }

 /* Quick synapse count in each temporary file */
 for (ifile=0;ifile<(np*np);ifile++){ 

      nsyn[ifile]=0;
      while(!feof(tmpsyn[ifile])){
	  fread(&syndummy,sizeof(struct SYN_TMPMKNODE), 1, tmpsyn[ifile]);
	        nsyn[ifile]=nsyn[ifile]+1;}
      nsyn[ifile]=nsyn[ifile]-1;

      fprintf(stdout,"nsyn[ifile]=%i ifile=%i\n",nsyn[ifile],ifile);
 
 }

 /* Rewind temporary synapse files */
 for (ifile=0;ifile<(np*np);ifile++){ 

     rewind(tmpsyn[ifile]);  
 
 }

 /* Read and sort temporary synapse files */
 /* Write final synapse and corrected neuron files */
 for (ifile=0;ifile<(np*np);ifile++){ 

     syn_tmp=(struct SYN_TMPMKNODE *)calloc(nsyn[ifile],sizeof(struct SYN_TMPMKNODE));
          if (syn_tmp==NULL){fprintf(stderr,"Can't allocate temporary synapse memory space\n");exit(1);}
     fread(syn_tmp,sizeof(struct SYN_TMPMKNODE),(size_t)nsyn[ifile],tmpsyn[ifile]);
     fclose(tmpsyn[ifile]);

     qsort( (void *)syn_tmp, (size_t)nsyn[ifile],sizeof(struct SYN_TMPMKNODE),sort_neu_num); /* Sort by presynaptic neuron index */
 
     neu_pre=calloc(Nneu,sizeof(int));
     is=0;syn_n=0;
     for(i=0;i<Nneu;i++){
	 neu_pre[i]=0;
         while(i==(syn_tmp+syn_n)->neuron_pre){
	     syn_n++;
	     neu_pre[i]++;}
	 qsort((void *)(syn_tmp+is),(size_t)neu_pre[i],sizeof(struct SYN_TMPMKNODE), sort_delay); /* sort by delay */
	 is=syn_n;
     }

     /* Write synapse file */
     for(is=0;is<nsyn[ifile];is++)
     {
	 synwrite.neuron=(syn_tmp+is)->neuron_post;
	 synwrite.weight=(syn_tmp+is)->weight;
	 synwrite.delay=(syn_tmp+is)->delay;
	 synwrite.preclass_num=(syn_tmp+is)->preclass_num;
         printf("post= %d,preclass= %d, weight= %d,delay= %d\n",synwrite.neuron,synwrite.preclass_num,synwrite.weight,synwrite.delay);
         idel=synwrite.delay*1000*DELTA_T;
	 
	 //if( (synwrite.delay<0) ){
	 //     printf("delay=%i idel=%i\n",synwrite.delay,idel);
	 //}

	 syntime[idel]+=1;
	 fwrite(&synwrite,sizeof(struct SYNAPSE),1,synfile[ifile]);
     }

     /* Correct and write neuron file */
     for(i=0;i<Nneu;i++){
	 fread(&neu,sizeof(struct NEURON),1,neufile[ifile]);
	  oldsyn_num=neu.syn_num;
	  neu.syn_num=neu_pre[i];
	  //fprintf(stdout,"oldsyn_num=%i neu_pre=%i neu.syn_num=%i ifile=%i\n",oldsyn_num,neu_pre[i],neu.syn_num,ifile);
          fwrite(&neu,sizeof(struct NEURON),1,neuwritefile[ifile]);
     }
     /* Rewrite archaic configuration file */

     strcpy(fname,argv[1]);
     sprintf(name,"%s.cfg.%d",fname,ifile);
     cfg=fopen(name,"w");
     fprintf(cfg,"%s\n","N,V,H,L,S,$");
     fprintf(cfg,"%s\n",fname);
     fprintf(cfg,"%d\n",lrodz);
     fprintf(cfg,"%s\n",fname);
     fprintf(cfg,"%lu\n",nsyn[ifile]);
     fprintf(cfg,"%s\n",fname);
     fprintf(cfg,"%d\n",Nneu);
     fprintf(cfg,"%d,%d\n",Nneu,histoint);
     fprintf(cfg,"0\n");
     fprintf(cfg,"0.000010"); /* Point to readjust lambda!!! */
     fclose(cfg);
 }

 /* Close and delete files */
 for (ifile=0;ifile<(np*np);ifile++){ 
     fclose(neufile[ifile]);
     fclose(neuwritefile[ifile]);
     fclose(synfile[ifile]);
 }


/* Write time histogram files synt and linkt */

  synt=fopen("synt.dat","w");
    if(!synt){fprintf(stderr, "Can't open synt.dat file\n");return 1;}

  linkt=fopen("linkt.dat","w");
    if(!linkt){fprintf(stderr, "Can't open linkt.dat file\n");return 1;}

  for (i=0;i<MAXTIME;i++){

      fprintf(synt,"%i  %i\n",i,syntime[i]);
      fprintf(linkt,"%i  %i\n",i,linktime[i]);

  }

  fclose(synt);
  fclose(linkt);

return 0;
} 

/* sort_neu_num - needed in qsort() */
 
int sort_neu_num(const void * a, const void * b)
{
NEURON_INDEX neu1, neu2;
struct SYN_TMPMKNODE *syn_tmp1=NULL, *syn_tmp2=NULL;
syn_tmp1=(struct SYN_TMPMKNODE *)a;
syn_tmp2=(struct SYN_TMPMKNODE *)b;
neu1=syn_tmp1->neuron_pre;
neu2=syn_tmp2->neuron_pre;
return(neu1-neu2);
}


/* sort_delay - needed in qsort */
 
int sort_delay(const void * a, const void *b){
char delay1[1],delay2[1];
struct SYN_TMPMKNODE *syn1=NULL, *syn2=NULL;
syn1=(struct SYN_TMPMKNODE *)a;
syn2=(struct SYN_TMPMKNODE *)b;
*delay1=syn1->delay;
*delay2=syn2->delay;
 return( strcmp(delay1,delay2));
}
