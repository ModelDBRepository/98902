/***************************************************************************
 *   Copyright (C) 2006 by Piotr Franaszczuk                               *
 *   pfranasz@jhmi.edu                                                      *
 *   Johns Hopkins University, Baltimore,MD                                *                                                                      *
 *                                                                         * 
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#define NO_OF_PROC 16
#define NO_OF_LINKS 1000000

#define RISC
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>

#include "lnet.h"


int main(int argc,char *argv[])
{


  NEURON_INDEX **link_out;

  int lin_no[NO_OF_PROC][NO_OF_PROC];
//  int check_out[NO_OF_PROC][NO_OF_PROC];
//  int check_link=0;
  int lout_no[NO_OF_PROC];
  int nproc,last;
  char fname[100];
  char name[NO_OF_PROC][200];
//  int in_sum=0,out_sum=0;


  FILE* fil_in, *fil_out, *fil_tab, *fil_dst[NO_OF_PROC];
  int ii,i,nin,nout,src_proc,dst_proc,src_neu,dst_neu,preclassnum,w;
  int axflag,old_src;
  double del;
  struct LINK tmp;

  if(argc <3 ){fprintf(stderr,"USAGE:\n\tcrl fname(no.ext) nproc\n");exit(1);}

  nproc=atoi(argv[2]);
  if(nproc>NO_OF_PROC || nproc<1)
  {fprintf(stderr,"nproc= %i> %i\n",nproc,NO_OF_PROC);exit(1);}

  link_out=(NEURON_INDEX**)malloc(sizeof(NEURON_INDEX*)*nproc);
  if(!link_out)
  { fprintf(stderr,"malloc link_out");exit(1);}
  for(i=0;i<nproc;i++)
  {
    link_out[i]=(NEURON_INDEX*)calloc(sizeof(NEURON_INDEX),NO_OF_LINKS);
    if(!link_out[i])
    { fprintf(stderr,"calloc link_out[%i]",i);exit(1);}
    bzero(lin_no[i],nproc*sizeof(int));
//    bzero(check_out[i],nproc*sizeof(int));
    sprintf(name[i],"%s.dst.%i",argv[1],i);
    fil_dst[i]=fopen(name[i],"w+");
    if(fil_dst[i]==NULL){fprintf(stderr,"fopen(%s) error\n",name[i]);exit(1);}
  }
  strcat(strcpy(fname,argv[1]),".link1");
  fil_in=fopen(fname,"r");
  if(fil_in==NULL){fprintf(stderr,"fopen(%s) error\n",fname);exit(1);}

  bzero(lout_no,nproc*sizeof(int));

  if((fil_tab=fopen(strcat(strcpy(fname,argv[1]),".conn.tab"),"w"))==NULL)
  {fprintf(stderr,"fopen(%s) error\n",fname);exit(1);}
  fprintf(fil_tab,"%i\n",nproc);

  old_src=0;
  for(;;)
  {  //assumes src proc are sorted in order
  // read next line
   fscanf(fil_in,"%i:%i %i:%i %i %lf %i\n",&src_proc,&src_neu,&dst_proc,&dst_neu,&preclassnum,&del,&w);
          
    if(ferror(fil_in))
    {fprintf(stderr,"fscanf error on input: %s.link1\n",argv[1]);exit(1);}
    
//    check_out[dst_proc][src_proc]++;
    if(lout_no[dst_proc]==NO_OF_LINKS)
    {
      fprintf(stderr,"no_of out_links > %i for  src:%i dst:%i\n", NO_OF_LINKS,src_proc,dst_proc);
      exit(2);
    }
    fprintf(fil_dst[dst_proc],"%i:%i %i:%i %i %5.2f %i\n",src_proc,src_neu,dst_proc,dst_neu,preclassnum,del,w);
    if(ferror(fil_dst[dst_proc]))
    {fprintf(stderr,"fwrite dst.%i error\n",dst_proc);exit(3);}
    lin_no[src_proc][dst_proc]++;

    //if (lin_no[src_proc][dst_proc]>(USHRT_MAX-1)){
    //	        {fprintf(stderr,"lin_no exceeding max, src_proc=%i, dst_proc=%i\n",src_proc,dst_proc);exit(3);}
    //	}

	if(feof(fil_in))
	{
		link_out[dst_proc][lout_no[dst_proc]++]=(NEURON_INDEX)src_neu;
		//if (lout_no[dst_proc]>(USHRT_MAX-1)){
		//        {fprintf(stderr,"lout_no exceeding max, dst_proc=%i\n",dst_proc);exit(3);}
		//}

		src_proc=old_src+1;
	}
    if(src_proc!=old_src)
//	fprintf(stderr,"%i new=%i old=%i\n",check_link,src_proc,old_src);
    
    if(src_proc>old_src)
    {  // write out data for old_src processor
      nout=0;
      for(i=0;i<nproc;i++)
      {
        if(lout_no[i]>0)nout++;
      }
      sprintf(fname,"%s.lout.%i",argv[1],old_src);
      fil_out=fopen(fname,"w");
    if(fil_out==NULL){fprintf(stderr,"fopen(%s) error\n",fname);exit(1);}
      fwrite(&nout, sizeof(int),1,fil_out);

      if(nout>0)
      {

        printf("proc %i nout %i\n",old_src,nout);

        for(i=0;i<nproc;i++)
        {
          if(i!=old_src && lout_no[i]>0)
          {
            printf("%i->%i:%i\n",old_src,i,lout_no[i]);
            fwrite(&i,sizeof(int),1,fil_out);
            fwrite(&lout_no[i],sizeof(int),1,fil_out);
            fwrite(link_out[i],sizeof(NEURON_INDEX),lout_no[i],fil_out);
            /*		  for(ii=0;ii<lout_no[j][i];ii++)
              {
                printf("%i,",link_out[j][i][ii]);
            }
            printf("\n"); */
          }
        }
        // internal loops in the same node  last
        if(lout_no[old_src]>0)
        {
          printf("%i->%i:%i\n",old_src,old_src,lout_no[old_src]);
          fwrite(&old_src,sizeof(int),1,fil_out);
          fwrite(&lout_no[old_src],sizeof(int),1,fil_out);
          fwrite(link_out[old_src],sizeof(NEURON_INDEX),lout_no[old_src],fil_out);
          /*		  for(ii=0;ii<lout_no[j][j];ii++)
          		  {
          		  printf("%i,",link_out[j][j][ii]);
          		  }
          		  printf("\n"); */
        }

        //	    printf("\n");

      }
      fclose(fil_out);
      nout=0;
      for(i=0;i<nproc;i++)
        if(i!=old_src && lout_no[i]>0)nout++;

      fprintf(fil_tab,"%i %i:",old_src,nout);
      for(i=0;i<nproc;i++)
      {
        if(i!=old_src && lout_no[i]>0)
          fprintf(fil_tab," %i",i);
      }
      fprintf(fil_tab,"\n");

      // new src_proc

/*      // check
      for(i=0;i<nproc;i++)
      {
        out_sum+=lout_no[i];
        if(check_out[i][old_src]!=lout_no[i])
          fprintf(stderr,"check: %i->%i,%i,%i\n",old_src,i,check_out[i][old_src],lout_no[i]);
      }
//      fprintf(stderr,"%i:out_sum=%i,check_link=%i\n",old_src,out_sum,check_link);
*/
      if(feof(fil_in))break;
      bzero(lout_no,nproc*sizeof(int));
      old_src=src_proc;

    }
  	link_out[dst_proc][lout_no[dst_proc]++]=(NEURON_INDEX)src_neu;
//	check_link++;
  }/* while(feof) */
  fclose(fil_in);
  fclose(fil_tab);
/*
  for(i=0;i<nproc;i++)
  {
    int j;
    for(j=0;j<nproc;j++)
    {
      if(check_out[j][i]!=lin_no[i][j])
        fprintf(stderr,"check: %i-%i,%i,%i\n",i,j,check_out[i][j],lin_no[i][j]);
    }
  }
*/
  for(i=0;i<nproc;i++)
  {
    free(link_out[i]);
    fflush(fil_dst[dst_proc]);
  }
  free(link_out);

  for(dst_proc=0;dst_proc<nproc;dst_proc++)
  {

    int first=1;
    nin=0;
    for(i=0;i<nproc;i++)
    {
//      in_sum+=lin_no[i][dst_proc];
      if(lin_no[i][dst_proc]>0)nin++;
    }
//    fprintf(stderr,"%i: in_sum=%i\n",dst_proc,in_sum);
    sprintf(fname,"%s.lin.%i",argv[1],dst_proc);
    fil_out=fopen(fname,"w");
  if(fil_out==NULL){fprintf(stderr,"fopen(%s) error\n",fname);exit(1);}
    fwrite(&nin,sizeof(int),1,fil_out);
    do
    {
      if(nin>0)
      {
        if(first)printf("proc %i nin %i\n",dst_proc,nin);
        rewind(fil_dst[dst_proc]);
        old_src=-1;
        while(!feof(fil_dst[dst_proc]))
        {
          fscanf(fil_dst[dst_proc],"%i:%i %i:%i %i %lf %i\n",&src_proc,&src_neu,&dst_proc,&dst_neu,&preclassnum,&del,&w);
          if(ferror(fil_dst[dst_proc]))
          {fprintf(stderr,"fread dst.%i error\n",dst_proc);exit(3);}
          if(src_proc!=dst_proc ^ first )continue;
          if(src_proc>old_src && lin_no[src_proc][dst_proc]>0)
          {
            printf("%i<-%i:%i\n",dst_proc,src_proc,lin_no[src_proc][dst_proc]);
            fwrite(&src_proc,sizeof(int),1,fil_out);
            fwrite(&lin_no[src_proc][dst_proc],sizeof(int),1,fil_out);
            old_src=src_proc;
          }
          tmp.neuron_post=(NEURON_INDEX)dst_neu;
          tmp.delay=rint(del/(1000.*DELTA_T));  /*input in msec */
          tmp.weight=w;  /* weight */
          tmp.first=tmp.last=EMPTY_LIST; /* was EMPTY_LIST */
          tmp.preclass_num=(NEURON_INDEX)preclassnum;
          fwrite(&tmp,sizeof(struct LINK),1,fil_out);
        }


        //	      printf("\n");
        first=!first;
      }

    }
    while(!first);
    fclose(fil_out);
    fclose(fil_dst[dst_proc]);
    remove(name[dst_proc]);

  }
//  printf("no_links=%i\n",check_link);
  return 0;

}


