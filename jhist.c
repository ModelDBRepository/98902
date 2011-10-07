/* 05/29/03 P. Kudela now multi-class supported */

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <math.h>
#include "lnet.h"


int histo_f,out_f,vout_f,noise_f,iext_f,link_f,stim_f;  /* flagi sterujace*/
char flags[100];
FILE *cfg;

void  set_flags()
{
	/*signed */char z;
	/*unsigned*/ char buf[100],*l;
	l=flags;
	histo_f=out_f=vout_f=noise_f=iext_f=link_f=stim_f=0;  /* flagi sterujace*/

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
	    //	    puts(buf);
	  
}


int main(int argc ,char *argv[])
{
  HIS *histi,*hi,*his;
  FILE **his_inp, *his_out;
  int iproc, i,j,nproc,npr,n,*his_nr,*shnr,his_sum,his_max,hnr,heof,no_kind,ix,iy,ns,l,k,nn;
  char fname[100],name[100];
  int *dum;
 


  if(argc<3)
    {fprintf(stderr,"USAGE:\n\tjhist <inpname> <outname> \n\t no ext in names\n");return 1;}
  

  strcpy(name,argv[1]);
  strcat(name,".node");
  if((cfg = fopen(name,"r"))==NULL)
    {fprintf(stderr,"ERROR:fopen(%s)\n",name);return 2;}
  if(fscanf(cfg,"%i\n",&no_kind)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",name);return 2;}
  if((shnr=(int*)calloc(no_kind,sizeof(int)))==NULL)
    {fprintf(stderr,"ERROR:calloc(shnr)\n");return 1;}
  his_max=0;
  for(i=0;i<no_kind;i++)
    {
    if(fscanf(cfg,"%i\n",shnr+i)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",name);return 2;}
    if(shnr[i]>his_max)his_max=shnr[i];
    }

   if(fscanf(cfg,"%i\n",&ns)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",name);return 2;}
   if(fscanf(cfg,"%i",&npr)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",name);return 2;}
  fclose(cfg);
  //npr=atoi(argv[3]);
  nproc=npr*npr;
    
  if((his_nr=(int*)calloc(nproc,sizeof(int)))==NULL)
    {fprintf(stderr,"ERROR:calloc(his_nr)\n");return 1;} 

  if((his_inp=(FILE**)calloc(nproc,sizeof(FILE*)))==NULL)
    {fprintf(stderr,"ERROR:calloc(his_inp)\n");return 1;} 
  his_sum=0;
  for(iproc=0;iproc<nproc;iproc++)
    {
      sprintf(fname,"%s.cfg.%d",argv[1],iproc);
      if((cfg = fopen(fname,"r"))==NULL)
	{fprintf(stderr,"ERROR:fopen(%s)\n",fname);return 2;} 

      set_flags();

      if(!histo_f)
	{
	  his_nr[iproc]=0;
	  fclose(cfg);
	  continue;
	}
      fscanf(cfg,"%s\n",fname);
      fscanf(cfg,"%i\n",&no_kind);
      fscanf(cfg,"%s\n",fname);
      fscanf(cfg,"%i\n",&n);
      fscanf(cfg,"%s\n",fname);
      fscanf(cfg,"%i\n",&n);
      fscanf(cfg,"%i",&his_nr[iproc]);

      his_sum+=his_nr[iproc];
      if(his_nr[iproc]>0)
	{
	  printf("%d %d %d\n", iproc,his_nr[iproc],his_sum);
	  sprintf(fname,"histo_%s.%d",argv[1],iproc);
	  if((his_inp[iproc] = fopen(fname,"r"))==NULL)
	    {fprintf(stderr,"ERROR:fopen(%s)\n",fname);return 2;} 
	}
      fclose(cfg);
    }
  if(his_sum==0)
    {printf("nothing to do\n");return 0;}

/*  strcpy(fname,argv[1]);
  strcat(fname,".his");
  if((cfg = fopen(fname,"r"))==NULL)
    {fprintf(stderr,"ERROR:fopen(%s)\n",fname);return 2;} 
  if(fscanf(cfg,"%i\n",&ns)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",fname);return 2;}
  if((shnr=(int*)calloc(no_kind,sizeof(int)))==NULL)
    {fprintf(stderr,"ERROR:calloc(shnr)\n");return 1;} 
  his_max=0;
  for(i=0;i<no_kind;i++)
    {
    if(fscanf(cfg,"%i\n",shnr+i)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",fname);return 2;}
    if(shnr[i]>his_max)his_max=shnr[i];
    }
  fclose(cfg);
*/

  if((histi=(HIS*)calloc(his_max*his_max*ns*npr,sizeof(HIS)))==NULL)
    {fprintf(stderr,"ERROR:calloc(histi)\n");return 1;} 

  sprintf(fname,"histo_%s",argv[2]);
  if((his_out = fopen(fname,"w"))==NULL)
    {fprintf(stderr,"ERROR:fopen(%s)\n",fname);return 2;}
  while(1)
    for(i=0;i<no_kind;i++)
      {	
	nn=shnr[i];
	if((hnr=nn*nn*ns)>0)
	  for(iy=0;iy<npr;iy++)
	    	for(j=0;j<ns;j++)
		  { 
		 for(ix=0,hi=histi,iproc=iy*npr;ix<npr;ix++,iproc++,hi+=hnr)	
		  {	 
		  if(fread(hi,sizeof(HIS),hnr,his_inp[iproc])!=hnr)
		    {
		      if(feof(his_inp[iproc]))
			return 0;
		      fprintf(stderr,"ERROR:fread(histo_%s.%d)\n",argv[1],iproc);
		      return 3;
		    }
		  //		  printf("read %i to %i proc %i\n",hnr,hi-histi,iproc);
		  }
		for(l=0,his=histi;l<nn;l++,his+=nn)	    
		  for(k=0,hi=his;k<ns*npr;k++,hi+=nn*nn)
		    {        
		      if(fwrite(hi,sizeof(HIS),nn,his_out)!=nn)
			{fprintf(stderr,"ERROR:fwrite(histo_%s)\n",argv[2]);return 4;}
		      //		      printf("writen %i from %i\n",nn,hi-histi);
		    }
	      }
      }

}
