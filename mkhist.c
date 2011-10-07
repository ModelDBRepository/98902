/* 05/29/03 P. Kudela now multi class supported */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "lnet.h"
int main(int argc, char* argv[])
{
  FILE *inp,*out,*nodefile;
  HIS* hist;
  int ne,ni,k,n,max_n,i,lrodz,*N,*Nh,ns,npr;
  char outname[100];
  if(argc<4)
    {
      fprintf(stderr,"\nUSAGE:\n\tmkhist <inp_name> <node_file> <which class)>\n");
      return 1;
    } 
  k=atoi(argv[3]);
  if(!(inp=fopen(argv[1],"r")))
    {
      fprintf(stderr,"Can't open input file: %s\n",argv[1]);
      return 2;
    }
  if(!(nodefile=fopen(argv[2],"r"))) 
    {
      fprintf(stderr,"Can't open node file: %s\n",argv[2]);
      return 2;
    }
  fscanf(nodefile,"%i\n",&lrodz);
  if(k>lrodz)
   {
     fprintf(stderr,"%i > no of kind\n",k);
     return 2;
   }
  N=(int*)malloc(lrodz*sizeof(int));
  for(i=0;i<lrodz;i++)
    {
    if(fscanf(nodefile,"%i\n",&N[i])!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",argv[2]);return 2;}
    }
  if(fscanf(nodefile,"%i\n",&ns)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",argv[2]);return 2;}
  if(fscanf(nodefile,"%i\n",&npr)!=1)
      {fprintf(stderr,"ERROR:fread(%s)\n",argv[2]);return 2;}
  fclose(nodefile);
  Nh=(int*)malloc(lrodz*sizeof(int));
  for(i=0;i<lrodz;i++)
  {
   Nh[i]=N[i]*ns*npr; 
   Nh[i]*=Nh[i];

   /* printf("lrodz=%i, Nh[i]=%i, ns=%i, npr=%i\n",lrodz,Nh[i],ns,npr); */

  }
  //ne=atoi(argv[2]);
  //ne=N[0]*ns*npr;
  //ne*=ne;
  //ni=atoi(argv[3]);
  //ni=N[1]*ns*npr;
  //ni*=ni;

  max_n=0;
  for(i=0;i<lrodz;i++)
   max_n=Nh[i]>max_n?Nh[i]:max_n;

  if(!(hist=calloc(max_n,sizeof(HIS))))
    { 
      fprintf(stderr,"error calloc(hist)\n");
      return 3;
    }
  n=0;
  while(1)
    {
      sprintf(outname,"HISTO.%i",n);
      if(!(out=fopen(outname,"w")))
         { fprintf(stderr,"Can't open output file %s\n",outname); return 2;}

      for(i=0;i<lrodz;i++)
        {
         if(fread(hist,sizeof(HIS),Nh[i],inp)!=Nh[i]) 
          {
           if(feof(inp))
           goto end;
           else
            {printf("Read error in input\n"); return 6;}
          }
         if(i==(k-1))
          {
           if(fwrite(hist,sizeof(HIS),Nh[i],out)!=Nh[i])          
            {fprintf(stderr,"Write  error in %s\n",outname);return 5;}
          }
        }
      fclose(out);
      n++;
    }
end:
  printf("%i histograms\n",n);
  return 0;
}
