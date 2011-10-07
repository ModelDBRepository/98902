#include <stdio.h>
#include <stdlib.h>
#define NO_OF_PROC 16
int main(int argc, char* argv[])
{
  char fname[200],name[NO_OF_PROC][200];
  FILE *fil_in, *fil_sort[NO_OF_PROC];
  int i,nproc;
  int src_proc,dst_proc,src_neu,dst_neu,preclassnum,w;
  double del;

  if(argc <3 ){fprintf(stderr,"USAGE:\n\tcrl fname(no.ext) nproc\n");exit(1);}

  nproc=atoi(argv[2]);
  if(nproc>NO_OF_PROC || nproc<1)
  {fprintf(stderr,"nproc= %i> %i\n",nproc,NO_OF_PROC);exit(1);}

  strcat(strcpy(fname,argv[1]),".link");
  fil_in=fopen(fname,"r");
  if(fil_in==NULL){fprintf(stderr,"fopen(%s) error\n",fname);exit(1);}
  for(i=0;i<nproc;i++)
  {
    sprintf(name[i],"%s.dst.%i",argv[1],i);
    fil_sort[i]=fopen(name[i],"w+");
    if(fil_sort[i]==NULL){fprintf(stderr,"fopen(%s) error\n",name[i]);exit(1);}
  }
  while(!feof(fil_in))
  {
    fscanf(fil_in,"%i:%i %i:%i %i %lf %i\n",&src_proc,&src_neu,&dst_proc,&dst_neu,&preclassnum,&del,&w);
    if(ferror(fil_in))
    {fprintf(stderr,"fscanf error on input: %s.link\n",argv[1]);exit(1);}
    fprintf(fil_sort[src_proc],"%i:%i %i:%i %i %lf %i\n",src_proc,src_neu,dst_proc,dst_neu,preclassnum,del,w);
  }
  fclose(fil_in);
  strcat(strcpy(fname,argv[1]),".link1");
  fil_in=fopen(fname,"w");
  if(fil_in==NULL){fprintf(stderr,"fopen(%s) error\n",fname);exit(1);}
  for(i=0;i<nproc;i++)
  {
    rewind(fil_sort[i]);
    while(!feof(fil_sort[i]))
    {
      fscanf(fil_sort[i],"%i:%i %i:%i %i %lf %i\n",&src_proc,&src_neu,&dst_proc,&dst_neu,&preclassnum,&del,&w);
	 fprintf(fil_in,"%i:%i %i:%i %i %5.2f %i\n",src_proc,src_neu,dst_proc,dst_neu,preclassnum,del,w);

  
    }
    fclose(fil_sort[i]);
    remove(name[i]);
  }
  fclose(fil_in);
  return 0;
}


