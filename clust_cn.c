/* Procedures for inter processor */
/* communication for simulations  */
/* P.J. Franaszczuk, JHU 2001     */  
//#define DEB 
#include "clust_cn.h"
extern unsigned short PORT;    // for version with command line parameter port
extern char *buf;
int nproc=0;                         /* number of processors to use */
int this_proc=-1;                  /* number of the processor for this host */
char *this_ip="000.000.000.000";   /* IP of this machine */ 
char **proc_ip;                    /* table of ip addresses of processors in dot notation */
char **proc_name;                  /* table of names of processors */
char *this_name=NULL;              /* name of this processor */
int *conn_sock=NULL;               /* array of  sockets descriptors */ 
int orig_sock=-1;                  /* socket for accept */
extern char *path;                 /* working path */
struct hostent *host;              /* temp pointer needed by gethostip*/
struct CONN cn_all;                /* all connections for this host */
struct CONN cn_inp;                /* input conn.sockets  fo host */
struct CONN cn_out;                /* out conn.sockets  for host */

int nfd;                           /* largest file descriptor to check in select */
fd_set fdr0;                       /* descriptor set for select */ 
int *buf_in;  //was short   

void close_all(int);  // defined in netclust.c

void close_socks()
{
  int i;
  if(conn_sock!=NULL)
    for(i=0;i<nproc;i++)
      if(conn_sock[i]>-1)
	{  /*short buf[2]={-1,0};
	      	   send(conn_sock[i],buf,2*sizeof(short),0);  */
	close(conn_sock[i]);
	conn_sock[i]=-1;
 	}
  if(orig_sock>-1){close(orig_sock); orig_sock=-1;}
  fflush(stderr);
  if(this_proc==0)
    {
      printf("\n Terminated \n");
      fflush(stdout);
    }

}
  

void init_path(char *argv0)
{
  char *c, buf[400];
  FILE* fil;
  int pid;
  this_name=NULL;
  if((c=gethostip(this_name))==NULL)
    { 
      char erbuf[100];
      sprintf(erbuf, "gethostip failed for this");
      print_err(ERROR,"init_path",erbuf,NULL,NULL);
    }
  this_ip=strdup(c); 
  this_name=strdup(host->h_name); 
  pid=getpid();
  fil=fopen(strcat(strcpy(buf,"netclust"),".pid"),"w");
  if(fil==NULL)
    print_err(FERROR,"init_path","pid",buf,fil);
  fprintf(fil,"%i",pid);
  fclose(fil);

}
    

void print_err(enum error_flag f, char* procnam, char * text, char* fname,FILE*fil )
{
  char buf[1000];
  switch(f)
    {
    case INFO: 
      fprintf(stderr,"%s\t%s: %s\n",this_name,procnam,text);fflush(stderr);
      return;
    case DEBUG:
#ifdef DEB
//      fprintf(stderr,"%s\tDEBUG:%s: %s\n",this_name,procnam,text);fflush(stderr);
#endif
      return;
    case WARNING:
      fprintf(stderr,"%s\tWARNING:%s: %s\n",this_name,procnam,text);fflush(stderr);
      return;
    case ERROR:
      fprintf(stderr,"%s\tERROR:%s: %s\n",this_name,procnam,text);
      break;
    case PERROR:
      sprintf(buf,"%s\tERROR:%s: %s",this_name,procnam,text);
      perror(buf);
      break;
    case FERROR:
      if(fil==NULL)
	fprintf(stderr,"%s\tERROR:%s: %s file: %s\n",this_name,procnam,text,fname);
      else if(feof(fil))
	fprintf(stderr,"%s\tEOF:%s: %s file: %s\n",this_name,procnam,text,fname);
      else if(ferror(fil))
	fprintf(stderr,"%s\tFERROR:%s: %s file: %s\n",this_name,procnam,text,fname);
      else
	fprintf(stderr,"%s\tERROR:%s: %s file: %s\n",this_name,procnam,text,fname);
      break;;
    }
  //  close_socks();exit(3);
  close_all(-1);
  return;
}
	
char* file_name(char *fname,char*name,char *ext)
{   char buf[400];
/* create name of file path is global */
 strcat(strcpy(buf,path),name);
 if(ext)sprintf(fname,"%s%s.%i",buf,ext,this_proc);
 else
   strcpy(fname,buf);
 return fname;
}
       
      
  

int compint(const void * i, const void * j)
{ return (*(int*)i-*(int*)j);}

int read_conn(char *name)
{
  /* reads file with connections CONNTAB_NAME
     format of the file:
     first line no of processors to use 
     next lines:  n  nc: c1 c2 c3                  
     where: n proc no(starting with 0)                       
              
     nc number of connections for this processor      
     c1, c2... c(nc) processors connected to this one 
     reads also nodes_table from file IPTAB_NAME       
  */
  FILE* fil;
  int i,j,l,*k;
  struct CONN *cn;
  char ** ip, **in, fname[250],buf[100], *c;
  char *procname="read_conn";
  struct CONN *tab_conn;/* list of connection struct. for all processors */

  /* read ip table */
  strcat(strcpy(buf,name),IPTAB_NAME);
  if((fil=fopen(file_name(fname,buf,NULL),"r"))==NULL)
    print_err(FERROR,procname,"fopen",fname,fil);

  if(fscanf(fil,"%i\n",&nproc)!=1)
    print_err(FERROR,procname,"fscanf(nproc)",fname,fil);   
  if(nproc<1)
    print_err(FERROR,procname,"nproc < 1",fname,fil);  
  proc_ip=ip=(char**)malloc(sizeof(char*)*nproc);
  if(proc_ip==NULL)print_err(ERROR,procname,"malloc(proc_ip)",NULL,NULL);
  proc_name=in=(char**)malloc(sizeof(char*)*nproc);
  if(proc_name==NULL)print_err(ERROR,procname,"malloc(proc_name)",NULL,NULL); 
    

  for(i=0;i<nproc;i++,ip++,in++)
    {
	
      if(fscanf(fil,"%s",buf)!=1 )                 
	{
	  char erbuf[100];
	  sprintf(erbuf, "fscanf:line %i",i);
	  print_err(FERROR,procname,erbuf,fname,fil);
	}    
      *in=strdup(buf);
      if( *in==NULL )print_err(ERROR,procname,"strdup(*in)",NULL,NULL);
      *in=strdup(buf);
      if((c=gethostip(*in))==NULL)
	{
	  char erbuf[100];
	  sprintf(erbuf, "gethostip failed for %s",*in);
	  print_err(ERROR,procname,erbuf,NULL,NULL);
      	    
	}
      *ip=strdup(c);
      if( *ip==NULL )print_err(ERROR,procname,"strdup(*ip)",NULL,NULL); 
         
          

      for(j=0;j<i;j++)
	if(!strcmp(proc_ip[j],*ip))
	  { 
	    char erbuf[100];
	    sprintf(erbuf, "repeated node %s in line %i and %i",*ip,j,i);
	    print_err(FERROR,procname,erbuf,fname,fil);
	  }    
    }
  fclose(fil);
	
  for(i=0;i<nproc;i++)
    if(!strcmp(this_ip,proc_ip[i]))break;
  if(i==nproc)
    {
      char erbuf[100];	
      sprintf(erbuf,"%s not in the table",this_name);
      print_err(ERROR,procname,erbuf,NULL,NULL);
    }
  this_proc=i; 
 
  /* read connections */
  strcat(strcpy(buf,name),CONNTAB_NAME);
  if((fil=fopen(file_name(fname,buf,NULL),"r"))==NULL)  
    print_err(FERROR,procname,"fopen",fname,fil);

  if(fscanf(fil,"%i\n",&i)!=1 || nproc!=i)
    print_err(FERROR,procname,"fscanf(nproc)",fname,fil);

  if(i<1)
    {
      print_err(WARNING,procname,"no interprocessor connections",NULL,NULL);
      cn_inp.nconn=cn_out.nconn=cn_all.nconn=0;
      cn_inp.lconn=cn_out.lconn=cn_all.lconn=NULL;
      return 0;
    }
                

  tab_conn = cn= (struct CONN *)malloc( sizeof( struct CONN )*nproc);
  if(cn==NULL) print_err(ERROR,procname,"malloc(tab_conn)",NULL,NULL);

  for(i=0;i<nproc;i++,cn++)
    {
      char col=0;
      if(fscanf(fil,"%i %i%c",&j,&cn->nconn,&col)!=3 
	 || cn->nconn < 0 
	 || cn->nconn > nproc 
	 || col!=':'
	 || i!=j)
	{
	  char erbuf[100];
	  sprintf(erbuf, "fscanf: proc no %i(%i)",i,j);
	  print_err(FERROR,procname,erbuf,fname,fil);
	}     
      else
	if(cn->nconn>0){
	  cn->lconn=k=(int*)malloc(sizeof(int)*cn->nconn);
	  if(k==NULL)print_err(ERROR,procname,"malloc(cn->lconn)",NULL,NULL);
	      
	  for(j=0;j<cn->nconn;j++,k++)
	    {
	      if(fscanf(fil,"%i",k)!=1 || *k<0 || *k==i)
		{ char erbuf[100];
		sprintf(erbuf, "read proc no %i(%i)",i,j);
		print_err(FERROR,procname,erbuf,fname,fil);
		}     
		    
	    }
	  fscanf(fil,"%[^\n]\n",buf); /* skip to end of line */
	  qsort(cn->lconn,cn->nconn,sizeof(int),compint);
	}
	else cn->lconn=NULL;
    }
  fclose(fil);

    
  cn_inp.lconn=(int*)malloc(sizeof(int)*nproc);  
  if(cn_inp.lconn==NULL)print_err(ERROR,procname,"malloc(cn_inp)",NULL,NULL);
 
  /* find inputs to this */
  for(i=0,l=0,cn=tab_conn;i<nproc;i++,cn++)
    {
      if(i!=this_proc)
	{
	  for(j=0;j<cn->nconn;j++)
	    if(cn->lconn[j]==this_proc)
	      {
		cn_inp.lconn[l++]=i;
		break;
	      }
	  free(cn->lconn);

	}
      else
	{
	  cn_out.nconn=cn->nconn;
	  cn_out.lconn=cn->lconn;

	}
    }
  if(l==0 && nproc>1)
    print_err(ERROR,procname,"no input connections to this",NULL,NULL);
      
  cn_inp.nconn=l;
  
  free(tab_conn);
  if(cn_out.nconn)
  { int ** b; //was short
  cn_out.buf=b=(int**)malloc(sizeof(int*)*cn_out.nconn); //was cn_out.buf=b=(short**)malloc(sizeof(short*)*cn_out.nconn)
    if(b==NULL)print_err(ERROR,procname,"malloc(cn_out.buf)",NULL,NULL);

    for(i=0;i<cn_out.nconn;i++,b++)
	if((*b=(int*)malloc(PACKETBUFFER+2))==NULL) //was if((*b=(short*)malloc(PACKETBUFFER+2))==NULL)
	print_err(ERROR,procname,"malloc(*cn_out.buf)",NULL,NULL);
    }

  if(cn_inp.nconn)
    { /* buffers for input */
	int ** b; //was short
	buf_in=(int*)malloc(PACKETBUFFER*nproc+2); //was buf_in=(short*)malloc(PACKETBUFFER*nproc+2);
      if(buf_in==NULL)print_err(ERROR,procname,"malloc(buf_in)",NULL,NULL);

      cn_inp.buf=b=(int**)malloc(sizeof(int*)*cn_inp.nconn); //was cn_inp.buf=b=(short**)malloc(sizeof(short*)*cn_inp.nconn)
      if(b==NULL)print_err(ERROR,procname,"malloc(cn_inp.buf)",NULL,NULL);

      for(i=0;i<cn_inp.nconn;i++,b++){
	  if((*b=(int*)malloc(PACKETBUFFER*nproc+2))==NULL) //was  if((*b=(short*)malloc(PACKETBUFFER*nproc+2))==NULL)
	  print_err(ERROR,procname,"malloc(*cn_inp.buf)",NULL,NULL);
	*b[0]=0;}
    }
  

  /* find all connections (inp & out) */

  cn_all.lconn=(int*)malloc(sizeof(int)*(nproc-1));  
  if(cn_all.lconn==NULL)print_err(ERROR,procname,"malloc(cn_all)",NULL,NULL);

  for(i=0,j=0,l=0;i<cn_inp.nconn && j<cn_out.nconn;)
    {
      if(cn_inp.lconn[i]<cn_out.lconn[j])
	cn_all.lconn[l++]=cn_inp.lconn[i++];

      else if(cn_inp.lconn[i]==cn_out.lconn[j])
	{ cn_all.lconn[l++]=cn_inp.lconn[i++];j++;}
      else
	cn_all.lconn[l++]=cn_out.lconn[j++];
    }
  while(i<cn_inp.nconn)cn_all.lconn[l++]=cn_inp.lconn[i++];
  while(j<cn_out.nconn)cn_all.lconn[l++]=cn_out.lconn[j++];
  cn_all.nconn=l;

  /*    
	fprintf(stderr,"\ncn_all:\n nconn=%i: ", cn_all.nconn);
	for(i=0;i<cn_all.nconn;i++)
	fprintf(stderr,"%i,",cn_all.lconn[i]);

	fprintf(stderr,"\ncn_inp:\n nconn=%i: ", cn_inp.nconn);
	for(i=0;i<cn_inp.nconn;i++)
	fprintf(stderr,"%i,",cn_inp.lconn[i]);

	fprintf(stderr,"\ncn_out:\n nconn=%i: ", cn_out.nconn);
	for(i=0;i<cn_out.nconn;i++)
	fprintf(stderr,"%i,",cn_out.lconn[i]);
    
	fprintf(stderr,"\n");
  */  
  return nproc;
	   
}



char * gethostip(char* name)
{
  /* return IP of named host as string
     and host structure in host 
     if name == NULL returns this host IP
     side effect midifies global struct host
  */
  char buf[100];
  size_t len=100;
  struct in_addr in;
  if(name==NULL)
    {
      if(gethostname(buf,len)<0)
	{
	  return NULL;
	}
      host=gethostbyname(buf);
    }
  else
    host=gethostbyname(name);
  
 
  if(host){
    memcpy(&in.s_addr,*host->h_addr_list, sizeof(in.s_addr));
    return inet_ntoa(in);
  }
  else
    {
      return NULL;
   
    }
}

void proc_conn()
{
  char *pname="proc_conn";
  int i,j,iconn,icc,nconn,mconn,*l,flag=1;                    /* flag for ioctl */
 
  struct sockaddr_in conn_adr, serv_adr;
  int conn_len=sizeof(conn_adr);
  fd_set fdset,fdacc, fdex;
  struct timeval tim = {10,0}; /* timeout for select*/ 
  int nrep = 20;  /* repetition for connect */
   
  
 
  if((conn_sock =(int*)malloc(sizeof(int)*nproc))==NULL)
    print_err(ERROR,pname,"memory:conn_sock",NULL,NULL);
  for(i=0;i<nproc;conn_sock[i++]=-1);

  if(this_proc < nproc-1)
    {
      memset(&serv_adr, 0, conn_len);	   /* Clear it out  */
      serv_adr.sin_family = AF_INET;
      serv_adr.sin_addr.s_addr = htonl(INADDR_ANY);             /* accept any */
      serv_adr.sin_port   = htons(PORT);  
  
      if ((orig_sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 	/* SOCKET  */
   	print_err(PERROR,pname,"orig_socket",NULL,NULL);

      if (ioctl(orig_sock, FIONBIO, &flag) < 0 )                /* ioctl */
	print_err(ERROR,pname,"orig_ioctl",NULL,NULL); 

      if (bind(orig_sock,(struct sockaddr*)&serv_adr,sizeof(serv_adr))<0)                    /* BIND */   
	print_err(PERROR,pname,"bind",NULL,NULL);

      if (listen(orig_sock,nproc-1)<0)                         /* listen  */   
	print_err(PERROR,pname,"bind",NULL,NULL);
      nfd=orig_sock+1;
    }


  
  memset(&conn_adr, 0, sizeof(conn_adr));	   /* Clear it out  */
  conn_adr.sin_family = AF_INET;
  conn_adr.sin_port   = htons( PORT );  
  

  
 
  FD_ZERO(&fdset);   /* set for connect */
  FD_ZERO(&fdacc);   /* set for access  */
  FD_ZERO(&fdex);   /* set for exceptions  */
 

  /* open sockets for processors to connect to (less < this; inp&out!! ) */
  i=0;
  while(i < cn_all.nconn && ((j=cn_all.lconn[i]) < this_proc)) 
    {
      
      if((conn_sock[j]=socket(AF_INET,SOCK_STREAM,0)) < 0 )
	print_err(PERROR,pname,"conn_socket",NULL,NULL); 

      if ((ioctl(conn_sock[j], FIONBIO, &flag)) < 0 ) 
	print_err(PERROR,pname,"conn_ioctl",NULL,NULL); 
      nfd=conn_sock[j]+1;     
      i++;
    } 

  iconn=i;              /* no of connections to make with connect*/
  nconn=cn_all.nconn-i;    /* no of connections to accept */ 
 
  for(i=0;i<nfd;i++)  /* exceptions */
    FD_SET(i,&fdex);

  
  if(nconn>0)
    {
      nconn=try_accept(orig_sock,nconn);
      if(nconn>0)
     	FD_SET(orig_sock,&fdacc);
    }


  mconn=iconn;
  for(i=0;i<mconn;i++)
    {
      if(try_connect(cn_all.lconn[i],&conn_adr))
	iconn--;
      else
	FD_SET(conn_sock[cn_all.lconn[i]],&fdset);
    }
  
 
  /* check if all connections made and no errors */
  icc=0;
  while(iconn>0 || nconn>0)
    {
      int j;
      fd_set fds,fdr,fde;
   
 
      memcpy(&fds,&fdset,sizeof(fdset));
      memcpy(&fdr,&fdacc,sizeof(fdacc));
      memcpy(&fde,&fdex,sizeof(fdacc));

      if((i=select(nfd,&fdr,&fds,&fde,&tim))<0)
	print_err(PERROR,pname,"select",NULL,NULL);
   
      if(i==0)
	{
	  char erbuf[100];
	  sprintf(erbuf," select timeout iconn=%i, nconn=%i",iconn,nconn);
	  print_err(ERROR,pname,erbuf,NULL,NULL);
	}

      for(j=0;j<nfd;j++)
	if(FD_ISSET(j,&fde))
	  { char erbuf[100];
	  sprintf( erbuf,"exception in fd %d",j);
	  print_err(PERROR,pname,erbuf,NULL,NULL);
	  }
      if(nconn && FD_ISSET(orig_sock,&fdr))
	{
	  nconn=try_accept(orig_sock,nconn);
	  i--;
	}  
      for(j=0;j<mconn && i>0;j++)
	{
	  char erbuf[100];
	  int er_so,er_len=sizeof(int);
	  int jj=cn_all.lconn[j];

	  if(!FD_ISSET(conn_sock[jj],&fds))continue;
	  
	  sprintf(erbuf,"connection to IP %s",proc_ip[jj]);
          if(getsockopt(conn_sock[jj],SOL_SOCKET,SO_ERROR,&er_so,&er_len))
	    {
	      strcat(erbuf,":getsockopt");
	      print_err(PERROR,pname,erbuf,NULL,NULL);
	    }
	  else if(er_so)
	    { 
	      if(er_so==ECONNREFUSED)
		{
		  close(conn_sock[jj]);
		  if((conn_sock[jj]=socket(AF_INET,SOCK_STREAM,0)) < 0 )
		    print_err(PERROR,pname,"reconn_socket",NULL,NULL); 

		  if ((ioctl(conn_sock[jj], FIONBIO, &flag)) < 0 ) 
		    print_err(PERROR,pname,"reconn_ioctl",NULL,NULL); 
		    
		  if(try_connect(jj,&conn_adr))
		    {
		      iconn--;
		      FD_CLR(conn_sock[jj],&fdset);
		      icc=0;
		    }
		  else
		    {
                      //		      if(icc++>tim.tv_sec)
                    if(icc++>nrep)
			{
			  strcat(erbuf,"\n\tconnect repetition limit");
			  print_err(ERROR,pname,erbuf,NULL,NULL);
			}
		    } 
		  
		  /*  strcat(erbuf,"getsockopt:ECONREFUSED");
		      print_err(DEBUG,pname,erbuf,NULL,NULL); */
		  continue;
		}
	      else
		{
		  errno=er_so;
		  strcat(erbuf,":getsockopt\n\tSO_ERROR");
		  print_err(PERROR,pname,erbuf,NULL,NULL);
		}
	    }
	  FD_CLR(conn_sock[jj],&fdset);
	  i--;iconn--;icc=0;
	  strcat(erbuf," successful");
	  print_err(DEBUG,pname,erbuf,NULL,NULL);
	}
 
      sleep(1);
    }
  if(orig_sock>-1) close(orig_sock);  /* no longer needed */ 

  /*    set mask for select for future reading */  

  FD_ZERO(&fdr0);
  nfd=-1;
  l=cn_inp.lconn;
  for(i=0;i<cn_inp.nconn;i++,l++)
    {
      FD_SET(conn_sock[*l],&fdr0);
      if(conn_sock[*l]>nfd)nfd=conn_sock[*l];
    }
  nfd++;
   
}
 
	 
int try_accept(int orig_sock, int nconn)
{  
  char *pname="try_accept";
  struct sockaddr_in conn_adr;
  int conn_len=sizeof(conn_adr),new_sock;
   
   
  memset(&conn_adr, 0, sizeof(conn_adr));
  print_err(DEBUG,pname,"entry",NULL,NULL); 
  while(nconn>0) {
  
    if((new_sock=accept(orig_sock,(struct sockaddr*)&conn_adr,&conn_len))<0)  /*accept */
      {
	if(errno!=EWOULDBLOCK) 
	  print_err(PERROR,pname,"accept",NULL,NULL);
	return nconn; 
      }
    else
      {  char erbuf[100];
      int j;
      char * ip = inet_ntoa(conn_adr.sin_addr);         /* get client ip  */
      for(j=0;j<nproc;j++)
	if(!strcmp(ip,proc_ip[j]))break;
      if(j==nproc)
	{
	  sprintf(erbuf,"accepted IP %s not in the table",ip);
	  print_err(ERROR,pname,erbuf,NULL,NULL);
	}
      else if(j<this_proc)
	{
	  sprintf(erbuf,"accepted IP %s < this",ip);
	  print_err(ERROR,pname,erbuf,NULL,NULL);
	}
      sprintf(erbuf,"accepted connection from %s",ip);
      print_err(DEBUG,pname,erbuf,NULL,NULL);
              
		 
      conn_sock[j]=new_sock;
      nconn--;
      }
	  
  }
  return 0;
 
}	 

int try_connect(int i, struct sockaddr_in *conn_adr)
{   
  char * pname="try_connect";
   

  if(!inet_aton(proc_ip[i],&(conn_adr->sin_addr)))     /* IP address */
    print_err(ERROR,pname,"inet_aton",NULL,NULL); 
      	  
  if(connect(conn_sock[i],(struct sockaddr *)conn_adr, sizeof(struct sockaddr)) < 0)
    { char erbuf[100];   
    if(errno!=EINPROGRESS)
      {	      
	sprintf(erbuf,"connect to %s",proc_ip[i]);
	print_err(PERROR,pname,erbuf,NULL,NULL); 
      }
    sprintf(erbuf,"trying to connect to %s",proc_ip[i]);
    print_err(DEBUG,pname,erbuf,NULL,NULL); 
    return 0;
    }
  else
    { char erbuf[100];
    sprintf(erbuf,"connected to %s",proc_ip[i]);
    print_err(DEBUG,pname,erbuf,NULL,NULL);
    return 1;
    } 
}
