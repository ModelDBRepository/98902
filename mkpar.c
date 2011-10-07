#define RISC
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
//#include <curses.h>
#include "lnet.h"

#define randomm(num) (int)((double)num*(double)rand()/((double)(RAND_MAX+1.0)))
#define randomize()     srand((unsigned)time(NULL))

#define randoml(num) (int)((double)num*drand48())
#define randomizel()     srand48((long int)time(NULL))


long filelength(int handle)
{
   FILE *fp;
   long int curpos, length;
   if((fp=fdopen(handle,"r+"))==NULL)
    return(-1);
   curpos = ftell(fp);
   fseek(fp, 0L, SEEK_END);
   length = ftell(fp);
   fseek(fp, curpos, SEEK_SET);
   return length;
}

long filelengthf(char *fil)
{
   FILE *fp;
   long int curpos, length;
   if((fp=fopen(fil,"r+"))==NULL)
    return(-1);
   curpos = ftell(fp);
   fseek(fp, 0L, SEEK_END);
   length = ftell(fp);
   fseek(fp, curpos, SEEK_SET);
   return length;
}

char konfig[5];


int main(int argc, char *argv[])
{

FILE *fp;
int fo,i,k,rez,update=0;
int kk;
long len;
struct PARAMETERS *par;
PARAMETER_INDEX N_class;
char opc[9], buf[100];
float   dt = 1000*DELTA_T;
//static VAR dt=.01;
static VAR onset_exc=0.5,decay_exc=3.0;
static VAR onset_inh=0.5,decay_inh=3.0;
VAR tau_o_exc[NEXCCLASS],tau_d_exc[NEXCCLASS];
VAR tau_o_inh[NINHCLASS],tau_d_inh[NINHCLASS];
VAR tmax_exc[NEXCCLASS],tmax_inh[NINHCLASS],Imax_exc[NEXCCLASS],Imax_inh[NINHCLASS];
VAR robol;
SPIKE_INTERVAL *tab_del=NULL;
VAR delay_ms=1.0;
char str[10]="";
char name[5],choose;
char *file_name;
VAR Amp=0.00655;

int get_cfg(char *);
int get_par(char *,int);

VAR improve_par(char *,VAR, char *);

file_name=calloc(1,sizeof(char));

if(argc<2)
 {
  printf("cfg file? : ");
  scanf("%s",file_name);
  getchar();
 }
else
  strcpy(file_name,argv[1]);

if(get_cfg(file_name)==-1)
 {printf("cannot read config from %s.cfg\n",file_name);return(-1);}
if((rez=get_par(file_name,1))<0)
 {
  printf("cannot read the number of neurons from %s.cfg\n",file_name);
  return(-1);
 }
if(rez>0)
 N_class=(PARAMETER_INDEX)rez;
 else
 {
  printf("Number of classes ?: ");
  scanf("%c",&N_class);
  N_class-='0';
  printf(" N_class = %d\n",N_class);
 }
if((fp=fopen(strcat(strcpy(name,file_name),".cfg"),"r"))==NULL)
 {
  printf("Error opening %5s file\n",name);
  return(-1);
 }

tab_del=malloc(sizeof(SPIKE_INTERVAL)*N_class);   /*tablica delayow*/
printf("number of classes is %d\n",N_class);
fseek(fp,0L,SEEK_SET);

fscanf(fp,"%8s",opc);
if(strncmp(&opc[0],"ONSET",5)==0)          /* nie potrzebne */

/* cfg file not needed anymore*/
fclose(fp);

if((fo=open(strcat(strcpy(name,file_name),".par"),O_RDWR|O_CREAT,S_IREAD|S_IWRITE|S_IRGRP|S_IROTH))==-1)
 {
  printf("Error opening %5s file\n",name);
   free(tab_del);
     return(-1);
 }
len = filelengthf(strcat(strcpy(name,file_name),".par"));

par=calloc(N_class,sizeof(struct PARAMETERS));

if( ( N_class >= 1 )&&( len > 0 ) )
 len-=N_class*(N_class*sizeof(SPIKE_INTERVAL));

if( ( N_class*sizeof(struct PARAMETERS))==len )
 {
  printf("update %s ?",name);
  fflush(stdin);
  choose=getchar();
  if((choose=='n')||(choose=='N'))
   {
    free(tab_del);
      free(par);
	close(fo);
	  return(0);
    }
  update=1;
  dt = improve_par("dt step ",dt,"6.4");
  for(k=0;k<N_class;k++ )    /* czytam istniejacy zbior */
   {
   rez=read(fo,par+k,sizeof(struct PARAMETERS));
    if(rez!=sizeof(struct PARAMETERS))
     {
      printf("Update: error reading %5s file\n",name);
	free(tab_del);
	  free(par);
	    close(fo);
	      return(-1);
     }

   rez=read(fo,tab_del,sizeof(SPIKE_INTERVAL)*N_class);
    if(rez!=(sizeof(SPIKE_INTERVAL)*N_class))
     {
      printf("Update: error reading delay tab from %s file\n",name);
	 free(tab_del);
	    free(par);
	      close(fo);
		return(-1);
     }
   }
 }
else 
 {
  if(len==0 )
   {
    for(k=0;k<N_class;k++ )
    {

    for (kk=0;kk<NEXCCLASS;kk++){	
	 tau_o_exc[kk]=onset_exc/dt;
         tau_d_exc[k]=decay_exc/dt;
         tmax_exc[kk]=(tau_o_exc[kk]*tau_d_exc[kk]*log(tau_o_exc[kk]/tau_d_exc[kk]))/(tau_o_exc[kk]-tau_d_exc[kk]);
         Imax_exc[kk]=exp(-tmax_exc[kk]/tau_d_exc[kk])-exp(tmax_exc[kk]/tau_o_exc[kk]);
         (par+k)->de_o[kk]=onset_exc;
         (par+k)->de_d[kk]=decay_exc;
    }

    for (kk=0;kk<NINHCLASS;kk++){
	 tau_o_inh[kk]=onset_inh/dt;	
         tau_d_inh[kk]=decay_inh/dt;
	 tmax_inh[kk]=(tau_o_inh[kk]*tau_d_inh[kk]*log(tau_o_inh[kk]/tau_d_inh[kk]))/(tau_o_inh[kk]-tau_d_inh[kk]);
	 Imax_inh[kk]=exp(-tmax_inh[kk]/tau_d_inh[kk])-exp(tmax_inh[kk]/tau_o_inh[kk]); 
         (par+k)->di_o[kk]=onset_inh;
         (par+k)->di_d[kk]=decay_inh;
    }


    (par+k)->V_L= -50.; 
    (par+k)->V_Na=55.;
    (par+k)->V_K=-72.;
#ifdef CALCIUM
    (par+k)->V_Ca=124.;
    (par+k)->G_Ca=1.0*dt;
    (par+k)->G_K_Ca=3.5*dt;
    (par+k)->G_a=12.5*dt;
    (par+k)->Kp=0.0002;
    (par+k)->Kd=0.5;
    (par+k)->Kc=2.;
    (par+k)->R=0.006*dt;
    (par+k)->tau_x=25.0/dt;
    (par+k)->tau_b=10./dt;
    (par+k)->Utot=100;
    (par+k)->Utot0=100;
    (par+k)->Utot1=100;
    (par+k)->Utot2=100;
    (par+k)->Utot3=100;
    (par+k)->Utot4=100;
    (par+k)->Utot5=100;
    (par+k)->dr=0.0005;
    (par+k)->Rsh0=0.017;
    (par+k)->Rsh1=(par+k)->Rsh0+(par+k)->dr;
    (par+k)->Rsh2=(par+k)->Rsh1+(par+k)->dr;
    (par+k)->Rsh3=(par+k)->Rsh2+(par+k)->dr;
    (par+k)->Rsh4=(par+k)->Rsh3+(par+k)->dr;
    (par+k)->Rsh5=(par+k)->Rsh4+(par+k)->dr;
    (par+k)->Rsh=(par+k)->Rsh5+(par+k)->dr;
    (par+k)->back=0.1*dt;
    (par+k)->forw=0.3*dt;
    (par+k)->PCa=0.15*dt;
    (par+k)->Cao=2;
    (par+k)->gpump=1.6*dt;
    (par+k)->Kpump=0.75;
    (par+k)->D=0.000000006;
#endif
    (par+k)->G_Na=120.*dt; 
    (par+k)->G_L=0.3*dt;
    (par+k)->G_K_s=15.0*dt;
    (par+k)->i_ext=15.0*dt;
    (par+k)->Esyn_i=-72.;
    (par+k)->Esyn_e=-10.;

    for (kk=0;kk<NEXCCLASS;kk++){

         (par+k)->Aep[kk]=(Amp/Imax_exc[kk])*dt;

    }

    for (kk=0;kk<NINHCLASS;kk++){

         (par+k)->Aip[kk]=(Amp/Imax_inh[kk])*dt;

    }

    (par+k)->An=120*Amp*dt;
    *(tab_del+k)=100;
    }

   } /*end if len==0*/
  else
   {
     printf("Size of %s = %ld and siezeof(struct PARAMETERS) = %d are different\n",name,len,sizeof(struct PARAMETERS));
     choose=getchar();
     fflush(stdin);
     if((choose=='n')||(choose=='N'))
      {
       free(tab_del);
	free(par);
	  close(fo);
	    return(-1);
      }
   }

 }
printf("dt time step is %5.3f ms\n",dt);

for(k=0;k<N_class;k++ )    /*petla ustawiania parametrow po N_class*/
{
 printf(" << class : %d\n",k+1);
 delay_ms=(VAR)(*(tab_del+k))*dt;
 printf(" \n");
 printf("synaptic delay = %6.2lf ms; ",delay_ms);
 printf(" \n");
 fflush(stdin);
 /* gets() works better than fgets() */
 gets(str);
 if((isdigit(str[0]))||(str[0]=='.'))
  {
  delay_ms=(VAR)atof(str);
  }
 (par+k)->del=tab_del;
 robol=delay_ms/dt;
 *((par+k)->del+k)=robol;
  
  for (kk=0;kk<NEXCCLASS;kk++){

       onset_exc=(par+k)->de_o[kk];
 
       printf("\nEPSP from class #(0- %d) %d\n",NEXCCLASS-1,kk);
       onset_exc=improve_par("onset EPSPs rise time       tau_o = ",onset_exc,"6.4");
 
       tau_o_exc[kk]=onset_exc/dt;

       (par+k)->de_o[kk]=onset_exc;

       decay_exc=(par+k)->de_d[kk];

       printf("\nEPSP from class #(0- %d) %d\n",NEXCCLASS-1,kk);
       decay_exc=improve_par("decay EPSPs fall time        tau_d = ",decay_exc,"6.4");

       tau_d_exc[kk]=decay_exc/dt;

       (par+k)->de_d[kk]=decay_exc;

  }

  for (kk=0;kk<NINHCLASS;kk++){

       onset_inh=(par+k)->di_o[kk];

       printf("\nIPSP from class #(%d - %d) %d\n",NEXCCLASS,NUMCLASS-1,kk+NEXCCLASS);
       onset_inh=improve_par("onset IPSPs decrement       tau_o = ",onset_inh,"6.4");

       tau_o_inh[kk]=onset_inh/dt;
       (par+k)->di_o[kk]=onset_inh;

       decay_inh=(par+k)->di_d[kk];

       printf("\nIPSP from class #(%d - %d) %d\n",NEXCCLASS,NUMCLASS-1,kk+NEXCCLASS);
       decay_inh=improve_par("decay IPSPs fall time       tau_d = ",decay_inh,"6.4");

       tau_d_inh[kk]=decay_inh/dt;
       (par+k)->di_d[kk]=decay_inh;

  }


if(update)
{

  for (kk=0;kk<NEXCCLASS;kk++){

       tmax_exc[kk]=(tau_o_exc[kk]*tau_d_exc[kk]*log(tau_o_exc[kk]/tau_d_exc[kk]))/(tau_o_exc[kk]-tau_d_exc[kk]);
       Imax_exc[kk]=exp(-tmax_exc[kk]/tau_d_exc[kk])-exp(-tmax_exc[kk]/tau_o_exc[kk]);

       /* (par+k)->Aep[kk]=(Amp/Imax_exc[kk])*dt; */

  }

  for (kk=0;kk<NINHCLASS;kk++){

       tmax_inh[kk]=(tau_o_inh[kk]*tau_d_inh[kk]*log(tau_o_inh[kk]/tau_d_inh[kk]))/(tau_o_inh[kk]-tau_d_inh[kk]);
       Imax_inh[kk]=exp(-tmax_inh[kk]/tau_d_inh[kk])-exp(-tmax_inh[kk]/tau_o_inh[kk]);

       /* (par+k)->Aip[kk]=(Amp/Imax_inh[kk])*dt; */

  }

}

(par+k)->V_L  =improve_par("equilibrium potential         V_L = ",(par+k)->V_L,"6.4");
(par+k)->V_Na =improve_par("                             V_Na = ",(par+k)->V_Na,"6.4");
(par+k)->V_K  =improve_par("                             V_K = ",(par+k)->V_K,"6.4");
#ifdef CALCIUM
(par+k)->V_Ca =improve_par("                             V_Ca = ",(par+k)->V_Ca,"6.4");
#endif
(par+k)->G_Na =improve_par("conductance                  G_Na = ",(par+k)->G_Na/dt,"6.4")*dt;
(par+k)->G_L  =improve_par("                              G_L = ",(par+k)->G_L/dt,"6.4")*dt;
(par+k)->G_K_s=improve_par("                              G_K = ",(par+k)->G_K_s/dt,"6.4")*dt;
#ifdef CALCIUM
(par+k)->G_Ca =improve_par("                             G_Ca = ",(par+k)->G_Ca/dt,"6.4")*dt;
(par+k)->G_K_Ca=improve_par("                          G_K(Ca) = ",(par+k)->G_K_Ca/dt,"6.4")*dt;
(par+k)->G_a  =improve_par("                              G_A = ",(par+k)->G_a/dt,"6.4")*dt;
(par+k)->Kp   =improve_par("parameter                      Kp = ",(par+k)->Kp,"6.4");
(par+k)->Kd   =improve_par("                               Kd = ",(par+k)->Kd,"6.4");
(par+k)->Kc   =improve_par("                               Kc = ",(par+k)->Kc,"6.4");
(par+k)->R    =improve_par("                                R = ",(par+k)->R/dt,"6.4")*dt;
(par+k)->tau_x=improve_par("                            tau_x = ",(par+k)->tau_x*dt,"6.4")/dt;
(par+k)->tau_b=improve_par("                            tau_b = ",(par+k)->tau_b*dt,"6.4")/dt;
(par+k)->Utot =improve_par("                             Utot = ",(par+k)->Utot,"6.4");
(par+k)->Utot0=improve_par("                            Utot0 = ",(par+k)->Utot0,"6.4");
(par+k)->Utot1=improve_par("                            Utot1 = ",(par+k)->Utot1,"6.4");
(par+k)->Utot2=improve_par("                            Utot2 = ",(par+k)->Utot2,"6.4");
(par+k)->Utot3=improve_par("                            Utot3 = ",(par+k)->Utot3,"6.4");
(par+k)->Utot4=improve_par("                            Utot4 = ",(par+k)->Utot4,"6.4");
(par+k)->Utot5=improve_par("                            Utot5 = ",(par+k)->Utot5,"6.4");
(par+k)->dr   =improve_par("                               dr = ",(par+k)->dr,"6.5");
(par+k)->Rsh0 =improve_par("                             Rsh0 = ",(par+k)->Rsh0,"6.4");
//(par+k)->Rsh =improve_par("                             Rsh  = ",(par+k)->Rsh,"6.4");
(par+k)->Rsh1=(par+k)->Rsh0+(par+k)->dr;
(par+k)->Rsh2=(par+k)->Rsh1+(par+k)->dr;
(par+k)->Rsh3=(par+k)->Rsh2+(par+k)->dr;
(par+k)->Rsh4=(par+k)->Rsh3+(par+k)->dr;
(par+k)->Rsh5=(par+k)->Rsh4+(par+k)->dr;
(par+k)->Rsh=(par+k)->Rsh5+(par+k)->dr;
/*
for(i=0;i<6;i++)
{
sprintf(buf,"                          Utot[%i] = ",i);
(par+k)->Utot[i]=improve_par(buf,(par+k)->Utot[i],"6.4");
}
for(i=0;i<6;i++)
{
sprintf(buf,"                           Rsh[%i] = ",i);    
(par+k)->Rsh[i]=improve_par(buf,(par+k)->Rsh[i],"6.4");
}
*/
(par+k)->back =improve_par("                             back = ",(par+k)->back/dt,"6.4")*dt;
(par+k)->forw =improve_par("                             forw = ",(par+k)->forw/dt,"6.4")*dt;
(par+k)->PCa  =improve_par("                              PCa = ",(par+k)->PCa/dt,"6.4")*dt;
(par+k)->Cao  =improve_par("                              Cao = ",(par+k)->Cao,"6.4");
(par+k)->gpump=improve_par("                            gpump = ",(par+k)->gpump/dt,"6.4")*dt;
(par+k)->Kpump=improve_par("                            Kpump = ",(par+k)->Kpump,"6.4");
(par+k)->D    =improve_par("                                D = ",(par+k)->D,"12.11");
#endif  /* endif CALCIUM */
(par+k)->i_ext = improve_par("external current amplitude  i_ext = ",(par+k)->i_ext/dt,"6.4")*dt;
(par+k)->Esyn_e=improve_par("EPSP equilibrium potential Esyn_e = ",(par+k)->Esyn_e,"6.4");
(par+k)->Esyn_i=improve_par("IPSP equilibrium potential Esyn_i = ",(par+k)->Esyn_i,"6.4");

for (kk=0;kk<NEXCCLASS;kk++){

      printf("\nEPSP from class #(0- %d) %d\n",NEXCCLASS-1,kk);
      (par+k)->Aep[kk]=improve_par("EPSP amplitude                Aep = ",(par+k)->Aep[kk]/dt,"6.4")*dt;
      printf(" \n");

}

for (kk=0;kk<NINHCLASS;kk++){

      printf("\nIPSP from class #(%d - %d) %d\n",NEXCCLASS,NUMCLASS-1,kk+NEXCCLASS);
      (par+k)->Aip[kk]=improve_par("IPSP amplitude                Aip = ",(par+k)->Aip[kk]/dt,"6.4")*dt;
      printf(" \n");

}

(par+k)->An=improve_par("\nnoise amplitude                An = ",(par+k)->An/dt,"6.4")*dt;

/*for(j=0;j<TAB_SIZE;j++)
{
i=(int)((j-V_INT_OFFSET)/10.0);
robol=(1.+exp(-2.*0.065*((VAR)i+31)));
(par+k)->G_Na_m[j]=(par+k)->G_Na/(robol*robol*robol);
(par+k)->tau[j]=(1./(0.08*(exp(0.055*((VAR)i+35.))+exp(-0.055*((VAR)i+35.)))))/dt;
#ifdef CALCIUM
(par+k)->W_inf[j]=1./(1.+exp(-2.*0.055*((VAR)i+35.)));
(par+k)->X_inf[j]=1./(1.+exp(-2.*2*((VAR)i+45.)));
(par+k)->A_inf[j]=1./(1.+exp(-2.*0.02*((VAR)i+20.)));
(par+k)->B_inf[j]=1./(1.+exp(-2.*-0.1*((VAR)i+70.)));
#endif 
}*/

printf("\n >>class %d parameters :",k+1);
printf("\nsynaptic delay:         del[%d]  %6.2f",k+1,(*((par)->del+k))*dt);

for (kk=0;kk<NEXCCLASS;kk++){

     printf("\nEPSP from class #(0- %d) %d\n",NEXCCLASS-1,kk);
     printf("\nRise Time EPSPs de_o    %6.2f",(par+k)->de_o[kk]);
     printf("\nFall Time EPSPs de_d    %6.2f",(par+k)->de_d[kk]);

}

for (kk=0;kk<NINHCLASS;kk++){

     printf("\nIPSP from class #(%d - %d) %d\n",NEXCCLASS,NUMCLASS-1,kk+NEXCCLASS);
     printf("\n            onset IPSPs di_o    %6.2f",(par+k)->di_o[kk]);
     printf("\n            decay IPSPs di_d    %6.2f",(par+k)->di_d[kk]);

}

printf("\nequilibrium potentials: V_L     %6.2f",(par+k)->V_L);
printf("\n                        V_Na    %6.2f",(par+k)->V_Na);
printf("\n                        V_K     %6.2f",(par+k)->V_K);
#ifdef CALCIUM
printf("\n                        V_Ca    %6.2f",(par+k)->V_Ca);
#endif
printf("\nmembrane conductances:  G_Na    %6.2f",(par+k)->G_Na/dt);
printf("\n                        G_L     %6.2f",(par+k)->G_L/dt);
printf("\n                        G_K_s   %6.2f",(par+k)->G_K_s/dt);
#ifdef CALCIUM
printf("\n                        G_Ca    %6.2f",(par+k)->G_Ca/dt);
printf("\n                        G_K_Ca  %6.2f",(par+k)->G_K_Ca/dt);
printf("\n                        G_a     %6.2f",(par+k)->G_a/dt);
printf("\n                        Kd      %6.2f",(par+k)->Kd);
printf("\n                        Kp      %6.5f",(par+k)->Kp);
printf("\n                        Kc      %6.2f",(par+k)->Kc);
printf("\n                        R       %6.5f",(par+k)->R/dt);
printf("\n                        tau_x   %6.2f",(par+k)->tau_x*dt);
printf("\n                        tau_b   %6.2f",(par+k)->tau_b*dt);
printf("\n                        Utot    %6.2f",(par+k)->Utot);
printf("\n                        back    %6.2f",(par+k)->back/dt);
printf("\n                        forw    %6.2f",(par+k)->forw/dt);
printf("\n                        PCa     %6.2f",(par+k)->PCa/dt);
printf("\n                        Cao     %6.2f",(par+k)->Cao);
printf("\n                      gpump     %6.2f",(par+k)->gpump/dt);
printf("\n                      Kpump     %6.2f",(par+k)->Kpump);
printf("\n                       Utot     %6.2f",(par+k)->Utot);
printf("\n                       Utot0    %6.2f",(par+k)->Utot0);
printf("\n                       Utot1    %6.2f",(par+k)->Utot1);
printf("\n                       Utot2    %6.2f",(par+k)->Utot2);
printf("\n                       Utot3    %6.2f",(par+k)->Utot3);
printf("\n                       Utot4    %6.2f",(par+k)->Utot4);
printf("\n                       Utot5    %6.2f",(par+k)->Utot5);
printf("\n                          dr    %6.5f",(par+k)->dr);
printf("\n                        Rsh0    %6.5f",(par+k)->Rsh0);
printf("\n                        Rsh     %6.5f",(par+k)->Rsh);
printf("\n                          D     %12.11f",(par+k)->D);

#endif
printf("\nI_ext current amlitude: i_ext   %6.2f",(par+k)->i_ext/dt);
printf("\nPSPs equilibrium pot.:  Esyn_e  %6.2f",(par+k)->Esyn_e);
printf("\n                        Esyn_i  %6.2f\n",(par+k)->Esyn_i);

for (kk=0;kk<NEXCCLASS;kk++){

     printf("\nEPSP from class #(0- %d) %d\n",NEXCCLASS-1,kk);
     printf("\nPSPs amplitudes:  EPSPs Aep     %6.4f",(par+k)->Aep[kk]/dt);
     printf(" \n");

}

for (kk=0;kk<NINHCLASS;kk++){

     printf("\nIPSP from class #(%d - %d) %d\n",NEXCCLASS,NUMCLASS-1,kk+NEXCCLASS);
     printf("\n                  IPSPs Aip     %6.4f",(par+k)->Aip[kk]/dt);
     printf(" \n");

}

printf("\n                  noise An      %6.4f",(par+k)->An/dt);
printf("\n\nall parameters OK?");
choose=getchar();
fflush(stdin);
if((choose=='N')||(choose=='n')){ --k; getchar();}
}
/*
koniec pentli ustawiania parametrow po N_class
*/

lseek(fo,0L,SEEK_SET);

for(k=0;k<N_class;k++)
{
rez=write(fo,par+k,sizeof(struct PARAMETERS));
if((rez<0)||((rez!=sizeof(struct PARAMETERS))&&(rez!=0)))
 {
  printf("Error writing %5s file\n",name);
    free(tab_del);
      free(par);
	return(-1);
 }

 rez=write(fo,(par+k)->del,sizeof(SPIKE_INTERVAL)*N_class);
  if((rez<0)||((rez!=(sizeof(SPIKE_INTERVAL)*N_class))&&(rez!=0)))
   {
    printf("Error writing delay tab into %s file\n",name);
      free(tab_del);
	free(par);
	  return(-1);
   }

}
free(tab_del);
free(par);
close(fo);
return 1;
}


VAR improve_par(char *str_out, VAR x, char *lng)
{
 char str_in[20];
 char format[10]; 
 printf("%s",str_out);
 sprintf(format,"%%%s",lng);
 strcat(format,"f; ");
 //printf("%6.4f",x);
 printf(format,x);
 fflush(stdin);
 /* using gets may be dangerous! 
    but here gets() works better than fgets() */ 
 gets(str_in); 
  if( (isdigit(str_in[0]))||(str_in[0]=='.')||(str_in[0]=='-') )
    return((VAR)atof(str_in));
   else 
    return(x);
}

int get_par(char * file_name ,int d) /*wyciaga parametr o nr.d (liczac od pojawienia sie pierwszej cyfry-*/
{             	                     /*  - w zbiorze z *.cfg)  */
FILE *fp;
int i=0,j,k,l=0;
char name[5];
if((fp=fopen(strcat(strcpy(name,file_name),".cfg"),"r"))==NULL)
 {
 printf("Error opening %5s file\n",name);
 return(-1);
 }
while(i==0)
 {
  k=getc(fp);
  i=isdigit(k);
  ++l;
  if(l>70)
   return(0);
 }
if(d<=6)
 {
  for(l=0;l<=d-1;l++)
  fscanf(fp,"%d",&j);
  fclose(fp);
  return(j);
 }
else
 {
  while(i!=',')
    {
     k=getc(fp);
     i=k;
    }
  for(l=0;l<=d-7;l++)
  fscanf(fp,"%d",&j);
  fclose(fp);
  return(j);
 }
}
int get_cfg(char * file_name)
{
 int i;
 FILE *fp;
 extern char konfig[5];
 char konf[5],tym[8],rez=' ';
 long curpos;
 char name[5];
 if((fp=fopen(strcat(strcpy(name,file_name),".cfg"),"r"))==NULL)
  {
   printf("Error opening %5s file\n",name);
   return(-1);
  }
 for(i=0;i<=4;i++)
  konfig[i]=0;
 fscanf(fp,"%8s",tym);
 if(strncmp(&tym[0],"ONE_KIND",8)!=0)fseek(fp,0L,SEEK_SET);
 curpos=ftell(fp);
 fscanf(fp,"%8s",tym);
 if(strncmp(&tym[0],"ONSET",5)!=0)fseek(fp,curpos,SEEK_SET);
 curpos=ftell(fp);
 if(strncmp(&tym[0],"SAME_DELAY",5)!=0)fseek(fp,curpos,SEEK_SET);
 curpos=ftell(fp);
 if(strncmp(&tym[0],"FULL",4)!=0)fseek(fp,curpos,SEEK_SET);
 curpos=ftell(fp);
 if(strncmp(&tym[0],"FAST",4)!=0)fseek(fp,curpos,SEEK_SET);
 rez=getc(fp);
 do
  {
   switch(rez)
    {
     case 'O':
      konf[0]='1';
      break;
     case 'N':
      konf[1]='1';
      break;
     case 'H':
      konf[2]='1';
      break;
     case 'I':
      konf[3]='1';
      break;
     case 'V':
      konf[4]='1';
      break;
     case EOF:
     return(-1);
    }
  }while((rez=getc(fp))!='$');
 for(i=0;i<=4;i++)
 *(konfig+i)=*(konf+i);
 if(fclose(fp)!=EOF)
   return(0);
  else
   return(-1);
}
