#define DELTA_T  0.00001
/* define trzeba zmieniac zaleznie od net.cfg */
#define OUT                0x1		/* bit flagi output z neuronu */
#define NOISE		   0x2          /* bit dodawania szumu */
#define HISTO              0x4          /* bit opcji histogramu */
#define IEXT               0x8          /* bit dodawania pradu zewnetrznego */
#define VOUT		   0x10		/* bit opcji outputu V */
#define BYL_SPIKE          0x20         /* bit uzywany w calkuj */
#define NEURON_KIND        0x40  
#define SPIKE_DET          0x80         /* for links */  
#define POSPOL             0X1          /* neuron affected by positive polarity stim */
#define NEGPOL             0X2          /* neuron affected by negative polarity stim */

#define O_BINARY 0 
/*#define ONE_KIND*/
#define ONSET
/*#define SAME_DELAY*/
/*#define NOISEOUT*/
/*#define FULL*/
/*#define FAST*/
#define CALCIUM
/*#define VOUT_DOS      do zamiana typu int w UNIX na format DOS*/

#define EPS 1.E-10

  /* definicje struktur i stalych w modelu sieci neuronow */

#include <limits.h>



typedef unsigned char  byte;


#define True  		1
#define False 		0
#define AP_tr        	10.       /* prog do sprawdzania czy jest spike */
#define V_INT_OFFSET	1024      /* wartosc integer dodawana do V dla otrzymania indeksu 0 */
#define TAB_SIZE        2048       /* rozmiar tablic w parametrach neuronow */

#define NEXCCLASS       4
#define NINHCLASS       3
#define NUMCLASS        7
#define NOISECLASS      0


#define EMPTY_LIST	   UCHAR_MAX /* was UCHAR_MAX */
#define ONE_SPIKE          (EMPTY_LIST-1)
#define MAX_INTERVAL	   USHRT_MAX
#define SPIKE_LIST_LENGTH  20

 #if SPIKE_LIST_LENGTH >= ONE_SPIKE
 #error SPIKE_LIST LENGTH must be < ONE_SPIKE
 #endif

#define MAX_NEURON  USHRT_MAX      /* change if type of NEURON_INDEX changes */
#define MAX_KIND    UCHAR_MAX      /* change if type of PARAMETER_INDEX changes */

typedef signed char  	WEIGHT;
typedef ushort  	DELAY;  //was byte, can try ushort
typedef byte  	PARAMETER_INDEX;
typedef ushort 	SYNAPSE_NUMBER;
typedef ushort 	NEURON_INDEX;
typedef byte	SPIKE_INDEX;
typedef ushort  SPIKE_INTERVAL;
typedef double   VAR;
typedef byte    HIS;

typedef int   LINK_WEIGHT;
typedef int   LINK_INTERVAL;

#define MAX_LINK_INTERVAL  INT_MAX
#define LINK_LIST_LENGTH   200  /* was 20 */

 #if LINK_LIST_LENGTH >= ONE_SPIKE
 #error LINK_LIST LENGTH must be < ONE_SPIKE
 #endif

struct LINK
{
        NEURON_INDEX    neuron_post;  /* nr neuronu post-synaptycznego */
        NEURON_INDEX    preclass_num; /* presynaptic class */
        LINK_WEIGHT     weight;       /* weight */        
        LINK_INTERVAL   delay;        /* link delay */
	LINK_INTERVAL   interval[ LINK_LIST_LENGTH ];
        SPIKE_INDEX     first,last;    /* was SPIKE_INDEX */       
          
};


struct SYNAPSE
{
	WEIGHT 		weight;

	DELAY  		delay;


	NEURON_INDEX    neuron;  /* nr neuronu postsynaptycznego */

        NEURON_INDEX    preclass_num; /* presynaptic class */

};


struct NEURON
{
	byte		flags;   /* miejsce na opcje bitowe np. output */
        byte            polarity;/* byte to store polarity stimulation information */
	PARAMETER_INDEX param;   /* indeks do tablicy parametrow klas neuronow */
	SYNAPSE_NUMBER	syn_num; /* liczba synaps na wyjsciu */
	VAR		Vi_o[NINHCLASS]; /* skladowa onset potenjalu poch. z inhibicji */
	VAR 		Ve_o[NEXCCLASS]; /* skladowa onset potencjalu poch. z ekscytacji */
	VAR		Vi_d[NINHCLASS]; /* skladowa decay potenjalu poch. z inhibicji */
	VAR 		Ve_d[NEXCCLASS]; /* skladowa decay potencjalu poch. z ekscytacji */

	VAR		W;  /* zmienna pomocnicza */
	VAR		V;  /* napiecie blony */
        VAR             C[7];
//  	VAR             C;  /* intercellular koncentracja jonow wapnia */
//        VAR             C0;
//        VAR             C1;
//        VAR             C2;
//        VAR             C3;
//        VAR             C4;
//        VAR             C5;
        VAR             X;  /* calcium activation variable */
        VAR             B;
        VAR             U[7];
//        VAR             U;
//        VAR             U0;
//        VAR             U1;
//        VAR             U2;
//        VAR             U3;
//        VAR             U4;
//        VAR             U5;
        SPIKE_INDEX     first,last;    /* wskazniki do cyklicznej listy spike'ow */
	SPIKE_INTERVAL  interval;/* czas od ostatniego spike'u */
	SPIKE_INTERVAL  sum_interval; /* suma intervalow z listy czyli czas od pierwszego */
	SPIKE_INTERVAL		spikes[ SPIKE_LIST_LENGTH ];
};


struct PARAMETERS
{
        SPIKE_INTERVAL *del;
	VAR     Aip[NINHCLASS];    /* amplituda potencjalu IPSP w t=0 (bez wagi) */
	VAR     Aep[NEXCCLASS];    /* amplituda potencjalu EPSP w t=0 (bez wagi) */
	VAR     An;     /* amplituda dodawanego szumu */
	VAR 	di_o[NINHCLASS];     /* decrement onset potencjalu IPSP */
	VAR	de_o[NEXCCLASS];     /* decrement onset potencjalu EPSP */
	VAR 	di_d[NINHCLASS];     /* decrement decay potencjalu IPSP */
	VAR	de_d[NEXCCLASS];     /* decrement decay potencjalu EPSP jonow */

		       /* reverse potentials */
	VAR	V_L;
	VAR	V_Na;
	VAR	V_K;
			/* przewodnictwa blony */
	VAR     G_Na;
        VAR 	G_L;
	VAR	G_K_s;  /* = G_K / (s**wp) */
	VAR     Esyn_i; /* potencjal odwrotny hamujacy */
	VAR     Esyn_e; /* potencjal odwrotny pobudzajacy */
	VAR     i_ext;  /* prad zewnetrzny */
        VAR     V_Ca;
        VAR     G_K_Ca; /* Potas aktywowany wapniem */ 
        VAR     G_Ca;   /* przewodnictwo wapnia */
        VAR     G_a;
        VAR     Kd;
        VAR     Kp;
        VAR     Kc;
        VAR     R; 
        VAR     tau_x;
        VAR     tau_b;
        VAR     forw;
        VAR     back;
        VAR     gpump;
        VAR     Kpump;
        VAR     Cao;
        VAR     Utot;
        VAR     PCa;

        VAR     D;       /* diffusion coeficient */
        VAR     Rsh0; /* array of shell radii */
        VAR     Rsh1;
        VAR     Rsh2;
        VAR     Rsh3;
        VAR     Rsh4;
        VAR     Rsh5;
        VAR     Rsh;
        VAR     dr;      /* shell thicknes  */
        VAR     Utot0; /* Utot array */
        VAR     Utot1;
        VAR     Utot2;
        VAR     Utot3;
        VAR     Utot4;
        VAR     Utot5;   
} ;


/* deklaracje funkcji */
void set_flags();
void read_params(char*);
void read_synapses(char*);
void read_neurons(char*,int);
void write_neurons(char*);
void read_links(char*,int);
void write_links(char*); 
int calkuj( struct NEURON  *neuron_nr);
void szum( struct NEURON *neuron_nr);
int save_stimulus(void);
void close_all(int sig);

void init();


int stimulus();

void link_update();
int sszalloc(int s);
int stszalloc(int s);
int stim_freq(struct NEURON *neuron_nr);

#define V_0 -60.02209
#define W_0 0.05994952
#define X_0 0.002450918
#define B_0 0.1196677

#define C_0 2.082966e-05
#define U_0 0.

//#define C0_0 0.000018
//#define C1_0 0.000018
//#define C2_0 0.000018
//#define C3_0 0.000018
//#define C4_0 0.000018
//#define C5_0 0.000018

//#define U0_0 0.
//#define U1_0 0.
//#define U2_0 0.
//#define U3_0 0.
//#define U4_0 0.
//#define U5_0 0.



#define W_inf(V) (1./(1.+exp(-2.*0.055*(V+35.))))
#define X_inf(V) (1./(1.+exp(-2.*0.2*(V+45.))))
#define A_inf(V) (1./(1.+exp(-2.*0.02*(V+20.))))
#define B_inf(V) (1./(1.+exp(-2.*-0.095*(V+70.))))

VAR G_Na_m(VAR V,struct PARAMETERS *par)
{
register VAR robol;
robol=(1.+exp(-2.*0.065*(V+31.)));
return( par->G_Na/(robol*robol*robol) );
}

#define tau(V) ((1./(0.08*2*cosh(0.055*(V+35.))))/(DELTA_T*1000))


//#define omc(C) (0.25+exp(80*(C-0.55))/(1+exp(80*(C-0.55)))-0.25*exp(80*(C-0.35))/(1+exp(80*(C-0.35))))
//#define om(C) omc(C) 
//#define om(C) ((omc(C) < 0.5)?omc(C):0.5)
//#define om(C)  (0.25+exp(80*(C-0.55))/(1+exp(80*(C-0.55)))-0.25*exp(80*(C-0.35))/(1+exp(80*(C-0.35)))-0.5*(exp(80*(i-0.55))/(1+exp(80*(i-0.55)))) ) 
//#define taul(C) ((1000+1/(0.0001+C*C*C/100))/(DELTA_T*1000))
//#define taul(C) ((1000+1/(0.0001+C*C*C/100))/0.01)

//#define B_NMDA(V) (1/(1+exp(-0.062*V)*0.28011))
//#define V_NMDA(C) (13.319652*log(2/C))
//#define V_NMDA(C)   130
//#define V_NMDA(C) 0

struct SYNAPSESD
{
	VAR	w;
	VAR	V_d;
	VAR	V_o;
	VAR	C;
};

struct LINKD
{
        VAR     w;
        VAR     V_d;
        VAR     V_o;
        VAR     C;
};
