/* definitions for interprocessor comm */
/* P.J. Franaszczuk, JHU 2001          */
#define DEB

#include <errno.h>
#include <stdio.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>


#define IPTAB_NAME    ".nodes.tab"
#define CONNTAB_NAME  ".conn.tab"

#define MAX_PACKET_SIZE 1448       /* if MTU is 1500 */
#define PACKETBUFFER    100*MAX_PACKET_SIZE
#define MAX_SPIKE_OUT   PACKETBUFFER/sizeof(short)-2
//#define PORT 1954   /* previous 1954 */



struct CONN
{
    int nconn;      /* no of connections for this processor */
    int *lconn;     /* list of connections  */
    int **buf;    /* pointer to head of list for buffers for each connection */ //was short
};


enum error_flag {INFO,WARNING,PERROR,ERROR,FERROR,DEBUG};
                  
int read_conn(char* argv1);
void print_err(enum error_flag f , char* pname, char *text,char*fname, FILE*fil );
char * gethostip(char *name);
void proc_conn();
void init_path(char* argv0);
char * file_name(char *fname,char* name,char* ext);
int try_accept(int orig_sock, int nconn);
int try_connect(int i, struct sockaddr_in *conn_adr);
void close_socks();

