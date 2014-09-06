/*
	MaM : Multiple alignment Manipulator
	
	Implemented by: Can ALKAN & Eray TUZUN
		
	[  calkan@gmail.com   ]
	[ eraytuzun@gmail.com ]

	Last Update: Feb 28, 2006

*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <malloc.h>
#include <sys/wait.h>

#define SEQ_LENGTH 150
#define MAXEXONCNT 5000
#define RUNPROGRAM 1000

enum parameters {PARAM_OUT, PARAM_NEXUS, PARAM_MEGA, PARAM_SEPERATE, PARAM_CONCAT};

#define TRUE  1
#define FALSE 0
#define CONFIG_FILE ".mam-config"

typedef int bool;
enum pOptions {CROSS_MATCH, REPEATMASKER, SIM4, TABLE, CONVERT, NONE};

extern char **seqs;
extern char **names;
extern int seqTot;
extern int maxLen;
extern char *ret;
extern int arraysize;
extern void slider(char *, char *, char *, char *, char *);

extern char tablefile[100];
extern bool MAX;
extern bool KEEP;
extern bool VERBOSE;
extern enum pOptions program;
extern bool SINGLE_OUT;
extern bool CONSENSUS_OUT;
extern bool CLUSTALW_OUT;
extern bool NEXUS_OUT;
extern bool MEGA_OUT;
extern bool FASTA_OUT;
extern bool PHYLIP_OUT;
extern bool HTML_OUT;
extern bool IDENTITY_OUT;
extern bool INCLUDE_CONSENSUS;
extern bool SLIDER_OUT;
extern bool CONSENSUS_GAPS;
extern char PC[2];
extern char slidewindow[10];
extern char slidewidth[10]; 


/* repeatfinder.c */

struct _pos{
  int begin;
  int end;
};
typedef struct _pos position;

position *pos;
position *expos;
position *expos2;
position *temparray;

position **filepos;


char PATH_RMASKER[SEQ_LENGTH*2];
char PATH_CMATCH[SEQ_LENGTH*2];
char PATH_SIM4[SEQ_LENGTH*2];
char PATH_GNUPLOT[SEQ_LENGTH*2];


extern FILE *gfopen(char *, char *);
void printoutput(position *,char *);
int max(int,int);
int min(int,int);
char overlap(position,position);
int mergearrays(position *,position *,position *,int,int);
void updatelocations(char *,position *);
int countspace(char *);
int prune(position [],position [],int);
void copyarray(position *,position *,int);
void extractaln(char *inputfile);
void tosstable(position [],position []);
void crossmatch(int,char *);
void repeatmasker(int,char *);
void simfour(int,char *);
int findstring(char *);
extern void SkipNlines(int,FILE *);
int parsetablearray(int index,position **);
void readtablefile(position **, char *);
void allToggles(bool);
void inputToggles();
int  outputToggles();
void getFileInput();
bool checkToggles(bool);
void defaults();
void parseParams(int, char **);
void help(char *, int);

/* seqtools.c */

char *computeCons(int, char **);
int readAln(char *);
int readClustal(FILE *);
int readNexus(FILE *);
int readMega(FILE *);
int readFasta(FILE *);
int readPhylip(FILE *);
//FILE *readSFile();
enum subseqSelection sMenu();
enum selection menu();
enum toggleSelection toggleMenu();
void subseqOut(int, char **, int [], int [], int, char *);
void singleSubseq(int, char **, int [], int [], int, char *);
void printOut(int, char **, char *);
void printNexus(int, char **, char *);
void printMega(int, char **, char *);
void printFasta(int, char **, char *);
void printPhylip(int, char **, char *);
void printHTML(int, char **, char *);
int myisalnum(char);
void processOutOptions(char [], char [], char [], int);
void convertFormat(char []);
void alnStats(char []);
int sequence_len(char *);
float gc_cont(char *);
void die(void);
void loadConfig(void);
void createConfig(char *);
void run_prog(char *, char *);
