/*
	MaM : Multiple alignment Manipulator
	
	Implemented by: Can ALKAN & Eray TUZUN
		
	[    calkan@gmail.com   ]
	[  eraytuzun@gmail.com  ]

	Last Update: Oct 18, 2005
	Summary: parsimony score
	Summary: config file, run_prog (Oct 18, 2005)

*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include "main.h"

#define pairwise 0
#define complete 1

int *flag;
float slide(int, char **, int);
float divergence(char *, char *, int);
void SkipNlines(int n,FILE *in);
int exonloc(int);
FILE *gfopen(char *, char *);
void plotout(char *,int,char *,char *,char *);

position *slidepos;
int posno=0;

int SLIDE_WIDTH;
int WINDOW_SIZE;


void slider(char *arg1, char *arg2, char *arg3, char *arg4, char *arg5){
  //argv[1]; alignment file
  //argv[2]; exonlocation file
  //argv[3]; P/C pairwise or complete deletion or parsimony score
  //argv[4]; Slide width 
  //argv[5]; Window size
  
  int iterator; //To get the begin&end locations of exonfile
  int i=0;
  int k=0;
  FILE *unique;
  FILE *repeats;
  FILE *exon; 
 
  float stddif; //error rate
  
  flag = (int *) malloc(maxLen * sizeof (int));

  repeats = gfopen("repeats","w"); //exonfile
  unique  = gfopen("unique","w"); //non-exon file
  
  
  SLIDE_WIDTH=atoi(arg4);
  WINDOW_SIZE=atoi(arg5);
  
  exon=gfopen(arg2,"r"); //to be changed
  slidepos = (position *) malloc(sizeof(position) * arraysize);  
  while(fscanf(exon, "%d", &iterator) >0){
    slidepos[posno].begin=iterator;
    fscanf(exon,"%d",&iterator);
    slidepos[posno++].end=iterator;
  }

  for (i=0;i<maxLen;i++) //Initialize
    flag[i]=1; 
  for (i=0;i<seqTot;i++){
    for (k=0; k<maxLen; k++){
      if ((seqs[i][k]=='-') && !(strcmp(arg3,"C"))) //If complete deletion requested
	flag[k] = 0;
      else if (!strcmp(arg3, "S")) // if parsimony score
	flag[k] = 2;
    }
  }

  for(i=0;i<strlen(seqs[0]);i+=SLIDE_WIDTH){
    printf("\rSliding Windows %f%%",(((float)(i+1)/(float)(strlen(seqs[0])))*100)); 
    stddif = slide(i, seqs, seqTot);
    if (exonloc(i)){
      fprintf(repeats, "%d %f\n", i,stddif);
      fprintf(unique, "%d %f\n", i,0.0);
    }
    else{ 
      fprintf(repeats, "%d %f\n", i,0.0);
      fprintf(unique, "%d %f\n", i,stddif);
    }
  }
  
  fclose(repeats);
  fclose(unique);
  fclose(exon);
  
  plotout("plotfile",0,arg1,arg2,arg3);
  run_prog(PATH_GNUPLOT, "-persist plotfile");
  /*
    system("gnuplot -persist plotfile");
    
   */
}//slider

int exonloc(int position) 
{
  int j;
  for (j=0;j<posno;j++)
    if ((position<slidepos[j].end) && (position>slidepos[j].begin))
      return 1;
  return 0;
} //exonloc

float slide(int startpos, char **seqs, int seqtot){ 
  int i, j;
  float r;
  float totalscore=0.0;
  float avgscore=0.0;
  int a,t,g,c,n;
  int max;
  int whichmax; //n:0, a:1 c:2 g:3 t:4
  int parscore=0;

  if (flag[0] == 2){ // parsimony score
    for (j=startpos;j<startpos+WINDOW_SIZE;j++){
      a = t = g = c = n = 0;
      for (i=0;i<seqtot;i++){
	if (tolower(seqs[i][j]) == '-') n++;
	else if (tolower(seqs[i][j]) == 'a') a++;
	else if (tolower(seqs[i][j]) == 'c') c++;
	else if (tolower(seqs[i][j]) == 'g') g++;
	else if (tolower(seqs[i][j]) == 't') t++;	
      }
      max = n; whichmax=0;
      if (a>max) {max=a; whichmax=1;}
      if (c>max) {max=c; whichmax=2;}
      if (g>max) {max=g; whichmax=3;}
      if (t>max) {max=t; whichmax=4;}
      switch(whichmax){
      case 0:
	parscore=a+c+g+t;
	break;
      case 1:
	parscore=n+c+g+t;
	break;
      case 2:
	parscore=a+n+g+t;
	break;
      case 3:
	parscore=a+c+n+t;
	break;
      case 4:
	parscore=a+c+g+n;
	break;
      }
      totalscore+=(float)parscore/(float)WINDOW_SIZE;
    }
    avgscore=(float)totalscore/(float)seqtot;
  }
  else{
    for (i=0;i<seqtot;i++){
      for (j=0;j<i;j++){
	r=divergence(seqs[i],seqs[j],startpos);
	totalscore+=r;
      } //inner for
    } //outer for
    avgscore=((float)totalscore/(float)((seqtot-1)*seqtot*2));
  }

  return avgscore;  
} // slide

float divergence(char *S, char *T, int startpos){ 
  //Flag=0 pairwise deletion
  //Flag=1 complete deletion if there is a '-' in one of the S[i] we discard the whole i column. 
  //Flag=2 parsimony score
  int i;
  int div=0;
  int length=0;
  float percentage;
  int endpos=strlen(seqs[0]);

  if ((startpos+WINDOW_SIZE)<endpos)
	endpos=startpos+WINDOW_SIZE;

  for (i=startpos;i<endpos;i++){
    if (flag[i] == 1) //Else it is a complete deletion dont count this position in calculations
      {
	// If S[i]=='-' or T[i]=='-' we discard the position.
	if (S[i] != T[i] && S[i]!='-' && T[i]!='-')
	  {
	    div++;
	    length++;
	  }
	if (S[i]==T[i] && S[i]!='-' && T[i]!='-')
	  length++;
      }
  }
  if (div==0)
    percentage=0.0;
  else
    percentage = ((float)div/(float)length);
  return percentage;
} // divergence

void SkipNlines(int n,FILE *in){
  char ch;
  int linecount=0;
  //Skips first n lines
  while(fscanf(in,"%c",&ch) > 0){
    
    if (ch=='\n')
      linecount++;
    if (linecount==n)
      break;
  }
} // skipNlines

FILE *gfopen(char *fname, char *mode){
  //Gracefully file open, gives an error message if it can't open else returns file pointer
  FILE *fp;
  if ((fp=fopen(fname,mode))==0){
    printf("Cannot open %s \n",fname);
    die();        
  }
  return fp;
} // gfopen


void plotout(char *plotfile,int printflag,char *alnname,char *tablename, char *porc){
  FILE *plot;
  char psfile[200];
  
  plot=gfopen(plotfile,"w"); //gnuplot file
  //Plot information comes here
  sprintf(psfile,"%s.%s",alnname,tablename);
  if (porc[0]=='S')
    fprintf(plot,"set title \"DivergenceRate Graph with \'%s\' cDNA file with parsimony score SlideWidth: \'%d\'  WindowWidth: \'%d\'  \" \n",tablename,SLIDE_WIDTH, WINDOW_SIZE);
  else
    fprintf(plot,"set title \"DivergenceRate Graph with \'%s\' cDNA file with \'%s\' deletion SlideWidth: \'%d\'  WindowWidth: \'%d\'  \" \n",tablename,porc,SLIDE_WIDTH, WINDOW_SIZE);
  fprintf(plot,"set xlabel \"LengthofAlignment\" \n");
  fprintf(plot,"set ylabel \"DivergenceRate\" \n");
  fprintf(plot,"plot \"unique\"  with impulses, \"repeats\" with boxes\n");
  fprintf(plot,"set terminal postscript color\n");
  fprintf(plot,"set output \"%s.ps\"\n",psfile);
  fprintf(plot,"replot");
  fclose(plot);
} // plotout



/*

      set palette
       set palette {
                  { gray | color }
                  { gamma <gamma> }
                  {   rgbformulae <r>,<g>,<b>
                    | defined { ( <gray1> <color1> {, <grayN> <colorN>}... ) }
                    | file '<filename>' {datafile-modifiers}
                    | functions <R>,<G>,<B>
                  }
                  { model { RGB | HSV | CMY | YIQ | XYZ } }

set palette model CMY
*/


