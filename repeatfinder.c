/*
	MaM : Multiple alignment Manipulator
	
	Implemented by: Can ALKAN & Eray TUZUN
		
	[   calkan@gmail.com   ]
	[  eraytuzun@gmail.com ]

	Last Update: March 20, 2006
	Summary: parsimony score (Mar 8, 2005)
	Summary: tablefile check moved from processOutOptions to readtablefile (Aug 10, 2005)
	Summary: Length information printed (Aug 25, 2005)
	Summary: PHYLIP support added (Oct 7, 2005)
	Summary: Convert option (Oct 10, 2005)
	Summary: -noupdate option (Oct 10, 2005)
	Summary: -V option (verbose) (Oct 10, 2005)
	Summary: -alnstats option (Oct 10, 2005)
	Summary: -rmasker_opts, -cmatch_opts, -sim4_opts (Oct 14, 2005)
	Summary: config file, run_prog (Oct 18, 2005)
	Summary: HTML support added (Feb 28, 2006)
	Summary: Consensus Gaps On/Off added (March 20, 2006)
	Summary: -noupdate option converted to -update=on/off
*/

char **seqs;
char **names;


#include "main.h"

int *length;  

int MaxExons=0;
int MaxFileNo=1;
char inFile[100];
char tablefile[100];
bool MAX = TRUE;
bool KEEP = TRUE;
enum pOptions program = REPEATMASKER;
bool SINGLE_OUT              = TRUE;
bool CONSENSUS_OUT           = FALSE;
bool CLUSTALW_OUT            = TRUE;
bool NEXUS_OUT               = FALSE;
bool MEGA_OUT                = FALSE;
bool FASTA_OUT               = FALSE;
bool PHYLIP_OUT              = FALSE;
bool HTML_OUT                = FALSE;
bool IDENTITY_OUT            = FALSE;
bool INCLUDE_CONSENSUS       = FALSE;
bool SLIDER_OUT              = FALSE;
bool NOUPDATE                = FALSE;
bool VERBOSE                 = FALSE;
bool ALNSTATS                = FALSE;
bool CONSENSUS_GAPS          = TRUE;
char PC[2];
char slidewindow[10];
char slidewidth[10];
int  seqTot;
int arraysize=2;

char RMASKER_OPTS[SEQ_LENGTH*2];
char CMATCH_OPTS[SEQ_LENGTH*2];
char SIM4_OPTS[SEQ_LENGTH*2];


int main(int argc, char **argv){
  char alignfile[50]; //Sequence1.aln
  char outputname[SEQ_LENGTH]; //format: argv[1].argv[2].out
  int i, cnt; //for loop variable
  int result=0;
  int maxLen;
  int diff=0;
  char consFile[SEQ_LENGTH];
  char *consensus;
  FILE *out;

  strcpy(RMASKER_OPTS, "-no_is -nolow");
  strcpy(CMATCH_OPTS, "");
  strcpy(SIM4_OPTS, "");

  loadConfig();


  if (argc == 1)
    allToggles(FALSE); // get the input file and show the menu
  else if (argc == 2 && argv[1][0]!='-'){
    strcpy(inFile, argv[1]);
    allToggles(TRUE);    // the input file is given, just show the menu
  }
  else
    parseParams(argc, argv); // parameter driven version


  if (ALNSTATS)
    alnStats(inFile);

  if (CONSENSUS_OUT){
    consensus = (char *)malloc((strlen(seqs[0])+1) * sizeof(char));
    strcpy(outputname, argv[1]);
    for (i=strlen(outputname)-1;i>=0;i--)
      if (outputname[i]=='.'){
	outputname[i]=0;
	break;
      }
    sprintf(consFile, "%s.consensus.fa", outputname);
    out = fopen(consFile, "w");
    fprintf(out, ">%s Consensus\n", argv[1]);
    strcpy(consensus, computeCons(seqTot, seqs));
    cnt = 0;
    while (cnt < strlen(consensus)){
      fprintf(out, "%c", consensus[cnt++]);
      if (cnt % 60 == 0)
	fprintf(out, "\n");
    } // while
    fprintf(out, "\n");
    fclose(out);
    free(consensus);
  } // if consensus_out      

  
  switch(program){
  case REPEATMASKER:
    if (PATH_RMASKER[0] == 0){
      fprintf(stderr, "RepeatMasker path is not defined in configuration file.\nCorrect it manually or remove $HOME/.mam-config file, install \"which\" system tool and run MaM again.\n");
      exit(0);
    }
    break;
  case CROSS_MATCH:
    if (PATH_CMATCH[0] == 0){
      fprintf(stderr, "cross_match path is not defined in configuration file.\nCorrect it manually or remove $HOME/.mam-config file, install \"which\" system tool and run MaM again.\n");
      exit(0);
    }
    break;
  case SIM4:
    if (PATH_SIM4[0] == 0){
      fprintf(stderr, "sim4 path is not defined in configuration file.\nCorrect it manually or remove $HOME/.mam-config file, install \"which\" system tool and run MaM again.\n");
      exit(0);
    }
    break;
  default:
    break;
  }

  if (SLIDER_OUT && PATH_GNUPLOT[0] == 0){
      fprintf(stderr, "gnuplot path is not defined in configuration file.\nCorrect it manually or remove $HOME/.mam-config file, install \"which\" system tool and run MaM again.\n");
      exit(0);
  }


  sprintf(outputname,"%s.%s.out",inFile, tablefile);

  if (program != CONVERT && program != NONE)
    extractaln(inFile); 
  
  maxLen = strlen(seqs[0]);

  length = (int *) malloc(sizeof(int) * (seqTot + 1));
  filepos = (position **) malloc(sizeof(position *) * (seqTot+1));

  for (i=1;i<=seqTot;i++){
    length[i]=0; 
    filepos[i] = (position *) malloc(sizeof(position) * (length[i]+1));   
  }
  if (program == TABLE)            // TABLE option
    readtablefile(filepos, tablefile);
  else if (program == CONVERT){ // convert format
    convertFormat(inFile);
    return 1;
  }
  
  expos = (position *) malloc(sizeof(position));  
  expos2 = (position *) malloc(sizeof(position));  
  temparray = (position *) malloc(sizeof(position));  
 
  for (i=1;i<=seqTot;i++)              
    {
      switch (program){
      case CROSS_MATCH:
	crossmatch(i, tablefile);         
	break;
      case REPEATMASKER:
	repeatmasker(i, tablefile);
	break;
      case SIM4:
	simfour(i,tablefile);   
	break;
      case TABLE:
	MaxExons=parsetablearray(i,filepos);
	break;
      case CONVERT:
	// no way it can come here
	break;
      case NONE:
	if (ALNSTATS)
	  printf("Alignment statistics are dumped to %s.stats\n", inFile); 
	if (CONSENSUS_OUT)
	  fprintf(stderr, "Consensus sequence is written to file %s.\n", consFile);
	exit(0);
	break;
      } // switch

      

      /*      printf("DEBUG: MAXEXONS: %d\n", MaxExons); */
      if (VERBOSE)
	printf("Updating the locations according to .aln file: %s\n", names[i-1]);
      else{
	fprintf(stderr, "\rProcessing %f%%", ((float)(i)/(float)(seqTot))*100);
	fflush(stderr);
      }
      sprintf(alignfile,"%s.aln",names[i-1]);
      //      printf("arraysize: %d\n", arraysize);
      
      diff=MaxExons-diff;
      if (diff>0)
      	arraysize=(arraysize+diff)*1.5;
      diff=MaxExons;

      expos = (position *) realloc(expos, sizeof(position) * arraysize);
      expos2 = (position *) realloc(expos2, sizeof(position) * arraysize);
      temparray = (position *) realloc(temparray, sizeof(position) * arraysize);
    
      if (expos==NULL||expos2==NULL||temparray==NULL){
	printf("null array bug arraysize:%d diff:%d\n", arraysize,diff);
	printf("name: %s, size %d\n", names[i-1], strlen(seqs[i-1]));
	return 0;
      }
      if (!NOUPDATE && program == TABLE){
	updatelocations(alignfile, pos);
	if (VERBOSE)
	  printf("Locations are updated\n");
      }
      
      MaxExons=prune(pos,expos2,MaxExons);

      result=mergearrays(expos,expos2,temparray,result,MaxExons);

      result=prune(temparray,expos,result);

    }//for
  result=prune(expos,expos2,result); 
  expos2[0].begin=result; 
  expos2[0].end=maxLen;             
  
  
  if (KEEP){
    printoutput(expos2,outputname); 
    processOutOptions(inFile, outputname, tablefile, result);
  }
  else
    {
      tosstable(expos2,expos);
      printoutput(expos,outputname); 
      processOutOptions(inFile, outputname, tablefile, result+1);
    }

  if (ALNSTATS)
    printf("Alignment statistics are dumped to %s.stats\n", inFile); 

  if (CONSENSUS_OUT)
    fprintf(stderr, "Consensus sequence written to file %s.\n", consFile);

  die();
  return 1;
}

void updatelocations(char *file1,position *itr)
{ 
  //File1 is the aligned file with spaces. This function updates the pos array according to number of spaces
  /* ADD POSITION CHECK, THERE MAY BE ERROR IN USER DEFINED TABLEFILE */

  FILE *in=gfopen(file1,"r");

  char ch;
  int noofgaps[maxLen];
  int i=1;
  int gapcount=0;


  SkipNlines(1,in);

  while(fscanf(in,"%c",&ch) > 0)
    {
      if (ch!='\n')
	{	
          if (ch=='-')
	    noofgaps[i]=++gapcount;
	  else
	    noofgaps[i++]=gapcount;
        }
    }

  fclose(in);
  for (i=1;i<=MaxExons;i++)
    {
      itr[i].begin=itr[i].begin+noofgaps[itr[i].begin];
      itr[i].end=itr[i].end+noofgaps[itr[i].end];
      
    }

}

void printoutput(position itr[],char *outputname)
{
  //Index: Max number of repeat locations
  //Writes the result to "output" 
  FILE *out=gfopen(outputname,"w");
  int i=1;
  int Index=itr[0].begin;
  //int Index=MaxExons;

  for (i=1;i<=Index;i++)
    {
      if (VERBOSE)
	printf("%d st exon's begin pos is %d end position is %d  \n",i,itr[i].begin,itr[i].end);
      fprintf(out,"%d ",itr[i].begin);
      fprintf(out,"%d\n",itr[i].end);      
    }
   
  fclose(out);
}
void extractaln(char *inputfile)
{
  //Input *.aln file
  //Output: n *.aln files S1.aln S2.aln ...Sn.aln

  FILE *out;
  FILE *outaln;
  int i, j, cnt;
  char fname[50];

  for (i = 0 ; i < seqTot; i++){
    sprintf(fname, "%s", names[i]);
    if ((out = fopen(fname, "w")) == NULL){
      printf("Unable to open file %s\n", fname);
      die();
    }
    
    // align file
    sprintf(fname, "%s.aln", names[i]);
    outaln = fopen(fname, "w");

    fprintf(out, ">%s\n", names[i]);
    fprintf(outaln, ">%s\n", names[i]); // align file

    cnt = 0; 
    for (j=0;j<maxLen;j++){
      if (seqs[i][j] != '-'){
	fprintf(out, "%c", seqs[i][j]);
      }
      
      cnt++;
      fprintf(outaln, "%c", seqs[i][j]);

      if (cnt % 60 == 0){
	fprintf(out, "\n");
	fprintf(outaln, "\n");
      }
    }

    fprintf(out, "\n");
    fclose(out);
    fprintf(outaln, "\n");
    fclose(outaln);
    if (VERBOSE){
      printf("File %s extracted\n", names[i]);
      printf("File %s.aln extracted\n", names[i]);
    }
  } // for
} // extractaln

void copyarray(position X[],position Y[],int m)
{
  int i;
  for (i=1;i<=m;i++)
    {
      Y[i].begin=X[i].begin;
      Y[i].end=X[i].end;
    }
} // copyarray

int prune(position A[],position B[], int m)
{
  int i=1;
  int j=1;
  int maxint;

  while (i<=m)
    {
      B[j].begin=A[i].begin;
      maxint=A[i].end;
      if (i!=m)
	while (((maxint)+2)>=(A[i+1].begin))
	  {
	    i++;
	    maxint=max(maxint,A[i].end);
	    if (i==m)
	      break; //Ugly fixxxxxxxxxxxx
	  }
      B[j++].end=maxint;
      i++;
    }
  return j-1;
} // prune

int mergearrays(position A[],position B[],position C[],int m, int n) 
{
  int i=1;
  int j=1;
  int k=1;
  while (i<=m && j<=n)
    {
      if (overlap(A[i],B[j])=='X')
	{
	  if (MAX)
	    {
	      C[k].begin=min(A[i].begin,B[j].begin);
	      C[k++].end=max(A[i++].end,B[j++].end);  
	    }
	  else  // if MIN
	    {
	      C[k].begin=max(A[i].begin,B[j].begin);
	      C[k++].end=min(A[i++].end,B[j++].end); 
	      if (C[k-1].begin==C[k-1].end)
		k--;    
	    }
	} 
      else
	{
	  if (overlap(A[i],B[j])=='A')
	    {
	      C[k].begin=A[i].begin; 
	      C[k++].end=A[i++].end; 
	    }
	  else
	    {
	      C[k].begin=B[j].begin;
	      C[k++].end=B[j++].end; 
	    }
	}
    }
  while (i<=m)
    {
      C[k].begin=A[i].begin;
      C[k++].end=A[i++].end; 
    }
  while(j<=n)
    {
      C[k].begin=B[j].begin;
      C[k++].end=B[j++].end; 
    }
  return k-1;	
} // mergearrays

char overlap(position x, position y)
{

  if ((x.end)<((y.begin)-2))
    return 'A';
  if ((y.end)<((x.begin)-2))
    return 'B';
  return 'X';
} // overlap

int min(int x,int y)
{
  if (x<y)
    return x;
  else
    return y;
}  // min

int max(int x,int y)
{
  if (x>y)
    return x;
  else
    return y;
} // max

int findstring(char *input)
{
  int i;
 
  for (i=1;i<=seqTot;i++)
    if (!strstr(input,names[i-1])==0)
      return i;
  return 0; //No match found
} // findstring


void tosstable(position postable[],position tosstable[])
{
//boundary checking added for first and last positions
  int i=1; //postable index
  int j=1; //tosstable index

  int noofentries=postable[0].begin;
  int maxlength=postable[0].end;
 
  if( postable[i].begin>2)
  {	
  	tosstable[j].begin=1; //0 or 1
  	tosstable[j++].end=postable[1].begin-1;
  }
	
  for (i=2;i<=noofentries;i++)
    {
      tosstable[j].begin=postable[i-1].end+1;
      tosstable[j++].end=postable[i].begin-1;
    }
 
  if (postable[i].end<(maxlength-2))
  { 
	 tosstable[j].begin=postable[i].end+1; //was i-1 before! bug
  	 tosstable[j].end=maxlength;
  }
  
  tosstable[0].begin=j-1; //Number of entries in the table
  tosstable[0].end=maxlength; //Max length
} //tosstable


void crossmatch(int i,char *exonname)
{
  FILE *in;
  char line[50]; //word
  int j; 
  int count=0;
  char outfile[50]; //name of the file
  char command[100];  //System command
  int linecnt=0;
  char ch;

  sprintf(outfile,"out%d",i);
  MaxExons=0;
  
  if (strcmp(exonname,"default"))
    /*
      sprintf(command,"cross_match %s %s>%s",names[i-1],exonname,outfile);
    */
    sprintf(command,"%s %s %s>%s",CMATCH_OPTS, names[i-1],exonname,outfile);
  
  else
    sprintf(command,"%s %s > %s",CMATCH_OPTS, names[i-1], outfile);
  /*
    sprintf(command,"cross_match %s > %s",names[i-1],outfile);
  */
  
  run_prog(PATH_CMATCH, command);
  /*
    system(command); //Call cross_match
  */
  
  if (VERBOSE)
    printf("Crossmatch finished.. Now saving the beginning and end locations\n");
   
  in=gfopen(outfile,"r");
  if (VERBOSE)
    printf("Opening file > %s\n",outfile);

  while(fscanf(in,"%c",&ch) > 0){
    if (ch == '\n' || ch=='\r')
      linecnt++;
  }
  
  rewind(in);
  
  pos = (position *) malloc(sizeof(position) * linecnt);
   
  //Preprocessing step finds the number of sequences and their names 
  while(fscanf(in,"%s",line) > 0)
    {
      if (!strstr(line,names[i-1])==0)
	{
	  count++;
	  if (count>4)
	    {
		  
	      fscanf(in,"%s",line);
	      if ((strcmp(line,"T")==0)||(strcmp(line,"C")==0)||(strcmp(line,"A")==0)||(strcmp(line,"G")==0))
		continue;
	      pos[++MaxExons].begin=atoi(line);
	      fscanf(in,"%s",line);
	      pos[MaxExons].end=atoi(line);
	      for(j=1;j<5;j++) 
		fscanf(in,"%s",line);	       
	    }//if2
	}//if1

    }//while
 
  fclose(in); //Close the out file
  //Cleaning stuff
   
  sprintf(command,"rm -f %s.log",names[i-1]);
  system(command);
  sprintf(command,"rm -f %s",outfile);
  system(command); 
} // crossmatch

void repeatmasker(int i,char *exonname)
{
  FILE *in;
  char line[250];//word
  char outfile[50];//name of the file
  char command[100]; //System command
  int linecnt=0;
  int j;
  char ch;

  if (VERBOSE)
    printf("Repeatmasker processing %s\n", names[i-1]);

  sprintf(outfile,"%s.out",names[i-1]);
  MaxExons=0;
  if (strcmp(exonname,"default")){  // with a cDNA file input
    if (VERBOSE)
      sprintf(command,"%s -lib %s %s", RMASKER_OPTS, exonname,names[i-1]);    
    /*
      sprintf(command,"repeatmasker  -no_is -nolow -lib %s %s",exonname,names[i-1]);    
     */
    else
      sprintf(command,"%s -lib %s %s >/dev/null",RMASKER_OPTS, exonname,names[i-1]);    
    /*
      sprintf(command,"repeatmasker  -no_is -nolow -lib %s %s >/dev/null",exonname,names[i-1]);    
     */
  }
  else{   // default run
    if (VERBOSE)
      sprintf(command,"%s %s",RMASKER_OPTS, names[i-1]); 
    /*
      sprintf(command,"repeatmasker  -no_is -nolow %s",names[i-1]); 
    */
    else
      sprintf(command,"%s %s >/dev/null", RMASKER_OPTS, names[i-1]); 
    /*
      sprintf(command,"repeatmasker  -no_is -nolow %s >/dev/null",names[i-1]); 
    */
  }
  
  run_prog(PATH_RMASKER, command);
  
  /*
    system(command); //Call repeat_masker
  */

  if (VERBOSE)
    printf("RepeatMasker finished.. Now saving the beginning and end locations\n");
   
  in=gfopen(outfile,"r");
  if (VERBOSE)
    printf("Opening file > %s\n",outfile);
 
  /* try first */

  fgets(line, 250, in);
  
  if (strstr(line, "There were no repetitive sequences detected")){
    fprintf(stderr, "RepeatMasker returned null result for repetitions.\n");
    sprintf(command,"rm -f %s.out ",names[i-1]);
    system(command); 
    sprintf(command,"rm -f %s.cat ",names[i-1]);
    system(command); 
    sprintf(command,"rm -f %s.stderr ",names[i-1]);
    system(command); 
    sprintf(command,"rm -f %s.tbl ",names[i-1]);
    system(command); 
    sprintf(command,"rm -f %s.masked ",names[i-1]);
    system(command); 
    sprintf(command,"rm -f %s.log ",names[i-1]);
    system(command);
    sprintf(command,"rm -rf RM_* ");
    system(command);
    die();
  }
  
  rewind(in);

  while(fscanf(in,"%c",&ch) > 0){
    if (ch == '\n' || ch=='\r')
      linecnt++;
  }

  rewind(in);
  
  pos = (position *) malloc(sizeof(position) * linecnt);

  //Preprocessing step finds the number of sequences and their names 

  /*
  while(fscanf(in,"%s",line) > 0)
    {
      if (!strstr(line,"There")==0)
	break;            //If no repetitive sequences found it starts with there
      if (!strstr(line,names[i-1])==0)
	{
	  printf("DEBUG ::  LINE : %s\n", line);
  	  fscanf(in,"%s",line);
	  pos[++MaxExons].begin=atoi(line);
	  fscanf(in,"%s",line);
	  pos[MaxExons].end=atoi(line);	       
	}//if1
    }//while
  */

  /*
    update on May 6, 2004:
    the above parsing may cause problems for certain file names
    the below version should be safer
  */

  /* get first line */
  fgets(line, 250, in);
  if (!strstr(line,"There")){
    /* If no repetitive sequences found it starts with there */
    /* otherwise, it should be safe */
    /* get second line */
    fgets(line, 250, in);
    /* pass third empty line */
    //fgets(line, 250, in);
    while (fscanf(in, "%s", line) > 0){
      for (j=0;j<4;j++){
	if (!(fscanf(in, "%s", line) > 0))
	  break ; /* first five entries are not important */
	
      }
      /*      printf("DEBUG ::  LINE : %s\n", line); */
      fscanf(in, "%s", line); /* here is begin location */
      /*      printf("DEBUG ::  LINE2 : %s\n", line); */
      pos[++MaxExons].begin=atoi(line);
      fscanf(in, "%s", line); /* here is end location */
      /*      printf("DEBUG ::  LINE3 : %s\n", line); */
      pos[MaxExons].end=atoi(line);	           
      /* forget the rest of the line */
      fgets(line, 250, in);
    }
  }
  
  /*
    printf("DEBUG: NAME : %s\n", names[i-1]);
    printf("DEBUG: MAXEXONS: %d\n", MaxExons);
  */

  //Cleaning stuff
  
  sprintf(command,"rm -f %s.out ",names[i-1]);
  system(command); 
  sprintf(command,"rm -f %s.cat ",names[i-1]);
  system(command); 
  sprintf(command,"rm -f %s.stderr ",names[i-1]);
  system(command); 
  sprintf(command,"rm -f %s.tbl ",names[i-1]);
  system(command); 
  sprintf(command,"rm -f %s.masked ",names[i-1]);
  system(command); 
  sprintf(command,"rm -f %s.log ",names[i-1]);
  system(command);
  sprintf(command,"rm -rf RM_* ");
  system(command);
  
  fclose(in); //Close the out file
} // repeatmasker

void simfour(int i,char *exonname)
{
  FILE *in;
  char outfile[50];//name of the file
  char command[5*SEQ_LENGTH]; //System command
  int iterator;
  int linecnt=0;
  char ch;

  sprintf(outfile,"out%d",i);
  MaxExons=0;

  sprintf(command, "%s %s %s A=5>%s", exonname,names[i-1], SIM4_OPTS, outfile); //?
  //sprintf(command, "sim4 %s %s A=5>%s", exonname,names[i-1], outfile);
  if (VERBOSE)
    fprintf(stderr, "\nCOMMAND: %s\n", command);
  run_prog(PATH_SIM4, command);
  //system(command);

  /*
    sprintf(command,"sim4 %s %s A=5 >%s",exonname,names[i-1],outfile); //?
    system(command); //Call sim4
  */
      
  if (VERBOSE)
    printf("Sim4 finished.. Now saving the beginning and end locations\n");
  
  if (VERBOSE)
    printf("Opening file > %s\n",outfile);
  in=gfopen(outfile,"r");

  while(fscanf(in,"%c",&ch) > 0){
    if (ch == '\n' || ch=='\r')
      linecnt++;
  }

  rewind(in);
  
  pos = (position *) malloc(sizeof(position) * linecnt);

  SkipNlines(2,in); //We arent interested in first two lines

  
  while(fscanf(in,"%d",&iterator) > 0)
    {
      pos[++MaxExons].begin=iterator;
      fscanf(in,"%d",&iterator);
      pos[MaxExons].end=iterator;
    }//while
  
  fclose(in); //Close the out file
  //Cleaning stuff
   
  sprintf(command,"rm -f %s ",outfile);
  system(command); 
} // sim4

void readtablefile(position **filepos, char *tablefile)
{
  FILE *in;
  char line[50];//word
  int index;
  int linecnt=1;
  char ch;
    
  in=gfopen(tablefile,"r");
  if (VERBOSE)
    printf("Opening file > %s\n",tablefile);
  
  while(fscanf(in,"%c",&ch) > 0){
    if (ch == '\n' || ch=='\r')
      linecnt++;
  }

  rewind(in);
  
  pos = (position *) malloc(sizeof(position) * linecnt);
                  
  //Preprocessing step finds the number of sequences and their names 


  while(fscanf(in,"%s",line) > 0)
    {
      if (!NOUPDATE)
	index = findstring(line);
      else
	index = 1; // doesn't matter anyway
      if (index>0)
	{
	  filepos[index] = (position *) realloc(filepos[index],(sizeof(position) * 
(linecnt)));  
	
	  fscanf(in,"%s",line);
	  length[index]=(length[index])+1;
           
	  filepos[index][length[index]].begin=atoi(line);
	  fscanf(in,"%s",line);
	  filepos[index][length[index]].end=atoi(line);
	
	  /*
	    update on Aug 10, 2005
	    moved from processOutOptions
	  */
	  
	  if (filepos[index][length[index]].begin >= filepos[index][length[index]].end){
	    printf("Error !! Be careful to enter CORRECT start & end positions\n");
	    printf("Start position should be SMALLER than the End position\n");
	    printf("Error in Line : %d:%d\n", filepos[index][length[index]].begin, filepos[index][length[index]].end);
	    break;
	  } // if
	  else if  (filepos[index][length[index]].begin< 1 || filepos[index][length[index]].begin > maxLen || filepos[index][length[index]].end < 1 || filepos[index][length[index]].end > maxLen){
	    printf("Error !! Be careful to enter CORRECT start & end positions\n");
	    printf("Start & End Positions should be between 1 and length of alingment: %d\n",maxLen);
	    printf("Error in Line : %d:%d\n", filepos[index][length[index]].begin, filepos[index][length[index]].end);
	    break;
	  }

	  /*
	  printf("debug %d %d\n", filepos[index][length[index]].begin,filepos[index][length[index]].end);
	  */
	}//if1
      else
	{
	  printf("Unexpected input in the table file exiting...\n");
	  die();
	}

    }//while

  
  //printf(" Length of the first tablefile is %d second file is %d\n",length[1],length[2]);
  fclose(in); //Close the out file

} // readtablefile

int parsetablearray(int index, position **filepos)
{
  int j;

  
  for (j=1;j<=length[index];j++)
    {
      pos[j].begin=filepos[index][j].begin;
      pos[j].end=filepos[index][j].end;
    }
  return length[index];
} // parsetablearray


void allToggles(bool inputFile){
  strcpy(tablefile, "default");
  if (!inputFile){
    printf("Enter Alignment File: ");
    fflush(stdout);
    fgets(inFile, 100, stdin);
    inFile[strlen(inFile)-1] = 0;
    if (inFile[0] == 0)
      exit(0);
  } // if
  seqTot = readAln(inFile);
  ret = (char *) malloc(maxLen+1);
  do {
    inputToggles();
  } while (!checkToggles(TRUE)); 
} // menu

void inputToggles(){
  int select=1;
  char choice[2];
  int pselect=0;
  char line[SEQ_LENGTH*2];
  while (select != 6 && select != RUNPROGRAM){
    printf("\n-=-=-=-=-=-=-=-=-=-= INPUT TOGGLES =-=-=-=-=-=-=-=-=-=-\n\n");
    printf("1 - Toggle MAX/MIN [%s]\n", MAX ? "MAX" : "MIN");
    printf("2 - Toggle KEEP/TOSS [%s]\n", KEEP ? "KEEP" : "TOSS");
    printf("3 - Select Program ");
    switch (program){
    case CROSS_MATCH:
      printf("[CROSSMATCH]\n");
      break;
    case REPEATMASKER:
      printf("[REPEATMASKER]\n");
      break;
    case SIM4:
      printf("[SIM4]\n");
      break;
    case TABLE:
      printf("[TABLE]\n");
      break;
    case CONVERT:
      printf("[CONVERT]\n");
      break;
    case NONE:
      printf("[NONE]\n");
      break;
    } // switch program 
    printf("4 - Select cDNA File [%s] \n", tablefile);
    if (program == REPEATMASKER)
      printf("5 - Set RepeatMasker Options [%s] \n", RMASKER_OPTS);
    else if (program == CROSS_MATCH)
      printf("5 - Set Cross_Match Options [%s] \n", CMATCH_OPTS);
    else if (program == SIM4)
      printf("5 - Set Sim4 Options [%s] \n", SIM4_OPTS);
    else if (program == TABLE)
      printf("5 - Toggle Update locations [%c] \n", NOUPDATE?'N':'Y');
    if (program == CONVERT || program == NONE){
      printf("5 - Output Toggles \n");
      printf("6 - Run MaM \n");
      printf("7 - Quit MaM\n");
    }
    else{
      printf("6 - Output Toggles \n");
      printf("7 - Run MaM \n");
      printf("8 - Quit MaM\n");
    }
    printf("\n\nSelection: ");
    fflush(stdout);
    scanf("%s", choice);
    choice[1] = 0;
    select = atoi(choice);
    switch (select){
    case 1:
      MAX = !MAX;
      break;
    case 2:
      KEEP = !KEEP;
      break;
    case 3:
      printf("\nSelect one of the programs : \n\n");
      do{
	printf(" 1 - Cross_Match\n");
	printf(" 2 - Repeatmasker\n");
	printf(" 3 - Sim4\n");
	printf(" 4 - Table\n");	
	printf(" 5 - Just Convert Format (don't edit anything)\n");	
	printf(" 6 - None (only -alnstats will work)\n");	
	printf(" 7 - Exit\n");
	printf("\n\nSelection: ");
	fflush(stdout);
	scanf("%s", choice);
	choice[1] = 0;
	pselect = atoi(choice);
      } while (!(pselect > 0 && pselect <= 7));
      if (pselect != 7){
	program = pselect - 1;
	if (!checkToggles(VERBOSE))
	  getFileInput();
      } // if pselect != 6
      break;
    case 4:
      getFileInput();
      break;
    case 5:
      if (program == CONVERT || program == NONE)
	select = outputToggles();
      else if (program == REPEATMASKER){
	printf("\nOverride RepeatMasker Options [%s]: ", RMASKER_OPTS);
	fflush(stdout);
	getchar();
	fgets(line, SEQ_LENGTH*2, stdin);
	line[strlen(line)-1] = 0;
	/*if (line[0] != 0)*/
	strcpy(RMASKER_OPTS, line); 
      }
      else if (program == CROSS_MATCH){	printf("\nOverride Cross_Match Options [%s]: ", CMATCH_OPTS);
	fflush(stdout);
	getchar();
	fgets(line, SEQ_LENGTH*2, stdin);
	line[strlen(line)-1] = 0;
	/*if (line[0] != 0) */
	strcpy(CMATCH_OPTS, line);	
      }
      else if (program == SIM4){
	printf("\nOverride Sim4 Options [%s]: ", SIM4_OPTS);
	fflush(stdout);
	getchar();
	fgets(line, SEQ_LENGTH*2, stdin);
	line[strlen(line)-1] = 0;
	/*if (line[0] != 0)*/
	strcpy(SIM4_OPTS, line);
      }
      else if (program == TABLE){
	/*
	printf("\nToggle Update locations [%c]: ", NOUPDATE?'N':'Y');
	fflush(stdout);
	getchar();
	fgets(line, SEQ_LENGTH*2, stdin);
	line[strlen(line)-1] = 0;
	*/
	NOUPDATE = !NOUPDATE;
	/*if (line[0] != 0)*/
	
      }
      break;
    case 6: 
      if (program == NONE || program == CONVERT)
	return;
      else
	select = outputToggles();
      break;
    case 7:
      if (program == NONE || program == CONVERT)
	exit(0);
      else
	return;
      break;
    case 8:
      if (program == NONE || program == CONVERT)
	break;
      else
	exit(0);
    default:
      break;
    } // switch select
  } // while select != 5
} // inputtoggles

int outputToggles(){
  int select=1;
  char choice[3];
  int sw, ww;
  while (select != 14){
    printf("\n-=-=-=-=-=-=-=-=-=-= OUTPUT TOGGLES =-=-=-=-=-=-=-=-=-=-\n\n");
    printf(" 1 - Toggle SINGLE/MULTIPLE Files [%s]\n", SINGLE_OUT ? "SINGLE" : "MULTIPLE");
    printf(" 2 - Toggle Slider [%s]\n", SLIDER_OUT ? "ON" : "OFF");
    printf(" 3 - Toggle Consensus Out [%s]\n", CONSENSUS_OUT ? "ON" : "OFF");
    printf(" 4 - Toggle Gaps in Consensus Out [%s]\n", CONSENSUS_GAPS ? "ON" : "OFF");
    printf(" 5 - Toggle CLUSTALW Format Output [%s]\n", CLUSTALW_OUT ? "ON" : "OFF");
    printf(" 6 - Toggle NEXUS Format Output [%s]\n", NEXUS_OUT ? "ON" : "OFF");
    printf(" 7 - Toggle MEGA Format Output [%s]\n", MEGA_OUT ? "ON" : "OFF");
    printf(" 8 - Toggle FASTA Format Output [%s]\n", FASTA_OUT ? "ON" : "OFF");
    printf(" 9 - Toggle PHYLIP Format Output [%s]\n", PHYLIP_OUT ? "ON" : "OFF");
    printf("10 - Toggle HTML Format Output [%s]\n", HTML_OUT ? "ON" : "OFF");
    printf("11 - Toggle Identity Dotted Representation Format Output [%s]\n", IDENTITY_OUT ? "ON" : "OFF");
    printf("12 - Toggle Alignment Statistics File Dump [%s]\n", ALNSTATS ? "ON" : "OFF");
    printf("13 - Toggle Including Consensus in Output [%s]\n", INCLUDE_CONSENSUS ? "ON" : "OFF");
    printf("14 - Input Toggles\n");
    printf("15 - Run MaM\n");
    printf("16 - Quit MaM\n");
    printf("\n\nSelection: ");
    fflush(stdout);
    scanf("%s", choice);
    select = atoi(choice);
    switch(select){
    case 1:
      SINGLE_OUT = !SINGLE_OUT;
      break;
    case 2:
      SLIDER_OUT = !SLIDER_OUT;
      if (SLIDER_OUT){
	strcpy(PC, "P");
	strcpy(slidewindow, "5");
	strcpy(slidewidth,  "5");
	printf("Pairwise, Complete Deletion or Parsimony Score [P/C/S]> ");
	scanf("%s", PC);
	PC[0] = toupper(PC[0]);
	while (!(PC[0] == 'P' || PC[0]=='C' || PC[0]=='S')){
	  printf("Enter 'P', 'C' or 'S'> ");
	  scanf("%s", PC);	
	  PC[0] = toupper(PC[0]);
	}  // while
	printf("Enter Slide Width > ");
	scanf("%d", &sw);
	printf("Enter Window Width > ");
	scanf("%d", &ww);
	sprintf(slidewidth, "%d", sw);
	sprintf(slidewindow, "%d", ww);
      } // if SLIDER
      break;
    case 3:
      CONSENSUS_OUT = !CONSENSUS_OUT;
      break;
    case 4:
      CONSENSUS_GAPS = !CONSENSUS_GAPS;
      break;
    case 5:
      CLUSTALW_OUT = !CLUSTALW_OUT;
      break;
    case 6:
      NEXUS_OUT = !NEXUS_OUT;
      break;
    case 7:
      MEGA_OUT = !MEGA_OUT;
      break;
    case 8:
      FASTA_OUT = !FASTA_OUT;
      break;
    case 9:
      PHYLIP_OUT = !PHYLIP_OUT;
      break;
    case 10:
      PHYLIP_OUT = !PHYLIP_OUT;
      break;
    case 11:
      IDENTITY_OUT = !IDENTITY_OUT;
      break;
    case 12:
      ALNSTATS = !ALNSTATS;
      break;
    case 13:
      INCLUDE_CONSENSUS = !INCLUDE_CONSENSUS;
      break;
    case 14:
      inputToggles();
      return 0;
      break;
    case 15:
      return RUNPROGRAM;
      break;
    case 16:
      exit(0);
      break;
    default:
      break;
    } // switch
  }// while
  return 0;
} // outputtoggles

void getFileInput(){
  FILE *fp;
  //char noupdate;

  switch (program){
  case CROSS_MATCH:
  case REPEATMASKER:
    printf("Enter cDNA file or \"default\": ");
    fflush(stdout);
    getchar();
    fgets(tablefile, 100, stdin);
    tablefile[strlen(tablefile)-1] = 0;
    if (tablefile[0] == 0 || !strcmp(tablefile, "default"))
      strcpy(tablefile, "default");
    else{
      if ((fp = fopen(tablefile, "r")) == NULL){
	printf("Cannot open %s. Reloading defaults [REPEATMASKER - default]\n", tablefile);
	program = REPEATMASKER;
	strcpy(tablefile, "default");
      } // if fopen
      else
	fclose(fp);
    }
    break;
  case SIM4:
    printf("Enter cDNA file : ");
    fflush(stdout);
    getchar();
    fgets(tablefile, 100, stdin);
    tablefile[strlen(tablefile)-1] = 0;
    if (!strcmp(tablefile, "default")){
      printf("SIM4 is not compatible with \"default\" option !!!\n");
      printf("Program option is converted to Repeatmasker\n");
      program = REPEATMASKER;
    } // if sim4check
    else{
      if ((fp = fopen(tablefile, "r")) == NULL){
	printf("Cannot open %s. Reloading defaults [REPEATMASKER - default]\n", tablefile);
	program = REPEATMASKER;
	strcpy(tablefile, "default");
      } // if fopen
      else
      fclose(fp);
    }
    break;
  default:
    printf("Enter tablefile : ");
    fflush(stdout);
    getchar();
    fgets(tablefile, 100, stdin);
    tablefile[strlen(tablefile)-1] = 0;
    if (!strcmp(tablefile, "default")){
      printf("Table option is not compatible with \"default\" option !!!\n");
      printf("Program option is converted to Repeatmasker\n");
      program = REPEATMASKER;
    } // if sim4check
    else{
      if ((fp = fopen(tablefile, "r")) == NULL){
	printf("Cannot open %s. Reloading defaults [REPEATMASKER - default]\n", tablefile);
	program = REPEATMASKER;
	strcpy(tablefile, "default");
      } // if fopen
      else
	fclose(fp);
    }
    /*
    printf("Update coordinates (coordinates from unaligned sequence)? [Y/n] :");
    fflush(stdout);
    noupdate=getchar();
    if (toupper(noupdate) == 'N')
      NOUPDATE=TRUE;
    else
      NOUPDATE=FALSE;
    */
    break;
  } // switch program 
} // getfileinput

bool checkToggles(bool verbose){
  bool ret=TRUE;
  if (program == TABLE && !strcmp(tablefile, "default")){
    if (verbose)
      printf("Table option is not compatible with \"default\" option !!!\n");
    ret = FALSE;
  }  // tablefile & default
  if (program != TABLE && NOUPDATE){
    if (verbose){
      printf("Warning: NOUPDATE flag won't work unless you select table option!!");
      printf("Omitting the NOUPDATE flag, and moving on...\n");
    }
  }
  if (program == SIM4 && !strcmp(tablefile, "default")){
    if (verbose)
      printf("SIM4 is not compatible with \"default\" option !!!\n");
    ret = FALSE;
  } // sim4 and default
  if (program == REPEATMASKER && strstr(RMASKER_OPTS, "-lib")){
    printf("Please DO NOT use -lib option when overriding default RepeatMasker option. Use the -exonfile switch instead.\n");
    ret = FALSE;
  }

  if (!ret && verbose)
    printf("Errors Found !!!\n");
  return ret;
} // check toggles


void parseParams(int argc, char **argv){

  int i, j;
  int cnt;
  char param[20];

  if (!strcmp(argv[1], "-h"))
    help(argv[0],0);

  if (!strcmp(argv[1], "-v")){
    printf("MaM version 1.4.2\n");
    printf("Latest update: 03/20/2006\n");
    exit(0);
  }
  
  if (!strcmp(argv[1], "-defaults"))
    defaults();
  
  if (argv[1][0] == '-')
    help(argv[0],0);

  strcpy(inFile, argv[1]); // Alignment File Name
  
  i = 2;
  cnt = 0;
  strcpy(tablefile, "default");
  strcpy(PC, "P");
  strcpy(slidewindow, "5");
  strcpy(slidewidth, "5");
  program = REPEATMASKER;

  for (i=2; i < argc; i++){
    if (!strcmp(argv[i], "-h"))
      help(argv[0],0);
    printf("%s\n", argv[i]);
    if (strstr(argv[i], "-clustal")){ // CLUSTAL OUT
      for (j=8;j<strlen(argv[i]);j++)
	param[j-8] = argv[i][j];
      param[j-8] = 0;
      if (!strcmp(param, "=on"))
	CLUSTALW_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	CLUSTALW_OUT = FALSE;
      else
	help(argv[0], i);
      
    } // CLUSTAL OUT

    else if (strstr(argv[i], "-nexus")){ // NEXUS OUT
      for (j=6;j<strlen(argv[i]);j++)
	param[j-6] = argv[i][j];
      param[j-6] = 0;
      if (!strcmp(param, "=on"))
	NEXUS_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	NEXUS_OUT = FALSE;
      else
	help(argv[0], i);
    } // NEXUS OUT
    
    else  if (strstr(argv[i], "-identity")){ // IDENTITY DOT REPR. OUT
      for (j=9;j<strlen(argv[i]);j++)
	param[j-9] = argv[i][j];
      param[j-9] = 0;
      if (!strcmp(param, "=on"))
	IDENTITY_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	IDENTITY_OUT = FALSE;
      else
	help(argv[0], i);
    } // IDENTITY DOT REPR. OUT

    else  if (strstr(argv[i], "-mega")){ // MEGA OUT
      for (j=5;j<strlen(argv[i]);j++)
	param[j-5] = argv[i][j];
      param[j-5] = 0;
      if (!strcmp(param, "=on"))
	MEGA_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	MEGA_OUT = FALSE;
      else
	help(argv[0], i);
    } // MEGA OUT

    else  if (strstr(argv[i], "-fasta")){ // FASTA OUT
      for (j=6;j<strlen(argv[i]);j++)
	param[j-6] = argv[i][j];
      param[j-6] = 0;
      if (!strcmp(param, "=on"))
	FASTA_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	FASTA_OUT = FALSE;
      else
	help(argv[0], i);
    } // FASTA OUT

    else  if (strstr(argv[i], "-phylip")){ // PHYLIP OUT
      for (j=7;j<strlen(argv[i]);j++)
	param[j-7] = argv[i][j];
      param[j-7] = 0;
      if (!strcmp(param, "=on"))
	PHYLIP_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	PHYLIP_OUT = FALSE;
      else
	help(argv[0], i);
    } // PHYLIP OUT

    else  if (strstr(argv[i], "-html")){ // HTML OUT
      for (j=5;j<strlen(argv[i]);j++)
	param[j-5] = argv[i][j];
      param[j-5] = 0;
      if (!strcmp(param, "=on"))
	HTML_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	HTML_OUT = FALSE;
      else
	help(argv[0], i);
    } // HTML OUT

    else  if (strstr(argv[i], "-alnstats")){ // ALNSTATS OUT
      for (j=9;j<strlen(argv[i]);j++)
	param[j-9] = argv[i][j];
      param[j-9] = 0;
      if (!strcmp(param, "=on"))
	ALNSTATS = TRUE;
      else if (!strcmp(param, "=off"))
	ALNSTATS = FALSE;
      else
	help(argv[0], i);
    } // ALNSTATS OUT

    else  if (strstr(argv[i], "-consensus")){ //  Toggle Consensus out
      for (j=10;j<strlen(argv[i]);j++)
	param[j-10] = argv[i][j];
      param[j-10] = 0;
      if (!strcmp(param, "=on"))
	CONSENSUS_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	CONSENSUS_OUT = FALSE;
      else
	help(argv[0], i);
    } //  Toggle Consensus out

    else  if (strstr(argv[i], "-update")){ //  Toggle NOUPDATE
      for (j=7;j<strlen(argv[i]);j++)
	param[j-7] = argv[i][j];
      param[j-7] = 0;
      if (!strcmp(param, "=on"))
	NOUPDATE = FALSE;
      else if (!strcmp(param, "=off"))
	NOUPDATE = TRUE;
    } //  Toggle Consensus out

    else  if (strstr(argv[i], "-include")){ //  Toggle Include Consensus in output
      for (j=8;j<strlen(argv[i]);j++)
	param[j-8] = argv[i][j];
      param[j-8] = 0;
      if (!strcmp(param, "=on"))
	INCLUDE_CONSENSUS = TRUE;
      else if (!strcmp(param, "=off"))
	INCLUDE_CONSENSUS = FALSE;
      else
	help(argv[0], i);
    } //  Toggle Include Consensus in output
    
    else  if (strstr(argv[i], "-slider")){ //  Toggle slider
      for (j=7;j<strlen(argv[i]);j++)
	param[j-7] = argv[i][j];
      param[j-7] = 0;
      printf("slider param: %s\n", param);
      if (!strcmp(param, "=on"))
	SLIDER_OUT = TRUE;
      else if (!strcmp(param, "=off"))
	SLIDER_OUT = FALSE;
      else
	help(argv[0], i);
    } //  Toggle slider

    else  if (strstr(argv[i], "-pc")){ //  pairwise/complete
      for (j=3;j<strlen(argv[i]);j++)
	param[j-3] = argv[i][j];
      param[j-3] = 0;
      if (!strcmp(param, "=p") || !strcmp(param, "=P"))
	strcpy(PC, "P");
      else if (!strcmp(param, "=c") || !strcmp(param, "=C"))
	strcpy(PC, "C");
      else if (!strcmp(param, "=s") || !strcmp(param, "=S"))
	strcpy(PC, "S");
      else
	help(argv[0], i);
    } //  pairwise/complete
    
    else  if (strstr(argv[i], "-sw=")){ //  slide width
      for (j=4;j<strlen(argv[i]);j++)
	param[j-4] = argv[i][j];
      param[j-4] = 0;
      if (param[0] == 0)
	help(argv[0], i);
      strcpy(slidewidth, param);
    } //  slide width

    else  if (strstr(argv[i], "-ww=")){ //  window width
      for (j=4;j<strlen(argv[i]);j++)
	param[j-4] = argv[i][j];
      param[j-4] = 0;
      if (param[0] == 0)
	help(argv[0], i);
      strcpy(slidewindow, param);
    } //  window width
    
    else  if (strstr(argv[i], "-program")){ //  program to use
      for (j=8;j<strlen(argv[i]);j++)
	param[j-8] = argv[i][j];
      param[j-8] = 0;
      if (!strcmp(param, "=crossmatch"))
	program = CROSS_MATCH;
      else if (!strcmp(param, "=repeatmasker"))
	program = REPEATMASKER;
      else if (!strcmp(param, "=sim4"))
	program = SIM4;
      else if (!strcmp(param, "=table"))
	program = TABLE;
      else if (!strcmp(param, "=convert"))
	program = CONVERT;
      else if (!strcmp(param, "=none"))
	program = NONE;
      else
	help(argv[0], i);      
    } //  program to use

    else  if (strstr(argv[i], "-keep")){ //  keep
      for (j=5;j<strlen(argv[i]);j++)
	param[j-5] = argv[i][j];
      param[j-5] = 0;
      if (!strcmp(param, "=on"))
	KEEP = TRUE;
      else if (!strcmp(param, "=off"))
	KEEP = FALSE;
      else
	help(argv[0], i);
    } //  toggle KEEP 

    else  if (strstr(argv[i], "-cgaps")){ //  cgaps
      for (j=6;j<strlen(argv[i]);j++)
	param[j-6] = argv[i][j];
      param[j-6] = 0;
      if (!strcmp(param, "=on"))
	CONSENSUS_GAPS = TRUE;
      else if (!strcmp(param, "=off"))
	CONSENSUS_GAPS = FALSE;
      else
	help(argv[0], i);
    } //  toggle CONSENSUS_GAPS


    else  if (strstr(argv[i], "-merge")){ //  max
      for (j=6;j<strlen(argv[i]);j++)
	param[j-6] = argv[i][j];
      param[j-6] = 0;
      if (!strcmp(param, "=max"))
	MAX = TRUE;
      else if (!strcmp(param, "=min"))
	MAX = FALSE;
      else
	help(argv[0], i);
    } //  toggle max
         
    else  if (strstr(argv[i], "-exonfile")){ //  input exon file
      for (j=10;j<strlen(argv[i]);j++)
	param[j-10] = argv[i][j];
      param[j-10] = 0;
      if (param[0] == 0)
	help(argv[0], i);      
      strcpy(tablefile, param);
    } //  input exon file

    else  if (strstr(argv[i], "-column")){ 
      for (j=7;j<strlen(argv[i]);j++)
	param[j-7] = argv[i][j];
      param[j-7] = 0;
      if (!strcmp(param, "=single"))
	SINGLE_OUT = TRUE;
      else if (!strcmp(param, "=multiple"))
	SINGLE_OUT = FALSE;
      else
	help(argv[0], i);
    } //  single/multiple file generation
          
    else  if (strstr(argv[i], "-rmasker_opts")){ 
      for (j=14;j<strlen(argv[i]);j++)
	param[j-14] = argv[i][j];
      param[j-14] = 0;
      strcpy(RMASKER_OPTS, param);
    } //  rmasker options 
          
    else  if (strstr(argv[i], "-cmatch_opts")){ 
      for (j=13;j<strlen(argv[i]);j++)
	param[j-13] = argv[i][j];
      param[j-13] = 0;
      strcpy(CMATCH_OPTS, param);
    } //  cmatch options 
          
    else  if (strstr(argv[i], "-sim4_opts")){ 
      for (j=11;j<strlen(argv[i]);j++)
	param[j-11] = argv[i][j];
      param[j-11] = 0;
      strcpy(SIM4_OPTS, param);
    } //  sim4 options 
          
    else if (!strcmp(argv[i], "-defaults"))
      defaults();
    else if (!strcmp(argv[i], "-V"))
      VERBOSE = TRUE;

    else
      help(argv[0], 0);    

  } // for i all params are read

  if (checkToggles(VERBOSE)){
    seqTot = readAln(inFile);
    ret = (char *) malloc(maxLen+1);
  }
  else
    help(argv[0], 0);
  
} // parseParams


void help(char *pName, int i){
  printf("\n\n%s%s%s%s",
	 "Multiple Alignment Manipulator \n", 
	 "Usage: ",pName, " [Alignment File] [option1][=on/=off] [option2=on/off] ... [option n=on/off]\n");
  printf("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
	 "\nOptions: \n",
	 "-exonfile     : Exon table file or cDNA file.\n\t\tEx: -exonfile=default , -exonfile=NPIP\n",
	 "-update       : Update the alignment/sequence coordinates On/Off.\n\t\tSet to off if you want to use the alignment coordinates in the tablefile.\n",
	 "-column       : Toggle Single/Multiple File Output\n\t\tEx -column=single  or  -column=multiple\n",
	 "-merge        : Toggle MAX/MIN\n\t\tEx: -merge=max  or  -merge=min\n",
	 "-keep         : Toggle KEEP/TOSS\n\t\tEx: -keep=on  or  keep=off\n",
	 "-program      : Select program. One of:\n\t\t-program=crossmatch\n\t\t-program=repeatmasker\n\t\t-program=sim4\n\t\t-program=table\n\t\t-program=convert\n\t\t-program=none\n",
	 "-slider       : Toggle Slider On/Off\n",
	 "-pc           : Select One Of:\n\t\tPairwise Deletion (-pc=p),\n\t\tComplete Deletion (-pc=c),\n\t\tParsimony Score   (-pc=s)\n",
	 "-sw           : Select Slide Width (-sw=10)\n",
	 "-ww           : Select Window Width (-ww=100)\n",
	 "-clustal      : Toggle Clustal Output Format On/Off\n",
	 "-nexus        : Toggle Nexus Output Format On/Off\n",
	 "-mega         : Toggle Mega Output Format On/Off\n",
	 "-fasta        : Toggle Fasta Output Format On/Off\n",
	 "-phylip       : Toggle Phylip Output Format On/Off\n",
	 "-html         : Toggle HTML Output Format On/Off\n",
	 "-identity     : Toggle Consensus Identity Output Format On/Off\n",
	 "-alnstats     : Toggle Alignment Statistics File Dump On/Off\n",
	 "-consensus    : Toggle Consensus Output On/Off\n",
	 "-cgaps        : Toggle Gaps in Consensus Output On/Off\n",
	 "-include      : Toggle Including Consensus in output On/Off\n",
	 "-rmasker_opts : Override Default RepeatMasker Options\n",
	 "-cmatch_opts  : Override Default Cross_Match Options\n",
	 "-sim4_opts    : Override Default Sim4 Options\n",
	 "-defaults     : See defaults\n",
	 "-v            : Version\n",
	 "-V            : Verbose\n",
	 "-h            : Help\n\n");
  if (i)
    printf("Error in arg %d \n",i);
  exit(0);
} // help 

void defaults(){
  printf("\n\nMaM Defaults:\n\n");
  printf("-exonfile     : \"default\"\n");
  printf("-update       : on\n");
  printf("-column       : single\n");
  printf("-merge        : max\n");
  printf("-keep         : on\n");
  printf("-program      : repeatmasker\n");
  printf("-slider       : off\n");
  printf("-pc           : P\n");
  printf("-sw           : 5\n");
  printf("-ww           : 5\n");
  printf("-clustal      : on\n");
  printf("-nexus        : off\n");
  printf("-mega         : off\n");
  printf("-fasta        : off\n");
  printf("-phylip       : off\n");
  printf("-html         : off\n");
  printf("-identity     : off\n");
  printf("-alnstats     : off\n");
  printf("-consensus    : off\n");
  printf("-cgaps        : on\n");
  printf("-include      : off\n");
  printf("-rmasker_opts : %s\n", RMASKER_OPTS);
  printf("-cmatch_opts  : %s\n", CMATCH_OPTS);
  printf("-sim4_opts    : %s\n", SIM4_OPTS);
  exit(0);
} // defaults


void die(void){

  int i;
  char command[300];
  printf("Deleting temporary files.. ");
  fflush(stdout);
  for (i=1;i<=seqTot;i++){              
    sprintf(command,"rm -f %s ",names[i-1]);
    system(command);
    sprintf(command,"rm -f %s.aln ",names[i-1]);
    system(command);
  } 
  printf("[OK]\n");
  exit(0);
}

void loadConfig(void){
  FILE *config;
  char cfname[SEQ_LENGTH];
  char line[2*SEQ_LENGTH];
  int i;

  PATH_RMASKER[0]=0;
  PATH_CMATCH[0]=0;
  PATH_SIM4[0]=0;
  PATH_GNUPLOT[0]=0;

  sprintf(cfname, "%s/%s", getenv("HOME"), CONFIG_FILE);

  config = fopen(cfname, "r");

  if (config == NULL){
    /* create new config file */
    createConfig(cfname);
    return;
  }

  line[0]=0;
  do{
    if (fgets(line, 2*SEQ_LENGTH, config) == NULL)
      break;
    line[strlen(line)-1] = 0;
    if (line[0] != 0 && line[0] != '#'){ // comment lines are omitted
      if (strstr(line, "REPEATMASKER")){
	for (i=strlen("REPEATMASKER="); i<strlen(line); i++)
	  PATH_RMASKER[i-strlen("REPEATMASKER=")] = line[i];
	PATH_RMASKER[i]=0;
      }
      else if (strstr(line, "CROSSMATCH")){
	for (i=strlen("CROSSMATCH="); i<strlen(line); i++)
	  PATH_CMATCH[i-strlen("CROSSMATCH=")] = line[i];
	PATH_CMATCH[i]=0;
      }
      else if (strstr(line, "SIM4")){
	for (i=strlen("SIM4="); i<strlen(line); i++)
	  PATH_SIM4[i-strlen("SIM4=")] = line[i];
	PATH_SIM4[i]=0;
      }
      else if (strstr(line, "GNUPLOT")){
	for (i=strlen("GNUPLOT="); i<strlen(line); i++)
	  PATH_GNUPLOT[i-strlen("GNUPLOT=")] = line[i];
	PATH_GNUPLOT[i]=0;
      }
    }
  }while (line[0] != 0);

  if (PATH_RMASKER[0]==0 && VERBOSE)
    fprintf(stderr, "Warning: RepeatMasker path is not in the config file.\n");
  else if (VERBOSE)
    fprintf(stderr, "RepeatMasker path: %s\n", PATH_RMASKER);

  if (PATH_CMATCH[0]==0 && VERBOSE)
    fprintf(stderr, "Warning: cross_match path is not in the config file.\n");
  else if (VERBOSE)
    fprintf(stderr, "cross_match path: %s\n", PATH_CMATCH);

  if (PATH_SIM4[0]==0 && VERBOSE)
    fprintf(stderr, "Warning: sim4 path is not in the config file.\n");
  else if (VERBOSE)
    fprintf(stderr, "sim4 path: %s\n", PATH_SIM4);

  if (PATH_GNUPLOT[0]==0 && VERBOSE)
    fprintf(stderr, "Warning: gnuplot path is not in the config file.\n");
  else if (VERBOSE)
    fprintf(stderr, "gnuplot path: %s\n", PATH_GNUPLOT);

}

void createConfig(char *cfname){
  FILE *config;
  FILE *pipe;
  char cmd[SEQ_LENGTH];
  
  cmd[0]=0;
  pipe = popen("which RepeatMasker 2>/dev/null", "r");
  if (pipe == NULL){
    fprintf(stderr, "\"which\" tool is not found in your system.\nCreate %s file manually or install \"which\" tool.\n", cfname);
    //exit(0);
  }
  
  fgets(cmd, SEQ_LENGTH, pipe);
  pclose(pipe);

  if (cmd[0]==0){

    pipe = popen("which repeatmasker 2>/dev/null", "r");
    if (pipe == NULL){
      fprintf(stderr, "\"which\" tool is not found in your system.\nCreate %s file manually or install \"which\" tool.\n", cfname);
      //exit(0);
    }
    
    fgets(cmd, SEQ_LENGTH, pipe);
    pclose(pipe);
    
    if (cmd[0]==0)
      fprintf(stderr, "\nRepeatMasker not found in PATH.\nPlease install it or manually configure the %s file.\n", cfname);
    else{
      /* else, the repeatmasker is found, set the path */
      cmd[strlen(cmd)-1] = 0;
      strcpy(PATH_RMASKER, cmd);
    } 
  }
  else{
  /* else, the repeatmasker is found, set the path */
    cmd[strlen(cmd)-1] = 0;
    strcpy(PATH_RMASKER, cmd);
  }

  /* cross_match */

  cmd[0]=0;
  pipe = popen("which cross_match 2>/dev/null", "r");
  if (pipe == NULL){
    fprintf(stderr, "\"which\" tool is not found in your system.\nCreate %s file manually or install \"which\" tool.\n", cfname);
    //exit(0);
  }

  fgets(cmd, SEQ_LENGTH, pipe);
  pclose(pipe);

  if (cmd[0]==0){
    fprintf(stderr, "\ncross_match not found in PATH.\nPlease install it or manually configure the %s file.\n", cfname);
  }
  else{
  /* else, the cross_match is found, set the path */
    cmd[strlen(cmd)-1] = 0;
    strcpy(PATH_CMATCH, cmd);
  }

  /* sim4 */

  cmd[0]=0;
  pipe = popen("which sim4 2>/dev/null", "r");
  if (pipe == NULL){
    fprintf(stderr, "\"which\" tool is not found in your system.\nCreate %s file manually or install \"which\" tool.\n", cfname);
    //exit(0);
  }

  fgets(cmd, SEQ_LENGTH, pipe);
  pclose(pipe);

  if (cmd[0]==0){
    fprintf(stderr, "\nsim4 not found in PATH.\nPlease install it or manually configure the %s file.\n", cfname);
  }
  
  else{
  /* else, the sim4 is found, set the path */
    cmd[strlen(cmd)-1] = 0;
    strcpy(PATH_SIM4, cmd);
  }
  
  /* gnuplot */

  cmd[0]=0;
  pipe = popen("which gnuplot 2>/dev/null", "r");
  if (pipe == NULL){
    fprintf(stderr, "\"which\" tool is not found in your system.\nCreate %s file manually or install \"which\" tool.\n", cfname);
    //exit(0);
  }

  fgets(cmd, SEQ_LENGTH, pipe);
  pclose(pipe);

  if (cmd[0]==0){
    fprintf(stderr, "\ngnuplot not found in PATH.\nPlease install it or manually configure the %s file.\n", cfname);
  }
  else{
  /* else, the gnuplot is found, set the path */
    cmd[strlen(cmd)-1] = 0;
    strcpy(PATH_GNUPLOT, cmd);
  }
  
  config = fopen(cfname, "w");
  fprintf(config, "REPEATMASKER=%s\n", PATH_RMASKER);
  fprintf(config, "CROSSMATCH=%s\n", PATH_CMATCH);
  fprintf(config, "SIM4=%s\n", PATH_SIM4);
  fprintf(config, "GNUPLOT=%s\n", PATH_GNUPLOT);
  fclose(config);
 
  fprintf(stderr,"\n\n");
  fprintf(stderr, "**********************************************************************\n");
  fprintf(stderr, "*           Config file $HOME/.mam-config is created.                *\n"); 
  fprintf(stderr, "*            Check the configuration file for any errors,            *\n");
  fprintf(stderr, "* then run MaM again to load the configuration and normal operation. *\n");
  fprintf(stderr, "**********************************************************************\n\n");
  exit(0);

}

void run_prog(char *prog, char *opts){
  char command[(strlen(prog) + strlen(opts) + 2)];
  int ret; 

  //command = (char *) malloc(sizeof(char) * (strlen(prog) + strlen(opts) + 2));

  sprintf(command, "%s %s 2>/dev/null", prog, opts);  

  ret = system(command);

  if (WIFSIGNALED(ret) &&
      (WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT || WTERMSIG(ret) == SIGABRT)){
    fprintf(stderr, "Program %s crashed. Please check your configuration files and overridden program options!\n", prog);
    die();
  }
  
  if (!strcmp(prog, PATH_RMASKER) && ret == 256){
    fprintf(stderr, "Program %s crashed. Please check your configuration files and overridden program options!\n", prog);
    die();
  }
  
  
  //free(command);
}
