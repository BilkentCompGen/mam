/*
	MaM : Multiple alignment Manipulator
	
	Implemented by: Can ALKAN & Eray TUZUN
		
	[   calkan@gmail.com    ]
	[  eraytuzun@gmail.com  ]

	Last Update: March 20, 2006
	Summary: merge=min location information check fixed
	Summary: Length information printed (Aug 25, 2005)
	Summary: PHYLIP output support added (Oct 07, 2005)
	Summary: Convert option (Oct 10, 2005)
	Summary: -alnstats option (Oct 10, 2005)
	Summary: PHYLIP input support added (Oct 14, 2005)
	Summary: HTML input support added (Feb 28, 2006)
	Summary: Consensus Gaps On/Off added (March 20, 2006)
*/

#include "main.h"

int numberSeven = 7;
int maxLen;
char *ret;

void processOutOptions(char alnFileName[], char subFileName[], char exonName[], int maxsubseq){
  int startpos[maxsubseq]; int endpos[maxsubseq];
  int subseqCnt;
  FILE *subFile;
  char slidearg1[50];
  //char slidearg2[50];
  int i;

  for (i = strlen(alnFileName)-1; i>=0; i--)
    if (alnFileName[i] == '.')
      break;
  if (i!=0)
    alnFileName[i] = 0;

  printf("[OK]\nTotal %d sequences\n",seqTot);
  printf("------------------------------\n");
  
  printf("\nSubsequence Extraction\n");
  subseqCnt = 0;
  subFile = fopen(subFileName, "r");
  if (subFile == NULL){
    printf("Error reading subsequences file\n");
    exit(0);
  }
  
  while (fscanf(subFile, "%d %d", &(startpos[subseqCnt]),&(endpos[subseqCnt])) > 0 ){   
    if (startpos[subseqCnt] >= endpos[subseqCnt]){
      printf("Warning: No concatenation is possible at this step\n");
      printf("Skipping this line, good luck in the next step!\n");
      /*
      printf("Error in Line : %d -- %d:%d\n",(subseqCnt+1), startpos[subseqCnt], endpos[subseqCnt]);
      */
      subseqCnt = 0;
      break;
    } // if
    /*
    else if (startpos[subseqCnt] < 1 || startpos[subseqCnt] > maxLen || endpos[subseqCnt] < 1 || endpos[subseqCnt] > maxLen){
      printf("Error !! Be careful to enter CORRECT start & end positions\n");
      printf("Start & End Positions should be between 1 and length of alingment: %d\n",maxLen);
      printf("Error in Line : %d -- %d:%d\n",(subseqCnt+1), startpos[subseqCnt], endpos[subseqCnt]);
      subseqCnt = -1;
      break;
      } */
    else 
      subseqCnt++;
  }
  fclose(subFile);
  if (subseqCnt > 0){
      printf("Total %d column%s detected\n", subseqCnt, (subseqCnt==1 ? " is":"s are"));
      if (SINGLE_OUT)
	singleSubseq(seqTot, seqs, startpos, endpos, subseqCnt, alnFileName);
      else
	subseqOut(seqTot, seqs, startpos, endpos, subseqCnt, alnFileName);
  } // if
  else{
    printf("No columns can be outputted! Check your table file or MERGE option.\n");
    printf("You can also check file: %s to see the output column locations\n", tablefile);
  }
  if (SLIDER_OUT){
    // also add the slideargs
    sprintf(slidearg1, "%s.aln", alnFileName);
    printf("Executing Slider  %s %s...\n", slidearg1, subFileName);
    slider(slidearg1, subFileName, PC, slidewidth, slidewindow);      
    //slidearg2 is changed to subFileName
  }   
} // main

char *computeCons(int seqtot, char **seqs){
  int i,j;
  int freqs[27];
  int max;
  int imax;
  
  for (i=0;i<strlen(seqs[0]);i++){
    for (j=0;j<27;j++)
      freqs[j] = 0;
    
    for (j=0;j<seqtot;j++){
      if (seqs[j][i] == '-'){
	if (CONSENSUS_GAPS)
	  freqs[26]++;
      }
      else if (seqs[j][i]>='A' || seqs[j][i]<='Z')
	freqs[seqs[j][i]-'A']++;
      else{
	fprintf(stderr, "Illegal character: %c\n", seqs[j][i]);
	die();
      }      
    }
    max=0;imax=0;
    for (j=0;j<27;j++)
      if (freqs[j]>=max){
	max = freqs[j]; imax=j;
      }
    if (imax==26)
      ret[i]='-';
    else
      ret[i] = imax+'A'; 
  }
  ret[i] = 0;
  return ret;      
}


char *computeConsBackup(int seqtot, char **seqs){
  /* i will rewrite this in a more beautiful way sometime */
  int i,j;
  int a,c,g,t,u,d,e,f,h,ii,k,l,m,n,p,q,r,s,v,w,y,cons,gap;

  for (i=0;i<strlen(seqs[0]);i++){
    a=0; c=0; g=0; t=0; u=0; gap=0;
    d=0; e=0; f=0; h=0; ii=0; 
    k=0; l=0; m=0; n=0; p=0; 
    q=0; r=0; s=0; v=0; w=0; y=0; 
    for (j=0;j<seqtot;j++){
      switch (seqs[j][i]){
      case 'A':
	a++;
	break;
      case 'C':
	c++;
	break;
      case 'G':
	g++;
	break;
      case 'T':
	t++;
	break;
      case 'U':
	u++;
	break;
      case 'D':
	d++;
	break;
      case 'E':
	e++;
	break;
      case 'F':
	f++;
	break;
      case 'H':
	h++;
	break;
      case 'I':
	ii++;
	break;
      case 'K':
	k++;
	break;
      case 'L':
	l++;
	break;
      case 'M':
	m++;
	break;
      case 'N':
	n++;
	break;
      case 'P':
	p++;
	break;
      case 'Q':
	q++;
	break;
      case 'R':
	r++;
	break;
      case 'S':
	s++;
	break;
      case 'V':
	v++;
	break;
      case 'W':
	w++;
	break;
      case 'Y':
	y++;
	break;
      case '-':
	if (CONSENSUS_GAPS)
	  gap++;
	break;
      default:
	break;
      }// switch
      cons = a;
      if (c > cons)
	cons = c;
      if (g > cons)
	cons = g;
      if (t > cons)
	cons = t;
      if (u > cons)
	cons = u;
      if (d > cons)
	cons = d;
      if (e > cons)
	cons = e;
      if (f > cons)
	cons = f;
      if (h > cons)
	cons = h;
      if (ii > cons)
	cons = ii;
      if (k > cons)
	cons = k;
      if (l > cons)
	cons = l;
      if (m > cons)
	cons = m;
      if (n > cons)
	cons = n;
      if (p > cons)
	cons = p;
      if (q > cons)
	cons = q;
      if (r > cons)
	cons = r;
      if (s > cons)
	cons = s;
      if (v > cons)
	cons = v;
      if (w > cons)
	cons = w;
      if (y > cons)
	cons = y;
      if (gap > cons)
	cons = gap;
      if (cons == a)
	ret[i] = 'A';
      else if (cons == c)
	ret[i] = 'C';
      else if (cons == g)
	ret[i] = 'G';
      else if (cons == t)
	ret[i] = 'T';
      else if (cons == u)
	ret[i] = 'U';
      else if (cons == d)
	ret[i] = 'D';
      else if (cons == e)
	ret[i] = 'E';
      else if (cons == f)
	ret[i] = 'F';
      else if (cons == h)
	ret[i] = 'H';
      else if (cons == ii)
	ret[i] = 'I';
      else if (cons == k)
	ret[i] = 'K';
      else if (cons == l)
	ret[i] = 'L';
      else if (cons == m)
	ret[i] = 'M';
      else if (cons == n)
	ret[i] = 'N';
      else if (cons == p)
	ret[i] = 'P';
      else if (cons == q)
	ret[i] = 'Q';
      else if (cons == r)
	ret[i] = 'R';
      else if (cons == s)
	ret[i] = 'S';
      else if (cons == v)
	ret[i] = 'V';
      else if (cons == w)
	ret[i] = 'W';
      else if (cons == y)
	ret[i] = 'Y';
      else
	ret[i] = '-';
    } // for j
  } // for i
  ret[i]=0;
  return ret;
} // computecons

void printOut(int seqTot, char **seqs, char *inFile){
  int i, j, k, cnt;
  char *consensus;
  FILE *out;
  char fname[50];
  int nameSize;
  int endpoint;
  cnt = 0;

  sprintf(fname,"%s.out",inFile);

  printf("Generating OUT File [%s]... ", fname);
  fflush(stdout);
  out = fopen(fname, "w");
  
  consensus = computeCons(seqTot, seqs);

  
  nameSize=strlen("Consensus");
  for (i=0; i<seqTot; i++)
    if (strlen(names[i]) > nameSize)
      nameSize = strlen(names[i]);

  while (cnt < strlen(seqs[0])){
    for (k=0; k<nameSize+4; k++)
      fprintf(out," ");
    endpoint=cnt+59;
    if (endpoint > strlen(seqs[0]))
      endpoint = strlen(seqs[0]) + 1;
    
    if (endpoint < 100)
	fprintf(out, "%d %57d\n",(cnt+1),endpoint);
    else if (endpoint < 1000)
      fprintf(out, "%d %56d\n",(cnt+1),endpoint);
    else if (endpoint < 10000)
      fprintf(out, "%d %55d\n",(cnt+1),endpoint);
    else
      fprintf(out, "%d %54d\n",(cnt+1),endpoint);

    fprintf(out, "Consensus");
    for (k=0; k<=(nameSize-strlen("Consensus")); k++)
      fprintf(out, " ");
    fprintf(out, "   ");   
    for (i=cnt; i<cnt+59; i++)
      if (i<strlen(consensus))
	fprintf(out, "%c", consensus[i]);
    fprintf(out,"\n");
    for (j=0; j<seqTot; j++){
      fprintf(out,"%s",names[j]);
      for (k=0; k<=(nameSize-strlen(names[j])); k++)
	fprintf(out, " ");
      fprintf(out, "   ");
      for (i=cnt; i<cnt+59; i++){
	if (i >= strlen(seqs[0]))
	  break;
	if (seqs[j][i] == consensus[i])
	  fprintf(out, ".");
	else
	  fprintf(out, "%c", seqs[j][i]);
      } // for i
      fprintf(out, "\n");
    } // for j
    fprintf(out, "\n");
    cnt = i;
    fflush(out);
  } // while
  printf("[OK]\n%c%s generated in OUT format\n", numberSeven, fname);
  fclose(out);
} // printOut



int readAln(char *inFile){
  FILE *alnFile;
  char tempSeq[SEQ_LENGTH];
  char ch;
  if ((alnFile = fopen(inFile, "r")) == NULL){
    printf("Input File %s Not Found !!!\n",inFile);
    exit(0);
  } // if
  
  do {
    fscanf(alnFile, "%c", &ch);
  } while (isspace(ch));
  
  if (ch == '>'){
    printf("Sequences in FASTA format\n");
    return readFasta(alnFile);  
  }

  fscanf(alnFile, "%s", tempSeq);
  
  if (!strcmp(tempSeq, "LUSTAL")){
    printf("Sequences in CLUSTAL W format\n");
    fgets(tempSeq, SEQ_LENGTH, alnFile);
    return readClustal(alnFile);
  }
  else if (!strcmp(tempSeq, "NEXUS")){
    printf("Sequences in NEXUS format\n");
    fgets(tempSeq, SEQ_LENGTH, alnFile);
    return readNexus(alnFile);  
  }
  else if (!strcmp(tempSeq, "mega")){
    printf("Sequences in MEGA format\n");
    fgets(tempSeq, SEQ_LENGTH, alnFile);
    return readMega(alnFile);  
  }
  else if (isdigit(tempSeq[0])){
    printf("Sequences in PHYLIP format\n");
    return readPhylip(alnFile);  
  }
  else {
    printf("Error in Input File: %s. Not a valid alignment file ! \n",inFile);
    exit(0);
  } // if 
  return 0;
} //readaln

int readPhylip(FILE *alnFile){
  char *tempSeq;
  char *tempSeq2;
  int seqcnt, seqlen;
  int i,k,l;

  rewind(alnFile);
  fscanf(alnFile, "%d %d", &seqcnt, &seqlen);
  tempSeq = (char *) malloc(SEQ_LENGTH);
  tempSeq2 = (char *) malloc(SEQ_LENGTH);

  seqs = (char **) malloc((seqcnt+1) * sizeof(char *));
  for (i=0; i<=seqcnt; i++)
    seqs[i] = (char *) malloc(seqlen+2);

  names = (char **) malloc((seqcnt+1) * sizeof(char *));

  for (i=0; i<=seqcnt; i++)
    names[i] = (char *) malloc(100);
  
  for (i=0; i<=seqcnt; i++){
    seqs[i][0] = 0;
    names[i][0] = 0;
  }

  /* first pass */
  for (i=0;i<seqcnt;i++){
    fscanf(alnFile, "%s", tempSeq);
    strcpy(names[i], tempSeq);
    fgets(tempSeq, SEQ_LENGTH, alnFile);
    l=0;
    for (k=0;k<strlen(tempSeq)-1;k++)
      if (isalpha(tempSeq[k]) || tempSeq[k]=='-')
	tempSeq2[l++]=tempSeq[k];
    tempSeq2[l]=0;
    strcpy(seqs[i], tempSeq2);
  }

  i = 0;
  tempSeq[0]=0;

  while(fgets(tempSeq, SEQ_LENGTH, alnFile) > 0){
    l=0;
    for (k=0;k<strlen(tempSeq);k++)
      if (isalpha(tempSeq[k]) || tempSeq[k]=='-')
	tempSeq2[l++]=tempSeq[k];
    tempSeq2[l]=0;
    if (isalpha(tempSeq2[0]) || tempSeq2[0]=='-'){
      strcat(seqs[i], tempSeq2);
      i++;
    }
    if (i==seqcnt) // go back
      i=0;
    tempSeq[0]=0;
  }

  free(tempSeq);
  free(tempSeq2);
  maxLen = seqlen;
  return seqcnt;
}

int readClustal(FILE *alnFile){
  int cnt;
  char ch;
  char *tempSeq;
  char firstName[80];
  char seqName[80];
  char aseq[250];
  int i;
  int seqcnt=0, seqlen=0;

  cnt = 0;
  tempSeq = (char *) malloc(SEQ_LENGTH);
  
  strcpy(firstName, "first");
  while (1){ // read sequences
    if (!(fscanf(alnFile, "%s", tempSeq) > 0 ))
      break;
    if (tempSeq[0] == '*'){
      ch = 0;
      while (ch!='\n') // skip ********
	fscanf(alnFile,"%c",&ch);
      if ((!fscanf(alnFile, "%s %s", seqName, aseq) > 0 ))
	break;
    }
    else{
      strcpy(seqName, tempSeq);
      fscanf(alnFile, "%s", aseq);
    }
    if (!strcmp(firstName,"first")){
      strcpy(firstName,seqName);
      seqlen += strlen(aseq);  
    }
    else if (!strcmp(firstName,seqName)){
      seqcnt = 0;
      seqlen += strlen(aseq);
    }
    seqcnt++;
  } // while fscanf
  seqs = (char **) malloc((seqcnt+1) * sizeof(char *));
  for (i=0; i<=seqcnt; i++)
    seqs[i] = (char *) malloc(seqlen+2);

  names = (char **) malloc((seqcnt+1) * sizeof(char *));

  for (i=0; i<=seqcnt; i++)
    names[i] = (char *) malloc(100);
  
  for (i=0; i<seqcnt; i++){
    seqs[i][0] = 0;
    names[i][0] = 0;
  }
  rewind(alnFile);
  strcpy(firstName, "first");
  ch=0;
  free(tempSeq);

  tempSeq = (char *) malloc(seqlen+2);
  while (ch!='\n') // skip header
    fscanf(alnFile,"%c",&ch);
  printf("Reading Sequences ... ");
  fflush(stdout);
  while ((fscanf(alnFile, "%s", tempSeq) > 0 )){
    if (tempSeq[0] == '*'){
      ch = 0; 
      fgets(tempSeq, SEQ_LENGTH, alnFile);
      if (!(fscanf(alnFile, "%s", seqName) > 0 ))
	break;
      fscanf(alnFile, "%s", aseq);
    }
    else{
      strcpy(seqName, tempSeq);
      fscanf(alnFile, "%s", aseq);
    }
    if (!strcmp(firstName,"first"))
      strcpy(firstName,seqName);
    else if (!strcmp(firstName,seqName))
      cnt = 0;
    if (cnt < seqcnt){
      strcat(seqs[cnt], aseq);
      if (names[cnt][0] == 0)
	strcpy(names[cnt], seqName);
    }
    cnt++;
  } // while fscanf

  maxLen = strlen(seqs[0]);

  free(tempSeq); 
  printf("[OK] %d sequences of length %d are read.\n",cnt,maxLen);
  return cnt; // return num of sequences
} // readclustal

int readNexus(FILE *alnFile){
  int cnt;
  char ch;
  char *tempSeq;
  char firstName[80];
  char seqName[80];
  char aseq[200];
  int i, j, k;
  int seqcnt=0, seqlen=0;
  
  tempSeq = (char *) malloc(SEQ_LENGTH);
  strcpy(tempSeq, "");
  strcpy(firstName, "first");

  do {
    fscanf(alnFile, "%s", tempSeq);
  } while (!strstr(tempSeq, "ntax") && !strstr(tempSeq, "nchar"));

  if (tempSeq[strlen(tempSeq)-1] == ';')
    tempSeq[strlen(tempSeq)-1] = 0;
  
  if (!strstr(tempSeq, "ntax") || !strstr(tempSeq,"nchar")){
    if (strstr(tempSeq, "nchar"))
      k = 0;
    else
      k = 1;
    if (!strcmp(tempSeq, "ntax") || !strcmp(tempSeq, "nchar")){ // ntax =number
      fscanf(alnFile, "%s", tempSeq); 
      if (!strcmp(tempSeq, "=")){ // ntax = number
	fscanf(alnFile, "%s", tempSeq);
	sscanf(tempSeq, "%d", &cnt);
      } 
      else{ // ntax =number 
	for (j=0; j<strlen(tempSeq)-1; j++)
	  tempSeq[j] = tempSeq[j+1];
	tempSeq[j] = 0;
	sscanf(tempSeq, "%d", &cnt);
      }
    } // if ntax
    else { // ntax=number
      for (j=0; j<strlen(tempSeq)-6+k; j++)
	tempSeq[j] = tempSeq[j+6-k];
      tempSeq[j] = 0;
      sscanf(tempSeq, "%d", &cnt);
    } // else
  } // if includes ntax
  if (k)
    seqcnt = cnt;
  else
    seqlen = cnt;
  do {
    fscanf(alnFile, "%s", tempSeq);
  } while (!strstr(tempSeq, "ntax") && !strstr(tempSeq, "nchar"));
  if (tempSeq[strlen(tempSeq)-1] == ';')
    tempSeq[strlen(tempSeq)-1] = 0;
  if (!strstr(tempSeq, "ntax") || !strstr(tempSeq,"nchar")){
    if (strstr(tempSeq,"nchar"))
      k = 0;
    else
      k = 1;
    if (!strcmp(tempSeq, "ntax") || !strcmp(tempSeq, "nchar")){ // ntax =number
      fscanf(alnFile, "%s", tempSeq); 
      if (tempSeq[strlen(tempSeq)-1] == ';')
	tempSeq[strlen(tempSeq)-1] = 0;
      if (!strcmp(tempSeq, "=")){ // ntax = number
	fscanf(alnFile, "%s", tempSeq);
	if (tempSeq[strlen(tempSeq)-1] == ';')
	  tempSeq[strlen(tempSeq)-1] = 0;
	sscanf(tempSeq, "%d", &cnt);
      } 
      else{ // ntax =number 
	for (j=0; j<strlen(tempSeq)-1; j++)
	  tempSeq[j] = tempSeq[j+1];
	tempSeq[j] = 0;
	sscanf(tempSeq, "%d", &cnt);
      }	
    } // if ntax
    else { // ntax=number
      for (j=0; j<strlen(tempSeq)-6+k; j++)
	tempSeq[j] = tempSeq[j+6-k];
      tempSeq[j] = 0;
      sscanf(tempSeq, "%d", &cnt);
    } // else
  } // if includes ntax
  if (k)
    seqcnt = cnt;
  else
    seqlen = cnt;

  seqs = (char **) malloc((seqcnt+1) * sizeof(char *));
  for (i=0; i<=seqcnt; i++)
    seqs[i] = (char *) malloc(seqlen+2);

  names = (char **) malloc((seqcnt+1) * sizeof(char *));
  for (i=0; i<=seqcnt; i++)
    names[i] = (char *) malloc(100);

  cnt = 0;

  for (i=0; i<seqcnt; i++){
    seqs[i][0] = 0;
    names[i][0] = 0;
  }

  while (strcmp(tempSeq, "matrix"))
    fscanf(alnFile, "%s", tempSeq);

  strcpy(firstName, "first");
  ch=0;
  free(tempSeq);
  tempSeq = (char *) malloc(seqlen+2);
  while (ch!='\n') // skip header
    fscanf(alnFile,"%c",&ch);
  printf("Reading Sequences ... ");
  fflush(stdout);
  while (1){ // read sequences
    if (!(fscanf(alnFile, "%s", tempSeq) > 0 ))
      break;
    if (!strcmp(tempSeq, ";"))
      break;
    strcpy(seqName, tempSeq);
    fscanf(alnFile, "%s", aseq);
    if (!strcmp(firstName,"first"))
      strcpy(firstName,seqName);
    else if (!strcmp(firstName,seqName))
      cnt = 0;
    strcat(seqs[cnt], aseq);
    if (names[cnt][0] == 0)
      strcpy(names[cnt], seqName);
    cnt++;
  } // while fscanf
  maxLen = strlen(seqs[0]);
  free(tempSeq); 
  for (i=0; i<seqcnt; i++){
    if (names[i][0] == '\''){
      for (j=0;j<strlen(names[i])-2;j++)
	names[i][j] = names[i][j+1];
      names[i][j] = 0;
    }
  } // for i
  printf("[OK] %d sequences read\n",seqcnt);
  return seqcnt; // return num of sequences
} // readnexus

int readMega(FILE *alnFile){
  int cnt;
  char ch;
  char *tempSeq;
  int i, j;
  int seqcnt=0, seqlen=0;
  cnt = 0;
  
  while (fscanf(alnFile, "%c", &ch)>0){
    if (ch == '#'){
      cnt++;
      if (cnt==1){
	do{
	  fscanf(alnFile, "%c", &ch);
	  seqlen++;
	} while (ch != '\r' && ch!='\n');
      } // if cnt==1
    } // if #
  } // while
  seqcnt = cnt;
  rewind(alnFile);

  tempSeq = (char *) malloc(seqlen+2);

  seqs = (char **) malloc((seqcnt+1) * sizeof(char *));
  for (i=0; i<=seqcnt; i++)
    seqs[i] = (char *) malloc(seqlen+2);

  names = (char **) malloc((seqcnt+1) * sizeof(char *));

  for (i=0; i<=seqcnt; i++)
    names[i] = (char *) malloc(100);

  cnt = 0;

  for (i=0; i<seqcnt; i++){
    seqs[i][0] = 0;
    names[i][0] = 0;
  }

  do{
    fscanf(alnFile, "%s", tempSeq);
  } while (strcmp(tempSeq, "#mega"));

  do{
    fscanf(alnFile, "%s", tempSeq);
  } while (tempSeq[0] != '#');
  // tempseq is the first sequence name
  strcpy(names[0], tempSeq);
  fscanf(alnFile, "%s", seqs[0]); // first sequence
  for (i=1; i<seqcnt; i++)
    fscanf(alnFile, "%s %s", names[i], seqs[i]); // name, sequence

  maxLen = strlen(seqs[0]);
  free(tempSeq); 
  for (i=0; i<seqcnt; i++){
    if (names[i][0] == '#'){
      for (j=0;j<strlen(names[i])-1;j++)
	names[i][j] = names[i][j+1];
      names[i][j] = 0;
    }
  } // for i

  printf("[OK] %d sequences read\n",seqcnt);

  return seqcnt;

} // readMega

int readFasta(FILE *alnFile){
  int cnt;
  char ch; 
  int i;
  int seqcnt=0, seqlen=0;

  cnt = 0; i=0;

  rewind(alnFile);
  while (fscanf(alnFile, "%c", &ch) > 0){
    if (ch == '>')
      cnt++;
    if (cnt == 1){
      while (fscanf(alnFile, "%c", &ch) > 0){
	if (ch!='>' && ch!='\r' && ch!='\n')
	  i++;
	if (ch == '>'){
	  cnt++;
	  break;
	}
      }
      seqlen = i;
    }
  }

  seqcnt = cnt;
  
  rewind(alnFile);
  printf("seqcnt: %d seqlen: %d\n", seqcnt, seqlen);

  seqs = (char **) malloc((seqcnt+1) * sizeof(char *));
  
  for (i=0; i<=seqcnt; i++)
    seqs[i] = (char *) malloc(seqlen+2);
  
  names = (char **) malloc((seqcnt+1) * sizeof(char *));

  for (i=0; i<=seqcnt; i++)
    names[i] = (char *) malloc(100);

  for (i=0; i<seqcnt; i++){
    seqs[i][0] = 0;
    names[i][0] = 0;
  }

  cnt = -1; 
  while (fscanf(alnFile, "%c", &ch) > 0){
    if (ch == '>'){
      cnt++;
      fgets(names[cnt], SEQ_LENGTH, alnFile);
      names[cnt][strlen(names[cnt])-1] = 0;
      if (VERBOSE)
	printf("seq-%d: %s\n", cnt, names[0]);
    }
    i = 0;
    if (cnt != 0)
      seqs[cnt][i++] = ch;
    do{
      if (!(fscanf(alnFile, "%c", &ch) > 0))
	break;
      if (ch!='>' && ch!='\r' && ch!='\n')
	seqs[cnt][i++] = ch;
    } while (ch != '>');
    seqs[cnt][i] = 0;
    if (ch == '>'){
      cnt++;
      if (cnt != seqcnt){
	fgets(names[cnt], SEQ_LENGTH, alnFile);
	names[cnt][strlen(names[cnt])-1] = 0;
	if (VERBOSE)
	  printf("seq-%d: %s\n", cnt, names[cnt]);
      }
    } // if
  } // while
	    
  maxLen = strlen(seqs[0]);
  printf("[OK] %d sequences read\n",seqcnt);
  return seqcnt;
} // readFasta

void subseqOut(int seqTot, char **seqs, int startpos[], int endpos[], int subseqCnt, char *inFile){
  char outfname[100];
  int i, j, cnt, loop, ii, jj;
  int nameSize=0;
  char **pSubseq;
  int spos[2];
  int epos[2];

  for (i=0; i<seqTot; i++)
    if (strlen(names[i]) > nameSize)
      nameSize = strlen(names[i]);
  
  for (loop=0;loop<subseqCnt;loop++){
    pSubseq = (char **) malloc(seqTot * sizeof(char *));
    for (i=0; i<seqTot; i++)
      pSubseq[i] = (char *) malloc(endpos[loop] - startpos[loop] + 1);

    cnt = startpos[loop]-1;
    sprintf(outfname, "%s_%d-%d", inFile, startpos[loop], endpos[loop]);

    while (cnt <= endpos[loop]-1){
      for (j=0, jj=0; j<seqTot; j++, jj++){

	ii = cnt - startpos[loop] + 1;
	for (i=cnt; i<cnt+59; i++){
	  if (i >= endpos[loop])
	    break;
	  pSubseq[jj][ii++] =  seqs[j][i];
	} // for i
      } // for j
      cnt = i;
    } // while
    
    pSubseq[jj-1][ii] = 0;
    spos[0] = 1;
    epos[0] = ii;
    singleSubseq(seqTot, pSubseq, spos, epos, 1, outfname);
    
    if (IDENTITY_OUT)
      printOut(seqTot, pSubseq, outfname);
    if (NEXUS_OUT)
      printNexus(seqTot, pSubseq, outfname);
    if (MEGA_OUT)
      printMega(seqTot, pSubseq, outfname);
    if (FASTA_OUT)
      printFasta(seqTot, pSubseq, outfname);
    if (PHYLIP_OUT)
      printPhylip(seqTot, pSubseq, outfname);
    if (HTML_OUT)
      printHTML(seqTot, pSubseq, outfname);
  } // for loop
} // subseqOut

void singleSubseq(int seqTot, char **seqs, int startpos[], int endpos[], int subseqCnt, char *inFile){
  char **pSubseq;
  int i, j, k, m, cnt;
  FILE *out;
  char outfname[50];
  int bool;
  char first;
  int nameSize=0;
  int start, end;
  char *consensus;

  for (i=0;i<subseqCnt;i++){ // sort
    for (j=0;j<subseqCnt;j++)
      if (i<j && startpos[i] > startpos[j]){
	start = startpos[i];
	end = endpos[i];
	startpos[i] = startpos[j];
	endpos[i] = endpos[j];
	startpos[j] = start;
	endpos[j] = end;	
      } // if
  } // for i
  
  for (i=0;i<subseqCnt-1;i++){ // check overlaps
    if (endpos[i] >= startpos[i+1]){
      printf("Subsequence Start & Positions Overlap !!! Check Input !!!\n");
      return;
    } // if
  } // for check
  
  if (INCLUDE_CONSENSUS)
    nameSize = strlen("Consensus");
  for (i=0; i<seqTot; i++)
    if (strlen(names[i]) > nameSize)
      nameSize = strlen(names[i]);  
  
  cnt = 0;
  pSubseq = (char **) malloc(seqTot * sizeof(char *));
  for (i=0; i<seqTot; i++)
    pSubseq[i] = (char *) malloc(maxLen+1); 
  for (i=0;i<subseqCnt;i++){
    for (j=0;j<seqTot;j++){
      m = cnt;
      for (k=startpos[i];k<=endpos[i];k++)
	pSubseq[j][m++] = seqs[j][k-1];
    } // for all seqs
    cnt = m;
  } // for subseqs
  for (j=0;j<seqTot;j++)
    pSubseq[j][m] = 0;

  strcpy(outfname, inFile);
  for (i=strlen(outfname)-1;i>=0;i++)
    if (outfname[i]=='.'){
      outfname[i] = 0;
      break;
    }
  strcat(outfname, ".sub.aln");
    
  if (CLUSTALW_OUT){

    printf("Generating Subsequences File %s ",outfname);
    fflush(stdout);
    
    out = fopen(outfname, "w");

    if (INCLUDE_CONSENSUS){
      consensus = (char *) malloc(strlen(pSubseq[0])+1);
      strcpy(consensus, computeCons(seqTot, pSubseq));
    }
    
    
    fprintf(out, "CLUSTAL W (1.83) multiple sequence alignment\n\n\n"); // header
    cnt = 0;
    while (cnt <= m-1){
      if (INCLUDE_CONSENSUS){
	fprintf(out, "Consensus");
	for (k=0; k<=(nameSize-strlen("Consensus")); k++)
	  fprintf(out, " ");
	fprintf(out, "   ");
	for (i=cnt; i<cnt+59; i++){
	  if (i >= m)
	    break;
	  fprintf(out, "%c", consensus[i]);
	}
	fprintf(out, "\n");
      }
      for (j=0; j<seqTot; j++){
	fprintf(out, "%s",names[j]);
	for (k=0; k<=(nameSize-strlen(names[j])); k++)
	  fprintf(out, " ");
	fprintf(out, "   ");
	for (i=cnt; i<cnt+59; i++){
	  if (i >= m)
	    break;
	  fprintf(out, "%c", pSubseq[j][i]);
	} // for i
	fprintf(out, "\n");
      } // for j
      for (k=0;k<nameSize;k++)
	fprintf(out," ");
      fprintf(out,"    ");
      
      for (j=cnt;j<i;j++){ // place perfect matches
	bool = 1;
	first = pSubseq[0][j];
	for (k=1; k <seqTot; k++)
	  if (pSubseq[k][j] != first)
	    bool = 0;
	if (bool)
	  fprintf(out,"*");
	else
	  fprintf(out," ");
      }      // for loop of placing perfect matches
      fprintf(out, "\n\n");
      cnt = i;
    } // while
    
    printf("[OK] \n");
    bool = 1;
    if (INCLUDE_CONSENSUS)
      free(consensus);
    fclose(out);
  } // if CLUSTALW_OUT
  
  if (IDENTITY_OUT)
    printOut(seqTot, pSubseq, outfname);
  if (NEXUS_OUT)
    printNexus(seqTot, pSubseq, outfname);
  if (MEGA_OUT)
    printMega(seqTot, pSubseq, outfname);
  if (FASTA_OUT)
    printFasta(seqTot, pSubseq, outfname);
  if (PHYLIP_OUT)
      printPhylip(seqTot, pSubseq, outfname);
  if (HTML_OUT)
      printPhylip(seqTot, pSubseq, outfname);

  for (i=0; i<seqTot; i++)
    free(pSubseq[i]);

  free(pSubseq);
  
  
} // singleSubseq 

void printHTML(int seqTot, char **seqs, char *inFile){
  int i, j, k, m, cnt;
  FILE *out;
  char outfname[50];
  bool willStar;
  char first;
  int nameSize=0;
  char *consensus;
  

  if (INCLUDE_CONSENSUS)
    nameSize = strlen("Consensus");
  for (i=0; i<seqTot; i++)
    if (strlen(names[i]) > nameSize)
      nameSize = strlen(names[i]);  
  
  strcpy(outfname, inFile);
  for (i=strlen(outfname)-1;i>=0;i++)
    if (outfname[i]=='.'){
      outfname[i] = 0;
      break;
    }

  strcat(outfname, ".sub.html");
    
  printf("Generating HTML File %s ",outfname);
  fflush(stdout);
  
  out = fopen(outfname, "w");
  
  consensus = (char *) malloc(strlen(seqs[0])+1);
  strcpy(consensus, computeCons(seqTot, seqs));
  
  fprintf(out,"<html><head><title>Multiple Sequence Alignment %s</title></head>\n", outfname);
  fprintf(out,"<body>\n<pre>\n");
  fprintf(out, "Multiple Sequence Alignment based on CLUSTAL W Format\n\n\n"); // header
  cnt = 0;
  m = strlen(seqs[0]);
  while (cnt < m){
    if (INCLUDE_CONSENSUS){
      fprintf(out, "Consensus");
      for (k=0; k<=(nameSize-strlen("Consensus")); k++)
	fprintf(out, " ");
      fprintf(out, "   ");
      for (i=cnt; i<cnt+49; i++){
	if (i >= m)
	  break;
	fprintf(out, "%c", consensus[i]);
      }
      fprintf(out, "\n");
    }
    for (j=0; j<seqTot; j++){
      fprintf(out, "%s",names[j]);
	for (k=0; k<=(nameSize-strlen(names[j])); k++)
	  fprintf(out, " ");
	fprintf(out, "   ");
	for (i=cnt; i<cnt+49; i++){
	  if (i >= m)
	    break;
	  willStar = 1;
	  first = seqs[0][i];
	  for (k=1; k <seqTot; k++)
	    if (seqs[k][i] != first)
	      willStar = 0;
	  if (willStar)
	    fprintf(out, "%c", seqs[j][i]);
	  else{
	    if (seqs[j][i] == consensus[i])
	      fprintf(out, "<span style=\"background-color:#00FFFF\">%c</span>", seqs[j][i]);
	    else
	      fprintf(out, "<span style=\"background-color:#FF33E9\">%c</span>", seqs[j][i]);
	  }
	} // for i
	fprintf(out, "\n");
    } // for j
    for (k=0;k<nameSize;k++)
      fprintf(out," ");
    fprintf(out,"    ");
    
    for (j=cnt;j<i;j++){ // place perfect matches
      willStar = 1;
      first = seqs[0][j];
      for (k=1; k <seqTot; k++)
	if (seqs[k][j] != first)
	  willStar = 0;
      if (willStar)
	fprintf(out,"*");
      else
	fprintf(out," ");
    }      // for loop of placing perfect matches
    fprintf(out, "\n\n");
    cnt = i;
  } // while
  
  fprintf(out,"</pre></body></html>\n");

  printf("[OK]\n%c%s generated in HTML format\n", numberSeven, outfname);
  willStar = 1;
  if (INCLUDE_CONSENSUS)
    free(consensus);
  fclose(out);
  
  free(consensus);
  
} // printHTML


void printNexus(int seqTot, char **seqs, char *inFile){
  int i, j, k, cnt;
  FILE *out;
  char fname[50];
  int nameSize;
  int len;
  char **names2;
  int quote[seqTot];
  char *consensus;

  bool doQuote = FALSE;
  cnt = 0;
  
  sprintf(fname,"%s.nexus",inFile);
  out = fopen(fname, "w");

  if (INCLUDE_CONSENSUS){
    consensus = (char *) malloc(strlen(seqs[0])+1);
    strcpy(consensus, computeCons(seqTot, seqs));
  }


  names2 = (char **) malloc(seqTot * sizeof(char *));
  for (i=0; i<seqTot; i++){
    names2[i] = (char *) malloc(100);
    cnt = 0;
    quote[i] = 0;
    for (j=0;j<strlen(names[i]);j++){
      if (!myisalnum(names[i][j])){
	quote[i] = 1;
	doQuote = TRUE;
      }
      if (names[i][j] != '\'')
	names2[i][cnt++] = names[i][j];
      else{
	names2[i][cnt++] = '\'';
	names2[i][cnt++] = '\'';
      }
    }
    names2[i][cnt] = 0;
  } 

  nameSize=0;
  if (INCLUDE_CONSENSUS)
    nameSize = strlen("Consensus");
  for (i=0; i<seqTot; i++){
    if (strlen(names[i]) > nameSize)
      nameSize = strlen(names[i]);
  }

  len = strlen(seqs[0]);

  printf("Generating NEXUS File ... ");
  fflush(stdout);
  
  fprintf(out,"#NEXUS");
  fprintf(out,"\n\nbegin data;\n");
  if (INCLUDE_CONSENSUS)
    fprintf(out,"\tdimensions ntax=%d nchar=%d;\n",(seqTot+1),len);
  else
    fprintf(out,"\tdimensions ntax=%d nchar=%d;\n",seqTot,len);
  fprintf(out,"\tformat datatype=DNA interleave gap=-;\n");
  fprintf(out,"\tmatrix\n"); 
  cnt = 0;
  while (cnt < len){
    if (INCLUDE_CONSENSUS){
      fprintf(out, "Consensus");
      for (k=0; k<=(nameSize-strlen("Consensus")); k++)
	fprintf(out, " ");
      fprintf(out, "   ");
      if (doQuote)
	fprintf(out, "  ");
      for (i=cnt; i<cnt+59; i++){
	if (i >= len)
	  break;
	fprintf(out, "%c", consensus[i]);
      }
      fprintf(out, "\n");
    }
    for (j=0; j<seqTot; j++){      
      if (quote[j])
	fprintf(out,"\'%s\'",names2[j]);      
      else
	fprintf(out,"%s",names2[j]);      
      
      for (k=0; k<=(nameSize-strlen(names2[j])); k++)
	fprintf(out, " ");
      fprintf(out, "   ");
      
      for (i=cnt; i<cnt+59; i++){
	if (i >= len)
	  break;
	  fprintf(out, "%c", seqs[j][i]);
      } // for i
      fprintf(out, "\n");
    } // for j
    fprintf(out, "\n");
    cnt = i;
    fflush(out);
  } // while
  fprintf(out,";\nend;\n");
  for (i=0; i<seqTot; i++)
    free(names2[i]);
  free(names2);

  if (INCLUDE_CONSENSUS)
    free(consensus);

  printf("[OK]\n%c%s generated in NEXUS format\n", numberSeven, fname);
  fclose(out);
} // printNexus

void printPhylip(int seqTot, char **seqs, char *inFile){
  int i, j, k, cnt;
  FILE *out;
  char fname[50];
  int nameSize;
  int len;
  char *consensus;

  printf("Generating PHYLIP File ... ");
  fflush(stdout);
  
  sprintf(fname,"%s.phy",inFile);
  out = fopen(fname, "w");

  len = strlen(seqs[0]);
  nameSize = 10; // phylip4 allows max name length 10
  /*
  for (i=0; i<seqTot; i++){
    if (strlen(names[i]) > nameSize)
      nameSize = strlen(names[i]);
  }
  */
  
  if (INCLUDE_CONSENSUS){
    consensus = (char *) malloc(strlen(seqs[0])+1);
    strcpy(consensus, computeCons(seqTot, seqs));
    fprintf(out,"\t%d\t%d\n",(seqTot+1),len);
    if (nameSize<9)
      nameSize=9;
  }
  else
    fprintf(out,"%d %d\n",seqTot,len);    

  cnt = 0;
  while (cnt < len){
    if (INCLUDE_CONSENSUS){
      if (cnt == 0)
	fprintf(out, "Consensus ");
      else
	fprintf(out, "          ");
      /*
      for (k=0; k<=(nameSize-strlen("Consensus")); k++)
      fprintf(out, " ");
      fprintf(out, "   ");
      */
      for (i=cnt; i<cnt+50; i++){
	if (i >= len)
	  break;
	if (i!=0 && i%10==0)
	  fprintf(out, " ");
	fprintf(out, "%c", consensus[i]);
      }
      fprintf(out, "\n");
    }
    for (j=0; j<seqTot; j++){
      if (cnt == 0){
	for (i=0;i<nameSize-1 && i<strlen(names[j]);i++)
	  fprintf(out, "%c",names[j][i]);
	for (k=1; k<=(nameSize-i); k++)
	  fprintf(out, " ");
      }
      else{
	for (k=1; k<nameSize; k++)
	  fprintf(out, " ");
      }
      //fprintf(out, "   ");
      for (i=cnt; i<cnt+50; i++){
	if (i >= len)
	  break;
	if (i!=0 && i%10==0)
	  fprintf(out, " ");
	fprintf(out, "%c", seqs[j][i]);
      } // for i
      fprintf(out, "\n");
    } // for j
    for (k=0;k<nameSize;k++)
      fprintf(out," ");
    //fprintf(out,"    ");
    
    fprintf(out, "\n");
    cnt = i;
  } // while
  
  printf("[OK] \n");

  if (INCLUDE_CONSENSUS)
    free(consensus);
  fclose(out);
  
  printf("[OK]\n%c%s generated in PHYLIP format\n", numberSeven, fname);
} // printPhylip


void printMega(int seqTot, char **seqs, char *inFile){
  int i, j, k, cnt;
  FILE *out;
  char fname[50];
  int nameSize;
  int len;
  char *consensus;
  
  cnt = 0;

  sprintf(fname,"%s.mega",inFile);

  printf("Generating MEGA File ... ");
  fflush(stdout);
  
  out = fopen(fname, "w");
  
  if (INCLUDE_CONSENSUS){
    consensus = (char *) malloc(strlen(seqs[0])+1);
    strcpy(consensus, computeCons(seqTot, seqs));
  }

  nameSize=0;
  if (INCLUDE_CONSENSUS)
    nameSize = strlen("Consensus");
  for (i=0; i<seqTot; i++)
    if (strlen(names[i]) > nameSize)
      nameSize = strlen(names[i]);
  len = strlen(seqs[0]);
  fprintf(out,"#mega\n");
  fprintf(out,"TITLE: Multiple Alignment Manipulator\n\n");

  if (INCLUDE_CONSENSUS){
    fprintf(out, "#Consensus");
    for (k=0; k<=(nameSize-strlen("Consensus")); k++)
      fprintf(out, " ");
    fprintf(out, "   ");
    fprintf(out, "%s\r", consensus);
  }
  
  for (j=0; j<seqTot; j++){
    fprintf(out,"#%s",names[j]);
    for (k=0; k<=(nameSize-strlen(names[j])); k++)
      fprintf(out, " ");
    fprintf(out, "   ");
    fprintf(out, "%s\r", seqs[j]);
  } // for j
  fprintf(out, "\n");
  fflush(out);
  if (INCLUDE_CONSENSUS)
    free(consensus);
  printf("[OK]\n%c%s generated in MEGA format\n", numberSeven, fname);
  fclose(out);
} // printMega

void printFasta(int seqTot, char **seqs, char *inFile){
  int i;
  int cnt;
  FILE *out;
  char fname[SEQ_LENGTH];
  char *consensus;
  
  sprintf(fname,"%s.fa",inFile);

  printf("Generating FASTA File [%s] ... ", fname);
  fflush(stdout);
  
  out = fopen(fname, "w");

  if (INCLUDE_CONSENSUS){
    //    consensus = (char *) malloc(strlen(seqs[0])+1);
    //strcpy(consensus, computeCons(seqTot, seqs));
    consensus = computeCons(seqTot, seqs);
  }

  if (INCLUDE_CONSENSUS){
    fprintf(out, ">%s Consensus\n", inFile);
    cnt = 0;
    while (cnt < strlen(consensus)){
      fprintf(out, "%c", consensus[cnt++]);
      if (cnt % 60 == 0)
	fprintf(out, "\n");
    } // while
    fprintf(out, "\n");
  }

  for (i = 0 ; i < seqTot; i++){
    if (seqTot == 1 && strstr(inFile, "consensus"))
      fprintf(out, ">%s Consensus\n", inFile);
    else
      fprintf(out, ">%s\n", names[i]);
    cnt = 0;
    while (cnt < strlen(seqs[i])){
      fprintf(out, "%c", seqs[i][cnt++]);
      if (cnt % 60 == 0)
	fprintf(out, "\n");
    } // while
    fprintf(out, "\n");
  } // for
  fprintf(out, "\n");
  fflush(out);
  fclose(out);
  printf("[OK]\n%c%s generated in FASTA format\n", numberSeven, fname);
} // printFasta


int myisalnum(char ch){
  switch (ch){
  case '(':
  case ')':
  case '[':
  case ']':
  case '{':
  case '}':
  case '/':
  case '\\':
  case ',':
  case ';':
  case ':':
  case '=':
  case '*':
  case '\"':
  case '\'':
  case '`':
  case '+':
  case '-':
  case '<':
  case '>':
    return 0;
    break;
  } // switch
  return 1;
} // myisalnum

void convertFormat(char alnFileName[]){
  int startpos[1];
  int endpos[1];
  startpos[0]=1;
  endpos[0]=maxLen;
  if (IDENTITY_OUT)
    printOut(seqTot, seqs, alnFileName);
  if (NEXUS_OUT)
    printNexus(seqTot, seqs, alnFileName);
  if (MEGA_OUT)
    printMega(seqTot, seqs, alnFileName);
  if (FASTA_OUT)
    printFasta(seqTot, seqs, alnFileName);
  if (CLUSTALW_OUT)
    singleSubseq(seqTot, seqs, startpos, endpos, 1, alnFileName);
  if (PHYLIP_OUT)
    printPhylip(seqTot, seqs, alnFileName);
  if (HTML_OUT)
    printHTML(seqTot, seqs, alnFileName);
}

void alnStats(char alnFileName[]){
  FILE *stats;
  char fname[200];
  int i,j;
  int max=0;
  int min=1000000;
  int len;
  int wmax=0, wmin=0;
  int totallen=0;
  int stars=0; // consensus columns
  bool same = TRUE;
  sprintf(fname, "%s.stats", alnFileName);
  stats = fopen(fname, "w");
  fprintf(stats, "ALIGNMENT STATISTICS FOR %s\n", alnFileName);
  fprintf(stats, "------------------------------------------------------\n\n");
  fprintf(stats, "Number of Sequences: %d\n", seqTot);
  fprintf(stats, "Length of Alignment: %d\n", maxLen);
  for (i=0;i<seqTot;i++){
    len = sequence_len(seqs[i]);
    if (len > max){
      max = len;
      wmax = i;
    }
    if (len < min){
      min = len;
      wmin = i;
    }
    totallen += len;
  }
  fprintf(stats, "Max Sequence Length: %d Sequence: %s\n", max, names[wmax]);
  fprintf(stats, "Min Sequence Length: %d Sequence: %s\n", min, names[wmin]);
  fprintf(stats, "Average Sequence Length: %f\n", ((float)(totallen)/(float)(seqTot)));

  for (j=0;j<maxLen;j++){
    same = TRUE; 
    for (i=1;i<seqTot;i++){
      if (seqs[i][j] != seqs[0][j]){
	same = FALSE;
	break;
      }
    }
    if (same)
      stars++;
  }
  fprintf(stats, "Number of Conserved Columns: %d, Percentage: %f\n", stars, ((float)stars/(float)maxLen));
  fprintf(stats, "\n------------------------------------------------------\n\n");
  fprintf(stats, "TABLE OF SEQUENCES\n");
  fprintf(stats, "------------------------------------------------------\n\n");
  fprintf(stats, "%20s%10s%7s\n", "SEQ NAME", "LENGTH", "GC %");
  for (i=0;i<seqTot;i++){
      fprintf(stats, "%20s%10d%7.2f\n", names[i], sequence_len(seqs[i]), gc_cont(seqs[i]));
  }

  if (seqTot == 2){
    fprintf(stats, "\n Mutations:\n");
    fprintf(stats, "Pos\t%s\t%s\n", names[0], names[1]);
    for (j=0;j<maxLen;j++){
      if (seqs[0][j] != seqs[1][j]){
	fprintf(stats, "%d\t%c\t%c\n", (j+1), seqs[0][j], seqs[1][j]);
      }
    }
  }
}

int sequence_len(char *s){
  int i;
  int len=0;
  for (i=0;i<strlen(s);i++)
    if (s[i] != '-')
      len++;
  return len;
}

float gc_cont(char *s){
  int i;
  int gc=0;
  for (i=0;i<strlen(s);i++)
    if (s[i] == 'G' || s[i] == 'C')
      gc++;
  return ((float)gc / (float)sequence_len(s));
}
