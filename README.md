
	MaM : Multiple alignment Manipulator
	
	Implemented by: Can ALKAN & Eray TUZUN
		
	[    calkan@gmail.com   ]
	[  eraytuzun@gmail.com  ]

	Last Update: March 20, 2006
	Version 1.4.2

Usage: mam [Alignment File] [option1][=on/=off] [option2=on/off] ... [option n=on/off]

Note:
	MaM saves its configuration in $HOME/.mam-config file. Run MaM

	once to create the default configuration file, and edit it if necessary.


Options:

-exonfile     : Exon table file or cDNA file.

		Ex: -exonfile=default , -exonfile=NPIP

-update       : Update the alignment/sequence coordinates On/Off.

                Set to off if you want to use the alignment coordinates in the tablefile.

		- Works only with -program=table

-column	      : Toggle Single/Multiple File Output

		Ex -column=single  or  -column=multiple

-merge	      : Toggle MAX/MIN

		Ex: -merge=max  or  -merge=min

-keep	      : Toggle KEEP/TOSS

		Ex: -keep=on  or  keep=off

-program      : Select program. One of:

		-program=crossmatch

		-program=repeatmasker

		-program=sim4

		-program=table

		-program=convert

		-program=none (to be able to run -alnstats and/or -consensus alone)

-slider	      : Toggle Slider On/Off

-pc	      : Select One of:

		Pairwise Deletion (-pc=p),

		Complete Deletion (-pc=c),

		Parsimony Score   (-pc=s)

-sw	      : Select Slide Width (-sw=10)

-ww	      : Select Window Width (-ww=100)

-clustal      : Toggle Clustal Output Format On/Off

-nexus	      : Toggle Nexus Output Format On/Off

-mega	      : Toggle Mega Output Format On/Off

-fasta        : Toggle Fasta Output Format On/Off

-phylip	      : Toggle Phylip Output Format On/Off

-html	      : Toggle HTML (two-color marking) Output Format On/Off

-identity     : Toggle Consensus Identity  Output Format On/Off

-alnstats     : Toggle Alignment Statistics File Dump On/Off

-consensus    : Toggle Consensus Output On/Off

-cgaps        : Toggle Gaps in Consensus Output On/Off (if Off, second most frequent in such a column will be outputted if the most frequent is gap)

-include      : Toggle Including Consensus in output On/Off

-rmasker_opts : Override Default RepeatMasker Options

-cmatch_opts  : Override Default Cross_Match Options

-sim4_opts    : Override Default Sim4 Options

-defaults     : See defaults

-v	      : Version

-V	      : Verbose

-h            : Help


The default options for MaM are:


-exonfile : "default"

-update	  : on

-column   : single

-merge    : max

-keep     : on

-program  : repeatmasker

-slider   : off

-pc       : P

-sw       : 5

-ww       : 5

-clustal  : on

-nexus    : off

-mega     : off

-fasta    : off

-phylip   : off

-HTML     : off

-identity : off

-consensus: off

-cgaps	  : on

-include  : off



Input/Output:



MaM supports Clustal, MEGA, FASTA, NEXUS, and PHYLIP input formats, and

Clustal, MEGA, NEXUS, FASTA, PHYLIP, HTML and consensus identity dot representation

output formats. 



For more explanation for the program, options, and some examples, please

visit: http://mam-bio.sourceforge.net


