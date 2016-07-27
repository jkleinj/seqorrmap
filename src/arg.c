/*==============================================================================
arg.c : parse command line arguments
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "arg.h"

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
	fprintf(stdout, "\n%s %s\n", PACKAGE, VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header()
{
    fprintf(stdout, "\nseqorrmap : map sequence correlations onto structure\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stdout, "\nCopyright (C) 2012-2013 Jens Kleinjung\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the COPYING file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
	fprintf(stdout, "\nno citation\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg, Argpdb *argpdb)
{
    arg->pdbInFileName = "";
    arg->corrInFileName = "";
    arg->gappedseqInFileName = "";
	arg->rLoCutoff = 0.1; /* lower radius cutoff in A */
	arg->rUpCutoff = FLT_MAX; /* upper radius cutoff in A */
	arg->cLoCutoff = -FLT_MAX; /* lower correlation cutoff */
	arg->cUpCutoff = FLT_MAX; /* upper correlation cutoff */
	arg->sLoCutoff = 0; /* lower sequence distance cutoff */
	arg->sUpCutoff = INT_MAX; /* upper sequence distance cutoff */
	arg->firstSeqResidue = 0; /* correlation map starts at this sequence residue */
	arg->firstStrResidue = 0; /* correlation map starts at this structure residue */
    arg->pdbOutFileName = "seqorrmap.pdb";
    arg->corrOutFileName = "seqorrmap_corr.dat";
    arg->tclOutFileName = "seqorrmap.vmd";
    arg->corrmatOutFileName = "seqorrmap.mat";
    arg->corrcompmatOutFileName = "seqorrmap_comp.mat";
    arg->logoOutFileName = "seqorrmap.out";
    arg->logeOutFileName = "seqorrmap.err";
	arg->select = (char **)safe_malloc(2 * sizeof(char *)); /* 2 selection lists */
	arg->select[0] = ""; /* selection list */
	arg->select[1] = "";
    arg->strict_aa = 0; /* from here: run mode flags */
    arg->strict_nr = 0;
    arg->silent = 0;
    arg->multiModel = 0;
	arg->gapped = 0;
	arg->systemChain = 0; /* systematic intra-program chain numbering */
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg, Argpdb *argpdb)
{
	if (strlen(arg->pdbInFileName) == 0)
		Error("Invalid PDB file name");
	if (strlen(arg->corrInFileName) == 0)
		Error("Invalid correlation file name");
	assert(arg->rLoCutoff >= 0. && "cutoff radius too small");
	assert(arg->rUpCutoff >= 0. && "cutoff radius too small");
	assert(arg->rUpCutoff > 0  && "sequence distance cutoff too small");
	assert(arg->rLoCutoff < arg->rUpCutoff && "lower RMSD cutoff should be smaller than upper cutoff");
	assert(arg->cLoCutoff < arg->cUpCutoff && "lower correlation cutoff should be smaller than upper cutoff");
	assert(arg->sLoCutoff < arg->sUpCutoff && "lower sequence distance cutoff should be smaller than upper cutoff");
}

/*____________________________________________________________________________*/
/** print command line arguments */
static void print_args(Arg *arg, Argpdb *argpdb)
{
    time_t now;
    time(&now);

	if (! arg->silent) {
		fprintf(stdout, "\tdate: %s", ctime(&now));

		fprintf(stdout, \
						"\nRun parameters\n"
						"\tpdbIn: %s\n"
						"\tcorr: %s\n"
						"\tgappedSeq: %s\n"
						"\trLoCutoff: %8.2e\n"
						"\trUpCutoff: %8.2e\n"
						"\tcLoCutoff: %8.2e\n"
						"\tcUpCutoff: %8.2e\n"
						"\tsLoCutoff: %d\n"
						"\tsUpCutoff: %d\n"
						"\tfirstSeqResidue: %d\n"
						"\tfirstStrResidue: %d\n"
						"\tselect[0]: %s\n"
						"\tselect[1]: %s\n"
						"\tpdbOut: %s\n"
						"\tcorrOut: %s\n"
						"\tcorrmatOut: %s\n"
						"\tcorrcompmatOut: %s\n"
						"\ttclOut: %s\n"
						"\tstrict_aa: %d\n"
						"\tstrict_nr: %d\n"
						"\tsystemChain: %d\n",
			arg->pdbInFileName, arg->corrInFileName, arg->gappedseqInFileName,
			arg->rLoCutoff, arg->rUpCutoff,
			arg->cLoCutoff, arg->cUpCutoff,
			arg->sLoCutoff, arg->sUpCutoff,
			arg->firstSeqResidue, arg->firstStrResidue,
			arg->select[0], arg->select[1],
			arg->pdbOutFileName, arg->corrOutFileName,
			arg->tclOutFileName,
			arg->corrmatOutFileName, arg->corrcompmatOutFileName,
			arg->strict_aa, arg->strict_nr,
			arg->systemChain);

		fflush(stdout);
	}
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg, Argpdb *argpdb)
{
	int c;
	const char usage[] = "\nseqorrmap [--pdb ...] [--corr ...] [OPTIONS ...]\n\
	   --pdbIn <PDB input>\t\t\t\t(mode: mandatory, type: char  , default: void)\n\
	   --corr <seqorr correlation input>\t\t(mode: mandatory, type: char  , default: void)\n\
	   --gappedSeq <gapped target sequence>\t\t(mode: optional , type: char  , default: void)\n\
	   --rLoCutoff <lower RMSD cutoff>\t\t(mode: optional , type: float , default: 0.1 Angstrom)\n\
	   --rUpCutoff <upper RMSD cutoff>\t\t(mode: optional , type: float , default: FLT_MAX Angstrom)\n\
	   --cLoCutoff <lower correlation cutoff>\t(mode: optional , type: float , default: -FLT_MAX)\n\
	   --cUpCutoff <upper correlation cutoff>\t(mode: optional , type: float , default: FLT_MAX)\n\
	   --sLoCutoff <lower sequence distance cutoff>\t(mode: optional , type: int   , default: 0)\n\
	   --sUpCutoff <upper sequence distance cutoff>\t(mode: optional , type: int   , default: INT_MAX)\n\
	   --firstSeqResidue <seq residue number>\t(mode: optional , type: int   , default: 0)\n\
	   --firstStrResidue <str residue number>\t(mode: optional , type: int   , default: 0)\n\
	   --strict_aa\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --strict_nr\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --silent\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --pdbOut <PDB output>\t\t\t(mode: optional , type: char  , default: seqorrmap.pdb)\n\
	   --corrOut <corr output>\t\t\t(mode: optional , type: char  , default: seqorrmap_corr.dat)\n\
	   --tclOut <TCL output>\t\t\t(mode: optional , type: char  , default: seqorrmap.vmd)\n\
	   --corrmatOut <matrix output>\t\t\t(mode: optional , type: char  , default: seqorrmap.mat)\n\
	   --corrcompmatOut <matrix output>\t\t(mode: optional , type: char  , default: seqorrmap_comp.mat)\n\
	   --select1 <first-residue selection list>\t(mode: optional , type: int   , default: OFF (INT_MIN))\n\
	   --select2 <second-residue selection list>\t(mode: optional , type: int   , default: OFF (INT_MIN))\n\
	   --systemChain\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --cite\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --version\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --help\n";

    set_defaults(arg, argpdb);

    if (argc < 3) {
		print_header();
        fprintf(stderr, "%s", usage);
		print_license();
        exit(1);
    }

    /** long option definition */
    static struct option long_options[] = {
        {"pdbIn", required_argument, 0, 1},
        {"corr", required_argument, 0, 2},
        {"gappedSeq", required_argument, 0, 3},
        {"rLoCutoff", required_argument, 0, 5},
        {"rUpCutoff", required_argument, 0, 6},
        {"cLoCutoff", required_argument, 0, 7},
        {"cUpCutoff", required_argument, 0, 8},
        {"sLoCutoff", required_argument, 0, 9},
        {"sUpCutoff", required_argument, 0, 10},
        {"firstSeqResidue", required_argument, 0, 11},
        {"firstStrResidue", required_argument, 0, 12},
        {"pdbOut", required_argument, 0, 13},
        {"corrOut", required_argument, 0, 14},
        {"tclOut", required_argument, 0, 15},
        {"corrmatOut", required_argument, 0, 16},
        {"corrcompmatOut", required_argument, 0, 17},
        {"select1", required_argument, 0, 18},
        {"select2", required_argument, 0, 19},
        {"strict_aa", no_argument, 0, 20},
        {"strict_nr", no_argument, 0, 21},
        {"systemChain", no_argument, 0, 22},
        {"silent", no_argument, 0, 30},
        {"cite", no_argument, 0, 31},
        {"version", no_argument, 0, 32},
        {"help", no_argument, 0, 33},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19: 20 21 22 30 31 32 33", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->pdbInFileName = optarg;
                break;
            case 2:
                arg->corrInFileName = optarg;
                break;
            case 3:
                arg->gappedseqInFileName = optarg;
                break;
            case 5:
                arg->rLoCutoff = atof(optarg);
                break;
            case 6:
                arg->rUpCutoff = atof(optarg);
                break;
            case 7:
                arg->cLoCutoff = atof(optarg);
                break;
            case 8:
                arg->cUpCutoff = atof(optarg);
                break;
            case 9:
                arg->sLoCutoff = atoi(optarg);
                break;
            case 10:
                arg->sUpCutoff = atoi(optarg);
                break;
            case 11:
                arg->firstSeqResidue = atoi(optarg);
                break;
            case 12:
                arg->firstStrResidue = atoi(optarg);
                break;
            case 13:
                arg->pdbOutFileName = optarg;
                break;
            case 14:
                arg->corrOutFileName = optarg;
                break;
            case 15:
                arg->tclOutFileName = optarg;
                break;
            case 16:
                arg->corrmatOutFileName = optarg;
                break;
            case 17:
                arg->corrcompmatOutFileName = optarg;
                break;
            case 18:
                arg->select[0] = optarg;
                break;
            case 19:
                arg->select[1] = optarg;
                break;
            case 20:
                arg->strict_aa = 1;
                break;
            case 21:
                arg->strict_nr = 1;
                break;
            case 22:
                arg->systemChain = 1;
                break;
            case 30:
                arg->silent = 1;
                break;
            case 31:
                print_citation();
                exit(0);
            case 32:
				print_version();
				print_license();
                exit(0);
            case 33:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(1);
        }
    }

	check_input(arg, argpdb);
	if (! arg->silent) {
		print_header();
		print_args(arg, argpdb);
	}

	if (strlen(arg->gappedseqInFileName) > 0)
		++ arg->gapped;

    return 0;
}

