/*==============================================================================
getseq.c : read FASTA sequence
Copyright (C) 2004 Jens Kleinjung and John Romein
Read the COPYING file for license information.
==============================================================================*/

#include "getseq.h"

/*____________________________________________________________________________*/
/* Reads a line starting with '>' and followed by the sequence name
	Leading white space is skipped. Reading proceeds until the end
	of the line is reached. The name is read into seq->name. */
__inline__ static int read_sequence_name(FILE *file, Seq *seq)
{
    int ch;
	unsigned int name_length = 0;
	unsigned int allocated = 64;

    while ((ch = getc(file)) != EOF && isspace(ch))
		;

    if (ch != '>')
		return 0;

    seq->name = safe_malloc(allocated);

	/*____________________________________________________________________________*/
    do {
		seq->name[name_length ++] = ch;

		if (name_length == allocated) {
			seq->name = safe_realloc(seq->name, allocated += 64);
		}
    } while ((ch = getc(file)) != EOF && ch != '\n' && isprint(ch));

    seq->name[name_length] = '\0';
    return 1;
}

/*____________________________________________________________________________*/
/* Reads the residues in a sequence, up to (but not including) the
	next sequence header (starting with '>'), or up to end of file.
	Residues may span multiple lines.  White space and gaps are skipped.
	Nonalpha characters are rejected, resulting in an error message.
	Alpha characters are converted to upper case.  The string is read
	into seq->residues and zero-terminated; the length is stored
	into seq->length. */
__inline__ static void read_sequence_residues(FILE *file, Seq *seq)
{
    int ch;
	unsigned int allocated = 64;

    seq->res = safe_malloc(allocated * sizeof(char));
    seq->dgap = safe_malloc(allocated * sizeof(int));
    seq->cgap = safe_malloc(allocated * sizeof(int));
    seq->length = 0;
	seq->nGap = 0;

    while ((ch = getc(file)) != EOF && ch != '>') {
		/* assign residues */
		if (isalpha(ch)) {
			seq->dgap[seq->length] = 0;
			seq->res[seq->length ++] = toupper(ch);
		/* count and assign gaps */
		} else if (ch != '.' || ch != '-') { 
			seq->dgap[seq->length] = 1;
			seq->cgap[seq->length] = ++ seq->nGap;
			seq->res[seq->length ++] = ch;
		/* illegal characters */
		} else if (!isspace(ch)) {
			fprintf(stderr, "Exiting: Illegal character '%c' in seq\n", ch);
			exit(1);
		/* integrity check */
		} else {
			assert(1 != 1 && "this should not happen");
		}
		/* assign more sequence memory if needed */
		if (seq->length == allocated) {
			allocated += 64;
			seq->res = safe_realloc(seq->res, allocated * sizeof(char));
			seq->dgap = safe_realloc(seq->dgap, allocated * sizeof(int));
			seq->cgap = safe_realloc(seq->cgap, allocated * sizeof(int));
		}
	}

	/*____________________________________________________________________________*/
    if (ch == '>')
		ungetc('>', file);

    assert(seq->length > 0);

    seq->res[seq->length] = '\0';
}

/*____________________________________________________________________________*/
/* read FASTA sequence */
void read_sequence(Arg *arg, Seq *seq)
{
	arg->gappedseqInFile = safe_open(arg->gappedseqInFileName, "r");

    if (read_sequence_name(arg->gappedseqInFile, seq))
		read_sequence_residues(arg->gappedseqInFile, seq);

	fclose(arg->gappedseqInFile);

    /*___________________________________________________________________________*/
    if (! arg->silent) fprintf(stdout, "\tsequence file: %s\n"
										"\tsequence file content:\n"
										"\t\tsequence name = %s\n"
										"\t\tsequence length = %d\n",
						arg->gappedseqInFileName, seq->name, seq->length);	
}

