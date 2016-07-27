/*==============================================================================
select.c : selection list
(C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "select.h"

/*____________________________________________________________________________*/
/** parse selection list */
/* new format: chain:residue-range */
/* A selection list of residues on the command line of the type
	"A:4-9 A:13 A:25-34" is converted into an array of 'start end' selections
	corresponding to selection[][0] and selection[][1] entries.
	A range value pairs like '4-9' yields selection[][0]=4 and selection[][1]=9,
	a single value like '13' yields selection[][0]=13 and selection[][1]=13.
*/
/* Starting at release 1.4.1,
	we use 2 selection lists, 'selection0' and 'selection1',
	because mapping might be asymmetrical between the elements of the pairlist,
	for example mapping 'A:4-9' to 'A:25-34' should not include mapping
	of 'A:4' to 'A:5'.
*/

void parse_selects(Arg *arg, char *select, Sel *sel)
{
	unsigned int i, j;
	unsigned int allocated = 64;
    const char delimiter0[2] = " "; /* delimiter types in selection string */
    const char delimiter1[2] = ":";
    const char delimiter2[2] = "-";
	Sel *psel = 0; /* pointer to selection list */
	char *token = 0; /* pointer to token */
	char *psave0 = 0; /* pointer to token position of delimiter0 */
	char *psave1 = 0;
	long int lcast; /* cast input string to long integer */

	/* for both selection lists */
	for (i = 0; i < 2; ++ i) {
		/* move pointer to current list */
		psel = &(sel[i]);
		/* if selection list exists */
		if (strlen(arg->select[i]) > 0) {
			/* initialise selections */
			psel->chainSelection = safe_malloc(allocated * sizeof(char [1]));
			psel->residueSelection = safe_malloc(allocated * sizeof(int [2]));
			psel->nSelection = 0;

			/* get the location of the first token */
			token = strtok_r(arg->select[i], delimiter0, &psave0);

			/* walk through all tokens delimited by delimiter0,
				here the 'space' separating selection segment(s) */
			do {
				/* walk through all tokens delimited by delimiter1,
					here the 'colon' separating chain identifier and residue number(s) */
				if (strstr(token, delimiter1)) {
					token = strtok(token, delimiter1); /* tokes points to chain identifier */
					assert(strlen(token) == 1 && "Chain ID should be a single character");
					psel->chainSelection[psel->nSelection][0] = token[0]; /* assign chain identifier */

					token = strtok(NULL, delimiter1); /* token points to residue number(s) */
				} else {
				/* in case no chain (and colon) is given, the current 'token'
					is passed unchanged to the next 'do' loop,
					resulting in an empty chain identifier */
					psel->chainSelection[psel->nSelection][0] = ' '; /* assign empty chain identifier */
				}
				do {
					/* split the tokens delimited by delimiter2,
						here the 'dash' separating start and end residues of the selected segment(s) */
					if (strstr(token, delimiter2)) {
						token = strtok(token, delimiter2); /* token points to start residue */
						assert(((lcast = strtol(token, NULL, 10)) != 0) && "Residue selection list should be digits");
						psel->residueSelection[psel->nSelection][0] = (int)lcast; /* assign start residue as integer value */

						token = strtok(NULL, delimiter2); /* token points to end residue */
						assert(((lcast = strtol(token, NULL, 10)) != 0) && "Residue selection list should be digits");
						psel->residueSelection[psel->nSelection][1] = (int)lcast; /* assign end residue as integer value */
					/* single residue tokens have identical start-end values */
					} else {
						/* token points to start residue, set start and end value to this */
						assert(((lcast = strtol(token, NULL, 10)) != 0) && "Check residue selection list");
						psel->residueSelection[psel->nSelection][0] = psel->residueSelection[psel->nSelection][1] = (int)lcast;
					}

					assert(psel->residueSelection[psel->nSelection][1] >= psel->residueSelection[psel->nSelection][0] && \
							"Residue selection pair should be in increasing order");

					++ psel->nSelection;

					/* allocate more memory if needed */
					if (psel->nSelection >= allocated) {
						allocated += 64;
						psel->residueSelection = safe_realloc(psel->residueSelection, allocated * sizeof(int [2]));
					}

				} while ((token = strtok_r(NULL, delimiter1, &psave1)));
			} while ((token = strtok_r(NULL, delimiter0, &psave0)));

			if (! arg->silent) {
				fprintf(stdout, "\t%d residue selection(s) in list %d\n",
							psel->nSelection, i+1);
				for (j = 0; j < psel->nSelection; ++ j) {
					fprintf(stdout, "\t\t%d/%d: '%c:%d-%d'\n",
						j+1, psel->nSelection, 
						psel->chainSelection[j][0],
						psel->residueSelection[j][0], psel->residueSelection[j][1]);
				}
			}
		} else {
			psel->nSelection = 0;
		}
	}
}

