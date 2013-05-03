#ifndef _CONVERT_H_
#define _CONVERT_H_

#include <stdio.h>
#include <string.h>

#define BUFF_SIZE	(4096)
#define VALUE_SIZE	(16)
#define KEY_SIZE	(24)

#define ERROR_WRONG_HEADER (1)

/** \brief convert hhr file to tsv file
 * Converts hhr (hhblits-result) file to tsv (tab separated values) file.
 * Following columns are included:
 * 1, query:		FASTA id to query sequence
 * 2, hit:			FASTA id to hit sequence
 * 3, p_value:		Probability
 * 4, e_value:		E-value
 * 5, score:		Score
 * 6, algn_cols:	Aligned columns
 * 7, identity:		Identity
 * 8, similarity:	Similarity
 * 9, sum_p:		Summed probability
 *
 * @param fhhr file pointer to hhr file (input)
 * @param ftsv file pointer to tsv file (output)
 */
int hhr2tsv(FILE *fhhr, FILE *ftsv);

/** \brief Write header of tsv file
 * @param ftsv File pointer ot tsv file (output)
 */
FILE * write_tsv_header(FILE *ftsv);

/** \brief Skip file curosr to next line
 *  @param file File pointer
 */
FILE * next_line(FILE *file);

/** \brief Skip white spaces (spaces, tabs and line breaks)
 * @param file File pointer to input file
 */
FILE * skip_white_spaces(FILE *file);

/** read string until space, tab, line break and EOF. 0 is appended at end of
 * the string.
 * @param file pointer to input file
 * @return character after the read word 
 */
int read_word(FILE *file, char* str, size_t len);

/** check whether a character indicates end of line
 * @param chr character to be checked
 */
int is_line_ended(int chr);

struct struct_entry {
	char pval[VALUE_SIZE];
	char eval[VALUE_SIZE];
	char score[VALUE_SIZE];
	char align_col[VALUE_SIZE];
	char identity[VALUE_SIZE];
	char similarity[VALUE_SIZE];
	char sumpval[VALUE_SIZE];
};
typedef struct struct_entry entry_t;

/** read a hhr alignment performance entry
 * @param fhhr file pointer to hhr file
 * @param entry pointer to an entry data structure
 */
void read_entry(FILE *fhhr, entry_t *entry);

/** print out an alignment performance entry
 * @param ftsv file pointer to tsv file
 * @param entry pointer to an entry
 */
void print_entry(FILE *ftsv, const entry_t *entry);

/** get key and value from a section in alignment performance entry
 * @param word section that contains key and value
 * @param key char pointer for key string
 * @param val char pointer for value string
 */
void get_key_and_value(const char* word, char* key, char* val);

#endif
