#include "convert.h"

int hhr2tsv(FILE *fhhr, FILE *ftsv)
{
	char buff[BUFF_SIZE] = {0},
	     buff2[BUFF_SIZE] = {0};
	char qid[BUFF_SIZE] = {0};
	const char *query = qid;
	entry_t entry;
	int len = 0, chr = 0;

	/* read hhr header */
	fread(buff, 1, 5, fhhr);
	if (strncmp(buff, "Query", 5) != 0) {
		return ERROR_WRONG_HEADER; 
	} else { 
		read_word(skip_white_spaces(fhhr),
				  qid,
				  BUFF_SIZE);
	}
	
	/* write tsv header */
	write_tsv_header(ftsv);
	
	/* search for hit entry */
	while ( !feof(fhhr)	) 
	{
		/* reset buff */
		buff[0] = 0;
		chr = read_word(fhhr, buff, BUFF_SIZE);
		/* got hit sequence id, print it */
		if (buff[0] == '>') {
			fprintf(ftsv, "%s\t", query);
			fprintf(ftsv, "%s\t", (const char*)(buff+1) );
		} 

		/* current line not ended, skip it */
		if (chr != '\n' && 
		    chr != '\0' && 
		    chr != EOF) 
		{
			next_line(fhhr);
		}

		/* read and print entry for tsv */
		if (buff[0] == '>')
		{
			read_entry(fhhr, &entry);
			print_entry(ftsv, &entry);
		}
	}

	return 0;
}

FILE* write_tsv_header(FILE *ftsv)
{
	fprintf(ftsv, 
	"query\thit\tp_value\te_value\tscore\talgn_cols\tidentity\tsimilarity\tsum_p\n");
	if (ferror(ftsv)) 
		perror("Error occurred by writing header");
	return ftsv;
}

FILE* next_line(FILE *file)
{
	int ch = 0;
	do {
		ch = fgetc(file);
	} while (ch != EOF && ch != '\n');
	return file;
}

FILE* skip_white_spaces(FILE *file)
{
	int ch = 0;
	do {
		ch = fgetc(file);
	} while(ch == ' ' || ch == '\t' || ch == '\n');
	fseek(file, -1, SEEK_CUR);
	return file;
}

int read_word(FILE *file, char* str, size_t len)
{
	int count = 0;
	int ch = 0;
	do {
		ch = fgetc(file);
		if (ch != EOF && 
			ch != ' ' && 
			ch != '\t' && 
			ch != '\n' && 
			ch != '\0' && 
			count < len) 
		{
			str[count++] = ch;
		}
	} while 
	(ch != EOF && 
	 ch != ' ' && 
	 ch != '\t' && 
	 ch != '\n' && 
	 ch != '\0' &&
	 count < len);
	str[count] = '\0';
	return ch;
}

void print_entry(FILE *ftsv, const entry_t *entry)
{
	fprintf(ftsv, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
			entry->pval,
			entry->eval,
			entry->score,
			entry->align_col,
			entry->identity,
			entry->similarity,
			entry->sumpval);
}

void read_entry(FILE *fhhr, entry_t *entry) 
{
	char buff[BUFF_SIZE] = {0},
	key[KEY_SIZE] = {0},
	val[VALUE_SIZE] = {0};
	int ter = 0, sections = 7;
	do {
		ter = read_word(fhhr, buff, BUFF_SIZE);
		if (buff[0]) 
		{
			get_key_and_value(buff, key, val);
			if ( strcmp(key, "Probab") == 0 ) {
				strcpy(entry->pval, val);
				sections--;
			} else if ( strcmp(key, "E-value") == 0 ) {
				strcpy(entry->eval, val);
				sections--;
			} else if ( strcmp(key, "Score") == 0 ) {
				strcpy(entry->score, val);
				sections--;
			} else if ( strcmp(key, "Aligned_cols") == 0 ) {
				strcpy(entry->align_col, val);
				sections--;
			} else if ( strcmp(key, "Identities") == 0 ) {
				strcpy(entry->identity, val);
				sections--;
			} else if ( strcmp(key, "Similarity") == 0 ) {
				strcpy(entry->similarity, val);
				sections--;
			} else if (strcmp(key, "Sum_probs") == 0 ) {
				strcpy(entry->sumpval, val);
				sections--;
			}
		}
		skip_white_spaces(fhhr);
	} while (ter != EOF && ter != '\n' && sections > 0);
}

void get_key_and_value(const char* word, char* key, char* val)
{
	int idx = 0, key_ended = 0, out_idx = 0;
	for (; word[idx] != '\0'; idx++) 
	{
		if (!key_ended && word[idx] != '=')
		{
			key[out_idx++] = word[idx];
		} 
		else if (word[idx] == '=') 
		{
			key_ended = 1;
			key[out_idx] = 0;
			out_idx = 0;
		} 
		else if (key_ended && word[idx] != '=')
		{
			val[out_idx++] = word[idx];
		}
	}
	val[out_idx] = 0;
}
