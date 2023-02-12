#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

typedef struct {
  char *id_start;
  int id_length;

  char *annot_start;
  int annot_length;

  char *def_start;
  int def_length;
} header_parts;

static size_t count_fasta(const char * data, size_t size) {
  size_t i;
  size_t  j;

  for(i=0,j = 0; j < size; j++, data++)
        if (*data == '>' && (data[-1] == '\n' || j == 0) ) i++;

  return i;
}

static header_parts locate_annotation(const char *header, size_t size) {
  /* Look for the beginning of the annotation or the end of the header */
  char *index;
  char *index_e;
  size_t  i,i_e;
  int  length;
  int in_json, in_quote;

  header_parts parts;

  parts.id_length = -1;
  parts.annot_length = -1;
  parts.def_length = -1;

  i = 0;
  index = (char*)header;

  // Looks at the first '>' (fasta) or '@' (fastq) indicating
  // the start of the header.

  while (!(*index == '>' || *index == '@')
           && i < size)
           i++,index++;

  if (i >= size) return parts;

  index++;
  length=0;

  parts.id_start = index;

  /* Look for the next white char to delimite the end of the Sequence id */

  while(*index != ' ' && *index!= '\r' && *index != '\n' && i < size)
    index++,i++,length++;

  if (i >= size) return parts;

  parts.id_length = length;

  while(*index == ' ' && *index!= '\r' && *index != '\n' && i < size)
    index++,i++;

  index_e = index;
  i_e = i;

    /* Look for the beginning of the annotation or the end of the header */

  while (*index!= '{' && *index!= '\r' && *index != '\n' && i < size)
    index++,i++;

  if (*index== '{') {
    length  = 0;
    parts.annot_start = index;
    in_json = 1;
    in_quote = 0;


    while (in_json > 0 && *index!= '\r' && *index != '\n' && i < size) {
      length++;
      index++;
      i++;
      if (*index=='"') {
        if (in_quote == 1) in_quote = 0; else in_quote=1;
      }

      if (in_quote == 0) {
        if (*index=='{') in_json++;
        if (*index=='}') in_json--;
      }
    }

    // We have correctly reach the end of the JSON
    if (in_quote == 0 && in_json == 0 && *index=='}') {
      index++;
      i++;
      length++;


      parts.annot_length = length;
    }

  }

  if  (parts.annot_length == -1) {
    index = index_e;
    i = i_e;
  }

  while(*index == ' ' && i < size)
    index++,i++;

  if (i >= size) return parts;

  /* Look for the next end of line to delimite Sequence header */
  parts.def_start = index;
  length=0;

  while(*index != 0 && *index!= '\r' && *index != '\n' && i < size)
    index++,i++,length++;

  if (i >= size) length--,index--;

  while((*index== 0 || *index== ' ' || *index== '\r' || *index == '\n') && index >= parts.def_start) {
    index--,length--;
  }

  length++;

  fprintf(stdout,"length : %d\n",length);

  if (length > 0) parts.def_length = length;

  return parts;
}

SEXP R_parse_fasta(SEXP raw) {
  char* pointer;
  char* seqbuffer;
  char* iseqbuffer;
  size_t seqbuffer_size;
  size_t size;
  size_t ssize;
  size_t nfasta;

  SEXP   seqids;
  SEXP   annotation;
  SEXP   headers;
  SEXP   sequences;
  SEXP   fastas;

  size_t i;
  size_t j;
  int   in_json,in_quote;

  char* index;
  char* index_json;
  char* index_e;
  char  save;

  seqbuffer = NULL;
  seqbuffer_size = 0;

  if (! IS_RAW(raw))
      error("the raw parameter must be of raw type");

	pointer = (char *) RAW(raw);
	size = LENGTH(raw);
	nfasta = count_fasta((const char*)pointer,size);

	seqids = PROTECT(NEW_STRING(nfasta));
	annotation = PROTECT(NEW_STRING(nfasta));
	headers = PROTECT(NEW_STRING(nfasta));
	sequences = PROTECT(NEW_STRING(nfasta));

	fastas = PROTECT(NEW_LIST(4));

	i = 0;
	j = 0;
	index = pointer;


	while (i < size) {

	  while (!(*index == '>'
              && (index[-1] == '\n' || i == 0))
            && i < size)
	    i++,index++;

	  if (i >= size) break;

    /* pointer pointes on the > char */
    index++;
    i++;
    index_e = index;

    /* Look for the next white char to delimite Sequence id */

    while(*index_e != ' ' && *index_e!= '\r' && *index_e != '\n' && i < size)
      index_e++,i++;

	  if (i >= size) break;

    save = *index_e;
    *index_e = 0;
    SET_STRING_ELT(seqids, j, mkChar(index));
    *index_e = save;

    while(*index_e == ' ' && i < size)
      index_e++,i++;

	  if (i >= size) break;

    index = index_e;

    /* Look for the beginning of the annotation or the end of the header */

    while (*index_e!= '{' && *index_e!= '\r' && *index_e != '\n' && i < size)
      index_e++,i++;

    index_json = index_e;

    if (*index_e== '{') in_json = 1; else in_json = 0;
    in_quote = 0;

    while (in_json > 0 && *index_e!= '\r' && *index_e != '\n' && i < size) {
      index_e++;
      i++;
      if (*index_e=='"') {
        if (in_quote == 1) in_quote = 0; else in_quote=1;
      }

      if (in_quote == 0) {
        if (*index_e=='{') in_json++;
        if (*index_e=='}') in_json--;
      }
    }

    if (in_quote == 0 && in_json == 0 && *index_e=='}') {
      index_e++;
      i++;
      save = *index_e;
      *index_e = 0;
        SET_STRING_ELT(annotation, j, mkChar(index_json));
      *index_e = save;
    } else
      SET_STRING_ELT(annotation, j, NA_STRING);

    while(*index_e == ' ' && i < size)
      index_e++,i++;

    if (i >= size) break;
    index = index_e;

    /* Look for the next end of line to delimite Sequence header */

    while(*index_e!= '\r' && *index_e != '\n' && i < size)
      index_e++,i++;

	  if (i >= size) break;

    save = *index_e;
    *index_e = 0;
    SET_STRING_ELT(headers, j, mkChar(index));
    *index_e = save;

    while ((*index_e== '\r' || *index_e == '\n') && i < size)
      index_e++,i++;

	  if (i >= size) break;

    index = index_e;

    while (!(*index_e== '>' && index_e[-1] == '\n') && i < size)
      index_e++,i++;

    ssize = index_e - index;

    if (ssize > seqbuffer_size) {
      seqbuffer_size = ssize;
      seqbuffer = R_alloc(seqbuffer_size,sizeof(char));
      if (seqbuffer == NULL)
        error("cannot assign requested memory");
    }

    iseqbuffer = seqbuffer;

    for (; index < index_e; index++) {
      if ((*index >= 'A' && *index <= 'Z') ||
          (*index >= 'a' && *index <= 'z')) {
        *iseqbuffer = *index;
        iseqbuffer++;
      }
    }
    *iseqbuffer=0;


    SET_STRING_ELT(sequences, j, mkChar(seqbuffer));

    j++;
  }

  SET_VECTOR_ELT(fastas,0,seqids);
	SET_VECTOR_ELT(fastas,1,annotation);
	SET_VECTOR_ELT(fastas,2,headers);
  SET_VECTOR_ELT(fastas,3,sequences);

	UNPROTECT(5);

	return fastas;
}

SEXP R_parse_obiJSONheader(SEXP title_lines) {
  SEXP   seqids;
  SEXP   annotations;
  SEXP   definitions;
  SEXP   headers;
  size_t nfasta;
  size_t i;
  header_parts parts;
  char save;

  headers = PROTECT(NEW_LIST(3));

  nfasta = LENGTH(title_lines);

  seqids = PROTECT(NEW_STRING(nfasta));
  annotations = PROTECT(NEW_STRING(nfasta));
  definitions = PROTECT(NEW_STRING(nfasta));

  for (i = 0 ; i < nfasta; i++) {
    SEXP h = STRING_ELT(title_lines,i);
    parts = locate_annotation(CHAR(h),LENGTH(h));

    if (parts.id_length > 0) {
      save = parts.id_start[parts.id_length];
      parts.id_start[parts.id_length] = 0;
      SET_STRING_ELT(seqids, i, mkChar(parts.id_start));
      parts.id_start[parts.id_length] = save;
    } else {
      SET_STRING_ELT(seqids, i, NA_STRING);
    }

    if (parts.annot_length > 0) {
      save = parts.annot_start[parts.annot_length];
      parts.annot_start[parts.annot_length] = 0;
      SET_STRING_ELT(annotations, i, mkChar(parts.annot_start));
      parts.annot_start[parts.annot_length] = save;
    } else {
      SET_STRING_ELT(annotations, i, NA_STRING);
    }

    if (parts.def_length > 0) {
      save = parts.def_start[parts.def_length];
      parts.def_start[parts.def_length] = 0;
      SET_STRING_ELT(definitions, i, mkChar(parts.def_start));
      parts.def_start[parts.def_length] = save;
    } else {
      SET_STRING_ELT(definitions, i, NA_STRING);
    }

  }

  SET_VECTOR_ELT(headers,0,seqids);
  SET_VECTOR_ELT(headers,1,annotations);
  SET_VECTOR_ELT(headers,2,definitions);

  UNPROTECT(4);

  return headers;
}
