/* The MIT License
   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>
   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
https://github.com/attractivechaos/klib
*/

#ifndef _biomcmc_AC_BMC2_KSEQ_H
#define _biomcmc_AC_BMC2_KSEQ_H

#include <zlib.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define BMC2_KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define BMC2_KS_SEP_TAB   1 // isspace() && !' '
#define BMC2_KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define BMC2_KS_SEP_MAX   2

#ifndef bmc2_klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define bmc2_klib_unused __attribute__ ((__unused__))
#else
#define bmc2_klib_unused
#endif
#endif /* bmc2_klib_unused */

#define __BMC2_KS_TYPE(type_t) \
	typedef struct __bmc2_kstream_t { \
		int begin, end; \
		int is_eof:2, bufsize:30; \
		type_t f; \
		unsigned char *buf; \
	} bmc2_kstream_t;

#define bmc2_ks_eof(bmc2_ks) ((bmc2_ks)->is_eof && (bmc2_ks)->begin >= (bmc2_ks)->end)
#define bmc2_ks_rewind(bmc2_ks) ((bmc2_ks)->is_eof = (bmc2_ks)->begin = (bmc2_ks)->end = 0)

#define __BMC2_KS_BASIC(SCOPE, type_t, __bufsize) \
	SCOPE bmc2_kstream_t *bmc2_ks_init(type_t f) \
	{ \
		bmc2_kstream_t *bmc2_ks = (bmc2_kstream_t*)calloc(1, sizeof(bmc2_kstream_t)); \
		bmc2_ks->f = f; bmc2_ks->bufsize = __bufsize; \
		bmc2_ks->buf = (unsigned char*)malloc(__bufsize); \
		return bmc2_ks; \
	} \
	SCOPE void bmc2_ks_destroy(bmc2_kstream_t *bmc2_ks) \
	{ \
		if (!bmc2_ks) return; \
		free(bmc2_ks->buf); \
		free(bmc2_ks); \
	}

#define __BMC2_KS_INLINED(__read) \
	static inline bmc2_klib_unused int bmc2_ks_getc(bmc2_kstream_t *bmc2_ks) \
	{ \
		if (bmc2_ks->is_eof && bmc2_ks->begin >= bmc2_ks->end) return -1; \
		if (bmc2_ks->begin >= bmc2_ks->end) { \
			bmc2_ks->begin = 0; \
			bmc2_ks->end = __read(bmc2_ks->f, bmc2_ks->buf, bmc2_ks->bufsize); \
			if (bmc2_ks->end < bmc2_ks->bufsize) bmc2_ks->is_eof = 1; \
			if (bmc2_ks->end == 0) return -1; \
		} \
		return (int)bmc2_ks->buf[bmc2_ks->begin++]; \
	} \
	static inline int bmc2_ks_getuntil(bmc2_kstream_t *bmc2_ks, int delimiter, bmc2_kstring_t *str, int *dret) \
	{ return bmc2_ks_getuntil2(bmc2_ks, delimiter, str, dret, 0); }

#ifndef BMC2_KSTRING_T
#define BMC2_KSTRING_T bmc2_kstring_t
typedef struct __bmc2_kstring_t {
	unsigned l, m;
	char *s;
} bmc2_kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define __BMC2_KS_GETUNTIL(SCOPE, __read) \
	SCOPE int bmc2_ks_getuntil2(bmc2_kstream_t *bmc2_ks, int delimiter, bmc2_kstring_t *str, int *dret, int append) \
	{ \
		if (dret) *dret = 0; \
		str->l = append? str->l : 0; \
		if (bmc2_ks->begin >= bmc2_ks->end && bmc2_ks->is_eof) return -1; \
		for (;;) { \
			int i; \
			if (bmc2_ks->begin >= bmc2_ks->end) { \
				if (!bmc2_ks->is_eof) { \
					bmc2_ks->begin = 0; \
					bmc2_ks->end = __read(bmc2_ks->f, bmc2_ks->buf, bmc2_ks->bufsize); \
					if (bmc2_ks->end < bmc2_ks->bufsize) bmc2_ks->is_eof = 1; \
					if (bmc2_ks->end == 0) break; \
				} else break; \
			} \
             if (delimiter == BMC2_KS_SEP_LINE)  { for (i = bmc2_ks->begin; i < bmc2_ks->end; ++i) if (bmc2_ks->buf[i] == '\n') break; \
			} else if (delimiter > BMC2_KS_SEP_MAX)    { for (i = bmc2_ks->begin; i < bmc2_ks->end; ++i) if (bmc2_ks->buf[i] == delimiter) break; \
			} else if (delimiter == BMC2_KS_SEP_SPACE) { for (i = bmc2_ks->begin; i < bmc2_ks->end; ++i) if (isspace(bmc2_ks->buf[i])) break; \
			} else if (delimiter == BMC2_KS_SEP_TAB)   { for (i = bmc2_ks->begin; i < bmc2_ks->end; ++i) if (isspace(bmc2_ks->buf[i]) && bmc2_ks->buf[i] != ' ') break; \
			} else i = 0; /* never come to here! */ \
			if (str->m - str->l < (size_t)(i - bmc2_ks->begin + 1)) { \
				str->m = str->l + (i - bmc2_ks->begin) + 1; \
				kroundup32(str->m); \
				str->s = (char*)realloc(str->s, str->m); \
			} \
			memcpy(str->s + str->l, bmc2_ks->buf + bmc2_ks->begin, i - bmc2_ks->begin); \
			str->l = str->l + (i - bmc2_ks->begin); \
			bmc2_ks->begin = i + 1; \
			if (i < bmc2_ks->end) { \
				if (dret) *dret = bmc2_ks->buf[i]; \
				break; \
			} \
		} \
		if (str->s == 0) { \
			str->m = 1; \
			str->s = (char*)calloc(1, 1); \
		} else if (delimiter == BMC2_KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; \
		str->s[str->l] = '\0'; \
		return str->l; \
	}

#define BMC2_KSTREAM_INIT2(SCOPE, type_t, __read, __bufsize) \
	__BMC2_KS_TYPE(type_t) \
	__BMC2_KS_BASIC(SCOPE, type_t, __bufsize) \
	__BMC2_KS_GETUNTIL(SCOPE, __read) \
	__BMC2_KS_INLINED(__read)

#define BMC2_KSTREAM_INIT(type_t, __read, __bufsize) BMC2_KSTREAM_INIT2(static, type_t, __read, __bufsize)

#define BMC2_KSTREAM_DECLARE(type_t, __read) \
	__BMC2_KS_TYPE(type_t) \
	extern int bmc2_ks_getuntil2(bmc2_kstream_t *bmc2_ks, int delimiter, bmc2_kstring_t *str, int *dret, int append); \
	extern bmc2_kstream_t *bmc2_ks_init(type_t f); \
	extern void bmc2_ks_destroy(bmc2_kstream_t *bmc2_ks); \
	__BMC2_KS_INLINED(__read)

/******************
 * FASTA/Q parser *
 ******************/

#define bmc2_kseq_rewind(bmc2_ks) ((bmc2_ks)->last_char = (bmc2_ks)->f->is_eof = (bmc2_ks)->f->begin = (bmc2_ks)->f->end = 0)

#define __BMC2_KSEQ_BASIC(SCOPE, type_t) \
	SCOPE bmc2_kseq_t *bmc2_kseq_init(type_t fd) \
	{ \
		bmc2_kseq_t *s = (bmc2_kseq_t*)calloc(1, sizeof(bmc2_kseq_t)); \
		s->f = bmc2_ks_init(fd); \
		return s; \
	} \
	SCOPE void bmc2_kseq_destroy(bmc2_kseq_t *bmc2_ks) \
	{ \
		if (!bmc2_ks) return; \
		free(bmc2_ks->name.s); free(bmc2_ks->comment.s); free(bmc2_ks->seq.s);	free(bmc2_ks->qual.s); \
		bmc2_ks_destroy(bmc2_ks->f); \
		free(bmc2_ks); \
	}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
#define __BMC2_KSEQ_READ(SCOPE) \
	SCOPE int bmc2_kseq_read(bmc2_kseq_t *seq) \
	{ \
		int c; \
		bmc2_kstream_t *bmc2_ks = seq->f; \
		if (seq->last_char == 0) { /* then jump to the next header line */ \
			while ((c = bmc2_ks_getc(bmc2_ks)) != -1 && c != '>' && c != '@'); \
			if (c == -1) return -1; /* end of file */ \
			seq->last_char = c; \
		} /* else: the first header char has been read in the previous call */ \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
		if (bmc2_ks_getuntil(bmc2_ks, BMC2_KS_SEP_LINE, &seq->name, &c) < 0) return -1; /* normal exit: EOF */ \
		if (c != '\n') bmc2_ks_getuntil(bmc2_ks, BMC2_KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ \
		if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ \
			seq->seq.m = 256; \
			seq->seq.s = (char*)malloc(seq->seq.m); \
		} \
		while ((c = bmc2_ks_getc(bmc2_ks)) != -1 && c != '>' && c != '+' && c != '@') { \
			if (c == '\n') continue; /* skip empty lines */ \
			seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
			bmc2_ks_getuntil2(bmc2_ks, BMC2_KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ \
		} \
		if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */ \
		if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
			seq->seq.m = seq->seq.l + 2; \
			kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
			seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
		} \
		seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ \
		if (c != '+') return seq->seq.l; /* FASTA */ \
		if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
			seq->qual.m = seq->seq.m; \
			seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
		} \
		while ((c = bmc2_ks_getc(bmc2_ks)) != -1 && c != '\n'); /* skip the rest of '+' line */ \
		if (c == -1) return -2; /* error: no quality string */ \
		while (bmc2_ks_getuntil2(bmc2_ks, BMC2_KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l); \
		seq->last_char = 0;	/* we have not come to the next header line */ \
		if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
		return seq->seq.l; \
	}

#define __BMC2_KSEQ_TYPE(type_t) \
	typedef struct { \
		bmc2_kstring_t name, comment, seq, qual; \
		int last_char; \
		bmc2_kstream_t *f; \
	} bmc2_kseq_t;

#define BMC2_KSEQ_INIT2(SCOPE, type_t, __read) \
	BMC2_KSTREAM_INIT2(SCOPE, type_t, __read, 16384) \
	__BMC2_KSEQ_TYPE(type_t) \
	__BMC2_KSEQ_BASIC(SCOPE, type_t) \
	__BMC2_KSEQ_READ(SCOPE)

#define BMC2_KSEQ_INIT(type_t, __read) BMC2_KSEQ_INIT2(static, type_t, __read)

#define BMC2_KSEQ_DECLARE(type_t) \
	__BMC2_KS_TYPE(type_t) \
	__BMC2_KSEQ_TYPE(type_t) \
	extern bmc2_kseq_t *bmc2_kseq_init(type_t fd); \
	void bmc2_kseq_destroy(bmc2_kseq_t *bmc2_ks); \
	int bmc2_kseq_read(bmc2_kseq_t *seq);

#endif
