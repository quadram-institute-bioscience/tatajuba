/*
 * CSCUTILS - A collection of various software routines uses in CSC projects
 * Copyright (C) 2015 Martin Koehler
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, see <http://www.gnu.org/licenses/>.
 *
 * https://github.com/mpimd-csc/flexiblas
 * https://github.com/mpimd-csc/flexiblas/blob/master/libcscutils/src/file/compress_io_xz.c
 */

#if _FILE_OFFSET_BITS == 64
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <lzma.h>

#define TRUE -1
#define FALSE 0
#define XZ_PRESET 8
#define XZ_BUFFERSIZE 4096

#define MAX(A,B) (((A)>(B))?(A):(B));
#define MIN(A,B) (((A)<(B))?(A):(B));

/* Internal structure to represent an uncompressed file. */
typedef struct {
  FILE *fp;
  char *path, mode;
  lzma_stream strm;
  uint8_t inbuf[XZ_BUFFERSIZE];
  uint8_t outbuf[XZ_BUFFERSIZE];
  uint8_t readbuf[XZ_BUFFERSIZE];
  size_t readpos;
  int eof;
  lzma_action action;
} xz_file_t;

xz_file_t * xz_open (const char *path, const char mode);
int xz_close (xz_file_t *f);
size_t xz_read (xz_file_t *f, uint8_t *ibuf, size_t len);
size_t xz_write(xz_file_t *f, char *cbuf, size_t len);


int init_encoder (lzma_stream *strm, uint32_t preset)
{
  lzma_ret ret = lzma_easy_encoder (strm, preset, LZMA_CHECK_CRC32);
  if (ret == LZMA_OK) return TRUE;    // Return successfully if the initialization went fine.
  const char *msg;  /* error */
  switch (ret) {
    case LZMA_MEM_ERROR: msg = "Memory allocation failed"; break;
    case LZMA_OPTIONS_ERROR: msg = "Specified preset is not supported"; break;
    case LZMA_UNSUPPORTED_CHECK: msg = "Specified integrity check is not supported"; break;
    default: msg = "Unknown error, possibly a bug"; break;
  }
  fprintf (stderr, "Error initializing the encoder: %s (error code %u)\n", msg, ret);
  return FALSE;
}

int init_decoder (lzma_stream *strm)
{
  lzma_ret ret = lzma_stream_decoder( strm, UINT64_MAX, 0x00 | LZMA_CONCATENATED );
  if (ret == LZMA_OK) return TRUE;  // Return successfully if the initialization went fine.
  const char *msg;  /* error */
  switch (ret) {
    case LZMA_MEM_ERROR: msg = "Memory allocation failed"; break;
    case LZMA_OPTIONS_ERROR: msg = "Unsupported decompressor flags"; break;
    default: msg = "Unknown error, possibly a bug"; break;
  }
  fprintf(stderr, "Error initializing the decoder: %s (error code %u)\n",msg, ret);
  return FALSE;
}

/* Local function which opens an uncompressed file. */
xz_file_t * xz_open (const char *path, const char mode) {
  xz_file_t *f;
  int err;

  if ((mode != 'w') && (mode != 'r')) {fprintf (stderr, "unrecognised mode %s\n", mode); return NULL; }
  /* Initlize the data structure  */
  f = (xz_file_t *) malloc(sizeof (xz_file_t));
  memset(&(f->strm), 0, sizeof(f->strm));
  f->path = strdup(path);
  f->mode = mode;

  /* Init LZMA  */
  if ( f->mode == 'w' ) {
    /* Open for writing */
    f->fp =fopen(f->path, "w");
    if ( !(f->fp)) {
      err = errno;
      fprintf(stderr, "opening file: %s failed, errno: %03d - %s\n", path, err, strerror(err));
      free(f);
      return NULL;
    }
    if ( !init_encoder(&(f->strm), XZ_PRESET) ) {
      fprintf (stderr, "Can not initialize lzma encoder\n");
      free(f);
      return NULL;
    }
    f->strm.next_in = NULL;
    f->strm.avail_in = 0;
    f->strm.next_out = f->outbuf;
    f->strm.avail_out = sizeof(f->outbuf);
    f->action = LZMA_RUN;
  } else {
    /* Open for READING */
    f->fp =fopen(f->path, "r");
    if ( !(f->fp)) {
      err = errno;
      fprintf (stderr, "opening file: %s failed, errno: %03d - %s\n", path, err, strerror(err));
      free(f);
      return NULL;
    }
    if ( !init_decoder(&(f->strm)) ) {
      fprintf(stderr,"Can not initialize lzma decoder\n");
      free(f);
      return NULL;
    }
    f->strm.next_in = NULL;
    f->strm.avail_in = 0;
    f->strm.next_out = f->outbuf;
    f->strm.avail_out = sizeof(f->outbuf);
    f->readpos = 0;
    f->eof = 0;
    f->action = LZMA_RUN;
  }
  return f;
}

/* Local function to close and free an uncompressed file */
int xz_close (xz_file_t *f) {
  lzma_ret ret;

  if (  f == NULL ) { fprintf(stderr, "Error: data == NULL\n"); return -1; }
  if ( f->mode == 'w') {
    f->action = LZMA_FINISH;
    ret = LZMA_OK;
    while ( ret != LZMA_STREAM_END ) {  /* Finish the encoding   */
      ret = lzma_code(&(f->strm), f->action);
      // If the output buffer is full or if the compression finished successfully, write the data from the output bufffer to the output file.
      if (f->strm.avail_out == 0 || ret == LZMA_STREAM_END) {
        // When lzma_code() has returned LZMA_STREAM_END, the output buffer is likely to be only partially
        // full. Calculate how much new data there is to be written to the output file.
        size_t write_size = sizeof(f->outbuf) - f->strm.avail_out;
        if (fwrite(f->outbuf, 1, write_size, f->fp) != write_size) { fprintf (stderr, "Write error: %s\n", strerror(errno)); return 0; }
        // Reset next_out and avail_out.
        f->strm.next_out = f->outbuf;
        f->strm.avail_out = sizeof(f->outbuf);
      }
      if ( ret != LZMA_STREAM_END && ret != LZMA_OK) { fprintf(stderr, "LZMA Encode error\n"); break; }
    }

    lzma_end(&(f->strm));
  } else {
    lzma_end(&(f->strm));
  }
  printf("close xz\n");
  fflush(f->fp);
  fclose(f->fp);
  f->fp = NULL;
  free(f->path);
  free(f);
  return 0;
}


/* Local function to read a block from an uncompressed file */
size_t xz_read (xz_file_t *f, uint8_t *ibuf, size_t len) {
  // lzma_action action = LZMA_RUN;
  lzma_ret ret;
  size_t pos = 0;
  size_t read_entries = 0;

  if ( f == NULL) return 0;
  if ( ibuf == NULL) return 0;
  if ( len == 0) return 0;

  if ( f->readpos == 0 && f->eof ) return 0;

  if ( f->readpos > 0 ) {
    size_t remain_to_copy = XZ_BUFFERSIZE - f->readpos;
    if ( remain_to_copy > len ) {
      memcpy(ibuf, &(f->readbuf[f->readpos]), sizeof(uint8_t) * len);
      f->readpos += len;
      return len;
    } else if ( remain_to_copy == len ) {
      memcpy(ibuf, &(f->readbuf[f->readpos]), sizeof(uint8_t) * len);
      f->readpos = 0;
      return len;
    } else {
      memcpy(ibuf, &(f->readbuf[f->readpos]), sizeof(uint8_t) * remain_to_copy);
      len = len - remain_to_copy;
      pos = remain_to_copy;
      read_entries = remain_to_copy;
      f->readpos = 0;
    }
  }

  f->strm.next_out = f->outbuf;
  f->strm.avail_out = sizeof(f->outbuf);

  while (! f->eof) {
    if (f->strm.avail_in == 0 && !feof(f->fp)) {
      f->strm.next_in = f->inbuf;
      f->strm.avail_in = fread(f->inbuf, 1, sizeof(f->inbuf),f->fp);

      if (ferror(f->fp)) { fprintf(stderr, "Read error: %s\n",strerror(errno)); return 0; }

      // Once the end of the input file has been reached, we need to tell lzma_code() that no more input
      // will be coming. As said before, this isn't required if the LZMA_CONATENATED flag isn't used when initializing the decoder.
      if (feof(f->fp)) f->action = LZMA_FINISH;
    }

    ret = lzma_code(&(f->strm), f->action);

    if (f->strm.avail_out == 0 || ret == LZMA_STREAM_END) {
      size_t write_size = sizeof(f->outbuf) - f->strm.avail_out;

      memcpy(f->readbuf, f->outbuf, sizeof(uint8_t) * write_size);

      if ( len > 0 ) {
        size_t towrite = MIN(len, write_size);
        memcpy(&ibuf[pos], f->readbuf, sizeof(uint8_t) * towrite);
        len = len - towrite;
        if ( len == 0 ) {
          f->readpos = towrite;
          read_entries += towrite;
          break;
        } else {
          pos += towrite;
          read_entries += towrite;
        }
      }

      f->strm.next_out = f->outbuf;
      f->strm.avail_out = sizeof(f->outbuf);
      if ( ret == LZMA_STREAM_END) f->eof = 1;
    }

    if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
      const char *msg;
      switch (ret) {
        case LZMA_MEM_ERROR: msg = "Memory allocation failed"; break;
        case LZMA_FORMAT_ERROR: msg = "The input is not in the .xz format"; break; // .xz magic bytes weren't found.
        case LZMA_OPTIONS_ERROR: msg = "Unsupported compression options"; break;
        case LZMA_DATA_ERROR: msg = "Compressed file is corrupt"; break;
        case LZMA_BUF_ERROR: msg = "Compressed file is truncated or otherwise corrupt"; break;
        case LZMA_MEMLIMIT_ERROR: msg = "The memory limit for decompression is too small"; break;
        default: msg = "Unknown error, possibly a bug"; break;
      }
      fprintf (stderr, "Decoder error: %s (error code %u)\n",msg, ret);
      return 0;
    }
  }
  if ( f->readpos == XZ_BUFFERSIZE) f->readpos = 0;
  return read_entries;
}

/* Local function to write a block to an uncompressed file. */
size_t xz_write(xz_file_t *f, char *cbuf, size_t len) {
  // lzma_action action = LZMA_RUN;
  lzma_ret ret;
  size_t pos = 0;

  if ( f == NULL) return -1;
  if ( cbuf == NULL) return -1;

  while (true) { // Fill the input buffer if it is empty.
    if (f->strm.avail_in == 0 && pos != len) {
      size_t tocopy = MIN(XZ_BUFFERSIZE, (len-pos));
      memcpy(f->inbuf, cbuf+pos, tocopy);
      f->strm.next_in = f->inbuf;
      f->strm.avail_in = tocopy;
      pos = pos + tocopy;
    } else if ( pos >= len && f->strm.avail_in == 0 ) return len;

    ret = lzma_code(&(f->strm), f->action);

    if (f->strm.avail_out == 0 || ret == LZMA_STREAM_END) {
      size_t write_size = sizeof(f->outbuf) - f->strm.avail_out;

      if (fwrite(f->outbuf, 1, write_size, f->fp) != write_size) { fprintf(stderr, "Write error: %s\n", strerror(errno)); return 0; }
      // Reset next_out and avail_out.
      f->strm.next_out = f->outbuf;
      f->strm.avail_out = sizeof(f->outbuf);
    }

    // Normally the return value of lzma_code() will be LZMA_OK until everything has been encoded.
    if (ret != LZMA_OK) {
      const char *msg = "" ;
      if (ret == LZMA_STREAM_END) return true;
      switch (ret) {
        case LZMA_MEM_ERROR: msg = "Memory allocation failed"; break;
        case LZMA_DATA_ERROR: msg = "File size limits exceeded"; break;
        default: break;
      }
      fprintf (stderr, "Encoder error: %s (error code %u)\n",msg, ret);
      return 0;
    }
  }
  return 0;
}

int main (int argc, char **argv)
{
  xz_file_t *f;
  char buf[256];
  size_t n;
  if (argc != 2) { fprintf (stderr, "I need one argument\n"); exit(1); }
  f = xz_open (argv[1], 'r');
  while ((n = xz_read (f, (uint8_t*)buf, 256)) > 0) printf ("<%s>%u\n", buf, n);
  xz_close (f);
  return 0;
}
