#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "Precomp.h"
#include "Alloc.h"
#include "LzFind.h"
#include "LzmaDec.h"
#include "LzmaEnc.h"
#include "7zFile.h"

#define COMP_LVL_MAX 9

#define IN_BUF_SIZE (1 << 16)
#define OUT_BUF_SIZE (1 << 16)

int encode_lzma(char** filenames_in, int nb_files, char* filename_out, int level_compression);
int decode_lzma(char* filename_in, int nb_files, char** filenames_out);

int PrintUserError(char *buffer);
int PrintErrorNumber(char *buffer, SRes val);
int PrintError(char *buffer, const char *message);
