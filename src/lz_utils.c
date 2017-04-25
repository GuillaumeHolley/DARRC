#include "lz_utils.h"
#include <stdint.h>
#include <inttypes.h>

int PrintError(char *buffer, const char *message)
{
  strcat(buffer, "\nError: ");
  strcat(buffer, message);
  strcat(buffer, "\n");
  return 1;
}

int PrintErrorNumber(char *buffer, SRes val)
{
  sprintf(buffer + strlen(buffer), "\nError code: %x\n", (unsigned)val);
  return 1;
}

int PrintUserError(char *buffer)
{
  return PrintError(buffer, "Incorrect command");
}

static SRes encode_lzma_internal(ISeqOutStream *outStream, ISeqInStream *inStream, UInt64 fileSize, char *rs, int level_compression)
{
    SRes res;
    CLzmaEncProps props;

    UNUSED_VAR(rs);

    CLzmaEncHandle enc = LzmaEnc_Create(&g_Alloc);
    if (enc == 0) return SZ_ERROR_MEM;

    LzmaEncProps_Init(&props);
    props.level = level_compression;
    res = LzmaEnc_SetProps(enc, &props);

    if (res == SZ_OK)
    {
        Byte header[LZMA_PROPS_SIZE + 8];
        size_t headerSize = LZMA_PROPS_SIZE;
        int i;

        res = LzmaEnc_WriteProperties(enc, header, &headerSize);

        for (i = 0; i < 8; i++) header[headerSize++] = (Byte)(fileSize >> (8 * i));

        if (outStream->Write(outStream, header, headerSize) != headerSize) res = SZ_ERROR_WRITE;
        else if (res == SZ_OK) res = LzmaEnc_Encode(enc, outStream, inStream, NULL, &g_Alloc, &g_Alloc);
    }

    LzmaEnc_Destroy(enc, &g_Alloc, &g_Alloc);
    return res;
}

int encode_lzma(char** filenames_in, int nb_files, char* filename_out, int level_compression)
{
    CFileSeqInStream inStream;
    CFileOutStream outStream;
    UInt64 fileSize;
    int res;

    char rs[800] = { 0 };

    FileSeqInStream_CreateVTable(&inStream);
    File_Construct(&inStream.file);

    FileOutStream_CreateVTable(&outStream);
    File_Construct(&outStream.file);

    size_t t4 = sizeof(UInt32);
    size_t t8 = sizeof(UInt64);
    if (t4 != 4 || t8 != 8) return PrintError(rs, "Incorrect UInt32 or UInt64");

    if (OutFile_Open(&outStream.file, filename_out) != 0) return PrintError(rs, "Can not open output file");
    //else PrintUserError(rs);

    for (int i = 0; i < nb_files; i++){

        if (InFile_Open(&inStream.file, filenames_in[i]) != 0) return PrintError(rs, "Can not open input file");

        File_GetLength(&inStream.file, &fileSize);
        //printf("fileSize is %" PRIu64 "\n", fileSize);
        res = encode_lzma_internal(&outStream.s, &inStream.s, fileSize, rs, level_compression);

        File_Close(&inStream.file);

        if (res != SZ_OK)
        {
            if (res == SZ_ERROR_MEM) return PrintError(rs, "Can not allocate memory");
            else if (res == SZ_ERROR_DATA) return PrintError(rs, "Data error");
            else if (res == SZ_ERROR_WRITE) return PrintError(rs, "Can not write output file");
            else if (res == SZ_ERROR_READ) return PrintError(rs, "Can not read input file");
            return PrintErrorNumber(rs, res);
        }
    }

    File_Close(&outStream.file);

    return 0;
}

static SRes decode_lzma_internal2(CLzmaDec *state, ISeqOutStream *outStream, ISeqInStream *inStream, UInt64 unpackSize)
{
    int thereIsSize = (unpackSize != (UInt64)(Int64)-1);
    Byte inBuf[IN_BUF_SIZE];
    Byte outBuf[OUT_BUF_SIZE];
    size_t inPos = 0, inSize = 0, outPos = 0;
    Int64 nb_bytes_back = 0;

    LzmaDec_Init(state);

    for (;;){

        if (inPos == inSize){
            inSize = IN_BUF_SIZE;
            RINOK(inStream->Read(inStream, inBuf, &inSize));
            inPos = 0;
        }

        SRes res;
        SizeT inProcessed = inSize - inPos;
        SizeT outProcessed = OUT_BUF_SIZE - outPos;
        ELzmaFinishMode finishMode = LZMA_FINISH_ANY;
        ELzmaStatus status;

        if (thereIsSize && outProcessed > unpackSize){
            outProcessed = (SizeT)unpackSize;
            finishMode = LZMA_FINISH_END;
        }

        res = LzmaDec_DecodeToBuf(state, outBuf + outPos, &outProcessed, inBuf + inPos, &inProcessed, finishMode, &status);
        inPos += inProcessed;
        outPos += outProcessed;
        unpackSize -= outProcessed;

        if (!unpackSize){
            nb_bytes_back = 0 - (inSize - inPos);
            RINOK(((ISeekInStream *) inStream)->Seek(inStream, &nb_bytes_back, SZ_SEEK_CUR));
        }

        if (outStream)
            if (outStream->Write(outStream, outBuf, outPos) != outPos) return SZ_ERROR_WRITE;

        outPos = 0;

        if ((res != SZ_OK) || (thereIsSize && unpackSize == 0)) return res;

        if (inProcessed == 0 && outProcessed == 0){
            if (thereIsSize || status != LZMA_STATUS_FINISHED_WITH_MARK) return SZ_ERROR_DATA;
            return res;
        }
    }
}

static SRes decode_lzma_internal(ISeqOutStream *outStream, ISeqInStream *inStream)
{
  UInt64 unpackSize;
  int i;
  SRes res = 0;

  CLzmaDec state;

  /* header: 5 bytes of LZMA properties and 8 bytes of uncompressed size */
  unsigned char header[LZMA_PROPS_SIZE + 8];

  /* Read and parse header */
  RINOK(SeqInStream_Read(inStream, header, sizeof(header)));

  unpackSize = 0;
  for (i = 0; i < 8; i++) unpackSize += (UInt64)header[LZMA_PROPS_SIZE + i] << (i * 8);

  //printf("unpackSize = %" PRIu64 "\n", unpackSize);

  LzmaDec_Construct(&state);
  RINOK(LzmaDec_Allocate(&state, header, LZMA_PROPS_SIZE, &g_Alloc));
  res = decode_lzma_internal2(&state, outStream, inStream, unpackSize);
  LzmaDec_Free(&state, &g_Alloc);

  return res;
}

int decode_lzma(char* filename_in, int nb_files, char** filenames_out)
{
    CFileInStream inStream;
    CFileOutStream outStream;
    int res;

    char rs[800] = { 0 };

    FileInStream_CreateVTable(&inStream);
    File_Construct(&inStream.file);

    FileOutStream_CreateVTable(&outStream);
    File_Construct(&outStream.file);

    size_t t4 = sizeof(UInt32);
    size_t t8 = sizeof(UInt64);
    if (t4 != 4 || t8 != 8) printf("Incorrect UInt32 or UInt64");

    if (InFile_Open(&inStream.file, filename_in) != 0) printf("Can not open input file");

    for (int i = 0; i < nb_files; i++){

        if (OutFile_Open(&outStream.file, filenames_out[i]) != 0) printf("Can not open output file");

        res = decode_lzma_internal(&outStream.s, (ISeqInStream *) &inStream.s);

        File_Close(&outStream.file);

        if (res != SZ_OK)
        {
            if (res == SZ_ERROR_MEM) printf("Can not allocate memory\n");
            else if (res == SZ_ERROR_DATA) printf("Data error\n");
            else if (res == SZ_ERROR_WRITE) printf("Can not write output file\n");
            else if (res == SZ_ERROR_READ) printf("Can not read input file\n");
            else if (res == SZ_ERROR_UNSUPPORTED) printf("Unsupported properties\n");
            else if (res == SZ_ERROR_INPUT_EOF) printf("It needs more bytes in input buffer (src)\n");

            PrintErrorNumber(rs, res);
        }
    }

    File_Close(&inStream.file);

    return 0;
}
