/*******************************************************************************
    COPYRIGHT 2002 CATERPILLAR INC.  ALL RIGHTS RESERVED

    This program is the property of Caterpillar Inc., includes Caterpillar's
    confidential and trade secret information, and is maintained as an
    unpublished copyrighted work.  It is not to be copied or used by others
    except under license from Caterpillar.  This program is also protected as an
    unpublished work in accordance with the copyright act of 1976.  In the event
    of either inadvertent or deliberate publication, Caterpillar Inc. intends to
    maintain copyright protection for this work under the relevant copyright
    laws pertaining to published works.  The inclusion of a copyright notice
    hereon is precautionary only, and does not imply publication or disclosure.

 
 Filename      : CubitFileIOWrapper.cpp

 Purpose       : Encapsulates file I/O operations for the Cubit file.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#include "CubitFileIOWrapper.hpp"
#include <cstring>

using namespace NCubitFile;

///////////////////////////////////////////////////////////////////////////////
// CIOWrapper -
///////////////////////////////////////////////////////////////////////////////

CIOWrapper::CIOWrapper(FILE* xpFile, UnsignedInt32 xintSourceEndian)
{
    if(!xpFile)
        throw CCubitFile::ePassedNullPointer;
    mpFile = xpFile;
    mintSwapEndian = (xintSourceEndian != CCubitFile::mintNativeEndian);
    mintBlockStart = mintBlockEnd = 0;
}

CIOWrapper::CIOWrapper(UnsignedInt32 swap_endian, FILE* xpFile )
{
    if(!xpFile)
        throw CCubitFile::ePassedNullPointer;
    mpFile = xpFile;
    mintSwapEndian = swap_endian;
    mintBlockStart = mintBlockEnd = 0;
}


CIOWrapper::CIOWrapper(FILE* xpFile, UnsignedInt32 xintAbsoluteOffset,
                       UnsignedInt32 xintRelativeOffset)
{
    if(!xpFile)
        throw CCubitFile::ePassedNullPointer;
    mpFile = xpFile;
    mintBlockStart = mintBlockEnd = 0;

    if(NCubitFile::SetLocation(mpFile, xintAbsoluteOffset + xintRelativeOffset, SEEK_SET))
        throw CCubitFile::eFileSeekError;
    UnsignedInt32 lintSourceEndian;
    if(fread(&lintSourceEndian, sizeof(UnsignedInt32), 1, mpFile) != 1)
        throw CCubitFile::eFileReadError;
    mintSwapEndian = (lintSourceEndian != CCubitFile::mintNativeEndian);
}

CIOWrapper::~CIOWrapper()
{
}

        
UnsignedInt32 CIOWrapper::BeginWriteBlock(UnsignedInt32 xintAbsoluteOffset)
// Begin the write of a continuous data block, return the location in the file
// where the block begins as a relative offset from the passed value (or
// absolute offset if the passed offset is 0).
{
    if(NCubitFile::SetLocation(mpFile, 0, SEEK_END))
        throw CCubitFile::eFileSeekError;
    UnsignedInt32 pos;
    mintBlockEnd = mintBlockStart = pos = NCubitFile::GetLocation(mpFile);
    return mintBlockStart - xintAbsoluteOffset;
}

void CIOWrapper::BeginRewriteBlock(UnsignedInt32 xintAbsoluteOffset,
                                   UnsignedInt32 xintRelativeOffset)
{
    mintBlockEnd = mintBlockStart = xintAbsoluteOffset + xintRelativeOffset;
    if(NCubitFile::SetLocation(mpFile, mintBlockStart, SEEK_SET))
        throw CCubitFile::eFileSeekError;
}

void CIOWrapper::Write(const UnsignedInt32* xpaintData, UnsignedInt32 xintCount)
{
    if(fwrite(xpaintData, sizeof(UnsignedInt32), xintCount, mpFile) != xintCount)
        throw CCubitFile::eFileWriteError;
    mintBlockEnd += (sizeof(UnsignedInt32) * xintCount);
}

void CIOWrapper::Write(const char* xpachrData, UnsignedInt32 xintCount,
                       UnsignedInt32 xint32bitPadded)
{
    if(xintCount) {
        if(fwrite(xpachrData, sizeof(char), xintCount, mpFile) != xintCount)
            throw CCubitFile::eFileWriteError;
    }
    UnsignedInt32 lintMod = 0;
    if(xint32bitPadded) {
        lintMod = (xintCount % sizeof(UnsignedInt32));
        if(lintMod) {
            lintMod = sizeof(UnsignedInt32) - lintMod;
            if(fwrite("\0\0\0\0", sizeof(char), lintMod, mpFile) != lintMod)
                throw CCubitFile::eFileWriteError;
        }
    }
    mintBlockEnd += (sizeof(char) * (xintCount + lintMod));
}

void CIOWrapper::Write(const double* xpadblData, UnsignedInt32 xintCount)
{
    if(fwrite(xpadblData, sizeof(double), xintCount, mpFile) != xintCount)
        throw CCubitFile::eFileWriteError;
    mintBlockEnd += (sizeof(double) * xintCount);
}

void CIOWrapper::Write(const char* xpachrData)
// Write a null terminated string to the data block, the strings length is 
// written first, followed by the data and padding to keep the file pointer on
// a word boundary (for ease of debugging).
{
    UnsignedInt32 lintLength = xpachrData ? std::strlen(xpachrData) : 0;
    UnsignedInt32 lintMod = (lintLength % sizeof(UnsignedInt32));
    if(fwrite(&lintLength, sizeof(UnsignedInt32), 1, mpFile) != 1)
        throw CCubitFile::eFileWriteError;
    if(lintLength) {
        if(fwrite(xpachrData, sizeof(char), lintLength, mpFile) != lintLength)
            throw CCubitFile::eFileWriteError;
        if(lintMod) {
            lintMod = sizeof(UnsignedInt32) - lintMod;
            if(fwrite("\0\0\0\0", sizeof(char),  lintMod, mpFile) != lintMod)
                throw CCubitFile::eFileWriteError;
        }
    }
    mintBlockEnd +=
        (sizeof(UnsignedInt32) + sizeof(char) * (lintLength + lintMod));
}

UnsignedInt32 CIOWrapper::EndWriteBlock()
// Completes the writing of a contiguous data block, checks for errors and
// returns the length of the block.
{
    UnsignedInt32 lintBlockEnd = NCubitFile::GetLocation(mpFile);
    if((UnsignedInt32)lintBlockEnd != mintBlockEnd)
        throw CCubitFile::eCorruptBlock;
    UnsignedInt32 lintLength = mintBlockEnd - mintBlockStart;
    mintBlockEnd = mintBlockStart = 0;
    return lintLength;
}


void CIOWrapper::BeginReadBlock(UnsignedInt32 xintAbsoluteOffset,
                                UnsignedInt32 xintRelativeOffset)
// Begin the read of a contiguous data block.
{
    if(NCubitFile::SetLocation(mpFile, xintAbsoluteOffset + xintRelativeOffset, SEEK_SET))
        throw CCubitFile::eFileSeekError;
}

void CIOWrapper::Read(UnsignedInt32* xpaintData, UnsignedInt32 xintCount)
{
    if(fread(xpaintData, sizeof(UnsignedInt32), xintCount, mpFile) != xintCount)
        throw CCubitFile::eFileReadError;
    if(mintSwapEndian)
        SwapEndian<UnsignedInt32>(xintCount, xpaintData);
}

void CIOWrapper::Read(char* xpachrData, UnsignedInt32 xintCount,
                      UnsignedInt32 xint32bitPadded)
{
    if(xintCount) {
        if(fread(xpachrData, sizeof(char), xintCount, mpFile) != xintCount)
            throw CCubitFile::eFileReadError;
    }
    if(xint32bitPadded) {
        char lachrPad[8]; // , *lpachrData = NULL;
        UnsignedInt32 lintMod = (xintCount % sizeof(UnsignedInt32));
        if(lintMod) {
            lintMod = sizeof(UnsignedInt32) - lintMod;
            if(fread(&lachrPad, sizeof(char), lintMod, mpFile) != lintMod)
                throw CCubitFile::eFileReadError;
        }
    }
}

void CIOWrapper::Read(double* xpadblData, UnsignedInt32 xintCount)
{
    if(fread(xpadblData, sizeof(double), xintCount, mpFile) != xintCount)
        throw CCubitFile::eFileReadError;
    if(mintSwapEndian)
        SwapEndian<double>(xintCount, xpadblData);
}

char* CIOWrapper::Read()
// Read a null terminated string from the data block, the string's length is 
// read first, followed by the data and then any padding.
{
    UnsignedInt32 lintLength, lintMod;
    char lachrPad[8], *lpachrData = NULL;
    if(fread(&lintLength, sizeof(UnsignedInt32), 1, mpFile) != 1)
        throw CCubitFile::eFileReadError;
    if(mintSwapEndian)
        SwapEndian<UnsignedInt32>(1, &lintLength);
    if(lintLength) {
        lpachrData = new char[lintLength + 1];
        if(!lpachrData)
            throw CCubitFile::eMemoryError;
        if(fread(lpachrData, sizeof(char), lintLength, mpFile) != lintLength) {
		    delete[] lpachrData;
            throw CCubitFile::eFileReadError;
		}	
        lpachrData[lintLength] = '\0';
        lintMod = (lintLength % sizeof(UnsignedInt32));
        if(lintMod) {
            lintMod = sizeof(UnsignedInt32) - lintMod;
            if(fread(&lachrPad, sizeof(char), lintMod, mpFile) != lintMod) {
				delete[] lpachrData;
                throw CCubitFile::eFileReadError;
			}
        }
    }
    return lpachrData;
}

void CIOWrapper::EndReadBlock()
{
    mintBlockStart = mintBlockEnd = 0;
}
    
UnsignedInt32 NCubitFile::GetLocation(FILE* f)
{
#ifdef _MSC_VER
  // normal ftell() returns long, which is a 32 bit signed integer.
  // we use this to increase our 2 GB limit to 4 GB.
  // To go past 4GB, we'd have to use a 64 bit integer instead of UnsignedInt32.
  __int64 offset = _ftelli64(f);
#else
  long offset = ftell(f);
#endif
  if(offset == -1L)
    throw CCubitFile::eFileTellError;
  return static_cast<UnsignedInt32>(offset);
}

UnsignedInt32 CIOWrapper::GetLocation()
{
  return NCubitFile::GetLocation(mpFile);
}

int NCubitFile::SetLocation(FILE* f, UnsignedInt32 offset, int whence)
{
#ifdef _MSC_VER
  return _fseeki64(f, offset, whence);
#else
  return fseek(f, offset, whence);
#endif
}
