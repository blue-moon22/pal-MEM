/* ================================================================= *
 *  file.h : Header file with supporting class definitions           *
 *                                                                   *
 *  pal-MEM: An efficient (MUMmer-like) tool to retrieve reads       *
 *          containing Maximum Exact Matches using hashing based     *
 *          algorithm                                                *
 *  Developer: Victoria Carr                                         *
 *                                                                   *
 *  Modified program from E-MEM                                      *
 *  https://github.com/lucian-ilie/E-MEM                             *
 *  Developer: Nilesh Khiste                                         *
 *                                                                   *
 *  This program is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public License as   *
 *  published by the Free Software Foundation, either version 3 of   *
 *  the License, or (at your option) any later version.              *
 *                                                                   *
 *  This program is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 *  GNU General Public License for more details.                     *
 *                                                                   *
 *  You should have received a copy of the GNU General Public        *
 *  License along with this program.                                 *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'LICENSE', which is part of this source code package.       *
 * ================================================================= */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <iterator>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <sys/stat.h>
//#include <chrono>

#include "pal-mem.h"
#include "file.h"

using namespace std;
using namespace boost;
//using namespace std::chrono;

/*
 * Function builds a kmer hash for a reference sequence.
 * Input: empty refHash
 * Output: populated refHash
 */
void buildRefHash(Knode* &refHash, uint64_t totalBits, seqFileReadInfo &RefFile)
{
    uint64_t j=0;
    uint64_t currKmerPos=0, currKmer=0;
    int32_t offset=0;
    int nextKmerPosition = commonData::minMemLen - commonData::kmerSize + 2;

    vector<mapObject>::iterator it;
    it = upper_bound(RefFile.blockOfNs.begin(), RefFile.blockOfNs.end(), currKmerPos, mapObject());
    while (currKmerPos<=totalBits)
    {
        if (currKmerPos + commonData::kmerSize - 2 > totalBits)
            break;

        if(RefFile.checkKmerForNs(currKmerPos, it)){
            currKmerPos+=nextKmerPosition; // Move L-K+2 bits = 50-28+1=23 char = 46 bits
            continue;
        }

        offset = currKmerPos%DATATYPE_WIDTH;
        j=currKmerPos/DATATYPE_WIDTH; // next loc in binReads

        currKmer = RefFile.binReads[j]; // Get the unsigned 64-bit
        currKmer <<= offset; // Get base by shifting left by offset

        if (offset > DATATYPE_WIDTH-commonData::kmerSize) // Kmer split in two integers
            // Append the remaining bits from the next bin to make the full kmer
            currKmer |= ((RefFile.binReads[j+1] & global_mask_left[(commonData::kmerSize-(DATATYPE_WIDTH-offset))/2 -1])>>(DATATYPE_WIDTH-offset));
        else
            currKmer &= global_mask_left[commonData::kmerSize/2 - 1];
        /* Add kmer to the hash table */
        refHash->addKmerNode(currKmer, currKmerPos);
        currKmerPos+=nextKmerPosition; // Move L-K+2 bits = 50-28+1=23 char = 46 bits
    }
}

/*
 * Function extends the kmer match in left/right direction for
 * possible MEMs.
 * Input: currRPos : current position of matching reference Kmer
 * Input: currRPos : current position of matching query Kmer
 * Input: totalRBits : total number of bits in reference
 * Input: totalQBits : total number of bits in query
 * Input: name : reference sequence string for output
 *
 */
void helperReportMem(uint64_t &currRPos, uint64_t &currQPos, uint64_t totalRBits, uint64_t totalQBits, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, mapObject &RefNpos, mapObject &QueryNpos, uint64_t &rQMEM, uint64_t &lRMEM, uint64_t &rRMEM, vector<seqData> &vecSeqInfo)
{
    /*
     * lRef and lQue are local variables for left extension of
     * reference and query sequence respectively. rRef and rQue
     * are their right counterparts.
     */
    uint64_t lRef=currRPos, lQue=currQPos; // Keeping lRef on currRPos-this makes offset computation simpler
    uint64_t offsetR=0,offsetQ=0;
    uint64_t rRef=currRPos+commonData::kmerSize, rQue=currQPos+commonData::kmerSize; // one character ahead of current match
    uint64_t extLength;
    uint64_t currR=0, currQ=0, lQtmp=0, rQtmp=0;
    int i=0,j=0,mismatch=0,binLength,shortestLenL,shortestLenR;
    uint64_t matchSize=0;
    seqData s;

    if (!(((QueryNpos.left==0x1)?true:QueryNpos.left<=lQue) && rQue<=QueryNpos.right))
        QueryFile.getKmerLeftnRightBoundForNs(lQue, QueryNpos);

    if (!(((RefNpos.left==0x1)?true:RefNpos.left<=lRef) && rRef<=RefNpos.right))
        RefFile.getKmerLeftnRightBoundForNs(lRef, RefNpos);

    if (RefNpos.right-((RefNpos.left==0x1)?0:RefNpos.left)+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    if (QueryNpos.right-((QueryNpos.left==0x1)?0:QueryNpos.left)+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    //match towards left
    while (lRef && lQue && ((QueryNpos.left==0x1)?true:QueryNpos.left<=lQue) && ((RefNpos.left==0x1)?true:RefNpos.left<=lRef))
    {
        if (!mismatch)
        {
            offsetR=(lRef)%DATATYPE_WIDTH;
            i=(lRef)/DATATYPE_WIDTH;
            offsetQ=(lQue)%DATATYPE_WIDTH;
            j=(lQue)/DATATYPE_WIDTH;

            if (offsetR > offsetQ)
                matchSize = offsetQ;
            else
                matchSize = offsetR;

            if (!matchSize)
                matchSize=2;

            if ((QueryNpos.left!=0x1) && (matchSize > lQue-QueryNpos.left))
                matchSize = lQue-QueryNpos.left;

            if ((RefNpos.left!=0x1) && (matchSize > lRef-RefNpos.left))
                matchSize = lRef-RefNpos.left;

            if (!matchSize)
                break;

            /*
             * There will never be case with offset=0 and i=0 because
             * i=0 happens only when lRef=0 and in that case we do not
             * enter this loop.
             */
            currR = RefFile.binReads[offsetR?i:i-1];
            if (offsetR)
                currR >>= DATATYPE_WIDTH-offsetR;
            currQ = QueryFile.binReads[offsetQ?j:j-1];
            if (offsetQ)
                currQ >>= DATATYPE_WIDTH-offsetQ;
        }

        if((currR & global_mask_right[matchSize/2 - 1]) != (currQ &  global_mask_right[matchSize/2 - 1])) {
            if (matchSize==2)
                break;

            mismatch=1;
            matchSize/=2;
            if (matchSize%2)
                matchSize+=1;
        }else {
            lRef-=matchSize;
            lQue-=matchSize;
            if (matchSize > static_cast<uint64_t>(commonData::maxMemLen)) {
                if (lRef?(lRef - RefNpos.left <= commonData::lenBuffer):!lRef) {
                    rRMEM = RefNpos.right;
                    lRMEM = RefNpos.left;
                }
                return;
            }

            if (mismatch) {
                if (matchSize==2)
                    break;
                currR >>= matchSize;
                currQ >>= matchSize;
            }
        }
    }

    /* Ignore reverse complements of the same ORF, i.e. where matching prefix/suffix of reference/query
     * Also excludes N mismatches
     */
    if (lRef?(lRef - RefNpos.left <= commonData::lenBuffer):!lRef) {
        rRMEM = RefNpos.right;
        lRMEM = RefNpos.left;
        return;
    }

    if (lQue?(lQue - QueryNpos.left <= commonData::lenBuffer):!lQue)
        return;

    if (totalRBits-lRef+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    if (totalQBits-lQue+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    //match towards right
    mismatch=0;
    while ((rRef <= totalRBits) && (rQue <= totalQBits) && (rRef <= RefNpos.right) && (rQue <= QueryNpos.right))
    {
        if (!mismatch)
        {
            offsetR=rRef%DATATYPE_WIDTH;
            i=rRef/DATATYPE_WIDTH;
            offsetQ=rQue%DATATYPE_WIDTH;
            j=rQue/DATATYPE_WIDTH;

            if (offsetR > offsetQ)
                matchSize = DATATYPE_WIDTH-offsetR;
            else
                matchSize = DATATYPE_WIDTH-offsetQ;

            if (rRef+matchSize > totalRBits)
                matchSize = totalRBits-rRef;

            if (rQue+matchSize > totalQBits)
                matchSize = totalQBits-rQue;

            if (rQue+matchSize > QueryNpos.right)
                matchSize = QueryNpos.right-rQue;

            if (rRef+matchSize > RefNpos.right)
                matchSize = RefNpos.right-rRef;

            if(!matchSize)
                matchSize=2;



            currR = RefFile.binReads[i];
            currR <<= offsetR;
            currQ = QueryFile.binReads[j];
            currQ <<= offsetQ;
        }

        if((currR & global_mask_left[matchSize/2 - 1]) != (currQ &  global_mask_left[matchSize/2 - 1])) {
            if (matchSize==2){
                rRef-=2;
                rQue-=2;
                break;
            }

            mismatch=1;
            matchSize/=2;
            if (matchSize%2)
                matchSize+=1;
        }else {
            if (mismatch) {
                if (matchSize==2) {
                    break;
                }
            }
            if ((rRef == totalRBits) || (rQue == totalQBits))
                return;

            currR <<= matchSize;
            currQ <<= matchSize;
            rRef+=matchSize;
            rQue+=matchSize;
            if (matchSize > static_cast<uint64_t>(commonData::maxMemLen)) {
                if ((RefNpos.right - rRef) <= commonData::lenBuffer) {
                    rRMEM = RefNpos.right;
                    lRMEM = RefNpos.left;
                }
                return;
            }
        }
    } // loop until extension reached to the right

    if ((RefNpos.right - rRef) <= commonData::lenBuffer) {
        rRMEM = RefNpos.right;
        lRMEM = RefNpos.left;
        return;
    }

    if ((QueryNpos.right - rQue) <= commonData::lenBuffer)
        return;

    /* Adjust rRef and rQue locations */

    if (rRef > RefNpos.right){
        rQue-=(rRef-RefNpos.right);
        rRef=RefNpos.right;
    }
    if (rQue > QueryNpos.right){
        rRef-=(rQue-QueryNpos.right);
        rQue=QueryNpos.right;
    }

    if (rRef > totalRBits){
        rQue-=(rRef-totalRBits);
        rRef=totalRBits;
    }
    if (rQue > totalQBits){
        rRef-=(rQue-totalQBits);
        rQue=totalQBits;
    }

    if (((RefNpos.right - rRef) <= commonData::lenBuffer) || (lRef?(lRef - RefNpos.left <= commonData::lenBuffer):!lRef)) {
        rRMEM = RefNpos.right;
        lRMEM = RefNpos.left;
        return;
    }

    if (((QueryNpos.right - rQue) <= commonData::lenBuffer) || (lQue?(lQue - QueryNpos.left <= commonData::lenBuffer):!lQue))
        return;

    if (rRef-lRef > static_cast<uint64_t>(commonData::maxMemLen) || rRef-lRef < static_cast<uint64_t>(commonData::minMemLen))
        return;

    lQtmp = ((QueryNpos.left == 1)?(QueryNpos.left + (QueryNpos.right - rQue) - 1):(QueryNpos.left + (QueryNpos.right - rQue)));
    rQtmp = ((QueryNpos.left == 1)?(QueryNpos.left + (QueryNpos.right - lQue) - 1):(QueryNpos.left + (QueryNpos.right - lQue)));
    arrayTmpFile.getInvertedRepeats(lQtmp, rQtmp, QueryFile, rRef, lRef, RefFile, vecSeqInfo);
    rQMEM = QueryNpos.right;
}

void reportMEM(Knode* &refHash, uint64_t totalBases, uint64_t totalQBases, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, vector<seqData> &vecSeqInfo)
{
    uint64_t totalQBits = CHARS2BITS(totalQBases); // convert char position to bit position of QueryFile.totalBases-1
    uint32_t copyBits=0;
    //auto start = high_resolution_clock::now();
    #pragma omp parallel num_threads(commonData::numThreads)
    {
        uint64_t currKmer=0, j=0;
        int32_t offset=0;
        uint32_t first=1;
        uint64_t rQMEM=0,rRMEM=0,lRMEM=0;
        int kmerWithNs=0;
        mapObject QueryNpos, RefNpos;
        vector<SeqPos> SeqPosVec;
        vector<mapObject>::iterator it;
        it = upper_bound(QueryFile.blockOfNs.begin(), QueryFile.blockOfNs.end(), 0, mapObject());

        #pragma omp single
        {
        /*
        * Number of copy bits during query kmer processing depends on kmer size.
        */
        if (DATATYPE_WIDTH-commonData::kmerSize > 32 )
            copyBits=32; //16 characters
        else if (DATATYPE_WIDTH-commonData::kmerSize > 16) //kmerSize = 40 so copyBits = 16
            copyBits=16; //8 characters
        else
            copyBits=8; //4 characters

        /* If copyBits more than 8, the for loop parallelisation will give
        * incorrect results - miss some Mems
        */
        if(commonData::numThreads > 1)
            copyBits=8; //4 characters // if more than one thread, then copyBits = 8
        }

        // parallelise search of query kmer in reference
        #pragma omp for
        for (uint64_t currKmerPos=0; currKmerPos<=totalQBits; currKmerPos+=2)
        {
            if ((rQMEM?currKmerPos >= (rQMEM + RANDOM_SEQ_SIZE*2):currKmerPos >= rQMEM)) {
                //QueryFile.getKmerLeftnRightBoundForNs(currKmerPos, QueryNpos);
                //cout << currKmerPos << endl;
                if ((currKmerPos + commonData::kmerSize - 2) > totalQBits)
                    continue;

                if(QueryFile.checkKmerForNs(currKmerPos, it)){
                    kmerWithNs=1;
                }

                j=currKmerPos/DATATYPE_WIDTH;// current location in binReads
                offset = currKmerPos%DATATYPE_WIDTH;
                if(first || !offset){
                    currKmer = QueryFile.binReads[j];
                    currKmer <<= offset;
                    if(offset > DATATYPE_WIDTH-commonData::kmerSize)
                      currKmer |= ((QueryFile.binReads[j+1] & global_mask_left[offset/2-1])>>(DATATYPE_WIDTH-offset));
                    first=0;
                }else
                    currKmer <<= 2;

                if(offset  && !(offset % copyBits))
                  currKmer |= ((QueryFile.binReads[j+1] & global_mask_left[offset/2-1])>>(DATATYPE_WIDTH-offset));

                if (kmerWithNs){
                    /* Do not process this Kmer, Ns in it */
                    kmerWithNs=0;
                    continue;
                }
                /* Find the K-mer in the refHash */
                uint64_t *dataPtr = NULL;
                if (refHash->findKmer(currKmer & global_mask_left[commonData::kmerSize / 2 - 1],
                                      dataPtr)) // dataPtr is position of the kmer in reference
                {
                    // We have a match
                    for (uint64_t n = 1; n <= dataPtr[0]; n++) { // currKmerPos is position of kmer in que
                        if (!(dataPtr[n] >= lRMEM && dataPtr[n] <= rRMEM) && (currKmerPos?currKmerPos > rQMEM:!currKmerPos))
                            helperReportMem(dataPtr[n], currKmerPos, CHARS2BITS(totalBases), CHARS2BITS(totalQBases),RefFile, QueryFile, arrayTmpFile, RefNpos, QueryNpos, rQMEM, lRMEM, rRMEM, vecSeqInfo);
                        else
                            break;
                    }
                }
            }
        }
    }
    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(stop - start);
    //cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
}

void searchQuery(Knode* &refHash, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, vector<seqData> &vecSeqInfo)
{
    QueryFile.resetCurrPos();
    reportMEM(refHash, RefFile.totalBases-1, QueryFile.totalBases-1, RefFile, QueryFile, arrayTmpFile, vecSeqInfo);
    // QueryFile.setCurrPos();
    // QueryFile.clearMapForNs();
    // QueryFile.clearTmpString();
}

void processReference(seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, vector<seqData> &vecSeqInfo)
{
    uint64_t numberOfKmers=0,n=0;
    int hashTableSizeIndex=0;
    Knode *refHash;

    numberOfKmers = ceil((RefFile.totalBases-commonData::kmerSize/2+1)/((commonData::minMemLen/2-commonData::kmerSize/2 + 1)) + 1);

    /* Set the size of the hash table to the numberofKmers. */
    for (n=0; n<450; ++n)
    {
        if (hashTableSize[n] > 1.75*numberOfKmers)
        {
            hashTableSizeIndex = n;
            break;
        }
    }

    Knode::currHashTabSize = hashTableSize[hashTableSizeIndex];  //Store the size of the hash table.
    if (hashTableSizeIndex)
        Knode::prevHashTabSize = hashTableSize[hashTableSizeIndex-1];
    else
        Knode::prevHashTabSize = 3;

    /* Create the refHash for K-mers. */
    cout << "Allocating hashtable..." << endl;
    refHash = new Knode[Knode::currHashTabSize];
    cout << "Building Ref hashtable..." << endl;
    buildRefHash(refHash, CHARS2BITS(RefFile.totalBases-1), RefFile);
    cout << "Searching query..." << endl;
    searchQuery(refHash, RefFile, QueryFile, arrayTmpFile, vecSeqInfo);

    delete [] refHash;
}

bool is_numeric(const string &str)
{
    return all_of(str.begin(), str.end(), ::isdigit);
}

void checkCommandLineOptions(uint32_t &options)
{
    if (!IS_FASTAU_DEF(options)){
        if (!IS_FASTA1_DEF(options) || !IS_FASTA2_DEF(options)){
            cout << "ERROR: both f1 forward and f2 reverse fasta files or fu unpaired fasta file must be passed!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!IS_DATABASE_DEF(options)){
        cout << "ERROR: query file must be passed!" << endl;
        exit(EXIT_FAILURE);
    }

    if (!IS_OUT_FILE_DEF(options)) {
        cout << "ERROR: output file must be passed!" << endl;
        exit(EXIT_FAILURE);
    }

    if (IS_NUM_THREADS_DEF(options)){
        if (commonData::numThreads <= 0){
            cout << "ERROR: -t cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if (IS_MIN_LENGTH_DEF(options)){
        if (commonData::minMemLen <= commonData::kmerSize){
            cout << "ERROR: -l cannot be less than or equal to the k-mer length 15!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if (IS_KMER_SIZE_DEF(options)){
        if (commonData::kmerSize <= 0){
            cout << "ERROR: -k cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }

        if (commonData::kmerSize > CHARS2BITS(28))
        {
            cout << "ERROR: kmer-size cannot be greater than 28" << endl;
            exit( EXIT_FAILURE );
        }
    }
}

void print_help_msg()
{
    cout <<  endl;
    cout << "pal-MEM Version 2.3.6, Oct. 24, 2022" << endl;
    cout << "Adapted from E-MEM Version 1.0.2, Dec. 12, 2017, by Nilesh Khiste and Lucian Ilie" << endl;
    cout <<  endl;
    cout << "pal-MEM outputs two fasta files and a tab-delimited file. One fasta file contains reads" << endl;
    cout << "with inverted repeats (IRs) and the other file contains reads without IRs." << endl;
    cout << "The tab-delimited file contains the original headers of read pairs containing the IRs." << endl;
    cout << endl;
    cout << "Usage: pal-mem [required] [optional]" << endl;
    cout << endl;
    cout << "Required:" << endl;
    cout << "-f1 <filename> -f2 <filename>\t" << "fasta files with forward and reverse paired-end reads" << endl;
    cout << "OR" << endl;
    cout << "-fu <filename>\t\t\t" << "fasta file with unpaired reads" << endl;
    cout << endl;
    cout << "-d <filename>\t\t\t" << "database fasta file" << endl;
    cout << "-o <output prefix>\t\t" << "prefix name of output files" << endl;
    cout << endl;
    cout << "Optional:" << endl;
    cout << "-l\t" << "set the minimum length of a match. Default: 24" << endl;
    cout << "-k\t" << "set the k-mer length. Default: 15" << endl;
    cout << "-t\t" << "number of threads. Default: 1" << endl;
    cout << "-h\t" << "show possible options" << endl;
}

int main (int argc, char *argv[])
{
    int32_t i=0, n=1;
    uint32_t options=0;
    seqFileReadInfo QueryFile, RefFile;
    outFileReadInfo OutFiles;
    string fasta1, fasta2, fastaU, database, outPrefix;

    // Check Arguments
    if (argc==1 || argc==2){
        print_help_msg();
        exit(EXIT_SUCCESS);
    }

    while(argv[n]) {

        if(boost::equals(argv[n], "-f1")){
            if (IS_FASTA1_DEF(options)){
                cout << "ERROR: f1 argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (IS_FASTAU_DEF(options)){
                cout << "ERROR: the f1 and f2 arguments must be used together or fu argument used alone!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FASTA1(options);
            fasta1 = argv[n+1];
            n+=2;
        }else if(boost::equals(argv[n], "-f2")){
            if (IS_FASTA2_DEF(options)){
                cout << "ERROR: f2 argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (IS_FASTAU_DEF(options)){
                cout << "ERROR: the f1 and f2 arguments must be used together or fu argument used alone!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FASTA2(options);
            fasta2 = argv[n+1];
            n+=2;
        }else if(boost::equals(argv[n], "-fu")) {
            if (IS_FASTAU_DEF(options)) {
                cout << "ERROR: fu argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (IS_FASTA1_DEF(options) || IS_FASTA2_DEF(options)) {
                cout << "ERROR: the f1 and f2 arguments must be used together or fu argument used alone!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FASTAU(options);
            fastaU = argv[n + 1];
            n += 2;
        }else if(boost::equals(argv[n], "-d")){
            if (IS_DATABASE_DEF(options)){
                cout << "ERROR: d argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_DATABASE(options);
            database = argv[n+1];
            n+=2;
        } else if(boost::equals(argv[n], "-o")){
            if (IS_OUT_FILE_DEF(options)) {
                cout << "ERROR: output prefix argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_OUT_FILE(options);
            outPrefix = argv[n+1];
            n+=2;
        }else if(boost::equals(argv[n],"-l")){
            if (IS_MIN_LENGTH_DEF(options)) {
                cout << "ERROR: Minimum length argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_MIN_LENGTH(options);
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -l option!" << endl;
                exit(EXIT_FAILURE);
            }
            commonData::minMemLen = 2*std::stoi(argv[n+1]);
            commonData::lenBuffer = commonData::minMemLen - 2;
            n+=2;
        }else if(boost::equals(argv[n],"-m")){
            if (IS_MAX_LENGTH_DEF(options)) {
                cout << "ERROR: Maximum length argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_MAX_LENGTH(options);
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -m option!" << endl;
                exit(EXIT_FAILURE);
            }
            commonData::maxMemLen = 2*std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-t")){
            if (IS_NUM_THREADS_DEF(options)) {
                cout << "ERROR: Number of threads argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -t option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_NUM_THREADS(options);
            commonData::numThreads = std::stoi(argv[n+1]);
            n+=2;
          }else if (boost::equals(argv[n],"-k")){
              if (IS_KMER_SIZE_DEF(options)) {
                  cout << "ERROR: Kmer size argument passed multiple times!" << endl;
                  exit(EXIT_FAILURE);
              }
              if (!argv[n+1] || !is_numeric(argv[n+1])){
                  cout << "ERROR: Invalid value for -k option!" << endl;
                  exit(EXIT_FAILURE);
              }
              SET_KMER_SIZE(options);
              commonData::kmerSize = 2*std::stoi(argv[n+1]);
              n+=2;
        }else if (argv[n][0] != '-'){
            cout << "ERROR: option must start with '-'!" << endl;
            exit(EXIT_FAILURE);
        }else if (boost::equals(argv[n],"-h")){
            print_help_msg();
            exit(EXIT_SUCCESS);
        }else {
            cout << "ERROR: Invalid option." << endl << flush;
            print_help_msg();
            exit( EXIT_FAILURE );
        }
    }

    checkCommandLineOptions(options);

    cout << "Opening tmp files..." << endl;
    sprintf(commonData::nucmer_path, "%s/%d_tmp", getenv("NUCMER_E_MEM_OUTPUT_DIRPATH")?getenv("NUCMER_E_MEM_OUTPUT_DIRPATH"):".",getpid());

    tmpFilesInfo arrayTmpFile(0);
    arrayTmpFile.openFiles(ios::out|ios::binary, 0, outPrefix + "_IR.tab");

    cout << "Opening fasta files..." << endl;
    if (IS_FASTA1_DEF(options) && IS_FASTA2_DEF(options)){
        vector<string> filenames = {fasta1, fasta2};
        QueryFile.setFiles(filenames);
    } else if (IS_FASTAU_DEF(options)) {
        vector<string> filenames = {fastaU};
        QueryFile.setFiles(filenames);
    }
    cout << "Opening query..." << endl;
    vector<string> database_file = {database};
    RefFile.setFiles(database_file);

    cout << "Generating reverse complement..." << endl;
    QueryFile.generateRevComplement();
    RefFile.setSize(QueryFile.getSize());
    RefFile.setNumSequences(QueryFile.getNumSequences());

    cout << "Creating seqData..." << endl;
    vector<seqData> querySeqInfo;
    querySeqInfo.reserve(QueryFile.getNumSequences());
    QueryFile.generateSeqPos(querySeqInfo);

    QueryFile.setReverseFile();

    arrayTmpFile.setNumMemsInFile(QueryFile.allocBinArray(0), QueryFile.getNumSequences());

    cout << "Allocate Ref bins..." << endl;
    RefFile.allocBinArray(1);

    RefFile.clearFileFlag();
    QueryFile.clearFileFlag();
    cout << "Encoding Query sequences..." << endl;

    if (QueryFile.readChunks())
    {
        for (i=0; i<commonData::d; i++) {
            cout << "Encoding chunk " << i+1 << " of Ref sequences..." << endl;
            if(RefFile.readChunks()){ // Encode sequence as 2-bits in RefFile object
                cout << "Getting inverted repeats..." << endl;
                processReference(RefFile, QueryFile, arrayTmpFile, querySeqInfo); // Build hashtable, query hashtable
                RefFile.setCurrPos();
                RefFile.clearMapForNs(); // clear block of Ns from memory

            }
            else
                break;
        }
    }

    cout << "Writing fasta files..." << endl;
    if (IS_FASTA1_DEF(options) && IS_FASTA2_DEF(options)){
        vector<string> filenames = {fasta1, fasta2};
        OutFiles.writeFiles(querySeqInfo, filenames, outPrefix, outPrefix);
    } else if (IS_FASTAU_DEF(options)) {
        vector<string> filenames = {fastaU};
        OutFiles.writeFiles(querySeqInfo, filenames, outPrefix, outPrefix);
    }

    QueryFile.closeFile();

    QueryFile.destroy();
    arrayTmpFile.removeTmp();
    arrayTmpFile.closeFile();
    RefFile.closeFile();

    RefFile.destroy();

    fflush(0);
    return 0;
}
