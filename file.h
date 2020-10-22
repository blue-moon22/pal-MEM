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
using namespace std;

#define CHARS2BITS(x) 		(2*(x)) //convert char position to bit position
#define DATATYPE_WIDTH          64 	// number of bits
#define RANDOM_SEQ_SIZE         10
#define NUM_TMP_FILES           24

class commonData {
  public:
    static int32_t minMemLen;
    static int32_t d;
    static int32_t numThreads;
    static int32_t kmerSize;
    static char nucmer_path[256];
    static int32_t lenBuffer;
};

int32_t commonData::minMemLen=48; // 2 bit representation
int32_t commonData::lenBuffer=46;
int32_t commonData::d=1;
int32_t commonData::numThreads=1;
int32_t commonData::kmerSize=30; //2 bit representation
char commonData::nucmer_path[256]={'\0'};\


class seqData {
  public:
    uint64_t start;
    uint64_t end;
    std::string seq;
    std::string header;
    int keep=1;
    seqData()
    {
        start=0;
        end=0;
    };

    bool operator ()(const seqData &obj1, const seqData &obj2)
    {
      return (obj2.start>obj1.end?true:false);
    }
};

class posData {
public:
    uint64_t L1Bound;
    uint64_t R1Bound;
    uint64_t L2Bound;
    uint64_t R2Bound;
    uint64_t palLength;
    std::string seq;

    posData() {}

    bool operator () (const posData &obj1, const posData &obj2)
    {
        if (obj1.L1Bound < obj2.L1Bound) {
            return true;
        } else if (obj1.L1Bound > obj2.L1Bound) {
            return false;
        } else {
            if (obj1.L2Bound < obj2.L2Bound)
                return true;
            else
                return false;
        }
    }
};

class mapObject {
    public:
        uint64_t left;
        uint64_t right;
        mapObject() {
          left=0;
          right=0;
        }

        mapObject(uint64_t x, uint64_t y) {
          left=x;
          right=y;
        }

        bool operator()(const uint64_t &x, const mapObject &y)
        {
          return x < y.left;
        }
};


class seqFileReadInfo {
    fstream file1, file2;
    int numSeqFiles;
    uint64_t size;
    string strTmp, strName;
    uint64_t binReadSize;
    uint64_t binReadsLocation;
    uint64_t currPos;
    uint64_t numSequences;

        string& randomStr()
        {
             static string str("NNNNNNNNNN");
             return str;
        }

        void processTmpString(uint64_t &sz, uint64_t &blockNCount)
        {
            string line = strTmp;
            strTmp.clear();
            totalBases=0;
            binReadsLocation=0;
            processInput(line, sz, blockNCount);
        }

      /*
       * Function converts a character sequence into an array of integers.
       * Input: character string
       * Output: array of integers, total number of bases
       */
    void processInput(string &str, uint64_t &sz, uint64_t &blockNCount)
    {
        int chooseLetter=0;

        /* Processing the sequences by encoding the base pairs into 2 bits. */
        for ( std::string::iterator it=str.begin(); it!=str.end(); ++it){
            if (totalBases == sz){
              strTmp += *it;
              continue;
            }else if (totalBases >= size) {
              strTmp += *it;
            }
            switch(*it)
                {
                    case 'A':
                    case 'a':
                    binReads[binReadsLocation] <<= 2; // shift left by 2 bits
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    case 'C':
                    case 'c':
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 1;
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    case 'G':
                    case 'g':
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 2;
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    case 'T':
                    case 't':
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 3;
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    default:
                    if(!blockNCount)
                        blockNCount=totalBases+1;
                        chooseLetter = rand() % 4;
                    if (chooseLetter == 0)
                        binReads[binReadsLocation] <<= 2;
                    else if (chooseLetter == 1)
                    {
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 1;
                    }
                    else if (chooseLetter == 2)
                    {
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 2;
                    }
                    else
                    {
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 3;
                    }
                }
            totalBases++;
            if ((totalBases%32)==0){
                binReadsLocation++;
            }
        }
    }

    public:
      uint64_t *binReads;
      uint64_t totalBases;
      std::vector <mapObject> blockOfNs;

      seqFileReadInfo() {
          size=0;
          currPos=0;
          binReadSize=0;
          binReadsLocation=0;
          numSequences=0;
          totalBases=0;
          numSeqFiles=0;
      }

      void setFiles(vector<string> &filenames)
      {
          numSeqFiles=filenames.size();
          if (numSeqFiles >= 1){
              file1.open(filenames[0], ios::in);
              if(!file1.is_open()) {
                  cout << "ERROR: unable to open " << filenames[0] << " file" << endl;
                  exit( EXIT_FAILURE );
              }
              if (numSeqFiles == 2){
                  file2.open(filenames[1], ios::in);
                  if(!file2.is_open()) {
                      cout << "ERROR: unable to open " << filenames[1] << " file" << endl;
                      exit( EXIT_FAILURE );
                  }
              }
          }
      }

      void setNumSequences(uint64_t numSeq) {
          numSequences = numSeq;
      }

      void setSize(uint64_t sizeSet) {
          size = sizeSet;
      }

      uint64_t &getNumSequences() {
          return numSequences;
      }

      uint64_t &getSize() {
          return size;
      }

      void openFile(string s, fstream &file){
          file.open(s, ios::in);
          if(!file.is_open()) {
              cout << "ERROR: unable to open "<< s << " file" << endl;
              exit( EXIT_FAILURE );
          }
      }

      void setReverseFile() {
          char buffer[256];
          memset(buffer,0,256);
          sprintf(buffer, "%s/revComp", commonData::nucmer_path);

          if (numSeqFiles >= 1){
              file1.close();
              if (numSeqFiles == 2){
                  file2.close();
              }
          }
          openFile(buffer, file1);
          numSeqFiles = 1;
      }

      void closeFile() {
          if (numSeqFiles >= 1){
              file1.close();
              if (numSeqFiles == 2){
                  file2.close();
              }
          }
      }

      void destroy() {
          currPos=0;
          binReadSize=0;
          binReadsLocation=0;
          totalBases=0;
          strName.clear();
          strTmp.clear();
          clearMapForNs();
          delete [] binReads;
      }

      void clearFileFlag()
      {
          if (numSeqFiles >= 1){
              file1.clear();
              file1.seekg(0, ios::beg);
              if (numSeqFiles == 2){
                  file2.clear();
                  file2.seekg(0, ios::beg);
              }
          }
      }

      uint64_t allocBinArray(int split)
      {
          if (split)
              size = size/commonData::d;
          binReadSize = floor((size+numSequences*RANDOM_SEQ_SIZE+commonData::d)/32+4);
          binReads = new uint64_t[binReadSize+1];
          memset(binReads, 0, sizeof(uint64_t)*(binReadSize+1));
          if (numSequences)
              blockOfNs.reserve(numSequences+10);
          return size;
      }

      void clearMapForNs()
      {
          blockOfNs.clear();
      }

      void clearTmpString()
      {
          strTmp.clear();
          strName.clear();
          clearMapForNs();
      }

      void getKmerLeftnRightBoundForNs(uint64_t &currKmerPos, mapObject &bounds)
      {
          uint64_t right=0;
          /*
           * Since we do all computation with bits, all our
           * positions are even. Here I return 1 (odd position),
           * an indication of no Ns towards left
           */

          if (!blockOfNs.size()){
              bounds.left=0x1;
              bounds.right=CHARS2BITS(totalBases-1);
              return;
          }

          vector<mapObject>::iterator it;
          it=upper_bound(blockOfNs.begin(), blockOfNs.end(), currKmerPos, mapObject());
          /* No N block beyond this point */
          if (it == blockOfNs.end())
              right = CHARS2BITS(totalBases-1);
          else
              right = (*it).left-2;

          /* This function never gets a position which is N */
          if (!currKmerPos || it==blockOfNs.begin()){
              bounds.left=0x1;
              bounds.right=right;
              return;
          }

          --it;

          bounds.left=(*it).right+2;
          bounds.right=right;
          return;
      }

      bool checkKmerForNs(uint64_t &currKmerPos, vector<mapObject>::iterator &it)
      {
          if (!blockOfNs.size())
              return false;

          while(it != blockOfNs.end())
          {
              if ((*it).left>currKmerPos)
                  break;
              else
                  ++it;
          }


          /* No N block beyond this point */
          if (it == blockOfNs.end()){
              --it;
              /* Current position within N block */
              if (((*it).left <=currKmerPos) && (currKmerPos <= (*it).right)){
                  ++it;
                  return true;
              }else{
                  ++it;
                  return false;
              }
          }

          if ((*it).left > (currKmerPos+commonData::kmerSize-2)){
              if (it != blockOfNs.begin()){
                  --it;
                  if ((*it).right < currKmerPos){
                      ++it;
                      return false;
                  }else {
                      ++it;
                      return true;
                  }
              }else
                  return false;
          }else {
              return true;
          }
      }


      void setCurrPos() {
          currPos+=size;;
      }

      uint64_t getCurrPos() {
          return currPos;
      }

      void resetCurrPos() {
          currPos=0;
      }

      bool readChunks()
      {
          string line;
          uint64_t blockNCount=0;
          int minSize = commonData::minMemLen/2-1;
          uint64_t sz=size+minSize;
          /* Process anything remaining from the last iteration */
          processTmpString(sz, blockNCount);

          if (numSeqFiles >= 1){
              while(getline( file1, line ).good()){
                  if(line[0] == '>' || (totalBases == sz)){
                      if( !strName.empty()){
                         if(line[0] != '>') {
                              processInput(line, sz, blockNCount);
                          }else {
                              processInput(randomStr(), sz, blockNCount);
                          }
                          if (totalBases == sz) {
                              if ((totalBases%32)!=0)
                              {
                                  uint64_t offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                                  binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                                  binReadsLocation++;
                                  binReads[binReadsLocation]=0;
                              }
                              if (blockNCount){
                                  blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                                  blockNCount=0;
                              }
                              return true;
                          }
                      }
                      if( !line.empty() ){
                          strName = line.substr(1);
                      }
                  } else if( !strName.empty() ){
                      processInput(line, sz, blockNCount);
                  }
              }
              if (numSeqFiles == 2) {
                  while(getline( file2, line ).good()){
                      if(line[0] == '>' || (totalBases == sz)){
                          if( !strName.empty()){
                             if(line[0] != '>') {
                                  processInput(line, sz, blockNCount);
                              }else {
                                  processInput(randomStr(), sz, blockNCount);
                              }
                              if (totalBases == sz) {
                                  if ((totalBases%32)!=0)
                                  {
                                      uint64_t offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                                      binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                                      binReadsLocation++;
                                      binReads[binReadsLocation]=0;
                                  }
                                  if (blockNCount){
                                      blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                                      blockNCount=0;
                                  }
                                  return true;
                              }
                          }
                          if( !line.empty() ){
                              strName = line.substr(1);
                          }
                      } else if( !strName.empty() ){
                          processInput(line, sz, blockNCount);
                      }
                  }
              }
          }
          if( !strName.empty() ){
              if ((totalBases%32)!=0)
              {
                  uint64_t offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                  binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                  binReadsLocation++;
                  binReads[binReadsLocation]=0;
              }
              if (blockNCount){
                  blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                  blockNCount=0;
              }
              if (!strTmp.size())
                  strName.clear();
              return true;
          }
          return false;
      }

      void flipCharacter(char &in, char &out)
      {
          switch(in)
          {
              case 'A':
              case 'a':
                  out='T';
                  break;
              case 'C':
              case 'c':
                  out='G';
                  break;
              case 'G':
              case 'g':
                  out='C';
                  break;
              case 'T':
              case 't':
                  out='A';
                  break;
              default:
                  out=in;
          }
      }

      void flipNswap(string &content)
      {
          string::iterator itBeg = content.begin();
          string::iterator itEnd = --content.end();
          char beg=0, end=0;
          uint64_t d=0;
          while ((d=distance(itBeg,itEnd)))
          {
              flipCharacter(*itBeg, end);
              flipCharacter(*itEnd, beg);
              (*itEnd)=end;
              (*itBeg)=beg;
              ++itBeg;
              --itEnd;
              if(d==1)
                  break;
          }
          if (!d)
              flipCharacter(*itEnd, *itEnd);
      }

      void writeReverseComplementString(string &name, string &content, fstream &file)
      {
          file << ">" << name << "\n";
          flipNswap(content);
          file << content ;
      }

      void generateSeqPos(vector<seqData> &vecSeqInfo) {
          seqData s;
          uint64_t i=0,j=0;
          string line;
          clearFileFlag();
          if (numSeqFiles >= 1){
              while(getline(file1, line).good() ){
                  if(line[0] == '>'){
                      if(!strName.empty()) {
                          s.start=CHARS2BITS(j);
                          s.end=CHARS2BITS(i-1);
                          s.seq.assign(strtok(const_cast<char *>(strName.c_str())," \t\n"));
                          vecSeqInfo.push_back(s);
                          s.seq.clear();
                          i+=RANDOM_SEQ_SIZE;
                          j=i;
                          strName.clear();
                      }
                      if(!line.empty())
                          strName=line.substr(1);
                  } else if( !strName.empty() ) {
                      i+=line.length();
                  }
              }
              if (numSeqFiles == 2){
                  while(getline(file2, line).good() ){
                      if(line[0] == '>'){
                          if(!strName.empty()) {
                              s.start=CHARS2BITS(j);
                              s.end=CHARS2BITS(i-1);
                              s.seq.assign(strtok(const_cast<char *>(strName.c_str())," \t\n"));
                              vecSeqInfo.push_back(s);
                              s.seq.clear();
                              i+=RANDOM_SEQ_SIZE;
                              j=i;
                              strName.clear();
                          }
                          if(!line.empty())
                              strName=line.substr(1);
                      } else if( !strName.empty() ) {
                          i+=line.length();
                      }
                  }
              }
          }
          if( !strName.empty() ) {
              i+=line.length();
              s.start=CHARS2BITS(j);
              s.end=CHARS2BITS(i-1);
              s.seq.assign(strtok(const_cast<char *>(strName.c_str())," \t\n"));
              vecSeqInfo.push_back(s);
              s.seq.clear();
              strName.clear();
          }
      }

      void generateRevComplement() {
          string line,content;
          fstream revFile;

          char buffer[256];
          memset(buffer,0,256);
          sprintf(buffer, "%s/revComp", commonData::nucmer_path);
          revFile.open(buffer, ios::out);
          if (!revFile.is_open())
          {
              cout << "ERROR: unable to open temporary reverse complement file" << endl;
              exit( EXIT_FAILURE );
          }

          clearFileFlag();

          if (numSeqFiles >= 1){
              while(getline( file1, line ).good()){
                  size += (line.length()+1);
                  if(line[0] == '>'){
                      if(!strName.empty()) {
                          numSequences++;
                          size += RANDOM_SEQ_SIZE;
                          writeReverseComplementString(strName, content, revFile);
                          content.clear();
                          strName.clear();
                      }
                      if(!line.empty())
                          strName=line.substr(1);
                  } else if( !strName.empty() ) {
                      content += "\n";
                      content += line;
                  }
              }
              // Catches the last sequence line
              if( !strName.empty() ) {
                  size += (line.length()+1);
                  content += "\n";
                  content += line;
                  numSequences++;
                  writeReverseComplementString(strName, content, revFile);
                  content.clear();
                  strName.clear();
              }
              if (numSeqFiles == 2){
                  while(getline( file2, line ).good()){
                      size += (line.length()+1);
                      if(line[0] == '>'){
                          if(!strName.empty()) {
                              numSequences++;
                              size += RANDOM_SEQ_SIZE;
                              writeReverseComplementString(strName, content, revFile);
                              content.clear();
                              strName.clear();
                          }
                          if(!line.empty())
                              strName=line.substr(1);
                      } else if( !strName.empty() ) {
                          content += "\n";
                          content += line;
                      }
                  }
                  // Catches the last sequence line
                  if( !strName.empty() ) {
                      size += (line.length()+1);
                      content += "\n";
                      content += line;
                      numSequences++;
                      writeReverseComplementString(strName, content, revFile);
                      content.clear();
                      strName.clear();
                  }
              }
          }
          revFile.close();
      }

      void convertToNucl(uint64_t bin, string &sequence) {
          if (bin == 0) {
              sequence += 'A';
          } else if (bin == 1) {
              sequence += 'C';
          } else if (bin == 2) {
              sequence += 'G';
          } else if (bin == 3) {
              sequence += 'T';
          }
      }

      void getFlankingSeq(string &sequence, uint64_t start, uint64_t end) {
          uint64_t currBin;
          uint32_t offset;

          for (uint64_t pos = (start / 2); pos != (end / 2); ++pos) {
            currBin = binReads[(pos * 2) / DATATYPE_WIDTH];
            offset = (pos * 2) % DATATYPE_WIDTH;
            currBin &= global_mask_right[(DATATYPE_WIDTH - offset) / 2 - 1];
            currBin >>= ((DATATYPE_WIDTH - 2) - offset);
            convertToNucl(currBin, sequence);
          }
      }
};

class outFileReadInfo {
    fstream outFile;

public:

    void writeString(string &content, fstream &file)
    {
        file << content << "\n";
    }

    void writeFiles(vector<seqData> &vecSeqInfo, vector<string> &filenames, string nonITRFilename, string ITRFilename) {
        seqData s;
        uint64_t lineNum = 0;
        string line;

        fstream palFile;
        palFile.open(ITRFilename, ios::out);
        if(!palFile.is_open()) {
            cout << "ERROR: unable to open " << ITRFilename << " file" << endl;
            exit( EXIT_FAILURE );
        }

        fstream outFile;
        outFile.open(nonITRFilename, ios::out);
        if(!outFile.is_open()) {
            cout << "ERROR: unable to open " << nonITRFilename << " file" << endl;
            exit( EXIT_FAILURE );
        }

        if (filenames.size() >= 1) {
            fstream file1;
            file1.open(filenames[0], ios::in);
            while (getline(file1, line).good()) {
                if (line[0] == '>') {
                    if (vecSeqInfo[lineNum].keep) {
                        line = strtok(const_cast<char *>(line.c_str()), " \t\n");
                        writeString(line, outFile);
                    } else {
                        writeString(vecSeqInfo[lineNum].header, palFile);
                    }
                } else {
                    if (vecSeqInfo[lineNum].keep)
                        writeString(line, outFile);
                    else
                        writeString(line, palFile);
                    ++lineNum;
                }
            }
            file1.close();
            if (filenames.size() == 2) {
                fstream file2;
                file2.open(filenames[0], ios::in);
                while (getline(file2, line).good()) {
                    if (line[0] == '>') {
                        if (vecSeqInfo[lineNum].keep) {
                            line = strtok(const_cast<char *>(line.c_str()), " \t\n");
                            writeString(line, outFile);
                        } else {
                            writeString(vecSeqInfo[lineNum].header, palFile);
                        }
                    } else {
                      if (vecSeqInfo[lineNum].keep)
                          writeString(line, outFile);
                      else
                          writeString(line, palFile);
                      ++lineNum;
                    }
                }
                file2.close();
            }
        }
        palFile.close();
        outFile.close();
    }
};

class SeqPos {
public:
    uint64_t lR;
    uint64_t rR;
    uint64_t lQ;
    uint64_t rQ;
    SeqPos() {
    }

    SeqPos(uint64_t lr, uint64_t rr, uint64_t lq, uint64_t rq)
    {
        lR=lr;
        rR=rr;
        lQ=lq;
        rQ=rq;
    }
};

class MemExt {
  public:
    uint64_t lR;
    uint64_t rR;
    uint64_t lQ;
    uint64_t rQ;
    MemExt() {
    }

    MemExt(uint64_t lr, uint64_t rr, uint64_t lq, uint64_t rq)
    {
         lR=lr;
         rR=rr;
         lQ=lq;
         rQ=rq;
    }

    bool operator () (const MemExt &obj1, const MemExt &obj2)
    {
      if (obj1.lQ<obj2.lQ)
          return true;
      else if (obj1.lQ>obj2.lQ)
          return false;
      else{
         if (obj1.lR<obj2.lR)
             return true;
         else
             return false;
      }
    }
};


class tmpFilesInfo {
    fstream *TmpFiles;
    fstream tabFile;
    vector <MemExt> MemExtVec;
    uint64_t numMemsInFile;

    bool checkMEMExt(uint64_t &lr, uint64_t &rr, uint64_t &lq, uint64_t &rq, seqFileReadInfo &QueryFile, seqFileReadInfo &RefFile) {
      if ((!lq && QueryFile.getCurrPos()) || rq == CHARS2BITS(QueryFile.totalBases-1)) {
         return true;
      }else if((!lr && RefFile.getCurrPos()) || rr == CHARS2BITS(RefFile.totalBases-1)) {
         return true;
      }
      return false;
    }

    void writeToFile(uint64_t lQ, uint64_t rQ, uint64_t lR, uint64_t rR) {
        MemExt m;
        m.lQ=lQ;
        m.lR=lR;
        m.rQ=rQ;
        m.rR=rR;
        TmpFiles[m.lQ/numMemsInFile].write((char *)&m, sizeof(MemExt));
    }

    void writeToVector(uint64_t lQ, uint64_t rQ, uint64_t lR, uint64_t rR) {
        MemExt m;
        m.lQ=lQ;
        m.lR=lR;
        m.rQ=rQ;
        m.rR=rR;
        MemExtVec.emplace_back(m);
    }


  public:

    tmpFilesInfo(int numFiles) {
        TmpFiles = new fstream[numFiles];
    }

    ~tmpFilesInfo() {
        delete [] TmpFiles;
    }

    void setNumMemsInFile(uint64_t size, uint64_t &numSequences) {
        numMemsInFile = ((2*(size*commonData::d+numSequences*RANDOM_SEQ_SIZE+commonData::d))/NUM_TMP_FILES);
    }

    static bool compare_reference (const MemExt &obj1, const MemExt &obj2)
    {
      if (obj1.lR<obj2.lR)
         return true;
      else if (obj1.lR>obj2.lR)
         return false;
      else{
         if (obj1.lQ<obj2.lQ)
             return true;
         else
             return false;
      }

    }

    static bool sortReverse (const MemExt &obj1, const MemExt &obj2)
    {
        if (obj1.lQ + obj1.lR > obj2.lQ + obj2.lR)
            return true;
        else if (obj1.lQ + obj1.lR < obj2.lQ + obj2.lR)
            return false;
        else{
            if (obj1.rQ + obj1.rR > obj2.rQ + obj2.rR)
                return true;
            else
                return false;
        }
    }

    static bool myUnique(const MemExt &obj1, const MemExt &obj2)
    {
      if((obj1.lQ==obj2.lQ) && (obj1.rQ==obj2.rQ) && (obj1.rR==obj2.rR) && (obj1.lR==obj2.lR))
          return true;
      else
          return false;
    }

    static bool uniqueIR(const posData &obj1, const posData &obj2)
    {
        if((obj1.L1Bound==obj2.L2Bound) && (obj1.R1Bound==obj2.R1Bound) && (obj1.L2Bound==obj2.L2Bound) && (obj1.R2Bound==obj2.R2Bound))
            return true;
        else
            return false;
    }

    static bool myUniqueQue(const MemExt &obj1, const MemExt &obj2)
    {
        if((obj1.lQ==obj2.lQ) && (obj1.rQ==obj2.rQ))
            return true;
        else
            return false;
    }

    void openFiles(ios_base::openmode mode, int numFiles, string tabFilename) {
        char buffer[256];
        memset(buffer,0,256);
        static int flag=0;

        sprintf(buffer, "%s", commonData::nucmer_path);
        if (!flag) {
            if(mkdir(buffer, S_IRWXU|S_IRGRP|S_IXGRP))
            {
                cout << "ERROR: unable to open temporary directory" << endl;
                exit( EXIT_FAILURE );
            }
            flag=1;
        }

        /* Last two files hold the sequence/pos mapping
         * for reference and query file respectively
         */
        for (int32_t i=0;i<numFiles;i++) {
            /* Temporary file to hold the mems */
            sprintf(buffer, "%s/%d", commonData::nucmer_path, i);
            TmpFiles[i].open(buffer, mode);
            if (!TmpFiles[i].is_open())
            {
                cout << "ERROR: unable to open temporary file" << endl;
                exit( EXIT_FAILURE );
            }
        }

        tabFile.open(tabFilename, ios::out);
        if(!tabFile.is_open()) {
            cout << "ERROR: unable to open " << tabFilename << " file" << endl;
            exit( EXIT_FAILURE );
        }
    }

    void closeFiles(int numFiles) {
        for (int32_t i=0;i<numFiles;i++){
           TmpFiles[i].close();
        }
    }

    fstream& getMapFile(int fIndex) {
        return TmpFiles[fIndex];
    }

    bool writeMemInTmpFiles(uint64_t &lRef, uint64_t &rRef, uint64_t &lQue, uint64_t &rQue, seqFileReadInfo &QueryFile, seqFileReadInfo &RefFile) {
       MemExt m;
       uint64_t currPosQ = CHARS2BITS(QueryFile.getCurrPos());
       uint64_t currPosR = CHARS2BITS(RefFile.getCurrPos());
       if (rRef-lRef+2 >= static_cast<uint64_t>(commonData::minMemLen)) {
           #pragma omp critical(writeFile)
           writeToFile(currPosQ+lQue, currPosQ+rQue, currPosR+lRef,  currPosR+rRef);
           return true;
       }else
           return false;
    }

    void mergeMemExtVector () {
        int flag=0;
        MemExt m;
        //if (commonData::d==1 && commonData::numThreads==1)
        //    return;

        if (MemExtVec.size() > 1) {
            sort(MemExtVec.begin(), MemExtVec.end(), MemExt());
            //Remove duplicates, this is happens due to removal of all checks.
            vector<MemExt>::iterator last=unique(MemExtVec.begin(), MemExtVec.end(), myUnique);
            do {
                flag=0;
                for (vector<MemExt>::iterator it=MemExtVec.begin(); it != (last-1); ++it) {
                    if (!(*it).lQ && !(*it).rQ && !(*it).lR && !(*it).rR )
                        continue;

                    vector<MemExt>::iterator dup = it;
                    ++dup;
                     for (; dup != last; ++dup) {
                        if (!(*dup).lQ && !(*dup).rQ && !(*dup).lR && !(*dup).rR )
                            continue;
                        if((*dup).lQ + static_cast<uint64_t>(commonData::minMemLen-2)-2 > (*it).rQ)
                            break;
                        if((*dup).lQ + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rQ) {
                            if((*dup).lR + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rR) {
                                flag=1;
                                (*it).rQ=(*dup).rQ;
                                (*it).rR=(*dup).rR;
                                (*dup).rQ=0;
                                (*dup).rR=0;
                                (*dup).lQ=0;
                                (*dup).lR=0;
                                break;
                            }
                        }
                    }
                    if (flag)
                        break;
                }
            } while (flag);

            last=unique(MemExtVec.begin(), last, myUnique);

            sort(MemExtVec.begin(), last, compare_reference);
            do {
                flag=0;
                for (vector<MemExt>::iterator it=MemExtVec.begin(); it != (last-1); ++it) {
                    if (!(*it).lQ && !(*it).rQ && !(*it).lR && !(*it).rR )
                        continue;
                    vector<MemExt>::iterator dup = it;
                    ++dup;
                    for (; dup != last; ++dup) {
                        if (!(*dup).lQ && !(*dup).rQ && !(*dup).lR && !(*dup).rR )
                            continue;
                        if((*dup).lR + static_cast<uint64_t>(commonData::minMemLen-2)-2 > (*it).rR)
                            break;
                        if((*dup).lR + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rR) {
                            if((*dup).lQ + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rQ) {
                                flag=1;
                                (*it).rQ=(*dup).rQ;
                                (*it).rR=(*dup).rR;
                                (*dup).rQ=0;
                                (*dup).rR=0;
                                (*dup).lQ=0;
                                (*dup).lR=0;
                                break;
                            }
                        }
                    }
                    if (flag)
                        break;
                }
            } while (flag);

            for (vector<MemExt>::iterator it=MemExtVec.begin(); it != last; ++it) {
                if ((*it).lQ || (*it).rQ || (*it).lR || (*it).rR ) {
                    writeToFile((*it).lQ, (*it).rQ, (*it).lR, (*it).rR);
                }
            }
        }
        MemExtVec.clear();
    }

    void extBin(seqFileReadInfo &RefFile, vector<uint64_t> &Bins, int &binLength, uint64_t &extLength, uint64_t &Bound) {
        int readBin=0, bin=0;
        uint64_t offset, totalOffset=0;

        Bins.resize(binLength);
        while (bin != binLength) {
            offset = (Bound) % DATATYPE_WIDTH;
            if (totalOffset % DATATYPE_WIDTH) {
                offset = DATATYPE_WIDTH - offset;
                Bins[bin] |= ((RefFile.binReads[(Bound) / DATATYPE_WIDTH + readBin] & global_mask_left[(totalOffset % DATATYPE_WIDTH)/2 - 1]) >> offset);
                ++bin;
            } else {
                Bins[bin] = RefFile.binReads[(Bound) / DATATYPE_WIDTH + readBin];
                Bins[bin] <<= offset;
                ++readBin;
                if (!offset)
                    ++bin;
            }
            totalOffset += offset;
        }
        if (extLength % DATATYPE_WIDTH)
            Bins[binLength-1] = Bins[binLength-1] & global_mask_left[(extLength % DATATYPE_WIDTH)/2 - 1];
    }

    void convertToNucl(uint64_t bin, string &sequence) {
        if (bin == 0) {
            sequence += 'A';
        } else if (bin == 1) {
            sequence += 'C';
        } else if (bin == 2) {
            sequence += 'G';
        } else if (bin == 3) {
            sequence += 'T';
        }
    }

    void getInvertedRepeats(uint64_t &lQue, uint64_t &rQue, seqFileReadInfo &QueryFile, uint64_t &rRef, uint64_t &lRef, seqFileReadInfo &RefFile, vector<seqData> &vecSeqInfo) {
        mapObject QueryNpos;
        seqData q, r;
        vector<seqData>::iterator seqitQ, seqitR;
        string currHeader;

        q.start=rQue;
        q.end=lQue;
        seqitQ = lower_bound(vecSeqInfo.begin(), vecSeqInfo.end(), q, seqData());

        r.start = rRef;
        r.end = lRef;
        seqitR = lower_bound(vecSeqInfo.begin(), vecSeqInfo.end(), r, seqData());

        tabFile << (*seqitQ).seq << "\t" << (*seqitR).seq << "\n";

        if ((*seqitQ).keep) {
            QueryFile.getKmerLeftnRightBoundForNs(lQue, QueryNpos);
            (*seqitQ).header = '>' + (*seqitQ).seq + "_LCoord_" + to_string(((lQue - (QueryNpos.left==1?QueryNpos.left=0:QueryNpos.left)) + 2)/2) + "_RCoord_" + to_string((rQue - (QueryNpos.left==1?QueryNpos.left=0:QueryNpos.left) + 2)/2);
            (*seqitQ).keep = 0;
        }
    }

    void removeTmp() {
        char buffer[256];
        memset(buffer,0,256);

        sprintf(buffer, "%s/revComp", commonData::nucmer_path);
        remove(buffer);

        int numFiles=NUM_TMP_FILES;
        for (int32_t i=0;i<numFiles;i++) {
            sprintf(buffer, "%s/%d", commonData::nucmer_path, i);
            remove(buffer);
        }

        sprintf(buffer, "%s", commonData::nucmer_path);
        remove(buffer);

    }

    void closeFile() {
        tabFile.close();
    }
};
