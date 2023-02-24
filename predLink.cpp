#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <unordered_map>

#define PACK_DENSITY 4

// 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions
#define MASK0   3 // 3 << 2 * 0
#define MASK1  12 // 3 << 2 * 1
#define MASK2  48 // 3 << 2 * 2
#define MASK3 192 // 3 << 2 * 3

using namespace std;

unordered_map<string, string> sc_alleles;
unordered_map<string, float>  sc_bhat;
unordered_map<string, bool>   sc_snplist;
unordered_map<string, int>    marker_set;
unordered_map <string, float*> geno_hash;

void decode_plink(char *output, const char *input, const int lengthInput){
  int i, k;
  char tmp, geno;
  int a1, a2;
  
  for(i=0;i<lengthInput;++i){
    tmp = input[i];
    k   = PACK_DENSITY * i;
    
    geno      = (tmp & MASK0);
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK1) >> 2; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK2) >> 4; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    
    geno      = (tmp & MASK3) >> 6; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
  }
}

int readGenotypes(string *IID,
                  float *posFoundSNPs, float *posAllSNPs,
                  float *dpqAll,float *dpqFound,
                  float *bhat, float *bhat_sq, bool *found, string bedfile, int n, int p,
                  char *packed, char *unpacked, int numBytes){
  ifstream influx;
  influx.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if(!influx){
    cerr << "[readGenotypes] Error reading file "<<bedfile<<endl;
  }
  influx.seekg(0, ifstream::end);
  influx.seekg(3, ifstream::beg);
  
  int     i,j;
  int    nEff;
  float x, mx;
  //float dpq_x;
  
  int k = -1;
  for(j=0;j<p;j++){
    influx.read((char*)packed, sizeof(char) * numBytes);
    decode_plink(unpacked, packed, numBytes);
    if(found[j]){
      nEff = 0;
      mx   = 0.;
      for(i=0;i<n;i++){
        x = (float) ((int) unpacked[i]);
        if(x!=3.){
          nEff++;
          mx += x;
        }
      }
      mx    = mx/nEff;
      //dpq_x = mx*(1.-mx*0.5);
      
      // cout<<dpqAll[j]<<" "<<dpq_x<<endl; // OK
      
      k++;
      posFoundSNPs[k] = posAllSNPs[j];
      dpqFound[k]     = dpqAll[j];
      bhat_sq[k]      = bhat[j] * bhat[j];
      
      for(i=0;i<n;i++){
        x = (float) ((int) unpacked[i]);
        if(x!=3.0){
          geno_hash[IID[i]][k] = (x-mx);
        }else{
          geno_hash[IID[i]][k] = 0.; // imputed to the mean
        }
      }
      
    }
  }
  influx.close();
  int nsnps = k + 1;
  return(nsnps);
}

float round_up(float value, int decimal_places) {
  const float multiplier = pow(10.0, decimal_places);
  return ceil(value * multiplier) / multiplier;
}

// Main funtion used on cov_files.txt with a sequential access
int main(int argc, char *argv[]){
  
  unordered_map<string, string> sc_alleles;
  unordered_map<string, float>  sc_bhat;
  unordered_map<string, float>  sc_freq;
  unordered_map<string, bool>   sc_snplist;
  unordered_map<string, int>    marker_set;

  // Input arguments
  string ibdFile      = "none";
  string outputprefix = "none";
  string Plinkprefix  = "none";
  string sumstatfile  = "none";
  string snplistfile  = "none";
  string chromosome   = "none";
  bool verbose        =   true;

  // Indices
  string sw;
  int i,j,k,l;
  
  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  int isnp           = -1;
  int ibeta          = -1;
  int iallele        = -1;
  int ifreq          = -1;
  
  float dcM         = 0.0;
  int stepMarkers   =   1;
  
  // Read arguments
  sw = argv[1];
  if (sw == "--help"){
    cerr<<"\t--ibd           : IBD probability output from Merlin."<<endl;
    cerr<<"\t--bfile         : binary PLINK format for genotypes."<<endl;
    cerr<<"\t--step          : Specifies step between markers to make prediction. Default is 1 (minimum is 1)."<<endl;
    cerr<<"\t--sumstat       : GWAS summary statistics file. E.g. --sumstat myfile col1[SNP] col2[Allele] col3[BETA] col4[Freq]."<<endl;
    cerr<<"\t--chr           : Specifies which chromosome to analyse (e.g. 1, 22, X, etc.)."<<endl;
    cerr<<"\t--out           : A prefix for the output file: [output-prefix].prdlink"<<endl;
    cerr<<"\t--silent        : Specify whether the different steps of the calculations should be displayed on screen."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  for(i = 1; i<argc;i++){
    sw = argv[i];
    if (sw == "--ibd"){
      ibdFile = argv[i + 1];
    }
    if (sw == "--step"){
      stepMarkers = atoi(argv[i + 1]);
    }
    if (sw == "--bfile"){
      Plinkprefix = argv[i + 1];
    }
    if (sw == "--sumstat"){
      sumstatfile = argv[i + 1];
      isnp        = atoi(argv[i + 2])-1;
      iallele     = atoi(argv[i + 3])-1;
      ibeta       = atoi(argv[i + 4])-1;
      ifreq       = atoi(argv[i + 5])-1;
    }
    if (sw == "--chr"){
      chromosome = argv[i + 1];
    }
    if (sw == "--out"){
      outputprefix = argv[i + 1];
    }
    if (sw == "--extract"){
      snplistfile = argv[i + 1];
    }
    if (sw == "--silent"){
      verbose = false;
    }
  }
  
  // Check input files
  if(ibdFile=="none"){
    cerr<<"\tIBD file is not specified. Use --ibd and check help [--help]."<<endl;
    exit(1);
  }
  if(stepMarkers<1){
    cerr<<"\tAn incorrect number of steps between markers has been specified. Number should be at least 1. Use --step and check help [--help]."<<endl;
    exit(1);
  }
  
  if(Plinkprefix=="none"){
    cerr<<"\tA prefix must be specified for files [prefix].bed, [prefix].bim and [prefix].fam. Use --bfile and check help [--help]."<<endl;
    exit(1);
  }
  if(sumstatfile=="none"){
    cerr<<"\tGWAS summary statistics file is not specified. Use --sumstat and check help [--help]."<<endl;
    exit(1);
  }
  if(isnp<0){
    cerr<<"\tSpecify positive column number for SNP ID. Use --sumstat and check help [--help]."<<endl;
    exit(1);
  }
  if(ibeta<0){
    cerr<<"\tSpecify positive column number for Effect size (BETA). Use --sumstat and check help [--help]"<<endl;
    exit(1);
  }
  if(iallele<0){
    cerr<<"\tSpecify positive column number for Reference allele. Use --sumstat and check help [--help]"<<endl;
    exit(1);
  }
  if(ifreq<0){
    cerr<<"\tSpecify positive column number for Allele frequency. Use --sumstat and check help [--help]"<<endl;
    exit(1);
  }
  if(chromosome=="none"){
    cerr<<"\tSpecify chromsome ID. Use --chr and check help [--help]"<<endl;
    exit(1);
  }
  
  string bedfile = Plinkprefix+".bed";
  string bimfile = Plinkprefix+".bim";
  string famfile = Plinkprefix+".fam";
  
  clock_t begin = clock();
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  cout <<"# Analysis starts: ";
  cout << (now->tm_year + 1900) << '-' 
       << (now->tm_mon + 1) << '-'
       <<  now->tm_mday << " at "
       <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
       << endl;
  
  if(verbose){
    cout<<"# INPUT files..."<<endl;
    cout<<"# BED: "<<Plinkprefix<<".bed.\n";
    cout<<"# BIM: "<<Plinkprefix<<".bim.\n";
    cout<<"# FAM: "<<Plinkprefix<<".fam.\n";
    cout<<"# IBD: "<<ibdFile<<".\n";
    cout<<"# STEP: "<<stepMarkers<<".\n";
    cout<<"# GWAS SUMMARY STATISTICS : "<<sumstatfile<<".\n";
    if(snplistfile!="none"){
      cout<<"# SNPs to extract: "<<snplistfile<<".\n";
    }
    cout<<"#\n";
  }

  // Few tools
  string line = "";
  string tok  = ""; 
  ifstream tmpStream;
  
  // Get number of SNPs
  if(verbose){
    cout<<"# [1] Counting the number of SNPs in genotype file...";
  }  
  int p = -1;
  tmpStream.open(bimfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    p++;
  }
  tmpStream.close();
  if(verbose){
    cout<<"[done]: found "<<p<<" SNPs."<<endl;
  }
  
  // Get sample size
  if(verbose){
    cout<<"# [2] Counting the number of samples in genotype file...";
  }    
  int n = -1; 
  tmpStream.open(famfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    n++;
  }
  tmpStream.close();  
  if(verbose){
    cout<<"[done]: found "<<n<<" samples."<<endl;
  }

  // Read ID
  if(verbose){
    cout<<"# [3] Reading individual IDs in genotype file...";
  }     
  int nColFamFile = 6;
  string *FID = new string[n];
  string *IID = new string[n];
  tmpStream.open(famfile.c_str());
  for(i=0;i<n;i++){
    for(k=0;k<nColFamFile;k++){
      tmpStream >> tok;
      if(k==0){
        FID[i] = tok;
      }
      if(k==1){
        IID[i] = tok;
      }    
    }
  }
  tmpStream.close();
  if(verbose){
    cout<<"[done]"<<endl;
  }
  
  
  if(verbose){
    cout<<"# [4] Reading "<<bimfile<<"...";
  }   
  int nColBimFile = 6;
  int dp = 2*p;
  string *snps      = new string[p];
  string *chrs      = new string[p];
  string *a1a2      = new string[dp];
  float *posAllSNPs = new float[p];
  
  tmpStream.open(bimfile.c_str());
  for(j=0;j<p;j++){
    for(k=0;k<nColBimFile;k++){
      tmpStream >> tok;
      if(k==0){
        chrs[j] = tok;
      }
      if(k==1){
        snps[j] = tok;
      }
      if(k==2){
        posAllSNPs[j] = atof(tok.c_str());
      }
      if(k==4){
        a1a2[j] = tok;
      }
      if(k==5){
        a1a2[j+p] = tok;
      }
    }
  }
  tmpStream.close();
  if(verbose){
    cout<<"[done]."<<endl;
  }
  
  // Reading the bhat
  if(verbose){
    cout<<"# [5] Reading GWAS summary statistics file: ";
  }
  int nSNPscore = -1;
  tmpStream.open(sumstatfile.c_str());
  getline(tmpStream,line); // header
  while(tmpStream){
    getline(tmpStream,line);
    nSNPscore++;
  }
  tmpStream.close();
  if(verbose){
    cout<<": found "<<nSNPscore<<" SNPs and ";
  }

  tmpStream.open(sumstatfile.c_str());
  string snp, allele, weight, freq;
  int nColSumstatFile = 0;
  getline(tmpStream,line);
  stringstream ss;
  ss << line;
  while( ss >> tok ){
    nColSumstatFile++;
  }
  if(verbose){
    cout<<nColSumstatFile<<" columns...";
  }
  
  for(j=0;j<nSNPscore;j++){
    for(k=0;k<nColSumstatFile;k++){
      tmpStream >> tok;
      if(k==isnp){
        snp = tok;
      }
      if(k==iallele){
        allele = tok;
      }
      if(k==ibeta){
        weight = tok;
      }
      if(k==ifreq){
        freq = tok;
      }
    }
    
    // check upper case
    if(allele=="a") allele="A";
    if(allele=="c") allele="C";
    if(allele=="g") allele="G";
    if(allele=="t") allele="T";
    
    sc_alleles.insert( { snp, allele });
    sc_bhat.insert( { snp, atof(weight.c_str()) });
    sc_freq.insert( { snp, atof(freq.c_str()) });
  }
  tmpStream.close();
  if(verbose){
    cout<<"[done]."<<endl;
  }

  if(snplistfile != "none"){
    if(verbose){
      cout<<"# *** Reading SNP list...";
    }
    int nSNPlist = -1;
    tmpStream.open(snplistfile.c_str());
    while(tmpStream){
      getline(tmpStream,line);
      nSNPlist++;
    }
    tmpStream.close();
    // Reading SNP list
    tmpStream.open(snplistfile.c_str());
    for(j=0;j<nSNPlist;j++){
      tmpStream >> snp;
      sc_snplist.insert( { snp, true });
    }
    tmpStream.close();
    if(verbose){
      cout<<"found "<<nSNPlist<<" SNPs. [done]."<<endl;
    }
  }

  if(verbose){
    cout<<"# [6] Checking overlap between .bim file and GWAS summary statistics file (and SNP list)...";
  }
  float *bhat   = new float[p];
  bool *found   = new bool[p];
  float *dpqAll = new float[p];
  for(j=0;j<p;j++){
    bhat[j]    = 0.0;
    dpqAll[j]  = 0.0;
    found[j]   = false;
  }
  int nfound = 0;
  if(snplistfile == "none"){
    for(j=0;j<p;j++){
      found[j] = ( (sc_alleles[snps[j]]==a1a2[j] or sc_alleles[snps[j]]==a1a2[j+p]) and chrs[j]==chromosome );
      if(found[j]){
        dpqAll[j] = 2. * sc_freq[snps[j]] * (1. - sc_freq[snps[j]]);
        if(sc_alleles[snps[j]]==a1a2[j]){
          bhat[j] = +sc_bhat[snps[j]];
        }else{
          bhat[j] = -sc_bhat[snps[j]];
        }
        nfound += found[j];
      }
    }
  }else{
    for(j=0;j<p;j++){
      found[j] = ( (sc_alleles[snps[j]]==a1a2[j] or sc_alleles[snps[j]]==a1a2[j+p]) and sc_snplist[snps[j]] and chrs[j]==chromosome );
      if(found[j]){
        dpqAll[j] = 2. * sc_freq[snps[j]] * (1. - sc_freq[snps[j]]);
        if(sc_alleles[snps[j]]==a1a2[j]){
          bhat[j] = +sc_bhat[snps[j]];
        }else{
          bhat[j] = -sc_bhat[snps[j]];
        }
        nfound += found[j];
      }
    }
  }

  if(verbose){
    cout<<"[done].\n#     \033[1;31mFound "<<nfound<<" SNPs in common with named alleles\033[0m.\n";
  }

  if(nfound==0){
    cerr<<"\n\t\033[1;31m**** No SNPs left after filter! ****\033[0m\n"<<endl;
    exit(1);
  }
  
  if(verbose){
    cout<<"# [7] Reading genotypes...";
  }
  for(i=0;i<n;i++){
    geno_hash[IID[i]] = new float[nfound];
  }
  
  float *posFoundSNPs = new float[nfound];
  float *dpqFound     = new float[nfound];
  float *bhat_sq = new float[nfound];
  
  int numBytes   = (int)ceil((float)n / PACK_DENSITY);
  char* packed   = new char[numBytes];
  char* unpacked = new char[numBytes * PACK_DENSITY];
  int     M      = readGenotypes(IID,
                                 posFoundSNPs,posAllSNPs,
                                 dpqAll,dpqFound,
                                 bhat,bhat_sq,found,bedfile,n,p,packed,unpacked,numBytes);
    
  if(verbose){
    cout<<"[done]. \033[1;31mIncluded "<<M<<" SNPs in analysis\033[0m.\n";
  }
  if(verbose){
    cout<<"# [8] Reading IBD file (first pass): ";
  }
  
  tmpStream.open(ibdFile.c_str());
  getline(tmpStream,line); // reading header
  
  int ncolIBDFile = 0;
  stringstream ss_ibd_header;
  ss_ibd_header << line;
  while(ss_ibd_header>>tok){
    ncolIBDFile++;
  }
  
  string marker_str;
  while(tmpStream and dcM==0.0){
    getline(tmpStream,line);
    stringstream ss_ibd_line_0;
    ss_ibd_line_0 << line;
    for(k=0;k<ncolIBDFile;k++){
      ss_ibd_line_0 >> marker_str;
      if(k==3 and line!=""){
        dcM  = atof(marker_str.c_str());
      }
    }
  }
  if(verbose){
    cout<<"\033[1;31mdetected step-size = "<<dcM<<" cM\033[0m.\n";
  }
  tmpStream.close();
  
  tmpStream.open(ibdFile.c_str());
  getline(tmpStream,line); // reading header
  float fl_marker   = 0.;
  int in_marker    =  0;
  int nLineIBD     = -1;

  while(tmpStream){
    getline(tmpStream,line);
    stringstream ss_ibd_line_1;
    ss_ibd_line_1 << line;
    for(k=0;k<ncolIBDFile;k++){
      ss_ibd_line_1 >> marker_str;
      if(k==3 and line!=""){
        fl_marker  = atof(marker_str.c_str()) / dcM;
        in_marker  = (int) fl_marker;
        marker_str = to_string(in_marker);
        if(marker_set.find(marker_str) == marker_set.end()){
          marker_set.insert({marker_str,1});
        }else{
          marker_set[marker_str]++;
        }
      }
    }
    nLineIBD++;
  }
  tmpStream.close();
  int nMarkers = marker_set.size();
  int npairs   = marker_set["0"];
  // Check if number of pairs is constant for all markers
  int nOK = 0;
  for(i=0;i<nMarkers;i++){
    int npi = marker_set[to_string(i)];
    nOK    += npi==npairs;
  }
  
  if(verbose){
    cout<<"#     Found "<<nLineIBD<<" lines (exluding header) and "<<ncolIBDFile<<" columns.\n";
    cout<<"#     Number of markers is "<<nMarkers<<" and number of pairs is "<<npairs<<"."<<endl;
    cout<<"#     \033[1;31mNumber of markers with "<<npairs<<" pairs = "<<nOK<<"/"<<nMarkers<<"\033[0m."<<endl;
  }

  if(verbose){
    cout<<"# [9] Reading IBD file (second pass)...";
  }
  
  float p1  = 0.;
  float p2  = 0.;

  string fid_str, id1_str, id2_str, mrk_str, p0_str, p1_str, p2_str;
  
  float *Pi       = new float[nLineIBD];
  bool *foundPair = new    bool[npairs];
  string *ID1     = new  string[npairs];
  string *ID2     = new  string[npairs];
  
  // initialise
  int pair_index = 0;
  for(pair_index=0;pair_index<npairs;pair_index++){
    ID1[pair_index]       = "none";
    ID2[pair_index]       = "none";
    foundPair[pair_index] = false;
  }

  int ibdfile_line_index = 0;
  int marker_index = 0;
  
  // re-open ibd file
  tmpStream.open(ibdFile.c_str());
  getline(tmpStream,line); // reading header
  
  // while(tmpStream){
  for(ibdfile_line_index=0;ibdfile_line_index<nLineIBD;ibdfile_line_index++){
    getline(tmpStream,line);
    if(line!=""){
      stringstream ss_ibd_line_2;
      ss_ibd_line_2 << line;
      
      ss_ibd_line_2 >> fid_str;
      ss_ibd_line_2 >> id1_str;
      ss_ibd_line_2 >> id2_str;
      ss_ibd_line_2 >> mrk_str;
      ss_ibd_line_2 >> p0_str;
      ss_ibd_line_2 >> p1_str;
      ss_ibd_line_2 >> p2_str;

      fl_marker     = atof(mrk_str.c_str()) / dcM;
      marker_index  = (int) fl_marker;
      p1            = atof(p1_str.c_str());
      p2            = atof(p2_str.c_str());
      
      // get pair_id
      // int remainder   = ibdfile_line_index % npairs;
      pair_index      = ibdfile_line_index / npairs;
      ID1[pair_index] = id1_str;
      ID2[pair_index] = id2_str;
      
      // check if present
      if( (geno_hash.find(id1_str) != geno_hash.end()) and (geno_hash.find(id2_str) != geno_hash.end()) ){
        foundPair[pair_index] = true;
      }else{
        foundPair[pair_index] = false;
      }
      
      // Calculate probabilities
      if(pair_index==npairs){
        
      }
      Pi[pair_index + npairs*marker_index] = 0.5*p1 + p2;
    }
  }
  tmpStream.close();
  
  int npairsEff = 0;
  for(pair_index=0;pair_index<npairs;pair_index++){
    npairsEff = npairsEff + foundPair[pair_index];
  }
  if(verbose){
    cout<<"[done].\n";
    cout<<"#     First cell filled is Pi[0,0] = "<<Pi[0]<<".\n";
    cout<<"#     Last cell filled is Pi["<<pair_index<<","<<marker_index<<"] = "<<Pi[pair_index + npairs*marker_index]<<".\n";
    cout<<"#     \033[1;31mNumber of pairs included in the analysis is "<<npairsEff<<"\033[0m.\n";
    cout<<"#     Examples of pairs: ("<<ID1[0]<<","<<ID2[0]<<"), ("<<ID1[1]<<","<<ID2[1]<<"), ("<<ID1[2]<<","<<ID2[2]<<"), etc.\n";
  }

  if(verbose){
    cout<<"# [9] Calculating co-inheritance probabilities. ";
  }

  // Prepare results files
  int size_prob = nMarkers * nfound; // w[k+nMarkers*j]
  float *probCoInheritance = new float[size_prob];
  float posMarker;
  float distance;
  int nLargeProb = 0;
  for(k=0;k<nMarkers;k++){
    posMarker  = ((float) k) * dcM;
    for(j=0;j<nfound;j++){
      distance = fabs(posMarker-posFoundSNPs[j]) / 100.; // in Morgan
      probCoInheritance[k + nMarkers*j] = exp(-4.0*distance); // -- change here on May 30th, 2020.
      if(probCoInheritance[k + nMarkers*j]>0.5){
        nLargeProb++;
      }
    }
  }
  if(verbose){
    cout<<"[done].\n";
    cout<<"#     Found "<<nLargeProb<<" co-inheritance probabilities larger than 0.5.\n";
  }

  // Check covariance of genotypes between pairs
  int np_j;
  float corr_x;
  float min_corr = +1.0;
  float max_corr = -1.0;
  float avg_corr =  0.0;
  for(j=0;j<nfound;j++){
    corr_x = 0.0;
    np_j   = 0;
    for(pair_index=0;pair_index<npairs;pair_index++){
      if(foundPair[pair_index]){
        corr_x += geno_hash[ID1[pair_index]][j] * geno_hash[ID2[pair_index]][j];
        np_j++;
      }
    }
    if(np_j != npairsEff){
      cerr<<"**** Error np_j = "<<np_j<<" while expecting = "<<npairsEff<<"\n";
      exit(1);
    }
    corr_x = corr_x / np_j;
    corr_x = corr_x / dpqFound[j];
    if(corr_x<min_corr) min_corr = corr_x;
    if(corr_x>max_corr) max_corr = corr_x;
    avg_corr += corr_x;
  }
  avg_corr = avg_corr / nfound;
  if(verbose){
    cout<<"#     Mean correlation of genotypes is "<<avg_corr<<" (min: "<<min_corr<<" - max: "<<max_corr<<").\n";
  }

  float *prdLink = new float[nMarkers];
  float *pL      = new float[nMarkers]; // old approach
  if(verbose){
    cout<<"# [X] Writing output file...\n";
  }
  string outfile = outputprefix+".prdlink.out";
  ofstream fileOut(outfile.c_str());
  fileOut<<"MARKER PRDLINK_QUICK PRDLINK_IBD"<<endl;
  for(k=0;k<nMarkers;k=k+stepMarkers){  
    if(verbose){
      // cout<<"#     \033[1;31mProcessing marker ["<<k<<"/"<<nMarkers-1<<"]\033[0m.\n";
    }
    pL[k]      = 0.0;
    prdLink[k] = 0.0;
    posMarker  = ((float) k) * dcM;
     
    for(j=0;j<nfound;j++){
      pL[k] += dpqFound[j] * bhat_sq[j] * probCoInheritance[k + nMarkers*j]; // note that bhat_sq[j] is already squared
      np_j = 0;
      for(l=0;l<npairs;l++){
        if(foundPair[l]){
          prdLink[k] += 2.0 * geno_hash[ID1[l]][j] * geno_hash[ID2[l]][j] * bhat_sq[j] * Pi[l + npairs*k] * probCoInheritance[k + nMarkers*j];
          np_j++;
        }
      }
      if(np_j != npairsEff){
        cerr<<"**** [Error] np_j = "<<np_j<<" while expecting = "<<npairsEff<<"\n";
        exit(1);
      }
    }
    prdLink[k] = 2.0 * prdLink[k] / npairsEff;
    fileOut<<posMarker<<" "<<pL[k]<<" "<<prdLink[k]<<endl;
  }
  fileOut.close();
  
  time_t t2 = time(0);   // get time now
  struct tm * now2 = localtime( & t2 );
  cout <<"# Analysis ends: ";
  cout << (now2->tm_year + 1900) << '-' 
       << (now2->tm_mon + 1) << '-'
       <<  now2->tm_mday << " at "
       <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
       << endl;
  clock_t end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;  
  if(verbose) cout<<"# Analysis took "<<elapsed_secs<<" second(s).\n";
  
  delete [] packed;
  delete [] unpacked;    
  delete [] FID;
  delete [] IID;
  delete [] bhat;
  delete [] snps;
  delete [] chrs;
  delete [] a1a2;
  delete [] found;
  delete [] dpqAll;
  delete [] dpqFound;
  delete [] posAllSNPs;
  delete [] posFoundSNPs;
  delete [] bhat_sq;
  
  delete [] ID1;
  delete [] ID2;
  delete [] Pi;
  delete [] probCoInheritance;
  
  delete [] prdLink;
  delete [] pL;
  delete [] foundPair;
  
  for(i=0;i<n;i++){
    // delete [] geno_hash[IID[i]];
  }

  return EXIT_SUCCESS;
}
