#ifndef _SOURCE_H
#define _SOURCE_H



#include "main.h"
#include "hyperparams.h"


using namespace std;



class Source
{
  string name;
  map<string, vector<double> *> frequencies;
  vector<map<string, vector<vector<double> > > > Pr_XXS;
  vector<map<string, vector<vector<double> > > > Pr_XOS;
  vector<map<string, vector<vector<double> > > > Pr_OS;
  vector<map<string, vector<vector<double> > > > Pr_XXS_copy;
  vector<map<string, vector<vector<double> > > > Pr_XOS_copy;
  vector<map<string, vector<vector<double> > > > Pr_OS_copy;

  vector<double> Pr_vec;
  vector<double> Pr_vec_copy;

  std::vector<std::string> locus_order;
  std::map<string, int> alleleCount;
  int k_stride;

  map<string, double> settings;

public:
  map<string, vector<double> *>* getLociFreqs(){return &frequencies;}
  vector<double>* getLocusFreqs(string locusIn){return frequencies.at(locusIn);}
  void setLocusFreqs(string locus, vector<double> * freqs);
  string getName(){return name;}
  // vector<map<string, vector<vector<double> > > > get_Pr_XXS() const {return Pr_XXS;}
  // vector<map<string, vector<vector<double> > > > get_Pr_XOS() const {return Pr_XOS;}
  // vector<map<string, vector<vector<double> > > > get_Pr_OS() const {return Pr_OS;}
  const vector<double>& get_Pr_vec() const {return Pr_vec;}

  void setName(string name){this->name = name;}
  void setPrVecs(map<string, vector<double> >*, vector<map<string, vector<double> > >*);
  void setPrVecs(const map<string, vector<double> > &, const vector<map<string, vector<double> > > &, const std::string &);
  void setLocusOrder(std::vector<std::string> locus_order){ this->locus_order = locus_order;}
  void setAlleleCount(std::map<string, int> alleleCount){ this->alleleCount = alleleCount;}
  void setKStride(int k_stride){this->k_stride = k_stride;}
  void initPrVecs();
  void copyProbabilities();
  void revertProbabilities();
  Source(map<string,double>& settingsIn){settings = settingsIn;}
}; // class Source



#endif // _SOURCE_H
