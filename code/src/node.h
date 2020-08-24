#ifndef _NODE_H
#define _NODE_H



#include "source.h"
#include "main.h"
#include "hyperparams.h"
#include "util.h"
#include <assert.h>
#include <cmath>
#include <string>

// extern double GENDIST_weight;



using namespace std;



class Node
{
  struct classcomp {
    bool operator() (const Node* lhs, const Node* rhs) const
    {return lhs->getCode() < rhs->getCode();}
  };

  vector<pair<Node *,int> > ancestors;
  int code;
  map<Node*,double,classcomp> edgeProbs;
  map<Node*,int,classcomp> genDists;
  map<Node*,int,classcomp> timeDists;
  map<Node*, std::tuple<double, double>,classcomp > spaceDists;
  map<string,vector<bool> *> loci;
  map<string,bool> lociData;
  double logLikGen;
  double logLikTime;
  double logLikSpace;
  double logLikTravelHist;
  vector<pair<Source *,int> > sources;
  int time;
  double lat;
  double lon;
  int history;
  string type;

  std::vector<std::string> locus_order;
  std::map<string, int> alleleCount;
  std::map<std::tuple<int, int, std::string>, double> locusCache;
  std::map<std::tuple<int, std::string, std::string>, double> srcLocusCache;
  int k_stride;

  map<string, double> settings;

  void calcLogLikGen(const vector<double> &, const Source &, const std::map<std::string, bool> &);
  void calcLogLikTime();
  void calcLogLikSpace(double);
  void calcLogLikTravelHist(double, double, double);

public:
  Node(){;}
  void assignAncestor(Node * nodeIn, int kIn){
    ancestors.push_back(pair<Node *,int>(nodeIn, kIn));}
  void assignCode(int codeIn){code = codeIn;}
  void assignSettings(map<string,double>& settingsIn){settings = settingsIn;}
  void assignSource(Source * sourceIn, int kIn){
    sources.push_back(pair<Source *,int>(sourceIn, kIn));}
  void calcEdgeProbs(double);
  int calcGenDist(Node *);
  void calcLogLik(const std::vector<double> &, const Source &, const std::map<std::string, bool> &, double, double, double, double);
  void calcPr_notcontrib_allele(double, double, double, double, int);
  int calcTimeDist(Node *);
  std::tuple<double,double> calcSpaceDist(Node *);
  double calcSpaceSigma(Node *, double, int);
  void clearCache();
  void clearLocusCache();
  void clearSrcLocusCache();
  void checkForNAs();
  void copyNode(Node *);
  bool empty();
  vector<pair<Node *,int> >* getAncestors(){return &ancestors;}
  int getCode() const {return code;}
  map<Node*,double,Node::classcomp>* getEdgeProbs(){return &edgeProbs;}
  map<Node*,int, Node::classcomp>* getGenDists(){return &genDists;}
  int getHistory(){return history;}
  map<string,bool>* getLociData(){return &lociData;}
  double getLogLik(){return (logLikGen + logLikTime + logLikSpace + logLikTravelHist);}
  double getLogLikGen(){return logLikGen;}
  double getLogLikTime(){return logLikTime;}
  double getLogLikSpace(){return logLikSpace;}
  double getLogLikTravelHist(){return logLikTravelHist;}
  map<string, vector<bool> *>* getLoci(){return &loci;}
  vector<bool>* getLocus(string locusIn){return loci[locusIn];}
  vector<pair<Source *,int> >* getSources(){return &sources;}
  std::map<std::tuple<int, int, std::string>, double> getLocusCache() {return locusCache;}
  std::map<std::tuple<int, std::string, std::string>, double> getSrcLocusCache() {return srcLocusCache;}
  int getTime(){return time;}
  double getLat(){return lat;}
  double getLon(){return lon;}
  map<Node*,int, Node::classcomp>* getTimeDists(){return &timeDists;}
  map<Node*,std::tuple<double, double>, Node::classcomp >* getSpaceDists(){return &spaceDists;}
  string getType(){return type;}
  bool hasAncSrc();
  bool isAncestor(Node *);
  bool isFounder();
  bool isOrphan();
  void loadTime(int);
  void removeAncestor(int);
  void removeAllAncestors();
  void removeSource();
  void resetAll();
  void resetLogLik(){logLikGen = 0.0, logLikTime = 0.0, logLikSpace = 0.0, logLikTravelHist = 0.0;}
  void setCode(int code){this->code = code;}
  void setGenDist(Node *, int);
  void set_k_anc(Node *, int);
  void set_k_src(Source *, int);
  void setKStride(int k_stride){this->k_stride = k_stride;}
  void setLocus(string locus, vector<bool>* alleles){loci[locus] = alleles;}
  void setLocusOrder(std::vector<std::string> locus_order) {this->locus_order = locus_order;}
  void setAlleleCount(std::map<string, int> alleleCount){ this->alleleCount = alleleCount;}
  void setTime(int node_time){this->time = node_time;}
  void setTimeDist(Node *, int);
  void setSpaceDist(Node *, std::tuple<double,double>);
  void setType(string typeIn){this->type = typeIn;}
  void setLat(double node_lat){this->lat = node_lat;}
  void setLon(double node_lon){this->lon = node_lon;}
  void setHistory(int node_history){this->history = node_history;}
  bool operator<(const Node&);
}; // const  class Node



#endif // _NODE_H
