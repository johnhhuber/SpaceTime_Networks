#include "network.h"

using namespace std;



double Network::add_anc(Node * nodeNow)
{
  double runif, logPropDiff = 0.0, probSum = 0.0, probSumMax = 0.0;
  map<Node *, double>::iterator mapItr;

  mapItr = nodeNow->getEdgeProbs()->begin();
  for(; mapItr != nodeNow->getEdgeProbs()->end(); mapItr++){
    if(!nodeNow->isAncestor(mapItr->first)){
      probSumMax += mapItr->second;
    }
  }

  runif = gsl_rng_uniform_pos(this->rGsl) * probSumMax;

  mapItr = nodeNow->getEdgeProbs()->begin();
  while(probSum < runif && mapItr != nodeNow->getEdgeProbs()->end()){
    if(nodeNow->isAncestor(mapItr->first)){
      mapItr++;
    }
    else{
      probSum += mapItr->second;
      mapItr++;
    }
  }
  --mapItr;

  if(nodeNow->isAncestor(mapItr->first)){
    return 100.0;
  }
  else{
    if(nodeNow->isFounder()){
      logPropDiff += log(settings["PROP_FOUNDER"] / double(nonLocalSources)) - log(settings["PROP_ADD"]);
    }
    if(nodeNow->isOrphan()){
      logPropDiff += log(settings["PROP_ORPHAN"]) - log(settings["PROP_ADD"]);
    }

    int kNew = 1;
    if(settings["PROP_K"]>0){
      kNew = getRandom_k();
      logPropDiff -= logProposal_k(kNew);
    }

    nodeNow->assignAncestor(mapItr->first, kNew);
    nodeNow->removeSource();

    return (logPropDiff - log(double(nodeNow->getAncestors()->size())) - log(mapItr->second / probSumMax));
  }
} // add_anc()



void Network::calcGenDists()
{
  int genDist;
  vector<Node *>::iterator nodeItr_i = nodes.begin(), nodeItr_j;

  for(; nodeItr_i != nodes.end(); nodeItr_i++){
    for(nodeItr_j = nodeItr_i + 1; nodeItr_j != nodes.end(); nodeItr_j++){
      genDist = (*nodeItr_i)->calcGenDist(*nodeItr_j);
      (*nodeItr_i)->setGenDist(*nodeItr_j, genDist);
      (*nodeItr_j)->setGenDist(*nodeItr_i, genDist);
    }
  }
} // calcGenDists()



void Network::initPrGObsOffspring(){
  std::vector<std::string>::iterator locItr;
  locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++){
    vector<double> v(4);
    Pr_g_obs_offspring.insert(pair<string,vector<double> >(*locItr, v));
  }
} // initPrGObsOffspring()



void Network::setPrGObsOffspring(){
  std::vector<std::string>::iterator locItr;
  locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++){
    Pr_g_obs_offspring[*locItr][0] = (1.0 - epsilon_pos[*locItr]);
    Pr_g_obs_offspring[*locItr][1] = (epsilon_neg[*locItr]);
    Pr_g_obs_offspring[*locItr][2] = (epsilon_pos[*locItr]);
    Pr_g_obs_offspring[*locItr][3] = (1.0 - epsilon_neg[*locItr]);
  }
} // setPrGObsOffspring()



void Network::setPrGObsOffspring(const std::string &locus){
    Pr_g_obs_offspring[locus][0] = 1.0 - epsilon_pos[locus];
    Pr_g_obs_offspring[locus][1] = epsilon_neg[locus];
    Pr_g_obs_offspring[locus][2] = epsilon_pos[locus];
    Pr_g_obs_offspring[locus][3] = 1.0 - epsilon_neg[locus];
} // setPrGObsOffspring()



void Network::setPrGInheritance() {
  std::vector<std::string>::iterator locItr;
  Pr_g_inheritance.clear();

  for(int k = 0; k < int(settings["MAX_K"]); k++){

    map<string, vector<double> > m;
    Pr_g_inheritance.push_back(m);

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++){

      vector<double> v(4);
      v[0] = pow(1.0 - mu / double(alleleCount[*locItr]), k + 1);
      v[1] = 1.0 - pow(lambda, k + 1);
      v[2] = 1.0 - pow(1.0 - mu / double(alleleCount[*locItr]), k + 1);
      v[3] = pow(lambda, k + 1);

      Pr_g_inheritance[k].insert(pair<string,vector<double> >(*locItr, v));

    }
  }
} // setPrGInheritance()



void Network::setPrGInheritance(const std::string &locus) {
  for(int k = 0; k < int(settings["MAX_K"]); k++){
    vector<double> v(4);
    v[0] = pow(1.0 - mu / double(alleleCount[locus]), k + 1);
    v[1] = 1.0 - pow(lambda, k + 1);
    v[2] = 1.0 - pow(1.0 - mu / double(alleleCount[locus]), k + 1);
    v[3] = pow(lambda, k + 1);
    Pr_g_inheritance[k][locus] = v;
  }
}



void Network::initPrGTrueParent(){
  std::vector<std::string>::iterator locItr;
  int num_alleles;

  Pr_g_true_parent.resize(int(settings["MAX_K"]));

  for(int k = 0; k < int(settings["MAX_K"]); k++){
    map<string, vector<vector<double> > *> probs_loci;

    Pr_g_true_parent[k] = probs_loci;

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++){

      vector<vector<double>> * probs_alleles = new vector<vector<double > >;
      Pr_g_true_parent[k].insert(pair<string,vector<vector<double> > *>(*locItr, probs_alleles));
      num_alleles = alleleCount[*locItr];

      for(size_t i = 0; i < num_alleles; i++){
        vector<double> v(4);
        Pr_g_true_parent[k][*locItr]->push_back(v);
      }

    }
  }
} // initPrGTrueParent()



void Network::setPrGTrueParent(){
  std::vector<std::string>::iterator locItr;
  vector<double>::iterator alleleItr;
  double freq_true;

  for(int k = 0; k < int(settings["MAX_K"]); k++){
    locItr = locus_order.begin();

    for(; locItr != locus_order.end(); locItr++){
      alleleItr = (*sources[LOCAL_SOURCE]->getLociFreqs())[*locItr]->begin();

      for(int a = 0; alleleItr != (*sources[LOCAL_SOURCE]->getLociFreqs())[*locItr]->end(); a++){
        vector<double> v = Pr_g_true_parent[k][*locItr]->at(a);

        freq_true = (*alleleItr - epsilon_pos[*locItr]) / (1.0 - epsilon_neg[*locItr] - epsilon_pos[*locItr]);

        v[0] = (1.0 - epsilon_pos[*locItr]) * (1.0 - freq_true) / (1.0 - *alleleItr);
        v[1] = epsilon_pos[*locItr] * (1.0 - freq_true) / (*alleleItr);
        v[2] = epsilon_neg[*locItr] * freq_true / (1.0 - *alleleItr);
        v[3] = (1.0 - epsilon_neg[*locItr]) * freq_true / (*alleleItr);


        (*Pr_g_true_parent[k][*locItr])[a] = v;
        alleleItr++;
      }
    }
  }
} // setPrGTrueParent()



void Network::setPrGTrueParent(const std::string &locus){
  vector<double>::iterator alleleItr;
  double freq_true;

  for(int k = 0; k < int(settings["MAX_K"]); k++){
    alleleItr = (*sources[LOCAL_SOURCE]->getLociFreqs())[locus]->begin();

    for(int a = 0; alleleItr != (*sources[LOCAL_SOURCE]->getLociFreqs())[locus]->end(); a++){
      vector<double> v = Pr_g_true_parent[k][locus]->at(a);

      freq_true = (*alleleItr - epsilon_pos[locus]) / (1.0 - epsilon_neg[locus] - epsilon_pos[locus]);

      v[0] = (1.0 - epsilon_pos[locus]) * (1.0 - freq_true) / (1.0 - *alleleItr);
      v[1] = epsilon_pos[locus] * (1.0 - freq_true) / (*alleleItr);
      v[2] = epsilon_neg[locus] * freq_true / (1.0 - *alleleItr);
      v[3] = (1.0 - epsilon_neg[locus]) * freq_true / (*alleleItr);
      (*Pr_g_true_parent[k][locus])[a] = v;
      alleleItr++;
    }
  }
} // setPrGTrueParent()



void Network::initPrVecs(){
  int s_i = 0;
  std::vector<std::string>::iterator locItr;

  // for(int k = 0; k < settings["MAX_K"]; k++){
  //   locItr = locus_order.begin();
  //
  //   for(; locItr != locus_order.end(); locItr++){
  //     string locusLabel = *locItr;
  //     s_i += 22 * alleleCount[locusLabel];
  //   }
  // }
  Pr_vec.resize(int(settings["MAX_K"]) * 22 * k_stride, 0);
  Pr_vec_copy.resize(int(settings["MAX_K"]) * 22 * k_stride, 0);

  map<string,Source *>::iterator sourceItr = sources.begin();
  for(; sourceItr != sources.end(); sourceItr++) {
    sourceItr->second->initPrVecs();
  }

} // initPrVecs()



void Network::setPrVecs(){
  std::vector<std::string>::iterator locItr;
  vector<double>::iterator alleleItr;
  string locusLabel;

  double local_Pr_g_obs_offspring_0;
  double local_Pr_g_obs_offspring_1;
  double local_Pr_g_obs_offspring_2;
  double local_Pr_g_obs_offspring_3;

  double local_Pr_g_inheritance_0;
  double local_Pr_g_inheritance_1;
  double local_Pr_g_inheritance_2;
  double local_Pr_g_inheritance_3;

  double local_Pr_g_true_parent_0;
  double local_Pr_g_true_parent_1;
  double local_Pr_g_true_parent_2;
  double local_Pr_g_true_parent_3;

  double inheritance_2_true_parent_0;
  double inheritance_0_true_parent_0;
  double inheritance_2_true_parent_1;
  double inheritance_0_true_parent_1;
  double inheritance_1_true_parent_3;
  double inheritance_1_true_parent_2;
  double offspring_1_inheritance_2_true_parent_0;
  double offspring_0_inheritance_1_true_parent_2;
  double offspring_0_inheritance_0_true_parent_0;
  double offspring_1_inheritance_2_true_parent_1;
  double offspring_0_inheritance_1_true_parent_3;
  double offspring_0_inheritance_0_true_parent_1;
  double offspring_3_inheritance_2_true_parent_0;
  double offspring_2_inheritance_1_true_parent_2;
  double offspring_2_inheritance_0_true_parent_0;
  double offspring_3_inheritance_2_true_parent_1;
  double offspring_2_inheritance_1_true_parent_3;
  double offpsring_2_inheritance_0_true_parent_1;

  int a_i = 0;

  for(int k = 0; k < int(settings["MAX_K"]); k++){

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++){
      locusLabel = *locItr;

      local_Pr_g_obs_offspring_0 = Pr_g_obs_offspring[locusLabel][0];
      local_Pr_g_obs_offspring_1 = Pr_g_obs_offspring[locusLabel][1];
      local_Pr_g_obs_offspring_2 = Pr_g_obs_offspring[locusLabel][2];
      local_Pr_g_obs_offspring_3 = Pr_g_obs_offspring[locusLabel][3];

      local_Pr_g_inheritance_0 = Pr_g_inheritance[k][locusLabel][0];
      local_Pr_g_inheritance_1 = Pr_g_inheritance[k][locusLabel][1];
      local_Pr_g_inheritance_2 = Pr_g_inheritance[k][locusLabel][2];
      local_Pr_g_inheritance_3 = Pr_g_inheritance[k][locusLabel][3];

      for(int a = 0; a < alleleCount[locusLabel]; a++){

        local_Pr_g_true_parent_0 = (*Pr_g_true_parent[k][locusLabel])[a][0];
        local_Pr_g_true_parent_1 = (*Pr_g_true_parent[k][locusLabel])[a][1];
        local_Pr_g_true_parent_2 = (*Pr_g_true_parent[k][locusLabel])[a][2];
        local_Pr_g_true_parent_3 = (*Pr_g_true_parent[k][locusLabel])[a][3];

        // XXXX = unk | unk | unk | unk

        // 0 = 00
        inheritance_2_true_parent_0 = local_Pr_g_inheritance_2 * local_Pr_g_true_parent_0;
        inheritance_0_true_parent_0 = local_Pr_g_inheritance_0 * local_Pr_g_true_parent_0;
        inheritance_2_true_parent_1 = local_Pr_g_inheritance_2 * local_Pr_g_true_parent_1;
        inheritance_0_true_parent_1 = local_Pr_g_inheritance_0 * local_Pr_g_true_parent_1;
        inheritance_1_true_parent_3 = local_Pr_g_inheritance_1 * local_Pr_g_true_parent_3;
        inheritance_1_true_parent_2 = local_Pr_g_inheritance_1 * local_Pr_g_true_parent_2;
        offspring_1_inheritance_2_true_parent_0 = local_Pr_g_obs_offspring_1 * inheritance_2_true_parent_0;
        offspring_0_inheritance_1_true_parent_2 = local_Pr_g_obs_offspring_0 * inheritance_1_true_parent_2;
        offspring_0_inheritance_0_true_parent_0 = local_Pr_g_obs_offspring_0 * inheritance_0_true_parent_0;
        offspring_1_inheritance_2_true_parent_1 = local_Pr_g_obs_offspring_1 * inheritance_2_true_parent_1;
        offspring_0_inheritance_1_true_parent_3 = local_Pr_g_obs_offspring_0 * inheritance_1_true_parent_3;
        offspring_0_inheritance_0_true_parent_1 = local_Pr_g_obs_offspring_0 * inheritance_0_true_parent_1;
        offspring_3_inheritance_2_true_parent_0 = local_Pr_g_obs_offspring_3 * inheritance_2_true_parent_0;
        offspring_2_inheritance_1_true_parent_2 = local_Pr_g_obs_offspring_2 * inheritance_1_true_parent_2;
        offspring_2_inheritance_0_true_parent_0 = local_Pr_g_obs_offspring_2 * inheritance_0_true_parent_0;
        offspring_3_inheritance_2_true_parent_1 = local_Pr_g_obs_offspring_3 * inheritance_2_true_parent_1;
        offspring_2_inheritance_1_true_parent_3 = local_Pr_g_obs_offspring_2 * inheritance_1_true_parent_3;
        offpsring_2_inheritance_0_true_parent_1 = local_Pr_g_obs_offspring_2 * inheritance_0_true_parent_1;

        Pr_vec[a_i] =
          log(local_Pr_g_obs_offspring_1 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_2 +
            offspring_1_inheritance_2_true_parent_0 +
            offspring_0_inheritance_1_true_parent_2 +
            offspring_0_inheritance_0_true_parent_0);

        // 1 = 01

        Pr_vec[a_i + 1] =
          log(local_Pr_g_obs_offspring_1 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_3 +
            offspring_1_inheritance_2_true_parent_1 +
            offspring_0_inheritance_1_true_parent_3 +
            offspring_0_inheritance_0_true_parent_1);

        // 2 = 10
        Pr_vec[a_i + 2] =
          log(local_Pr_g_obs_offspring_3 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_2 +
            offspring_3_inheritance_2_true_parent_0 +
            offspring_2_inheritance_1_true_parent_2 +
            offspring_2_inheritance_0_true_parent_0);

        // 3 = 11
        Pr_vec[a_i + 3] =
          log(local_Pr_g_obs_offspring_3 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_3 +
            offspring_3_inheritance_2_true_parent_1 +
            offspring_2_inheritance_1_true_parent_3 +
            offpsring_2_inheritance_0_true_parent_1);

        // XXOX = unk | unk | 0 | unk

        // 0 = 00
        Pr_vec[a_i + 4] =
          log(offspring_1_inheritance_2_true_parent_0 + offspring_0_inheritance_0_true_parent_0);

        // 1 = 01
        Pr_vec[a_i + 5] =
          log(offspring_1_inheritance_2_true_parent_1 + offspring_0_inheritance_0_true_parent_1);

        // 2 = 10
        Pr_vec[a_i + 6] =
          log(offspring_3_inheritance_2_true_parent_0 + offspring_2_inheritance_0_true_parent_0);

        // 3 = 11
        Pr_vec[a_i + 7] =
          log(offspring_3_inheritance_2_true_parent_1 + offpsring_2_inheritance_0_true_parent_1);


        // XOXX = unk | 0 | unk | unk

        // 0 = 00
        Pr_vec[a_i + 8] =
          log(offspring_0_inheritance_1_true_parent_2 + offspring_0_inheritance_0_true_parent_0);

        // 1 = 01
        Pr_vec[a_i + 9] =
          log(offspring_0_inheritance_1_true_parent_3 + offspring_0_inheritance_0_true_parent_1);

        // 2 = 10
        Pr_vec[a_i + 10] =
          log(offspring_2_inheritance_1_true_parent_2 + offspring_2_inheritance_0_true_parent_0);

        // 3 = 11
        Pr_vec[a_i + 11] =
          log(offspring_2_inheritance_1_true_parent_3 + offpsring_2_inheritance_0_true_parent_1);


        // XOOX = unk | 0 | 0 | unk

        // 0 = 00
        Pr_vec[a_i + 12] =
          log(offspring_0_inheritance_0_true_parent_0);

        // 1 = 01
        Pr_vec[a_i + 13] =
          log(offspring_0_inheritance_0_true_parent_1);

        // 2 = 10
        Pr_vec[a_i + 14] =
          log(offspring_2_inheritance_0_true_parent_0);

        // 3 = 11
        Pr_vec[a_i + 15] =
          log(offpsring_2_inheritance_0_true_parent_1);


        // XOX = unk | 0 | unk

        // 0 = 00
        Pr_vec[a_i + 16] =
          log(inheritance_2_true_parent_0 + inheritance_0_true_parent_0);

        // 1 = 01
        Pr_vec[a_i + 17] =
          log(inheritance_2_true_parent_1 + inheritance_0_true_parent_1);


        // OXX = 0 | unk | unk

        // 0 = 00
        Pr_vec[a_i + 18] =
          log(inheritance_1_true_parent_2 + inheritance_0_true_parent_0);

        // 1 = 01
        Pr_vec[a_i + 19] =
          log(inheritance_1_true_parent_3 + inheritance_0_true_parent_1);


        // OOX = 0 | 0 | unk

        // 0 = 00
        Pr_vec[a_i + 20] =
          log(inheritance_0_true_parent_0);
        // 1 = 01
        Pr_vec[a_i + 21] =
          log(inheritance_0_true_parent_1);


        a_i += 22;
      }
    }
  }
} // setPrVecs()



void Network::setPrVecs(const std::string &locusLabel){
  std::vector<std::string>::iterator locItr;
  vector<double>::iterator alleleItr;

  double local_Pr_g_obs_offspring_0;
  double local_Pr_g_obs_offspring_1;
  double local_Pr_g_obs_offspring_2;
  double local_Pr_g_obs_offspring_3;

  double local_Pr_g_inheritance_0;
  double local_Pr_g_inheritance_1;
  double local_Pr_g_inheritance_2;
  double local_Pr_g_inheritance_3;

  double local_Pr_g_true_parent_0;
  double local_Pr_g_true_parent_1;
  double local_Pr_g_true_parent_2;
  double local_Pr_g_true_parent_3;

  double inheritance_2_true_parent_0;
  double inheritance_0_true_parent_0;
  double inheritance_2_true_parent_1;
  double inheritance_0_true_parent_1;
  double inheritance_1_true_parent_3;
  double inheritance_1_true_parent_2;
  double offspring_1_inheritance_2_true_parent_0;
  double offspring_0_inheritance_1_true_parent_2;
  double offspring_0_inheritance_0_true_parent_0;
  double offspring_1_inheritance_2_true_parent_1;
  double offspring_0_inheritance_1_true_parent_3;
  double offspring_0_inheritance_0_true_parent_1;
  double offspring_3_inheritance_2_true_parent_0;
  double offspring_2_inheritance_1_true_parent_2;
  double offspring_2_inheritance_0_true_parent_0;
  double offspring_3_inheritance_2_true_parent_1;
  double offspring_2_inheritance_1_true_parent_3;
  double offpsring_2_inheritance_0_true_parent_1;

  local_Pr_g_obs_offspring_0 = Pr_g_obs_offspring[locusLabel][0];
  local_Pr_g_obs_offspring_1 = Pr_g_obs_offspring[locusLabel][1];
  local_Pr_g_obs_offspring_2 = Pr_g_obs_offspring[locusLabel][2];
  local_Pr_g_obs_offspring_3 = Pr_g_obs_offspring[locusLabel][3];

  int a_i = 0;

  for(int k = 0; k < int(settings["MAX_K"]); k++){

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++){
      if(*locItr != locusLabel) {
        a_i += 22 * alleleCount[*locItr];
      } else {

        local_Pr_g_inheritance_0 = Pr_g_inheritance[k][locusLabel][0];
        local_Pr_g_inheritance_1 = Pr_g_inheritance[k][locusLabel][1];
        local_Pr_g_inheritance_2 = Pr_g_inheritance[k][locusLabel][2];
        local_Pr_g_inheritance_3 = Pr_g_inheritance[k][locusLabel][3];

        for(int a = 0; a < alleleCount[locusLabel]; a++){

          local_Pr_g_true_parent_0 = (*Pr_g_true_parent[k][locusLabel])[a][0];
          local_Pr_g_true_parent_1 = (*Pr_g_true_parent[k][locusLabel])[a][1];
          local_Pr_g_true_parent_2 = (*Pr_g_true_parent[k][locusLabel])[a][2];
          local_Pr_g_true_parent_3 = (*Pr_g_true_parent[k][locusLabel])[a][3];

          // XXXX = unk | unk | unk | unk

          // 0 = 00
          inheritance_2_true_parent_0 = local_Pr_g_inheritance_2 * local_Pr_g_true_parent_0;
          inheritance_0_true_parent_0 = local_Pr_g_inheritance_0 * local_Pr_g_true_parent_0;
          inheritance_2_true_parent_1 = local_Pr_g_inheritance_2 * local_Pr_g_true_parent_1;
          inheritance_0_true_parent_1 = local_Pr_g_inheritance_0 * local_Pr_g_true_parent_1;
          inheritance_1_true_parent_3 = local_Pr_g_inheritance_1 * local_Pr_g_true_parent_3;
          inheritance_1_true_parent_2 = local_Pr_g_inheritance_1 * local_Pr_g_true_parent_2;
          offspring_1_inheritance_2_true_parent_0 = local_Pr_g_obs_offspring_1 * inheritance_2_true_parent_0;
          offspring_0_inheritance_1_true_parent_2 = local_Pr_g_obs_offspring_0 * inheritance_1_true_parent_2;
          offspring_0_inheritance_0_true_parent_0 = local_Pr_g_obs_offspring_0 * inheritance_0_true_parent_0;
          offspring_1_inheritance_2_true_parent_1 = local_Pr_g_obs_offspring_1 * inheritance_2_true_parent_1;
          offspring_0_inheritance_1_true_parent_3 = local_Pr_g_obs_offspring_0 * inheritance_1_true_parent_3;
          offspring_0_inheritance_0_true_parent_1 = local_Pr_g_obs_offspring_0 * inheritance_0_true_parent_1;
          offspring_3_inheritance_2_true_parent_0 = local_Pr_g_obs_offspring_3 * inheritance_2_true_parent_0;
          offspring_2_inheritance_1_true_parent_2 = local_Pr_g_obs_offspring_2 * inheritance_1_true_parent_2;
          offspring_2_inheritance_0_true_parent_0 = local_Pr_g_obs_offspring_2 * inheritance_0_true_parent_0;
          offspring_3_inheritance_2_true_parent_1 = local_Pr_g_obs_offspring_3 * inheritance_2_true_parent_1;
          offspring_2_inheritance_1_true_parent_3 = local_Pr_g_obs_offspring_2 * inheritance_1_true_parent_3;
          offpsring_2_inheritance_0_true_parent_1 = local_Pr_g_obs_offspring_2 * inheritance_0_true_parent_1;

          Pr_vec[a_i] =
            log(local_Pr_g_obs_offspring_1 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_2 +
              offspring_1_inheritance_2_true_parent_0 +
              offspring_0_inheritance_1_true_parent_2 +
              offspring_0_inheritance_0_true_parent_0);

          // 1 = 01

          Pr_vec[a_i + 1] =
            log(local_Pr_g_obs_offspring_1 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_3 +
              offspring_1_inheritance_2_true_parent_1 +
              offspring_0_inheritance_1_true_parent_3 +
              offspring_0_inheritance_0_true_parent_1);

          // 2 = 10
          Pr_vec[a_i + 2] =
            log(local_Pr_g_obs_offspring_3 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_2 +
              offspring_3_inheritance_2_true_parent_0 +
              offspring_2_inheritance_1_true_parent_2 +
              offspring_2_inheritance_0_true_parent_0);

          // 3 = 11
          Pr_vec[a_i + 3] =
            log(local_Pr_g_obs_offspring_3 * local_Pr_g_inheritance_3 * local_Pr_g_true_parent_3 +
              offspring_3_inheritance_2_true_parent_1 +
              offspring_2_inheritance_1_true_parent_3 +
              offpsring_2_inheritance_0_true_parent_1);

          // XXOX = unk | unk | 0 | unk

          // 0 = 00
          Pr_vec[a_i + 4] =
            log(offspring_1_inheritance_2_true_parent_0 + offspring_0_inheritance_0_true_parent_0);

          // 1 = 01
          Pr_vec[a_i + 5] =
            log(offspring_1_inheritance_2_true_parent_1 + offspring_0_inheritance_0_true_parent_1);

          // 2 = 10
          Pr_vec[a_i + 6] =
            log(offspring_3_inheritance_2_true_parent_0 + offspring_2_inheritance_0_true_parent_0);

          // 3 = 11
          Pr_vec[a_i + 7] =
            log(offspring_3_inheritance_2_true_parent_1 + offpsring_2_inheritance_0_true_parent_1);


          // XOXX = unk | 0 | unk | unk

          // 0 = 00
          Pr_vec[a_i + 8] =
            log(offspring_0_inheritance_1_true_parent_2 + offspring_0_inheritance_0_true_parent_0);

          // 1 = 01
          Pr_vec[a_i + 9] =
            log(offspring_0_inheritance_1_true_parent_3 + offspring_0_inheritance_0_true_parent_1);

          // 2 = 10
          Pr_vec[a_i + 10] =
            log(offspring_2_inheritance_1_true_parent_2 + offspring_2_inheritance_0_true_parent_0);

          // 3 = 11
          Pr_vec[a_i + 11] =
            log(offspring_2_inheritance_1_true_parent_3 + offpsring_2_inheritance_0_true_parent_1);


          // XOOX = unk | 0 | 0 | unk

          // 0 = 00
          Pr_vec[a_i + 12] =
            log(offspring_0_inheritance_0_true_parent_0);

          // 1 = 01
          Pr_vec[a_i + 13] =
            log(offspring_0_inheritance_0_true_parent_1);

          // 2 = 10
          Pr_vec[a_i + 14] =
            log(offspring_2_inheritance_0_true_parent_0);

          // 3 = 11
          Pr_vec[a_i + 15] =
            log(offpsring_2_inheritance_0_true_parent_1);


          // XOX = unk | 0 | unk

          // 0 = 00
          Pr_vec[a_i + 16] =
            log(inheritance_2_true_parent_0 + inheritance_0_true_parent_0);

          // 1 = 01
          Pr_vec[a_i + 17] =
            log(inheritance_2_true_parent_1 + inheritance_0_true_parent_1);


          // OXX = 0 | unk | unk

          // 0 = 00
          Pr_vec[a_i + 18] =
            log(inheritance_1_true_parent_2 + inheritance_0_true_parent_0);

          // 1 = 01
          Pr_vec[a_i + 19] =
            log(inheritance_1_true_parent_3 + inheritance_0_true_parent_1);


          // OOX = 0 | 0 | unk

          // 0 = 00
          Pr_vec[a_i + 20] =
            log(inheritance_0_true_parent_0);
          // 1 = 01
          Pr_vec[a_i + 21] =
            log(inheritance_0_true_parent_1);

          a_i += 22;
        }
      }
    }
  }
} // setPrVecs()



void Network::calcEdgeProbs()
{
  calcGenDists();
  calcTimeDists();
  calcSpaceDists();

  for(vector<Node *>::iterator nodeItr = nodes.begin();
      nodeItr != nodes.end();
      nodeItr++)
  {
    (*nodeItr)->calcEdgeProbs(beta);
  }
} // calcEdgeProbs()



double Network::calcLogLik()
{
  vector<Node *>::iterator itrNode;

  double logLikOut = 0.0;

  for(itrNode = nodes.begin(); itrNode != nodes.end(); itrNode++){
    (*itrNode)->resetLogLik();
    (*itrNode)->calcLogLik(Pr_vec, *sources[LOCAL_SOURCE], validLoci, diffusion_coef, tau_s, tau_l, prop_reported_travel);
    logLikOut += (*itrNode)->getLogLik();
  }

  if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
    validateLoci();

  return logLikOut;
} // calcLogLik()



void Network::calcSourcesMinMax()
{
  map<string,Source *>::iterator sourceItr;
  std::vector<std::string>::iterator locItr;
  vector<double>::iterator alleleItr;
  double prevMax_tmp, prevMin_tmp;

  for(locItr = locus_order.begin(); locItr != locus_order.end(); locItr++){
    prevMax_tmp = 0.0;
    prevMin_tmp = 1.0;

    for(sourceItr = sources.begin(); sourceItr != sources.end(); sourceItr++){
      alleleItr = sourceItr->second->getLocusFreqs(*locItr)->begin();
      for(; alleleItr != sourceItr->second->getLocusFreqs(*locItr)->end(); alleleItr++){
        if(*alleleItr > prevMax_tmp) {
          prevMax_tmp = *alleleItr;
        }
        if(*alleleItr < prevMin_tmp) {
          prevMin_tmp = *alleleItr;
        }
      }
    }
    prevMax.insert(pair<string, double>(*locItr, prevMax_tmp));
    prevMin.insert(pair<string, double>(*locItr, prevMin_tmp));
  }

} // calcSourcesMinMax()



void Network::calcTimeDists()
{
  int timeDist;
  vector<Node *>::iterator nodeItr_i = nodes.begin(), nodeItr_j;

  for(; nodeItr_i != nodes.end(); nodeItr_i++){
    for(nodeItr_j = nodeItr_i + 1; nodeItr_j != nodes.end(); nodeItr_j++){
      timeDist = (*nodeItr_i)->calcTimeDist(*nodeItr_j);
      (*nodeItr_i)->setTimeDist(*nodeItr_j, timeDist);
      (*nodeItr_j)->setTimeDist(*nodeItr_i, timeDist);
    }
  }
} // calcTimeDists()



void Network::calcSpaceDists()
{
  vector<Node *>::iterator nodeItr_i = nodes.begin(), nodeItr_j;

  for(; nodeItr_i != nodes.end(); nodeItr_i++)
  {
    for(nodeItr_j = nodeItr_i + 1; nodeItr_j != nodes.end(); nodeItr_j++)
    {
      std::tuple<double, double> spaceDist = (*nodeItr_i)->calcSpaceDist(*nodeItr_j);
      (*nodeItr_i)->setSpaceDist(*nodeItr_j, spaceDist);
      (*nodeItr_j)->setSpaceDist(*nodeItr_i, spaceDist);
    }
  }
} //calcSpaceDists()



map<string, double> Network::change_epsilon_neg(const std::string &locus)
{
  map<string, double> epsilon_neg_out = epsilon_neg;

  double new_epsilon_neg = gsl_rng_uniform(this->rGsl) * (1.0 - prevMax[locus]);
  assert(new_epsilon_neg < (1.0 - prevMax[locus]));
  epsilon_neg_out[locus] = new_epsilon_neg;

  return epsilon_neg_out;
} // change_epsilon_neg()



map<string, double> Network::change_epsilon_pos(const std::string &locus)
{
  map<string, double> epsilon_pos_out = epsilon_pos;

  double new_epsilon_pos = gsl_rng_uniform(this->rGsl) * prevMin[locus];
  assert(new_epsilon_pos < prevMin[locus]);
  epsilon_pos_out[locus] = new_epsilon_pos;

  return epsilon_pos_out;
} // change_epsilon_pos()



double Network::change_k_anc(Node * nodeIn)
{
  double logPropDiff = 0.0;
  vector<pair<Node *,int> >* ancestors =
    nodeIn->getAncestors();

  int ancIn = gsl_rng_uniform_int(this->rGsl, ancestors->size());
  int kIn = (*ancestors)[ancIn].second;

  // get new k and logProposal for change to new k
  double probSumMax = 0;
  for(int i=1; i<=settings["MAX_K"]; i++){
    if(i != kIn){
      probSumMax += exp(logProposal_k(i));
    }
  }
  int kOut = 0;
  double probSumK = 0.0;
  double runif = gsl_rng_uniform_pos(this->rGsl) * probSumMax;
  do{
    kOut++;
    if(kOut != kIn){
      probSumK += exp(logProposal_k(kOut));
    }
  }while(probSumK < runif && kOut < settings["MAX_K"]);
  logPropDiff -= logProposal_k(kOut)-log(probSumMax);

  // get logProposal for change to old k
  probSumMax = 0;
  for(int i=1; i<=settings["MAX_K"]; i++){
    if(i != kOut){
      probSumMax += exp(logProposal_k(i));
    }
  }
  logPropDiff += logProposal_k(kIn)-log(probSumMax);

  nodeIn->set_k_anc((*ancestors)[ancIn].first, kOut);

  return(logPropDiff);
} // change_k_anc()



void Network::change_k_src(Node * nodeIn)
{
  vector<pair<Source *,int> >* sources =
    nodeIn->getSources();

  int srcIn = gsl_rng_uniform_int(this->rGsl, sources->size());
  int kIn = (*sources)[srcIn].second, kOut;

  do{
    kOut = gsl_rng_uniform_int(this->rGsl, int(settings["MAX_K"])) + 1;
  }while(kIn == kOut);

  nodeIn->set_k_src((*sources)[srcIn].first, kOut);
} // change_k_src()



double Network::change_lambda()
{
  double lambda_out;
  if(lambda >= (1.0 - settings["LAMBDA_width"])) {
    lambda_out = runif_interval(1.0 - settings["LAMBDA_width"], settings["LAMBDA_width"]);
  } else if(lambda <= (settings["LAMBDA_min"] + settings["LAMBDA_width"])) {
    lambda_out = runif_interval(settings["LAMBDA_min"] + settings["LAMBDA_width"], settings["LAMBDA_width"]);
  } else {
    lambda_out = runif_interval(lambda, settings["LAMBDA_width"]);
  }

  return lambda_out;
} // change_lambda()



double Network::change_mu()
{
  double mu_out;
  if(mu <= settings["MU_width"]){
    mu_out = runif_interval(settings["MU_width"], settings["MU_width"]);
  }
  else if(mu >= (settings["MU_max"] - settings["MU_width"])){
    mu_out = runif_interval(settings["MU_max"] - settings["MU_width"], settings["MU_width"]);
  }
  else{
    mu_out = runif_interval(mu, settings["MU_width"]);
  }

  return mu_out;
} // change_mu()



double Network::change_tau_s()
{
  double tau_s_out;

  do
  {
    //tau_s_out = tau_s + gsl_ran_gaussian(rGsl, 0.000001);
    tau_s_out = tau_s + gsl_ran_gaussian(rGsl, 0.25);
  } while(tau_s_out < 0.0 || tau_s_out > 1.0);

  return tau_s_out;
}// Network::change_tau_s()



double Network::change_tau_l()
{
  double tau_l_out;

  do
  {
    //tau_l_out = tau_l + gsl_ran_gaussian(rGsl, 0.00001);
    tau_l_out = tau_l + gsl_ran_gaussian(rGsl, 0.25);
  } while(tau_l_out < 0.0 || tau_l_out > 1.0);

  return tau_l_out;
}// Network::change_tau_l()


double Network::change_diffusion_coef()
{
  double diffusion_coef_out;
  do
  {
    diffusion_coef_out = diffusion_coef + gsl_ran_gaussian(rGsl, 2.5);
  } while(diffusion_coef_out <= settings["DIFFUSION_min"]);

  return diffusion_coef_out;
} // change_diffusion_coef()


Source* Network::chooseSource(){
  int whichSource = gsl_rng_uniform_int(this->rGsl, nonLocalSources);

  Source* src;

  map<string,Source *>::iterator mapItr = sources.begin();
  for(int i = 0; mapItr != sources.end(); mapItr++){
    if(mapItr->first != LOCAL_SOURCE){
      if(i == whichSource){
        src = mapItr->second;
      }
      i++;
    }
  }
  return src;
} // chooseSource()



// bool Network::findNode(Node * nodeIn){
//   vector<Node *>::iterator vecItr;
//   for(vecItr = nodes.begin(); vecItr != nodes.end(); vecItr++)
//     if(*vecItr == nodeIn)
//       return 1;
//
//   return 0;
// } // findNode()

void Network::hotstartLoadParameters() {
  ifstream file(scalars_destination);
  int length = 0;
  std::string parameters;
  std::string header;
  char delim = ',';
  char epsilon_delim = '_';
  std::vector<std::string> parameter_vec;
  std::vector<std::string> header_vec;
  std::vector<std::string> parameter_entry;
  char c = '\0';

  if(!file.is_open()) {
    cout << scalars_file_path;
    perror( "Error While Opening Parameters File");
  }

  // Header
  std::getline(file, header);
  header_vec = util::split(header, delim);

  // Get last line
  file.seekg(0, file.end);
  length = file.tellg();
  for(int i = length - 2; i > 0; i--) {
    file.seekg(i);
    c = file.get();
    if(c == '\r' || c == '\n') {
      break;
    }
  }

  // Load Parameters
  std::getline(file, parameters);
  parameter_vec = util::split(parameters, delim);

  double parameter_value = 0;

  for (std::vector<std::string>::size_type i = 0; i != header_vec.size(); i++) {
    std::string parameter_label = header_vec[i];
    cout << header_vec[i] << endl;
    cout << parameter_vec[i] << endl;

    if(parameter_vec[i] == "nan") {
      parameter_value = std::numeric_limits<double>::quiet_NaN();
    } else {
      parameter_value = std::stod(parameter_vec[i], nullptr);
    }

    if (util::toLower(parameter_label) == "lambda") {
      lambda = parameter_value;
    } else if (util::toLower(parameter_label) == "mu") {
      mu = parameter_value;
    }
    else if (util::toLower(parameter_label) == "diffusion_coef") {
      diffusion_coef = parameter_value;
    }

    else if (util::toLower(parameter_label) == "tau_s") {
      tau_s = parameter_value;
    }

    else if (util::toLower(parameter_label) == "tau_l") {
      tau_l = parameter_value;
    }
    else {
      parameter_entry = util::split(parameter_label, epsilon_delim);

      if (parameter_entry.size() == 3) {
        std::string locus = parameter_entry[2];
        std::string epsilon_type = util::toLower(parameter_entry[1]);

        if (epsilon_type == "pos") {
          epsilon_pos[locus] = parameter_value;
        } else if (epsilon_type == "neg") {
          epsilon_neg[locus] = parameter_value;
        }
      }
    }
  }

  file.close();
}

void Network::hotstartLoadNetwork() {
  ifstream file(network_destination);
  char network_delim = ';';
  char data_delim = '-';
  string network_entry;
  string from_string;
  string to_string;
  string k_string;

  Source * src;
  Node * from;
  Node * to;
  int k;
  vector<string> network_data;
  vector<string> edge_data;

  int length = 0;
  char c = '\0';

  if(!file.is_open()) {
    perror("Error trying to open network file");
  }

  file.seekg(0, file.end);
  length = file.tellg();
  for(int i = length - 2; i > 0; i--) {
    file.seekg(i);
    c = file.get();
    if(c == '\r' || c == '\n') {
      break;
    }
  }

  std::getline(file, network_entry);
  network_data = util::split(network_entry, network_delim);
  for(vector<string>::iterator data = network_data.begin(); data != network_data.end(); data++) {
    //data is a pointer to a string element
    edge_data = util::split(*data, data_delim);
    from_string = edge_data[0];
    k_string = edge_data[1];
    to_string = edge_data[2];

    k = atoi(k_string.c_str());
    to = nodes[atoi(to_string.c_str()) - 1];

    map<string, Source *>::iterator from_src = sources.find(from_string);
    if(from_src != sources.end()) {
      src = from_src->second;
      to->assignSource(src, k);
    } else {
      from = nodes[atoi(from_string.c_str()) - 1];
      to->assignAncestor(from, k);
    };
  };
  file.close();
}

void Network::initialize_epsilon_neg()
{
  std::vector<std::string>::iterator locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++){
    epsilon_neg.insert(pair<string, double>(*locItr, gsl_rng_uniform(this->rGsl) * (1.0 - prevMax[*locItr])));
  }
} // initialize_epsilon_neg()



void Network::initialize_epsilon_pos()
{
  std::vector<std::string>::iterator locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++){
    epsilon_pos.insert(pair<string, double>(*locItr, gsl_rng_uniform(this->rGsl) * prevMin[*locItr]));
  }
} // initialize_epsilon_pos()


// Definition of the discrete probability distribution of k
void Network::initialize_k_proposal()
{
  log_proposal_k.clear();
  if(settings["PROP_K"]>0){
    if(settings["K_p"]==1){
      perror("Error: Invalid combination of values in setting file 'K_p=1' and 'PROP_K>0'.");
      throw;
    }
    else if(settings["K_p"]==0){
      for(unsigned int i = 1; i <= settings["MAX_K"]; i++){
        log_proposal_k.push_back(log(1/settings["MAX_K"]));
      }
    }
    else{
      double logPropSumMax = log(gsl_cdf_geometric_P((unsigned int)settings["MAX_K"], settings["K_p"]));
      for(unsigned int i = 1; i <= settings["MAX_K"]; i++){
        log_proposal_k.push_back(log(gsl_ran_geometric_pdf(i, settings["K_p"])) - logPropSumMax);
      }
    }
  }else{ // initialize log_proposal_k if settings["PROP_K"]==0, but log_proposal_k should never be used with setting["PROP_K"]==0
    log_proposal_k.push_back(log(1));
    for(unsigned int i = 2; i <= settings["MAX_K"]; i++){
      log_proposal_k.push_back(log(0));
    }
  }
} // initialize_k_proposal()


void Network::initializeParams()
{
  initialize_k_proposal();

  initialize_epsilon_neg();
  epsilon_neg_accept = 0;
  epsilon_neg_reject = 0;

  initialize_epsilon_pos();
  epsilon_pos_accept = 0;
  epsilon_pos_reject = 0;

  lambda = 0.5;
  lambda_accept = 0;
  lambda_reject = 0;

  mu = 0.05;
  mu_accept = 0;
  mu_reject = 0;

  //diffusion_coef = settings["DIFFUSION_COEF"];
  diffusion_coef = 2 * settings["DIFFUSION_min"];
  diffusion_coef_accept = 0;
  diffusion_coef_reject = 0;

  //tau_s = 1.0 - exp(-20.0);
  if(settings["BELIEVE_TRAVELHIST"] > 0.0)
  {
    tau_s = 1.0;
  }else{
    tau_s = 0.80;
  }
  tau_s_accept = 0;
  tau_s_reject = 0;

  //tau_l = exp(-20.0);
  if(settings["BELIEVE_TRAVELHIST"] > 0.0)
  {
    tau_l = 0.0;
  }else{
    tau_l = 0.20;
  }
  tau_l_accept = 0;
  tau_l_reject = 0;

  change_k_accept = 0;
  change_k_reject = 0;

  switch_anc_anc_accept = 0;
  switch_anc_anc_reject = 0;

  add_anc_accept = 0;
  add_anc_reject = 0;

  remove_anc_accept = 0;
  remove_anc_reject = 0;

  make_founder_accept = 0;
  make_founder_reject = 0;

  make_orphan_accept = 0;
  make_orphan_reject = 0;
} // initializeParams()



void Network::initialize()
{
  loadSourcesFile();
  writeHeader();
  calcSourcesMinMax();
  loadNodesFile();
  randomTree();
  loadProposalProbs();
  calcEdgeProbs();
  initializeParams();
  initPrVecs();
  initPrGObsOffspring();
  initPrGTrueParent();
  invalidateLoci();
  recalculateProbs();
} // initialize()

void Network::hotstartInitialize()
{
  loadSourcesFile();
  calcSourcesMinMax();
  loadNodesFile();
  hotstartLoadNetwork();
  loadProposalProbs();
  calcEdgeProbs();
  initializeParams();
  initPrVecs();
  initPrGObsOffspring();
  initPrGTrueParent();
  invalidateLoci();
  hotstartLoadParameters();
  recalculateProbs();
}



bool Network::isCyclic(Node * origin, Node * target)
{
  vector<pair<Node *, int> > *anc;
  vector<pair<Node *, int> >::iterator nodeItr;

  anc = target->getAncestors();
  nodeItr = anc->begin();

  if(anc->empty()) {
      return false;
  }

  for(; nodeItr != anc->end(); nodeItr++) {
    if(origin->getCode() == nodeItr->first->getCode()){
      return true;
    }
  }

  nodeItr = anc->begin();
  for(; nodeItr != anc->end(); nodeItr++) {
    if(isCyclic(origin, nodeItr->first)) {
      return true;
    }
  }
  return false;
} // isCyclic()



bool Network::isCyclic(Node * origin)
{
  vector<pair<Node *, int> > *anc;
  vector<pair<Node *, int> >::iterator nodeItr;

  anc = origin->getAncestors();
  nodeItr = anc->begin();

  if(anc->empty()) {
    return false;
  }

  for(; nodeItr != anc->end(); nodeItr++) {
    if(origin->getCode() == nodeItr->first->getCode()){
      return true;
    }
  }

  nodeItr = anc->begin();
  for(; nodeItr != anc->end(); nodeItr++) {
    if(isCyclic(origin, nodeItr->first)) {
      return true;
    }
  }
  return false;
} // isCyclic()



void Network::loadNetworkFile() {
  ifstream file(network_file_path);
  char network_delim = ';';
  char data_delim = '-';
  string network_entry;
  string from_string;
  string to_string;
  string k_string;

  Source * src;
  Node * from;
  Node * to;
  int k;
  vector<string> network_data;
  vector<string> edge_data;

  if(!file.is_open()) {
    perror("Error trying to open network file");
  }
  while(getline(file, network_entry)) {
    //split line into vector of strings
    network_data = util::split(network_entry, network_delim);
    for(
      vector<string>::iterator data = network_data.begin();
      data != network_data.end();
      data++) {
        //data is a pointer to a string element
        edge_data = util::split(*data, data_delim);
        from_string = edge_data[0];
        k_string = edge_data[1];
        to_string = edge_data[2];

        k = atoi(k_string.c_str());
        to = nodes[atoi(to_string.c_str()) - 1];

        map<string, Source *>::iterator from_src = sources.find(from_string);
        if(from_src != sources.end()) {
          src = from_src->second;
          to->assignSource(src, k);
        } else {
          from = nodes[atoi(from_string.c_str()) - 1];
          to->assignAncestor(from, k);
        };
      };
    };
  file.close();
} // loadNetworkFile()

void Network::loadScalarsFile() {
  ifstream file(scalars_file_path);
  char scalar_delim = ',';
  string line;
  vector<string> scalar_data;

  if(!file.is_open()) {
    perror("Error While Opening Scalars File");
  }

  getline(file, line);
  scalar_data = util::split(line, scalar_delim);
  //epsilon_neg = stod(scalar_data[1], nullptr);
  //epsilon_pos = stod(scalar_data[3], nullptr);
  lambda = stod(scalar_data[5], nullptr);
  mu = stod(scalar_data[7], nullptr);

  file.close();
} // loadScalarsFile()

void Network::loadNodesFile(){
  ifstream file(node_file_path);
  int code = 1;
  string line;
  char delim = ',';
  vector<string> header;
  vector<string> node_data;
  bool firstline = true;

  num_reported_travel = 0;
  num_travel_history = 0;

  if (!file.is_open()) {
    perror("Error While Opening Nodes File");
  }
  while(getline(file, line)) {
    if(firstline){
      header = util::split(line, delim);
      firstline = false;
    } else {
      node_data = util::split(line, delim);
      loadNode(&node_data, &header, code);
      code++;
    }
  }
  file.close();

  prop_reported_travel = (1.0 * num_reported_travel) / (1.0 * num_travel_history);
} // loadNodesFile()


void Network::loadParametersFile(){
  ifstream file(parameters_file_path);
  std::string line;
  char delim = ',';
  char epsilon_delim = '_';
  std::vector<std::string> parameter_entry;

  if(!file.is_open()) {
    cout << parameters_file_path;
    perror(" Error While Opening Parameters File");
  }

  while(getline(file, line)) {
    parameter_entry = util::split(line, delim);

    if(parameter_entry.size() != 2) {
      throw invalid_argument("Parameter file malformed");
    }

    std::string parameter_label = parameter_entry[0];
    double parameter_value = std::stod(parameter_entry[1], nullptr);

    if (util::toLower(parameter_label) == "lambda") {
      lambda = parameter_value;
    } else if (util::toLower(parameter_label) == "mu") {
      mu = parameter_value;
    }
    else if (util::toLower(parameter_label) == "diffusion_coef") {
      diffusion_coef = parameter_value;
    }

    else if (util::toLower(parameter_label) == "tau_s") {
      tau_s = parameter_value;
    }

    else if (util::toLower(parameter_label) == "tau_l") {
      tau_l = parameter_value;
    }
    else {
      parameter_entry = util::split(parameter_label, epsilon_delim);

      if (parameter_entry.size() != 3) {
        throw invalid_argument("Parameter file malformed");
      }

      std::string locus = parameter_entry[2];
      std::string epsilon_type = util::toLower(parameter_entry[1]);

      if (epsilon_type == "pos") {
        epsilon_pos[locus] = parameter_value;
      } else if (epsilon_type == "neg") {
        epsilon_neg[locus] = parameter_value;
      } else {
        throw invalid_argument("Parameter file malformed");
      }
    }
  }
}



void Network::loadNode(vector<string> * node_data, vector<string> * header, int code){
  int index = 0;
  Node * node = new Node();
  node->assignCode(code);
  node->assignSettings(settings);
  if(node_data->size() != header->size()){
    throw invalid_argument("Node Data Size Inconsistent");
  }

  for(auto data = node_data->begin(); data != node_data->end(); data++){
    string label;
    label = header->at(index);

    if(label == "time" or label == "Time"){
      int node_time = stoi(*data, nullptr);
      node->setTime(node_time);
    }
    else if(label == "type" or label == "Type"){
      node->setType(*data);
    }
    else if(label == "lat" or label == "Lat"){
      if(*data == "NA")
      {
        //node->setLat(std::numeric_limits<double>::quiet_NaN());
        node->setLat(MISSING_VALUE);
        // node->setLat(0.0 / 0.0);
      }
      else
      {
        double node_lat = stod(*data, nullptr);
        node->setLat(node_lat);
      }
    }
    else if(label == "lon" or label == "Lon"){
      if(*data == "NA")
      {
        // node->setLon(std::numeric_limits<double>::quiet_NaN());
        node->setLon(MISSING_VALUE);
        // node->setLon(0.0 / 0.0);
      }
      else
      {
        double node_lon = stod(*data, nullptr);
        node->setLon(node_lon);
      }
    }
    else if(label == "history" or label == "History"){
      if(*data == "NA")
      {
        node->setHistory(MISSING_VALUE);
      }else{
        int node_history = stoi(*data, nullptr);
        if(node_history == 1)
        {
          num_reported_travel++;
        }
        node->setHistory(node_history);
        num_travel_history++;
      }
    }

    else{
      vector<bool> * alleles = new vector<bool>;
      for(std::string::iterator allele = data->begin(); allele != data->end(); allele++){
        if(*allele == '1'){
          alleles->push_back(1);
        }
        else if(*allele == '0'){
          alleles->push_back(0);
        }
      }
      node->setLocus(label, alleles);
    }
    index++;
  }
  node->checkForNAs();
  node->setLocusOrder(locus_order);
  node->setAlleleCount(alleleCount);
  node->setKStride(k_stride);
  node->resetLogLik();
  nodes.push_back(node);
} // loadNode()



void Network::loadProposalProbs()
{
  proposal_probs.push_back(settings["PROP_EPSILON_NEG"]);
  proposal_probs.push_back(settings["PROP_EPSILON_POS"]);
  proposal_probs.push_back(settings["PROP_LAMBDA"]);
  proposal_probs.push_back(settings["PROP_MU"]);
  proposal_probs.push_back(settings["PROP_K"]);
  proposal_probs.push_back(settings["PROP_SWITCH"]);
  proposal_probs.push_back(settings["PROP_ADD"]);
  proposal_probs.push_back(settings["PROP_FOUNDER"]);
  proposal_probs.push_back(settings["PROP_ORPHAN"]);
  proposal_probs.push_back(settings["PROP_TAU_S"]);
  proposal_probs.push_back(settings["PROP_TAU_L"]);
  proposal_probs.push_back(settings["PROP_DIFFUSION_COEF"]);

} // loadProposalProbs()



void Network::loadSourcesFile() {
  ifstream file(sources_file_path);
  char delim = ',';
  string line;
  vector<string> source_data;
  max_alleles = 0;

  std::set<std::string> locus_set;

  if (!file.is_open()) {
    perror("Error While Opening File");
  }
  while(getline(file, line)) {
    source_data = util::split(line, delim);
    if(source_data.size() > 2){
      vector<double> * freqs = new vector<double>;
      string source_name = source_data[0];
      string locus = source_data[1];
      for(
        vector<string>::iterator freq_str = source_data.begin() + 2;
        freq_str != source_data.end() && freq_str->length() != 0;
        freq_str++){
          double f = stod(*freq_str, nullptr);
          freqs->push_back(f);
      }
      alleleCount[locus] = freqs->size();
      locus_set.insert(locus);
      loadSource(source_name, locus, freqs);
      if(freqs->size() > max_alleles){
        max_alleles = freqs->size();
      }
    }
  }
  file.close();

  locus_order = std::vector<string>(locus_set.begin(), locus_set.end());
  std::vector<string>::iterator locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++){
    k_stride += alleleCount[*locItr];
  }

  map<string, Source *>::iterator src_itr;
  Source * src;

  for(src_itr = sources.begin(); src_itr != sources.end(); src_itr++) {
    src = src_itr->second;
    src->setLocusOrder(locus_order);
    src->setAlleleCount(alleleCount);
    src->setKStride(k_stride);
  }

} // loadSourcesFile



void Network::loadSource(string label, string locus, vector<double> * frequencies) {
  Source * source;
  // string * source_label = new string(label);
  // string * locus_label = new string(locus);
  if(sources.count(label) > 0) {
    source = sources[label];
  }
  else {
    source = new Source(settings);
    sources[label] = source;
    if(label != LOCAL_SOURCE){
     nonLocalSources++;
    }
    numSources++;
  }
  source->setName(label);
  source->setLocusFreqs(locus, frequencies);
} // loadSource



double Network::logPrior_diffusion(double diffusion_in)
{
  return 0.0; // log(gsl_ran_gamma_pdf(diffusion_in - settings["DIFFUSION_min"], settings["DIFFUSION_a"], settings["DIFFUSION_b"]));
} // logPrior_lambda()



double Network::logPrior_lambda(double lambda_in)
{
  return log(gsl_ran_beta_pdf(lambda_in, settings["LAMBDA_a"], settings["LAMBDA_b"]));
} // logPrior_lambda()



double Network::logPrior_mu(double mu_in)
{
  return log(gsl_ran_beta_pdf(mu_in, settings["MU_a"], settings["MU_b"]));
} // logPrior_mu()



double Network::logPrior_tau_s(double tau_s_in)
{
  return log(gsl_ran_beta_pdf(tau_s_in, settings["TAU_S_a"], settings["TAU_S_b"]));
}// Network::logPrior_tau_s()



double Network::logPrior_tau_l(double tau_l_in)
{
  return log(gsl_ran_beta_pdf(tau_l_in, settings["TAU_L_a"], settings["TAU_L_b"]));
}// Network::logPrior_tau_l()



int Network::getRandom_k()
{
  int k = 0;
  double probSumK = 0.0;
  double runif = gsl_rng_uniform_pos(this->rGsl);
  do{
    k++;
    probSumK += exp(logProposal_k(k));
  } while(probSumK < runif && k < settings["MAX_K"]);
  return(k);
}


double Network::logProposal_k(int k)
{
  return(log_proposal_k[k-1]);
}


double Network::logProposal_diffusion(double diffusion_old, double diffusion_new, double proposal_width)
{
	return log(gsl_ran_gaussian_pdf(diffusion_new - diffusion_old, proposal_width) / (1.0 - gsl_cdf_gaussian_P(-diffusion_old, proposal_width)));
}



double Network::logProposal_tau_l(double tau_l_old, double tau_l_new, double proposal_width)
{
	return log(gsl_ran_gaussian_pdf(tau_l_new - tau_l_old, proposal_width) / (gsl_cdf_gaussian_P(1.0 - tau_l_old, proposal_width) - gsl_cdf_gaussian_P(-tau_l_old, proposal_width)));
}



double Network::logProposal_tau_s(double tau_s_old, double tau_s_new, double proposal_width)
{
	return log(gsl_ran_gaussian_pdf(tau_s_new - tau_s_old, proposal_width) / (gsl_cdf_gaussian_P(1.0 - tau_s_old, proposal_width) - gsl_cdf_gaussian_P(-tau_s_old, proposal_width)));
}



double Network::make_founder(Node * nodeIn)
{
  double logPropDiff = 0.0, increment = 0.0;
  int kOld = 1;

  if(nodeIn->isOrphan()){
    logPropDiff += log(settings["PROP_ORPHAN"]) - log(settings["PROP_FOUNDER"] / double(nonLocalSources));
    nodeIn->removeSource();
    nodeIn->assignSource(chooseSource(), 1);
  }
  else{
    for(
      vector<pair<Node *,int> >::iterator itrAnc = nodeIn->getAncestors()->begin();
      itrAnc != nodeIn->getAncestors()->end();
      itrAnc++)
    {
      kOld = (*itrAnc).second;
      logPropDiff += log(nodeIn->getEdgeProbs()->find((*itrAnc).first)->second) + logProposal_k(kOld);
      logPropDiff += log(double(nodeIn->getAncestors()->size()) - (increment++));
    }

    logPropDiff += double(nodeIn->getAncestors()->size()) * log(settings["PROP_ADD"]) - log(settings["PROP_FOUNDER"] / double(nonLocalSources));

    nodeIn->removeAllAncestors();
    if(nodeIn->getAncestors()->empty()){
      nodeIn->assignSource(chooseSource(), 1);
    }
  }

  return logPropDiff;
} // make_founder()



double Network::make_orphan(Node * nodeIn)
{
  double logPropDiff = 0.0, increment = 0.0;

  if(nodeIn->isFounder()){
    logPropDiff += log(settings["PROP_FOUNDER"] / double(nonLocalSources)) - log(settings["PROP_ORPHAN"]);
    nodeIn->removeSource();
    nodeIn->assignSource(sources[LOCAL_SOURCE], 1);
  }
  else{
    for(
      vector<pair<Node *,int> >::iterator itrAnc = nodeIn->getAncestors()->begin();
      itrAnc != nodeIn->getAncestors()->end();
      itrAnc++)
    {
      logPropDiff += log(nodeIn->getEdgeProbs()->find(nodeIn->getAncestors()->back().first)->second);
      logPropDiff += log(double(nodeIn->getAncestors()->size()) - (increment++));
    }

    logPropDiff += double(nodeIn->getAncestors()->size()) * log(settings["PROP_ADD"]) - log(settings["PROP_ORPHAN"]);

    nodeIn->removeAllAncestors();
    if(nodeIn->getAncestors()->empty()){
      nodeIn->assignSource(sources[LOCAL_SOURCE], 1);
    }
  }

  return logPropDiff;
} // make_orphan()



void Network::randomTree()
{
  // declare variables
  vector<Node *> unassignedNodes = nodes;
  vector<Node *> ancestors;
  vector<Node *>::iterator itrAnc;
  int a, d, v, V;
  for(itrAnc = unassignedNodes.begin(); itrAnc != unassignedNodes.end(); itrAnc++){;}

  // assign the root ancestor
  a = gsl_rng_uniform_int(this->rGsl, unassignedNodes.size());
  unassignedNodes[a]->assignSource(chooseSource(), 1);

  ancestors.push_back(unassignedNodes[a]);
  unassignedNodes.erase(unassignedNodes.begin() + a);
  // while there are still unassigned nodes
  while(!unassignedNodes.empty()){
    // randomly pick from the pool of ancestors
    a = gsl_rng_uniform_int(this->rGsl, ancestors.size());

    // for a Poisson-distributed number of descendants
    V = gsl_ran_poisson(this->rGsl, settings["R0_initial"]) + 1;
    for(v = 0; v < V; v++){
      if(!unassignedNodes.empty()){
        itrAnc = ancestors.begin() + a;
        // select an unassigned node to be made a descendant
        d = gsl_rng_uniform_int(this->rGsl, unassignedNodes.size());

        // assign the current ancestor to the descendant
        unassignedNodes[d]->assignAncestor(*itrAnc, 1);

        // place the descendant in the pool of ancestors
        ancestors.push_back(unassignedNodes[d]);

        // remove the descendant from the unassigned pool
        unassignedNodes.erase(unassignedNodes.begin() + d);
      }
    }

    // remove the current ancestor from the pool of ancestors
    ancestors.erase(ancestors.begin() + a);
  }
} // randomTree()



void Network::copyProbabilities() {
  Pr_g_obs_offspring_copy = Pr_g_obs_offspring;
  Pr_g_inheritance_copy = Pr_g_inheritance;
  Pr_g_true_parent_copy = Pr_g_true_parent;
  // Pr_XXXX_copy = Pr_XXXX;
  // Pr_XXOX_copy = Pr_XXOX;
  // Pr_XOXX_copy = Pr_XOXX;
  // Pr_XOOX_copy = Pr_XOOX;
  // Pr_XOX_copy = Pr_XOX;
  // Pr_OXX_copy = Pr_OXX;
  // Pr_OOX_copy = Pr_OOX;
  Pr_vec_copy = Pr_vec;

  map<string,Source *>::iterator sourceItr = sources.begin();
  for(; sourceItr != sources.end(); sourceItr++) {
    sourceItr->second->copyProbabilities();
  }
} // copyProbabilities()



void Network::revertProbabilities() {
  Pr_g_obs_offspring = Pr_g_obs_offspring_copy;
  Pr_g_inheritance = Pr_g_inheritance_copy;
  Pr_g_true_parent = Pr_g_true_parent_copy;
  // Pr_XXXX = Pr_XXXX_copy;
  // Pr_XXOX = Pr_XXOX_copy;
  // Pr_XOXX = Pr_XOXX_copy;
  // Pr_XOOX = Pr_XOOX_copy;
  // Pr_XOX = Pr_XOX_copy;
  // Pr_OXX = Pr_OXX_copy;
  // Pr_OOX = Pr_OOX_copy;
  Pr_vec = Pr_vec_copy;

  map<string,Source *>::iterator sourceItr = sources.begin();
  for(; sourceItr != sources.end(); sourceItr++) {
    sourceItr->second->revertProbabilities();
  }
} // revertProbabilities()



void Network::recalculateProbs(){
  setPrGObsOffspring();
  setPrGInheritance();
  setPrGTrueParent();
  setPrVecs();
  map<string,Source *>::iterator sourceItr = sources.begin();
  for(; sourceItr != sources.end(); sourceItr++) {
    sourceItr->second->setPrVecs(&Pr_g_obs_offspring, &Pr_g_inheritance);
  }
} // recalculateProbs()



void Network::recalculateProbs(const std::string &locus){
  setPrGObsOffspring(locus);
  setPrGInheritance(locus);
  setPrGTrueParent(locus);
  setPrVecs(locus);

  map<string,Source *>::iterator sourceItr = sources.begin();
  for(; sourceItr != sources.end(); sourceItr++) {
    sourceItr->second->setPrVecs(Pr_g_obs_offspring, Pr_g_inheritance, locus);
  }
} // recalculateProbs()



double Network::remove_anc(Node * nodeIn){
  double logPropDiff =
    log(nodeIn->getEdgeProbs()->find(nodeIn->getAncestors()->back().first)->second) +
    log(double(nodeIn->getAncestors()->size())) +
    log(settings["PROP_ADD"]) - log(settings["PROP_REMOVE"]);

  double probSumMax =
    nodeIn->getEdgeProbs()->find(nodeIn->getAncestors()->back().first)->second;

  int toRemove = gsl_rng_uniform_int(this->rGsl, nodeIn->getAncestors()->size());
  nodeIn->removeAncestor(toRemove);
  if(nodeIn->getAncestors()->empty())
    nodeIn->assignSource(chooseSource(), 1);

  map<Node *, double>::iterator mapItr = nodeIn->getEdgeProbs()->begin();
  for(; mapItr != nodeIn->getEdgeProbs()->end(); mapItr++){
    if(!nodeIn->isAncestor(mapItr->first)){
      probSumMax += mapItr->second;
    }
  }

  logPropDiff -= log(probSumMax);

  return logPropDiff;
} // remove_anc()



double Network::runif_interval(double start, double width)
{
  double unif_min = fmax(start - width, 0.0);
  double unif_max = fmin(start + width, 1.0);
  double random_point = gsl_rng_uniform(this->rGsl);

  return(unif_min + (unif_max - unif_min) * random_point);
} // runif_interval()



void Network::setBeta(double beta) {
  this->beta = beta;
} // setBeta()



void Network::simulate()
{
  logLik = calcLogLik();

  Node nodeSave;

  for(int i = 0; i < settings["SWAP_FREQUENCY"]; i++){
    update(&nodeSave);
  }

  if(beta == 1.0 && WRITE_STATUS == 1){
    writeScalars();
    writeNetwork();
  }
} // simulate()



double Network::switch_anc_anc(Node * nodeNow)
{
  double probSum = 0.0, probSumMax = 0.0, runif = 0.0, logPropDiff = 0.0;
  map<Node *, double>::iterator mapItr = nodeNow->getEdgeProbs()->begin();

  // get the probability associated with adding a new node as an ancestor
  // probSum = probSumMax = 0.0;
  mapItr = nodeNow->getEdgeProbs()->begin();

  for(; mapItr != nodeNow->getEdgeProbs()->end(); mapItr++){
    if(!nodeNow->isAncestor(mapItr->first)){
      probSumMax += mapItr->second;
    }
  }

  // pick a random node to add as an ancestor and add it
  runif = gsl_rng_uniform_pos(this->rGsl) * probSumMax;
  mapItr = nodeNow->getEdgeProbs()->begin();

  while(probSum < runif && mapItr != nodeNow->getEdgeProbs()->end()){
    if(nodeNow->isAncestor(mapItr->first)){
      mapItr++;
    }
    else{
      probSum += mapItr->second;
      mapItr++;
    }
  }
  mapItr--;

  int kNew = 1;
  if(settings["PROP_K"]>0){
    kNew = getRandom_k();
  }

  nodeNow->assignAncestor(mapItr->first, kNew);
  logPropDiff -= log(mapItr->second / probSumMax) + logProposal_k(kNew);

  // pick an ancestor to remove and remove it
  int toRemove = gsl_rng_uniform_int(this->rGsl, nodeNow->getAncestors()->size());
  Node * nodeRemoved = nodeNow->getAncestors()->at(toRemove).first;
  int kOld = nodeNow->getAncestors()->at(toRemove).second;
  nodeNow->removeAncestor(toRemove);

  // get the probability associated with adding a new node as an ancestor
  probSumMax = 0.0;
  mapItr = nodeNow->getEdgeProbs()->begin();
  for(; mapItr != nodeNow->getEdgeProbs()->end(); mapItr++){
    if(!nodeNow->isAncestor(mapItr->first)){
      probSumMax += mapItr->second;
    }
  }
  logPropDiff += log((*nodeNow->getEdgeProbs())[nodeRemoved] / probSumMax) + logProposal_k(kOld);

  return logPropDiff;
} // switch_anc_anc()



void Network::update(Node * nodeSave)
{
  bool updated;
  int numUpdates;
  int whichUpdate;
  double probSum, runif;

  updated = 0;
  numUpdates = 12;
  while(!updated){
    nodeSave->resetAll();

    whichUpdate = 0;
    probSum = 0.0;
    runif = gsl_rng_uniform_pos(rGsl);

    while(probSum < runif && whichUpdate < numUpdates){
      probSum += proposal_probs[whichUpdate];
      whichUpdate++;
    }
    --whichUpdate;

    switch(whichUpdate){
      case 0:
        update_epsilon_neg();
        updated = 1;
        break;
      case 1:
        update_epsilon_pos();
        updated = 1;
        break;
      case 2:
        update_lambda();
        updated = 1;
        break;
      case 3:
        update_mu();
        updated = 1;
        break;
      case 4:
        update_k(nodeSave);
        updated = 1;
        break;
      case 5:
        updated = update_switch_anc_anc(nodeSave);
        break;
      case 6:
        updated = update_add_anc(nodeSave);
        break;
      case 7:
        updated = update_make_founder(nodeSave);
        break;
      case 8:
        updated = update_make_orphan(nodeSave);
        break;
      case 9:
        update_tau_s();
        updated = 1;
        break;
      case 10:
        update_tau_l();
        updated = 1;
        break;
      case 11:
        update_diffusion_coef();
        updated = 1;
        break;
    }
  }

} // update()



bool Network::update_add_anc(Node * nodeSave)
{
  bool updated = 0;
  double logPropDiff, logLikNew;
  int randNode = gsl_rng_uniform_int(this->rGsl, nodes.size());

  if(nodes[randNode]->isFounder() && (settings["BELIEVE_TRAVELHIST"] > 0.0) && (nodes[randNode]->getHistory() == 1)){
    return updated;
  }

  if(!nodes[randNode]->isFounder() && !(settings["SUPERINFECTION_TOGGLE"] > 0.0)){
    return updated;
  }

  if(nodes[randNode]->getAncestors()->size() >= 1){
    return updated;
  }

  nodeSave->copyNode(nodes[randNode]);

  logPropDiff = add_anc(nodes[randNode]);

  if(logPropDiff < 100){
    if(!isCyclic(nodes[randNode])){
      logLikNew = updateLogLik(nodes[randNode]);

      if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik) + logPropDiff), 1.0)){
        nodes[randNode]->copyNode(nodeSave);
        add_anc_reject++;
      }
      else{
        logLik = logLikNew;
        add_anc_accept++;
      }
      updated = 1;
    }
    else{
      nodes[randNode]->copyNode(nodeSave);
    }
  }

  return updated;
} // update_add_anc()



void Network::validateLoci(){
  std::vector<std::string>::iterator locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++){
    validLoci[*locItr] = true;
  }
}



void Network::invalidateLoci(){
  std::vector<std::string>::iterator locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++){
    validLoci[*locItr] = false;
  }
}



void Network::validateLocus(const std::string &locus){
  validLoci[locus] = true;
}



void Network::invalidateLocus(const std::string &locus){
  validLoci[locus] = false;
}



void Network::update_epsilon_neg()
{
  int whichLocus = gsl_rng_uniform_int(this->rGsl, locus_order.size());
  std::string locus = locus_order[whichLocus];
  map<string, double> epsilon_neg_old = epsilon_neg;
  map<string, double> epsilon_neg_new = change_epsilon_neg(locus);
  double logLikNew;
  double logPriorDiff;
  double logLikOld;

  epsilon_neg = epsilon_neg_new;
  copyProbabilities();
  recalculateProbs(locus);
  logLikOld = logLik;
  invalidateLocus(locus);
  logLikNew = calcLogLik();
  logPriorDiff = 0.0;

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik + logPriorDiff)), 1.0)){
    epsilon_neg = epsilon_neg_old;
    revertProbabilities();
    invalidateLocus(locus);
    logLik = logLikOld;
    epsilon_neg_reject++;
  }
  else{
    logLik = logLikNew;
    epsilon_neg_accept++;
  }
} // update_epsilon_neg()



void Network::update_epsilon_pos()
{
  int whichLocus = gsl_rng_uniform_int(this->rGsl, locus_order.size());
  std::string locus = locus_order[whichLocus];
  map<string, double> epsilon_pos_old = epsilon_pos;
  map<string, double> epsilon_pos_new = change_epsilon_pos(locus);
  double logLikNew;
  double logPriorDiff;
  double logLikOld;

  epsilon_pos = epsilon_pos_new;
  copyProbabilities();
  recalculateProbs(locus);
  logLikOld = logLik;
  invalidateLocus(locus);
  logLikNew = calcLogLik();
  logPriorDiff = 0.0;

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik + logPriorDiff)), 1.0)){
    epsilon_pos = epsilon_pos_old;
    revertProbabilities();
    invalidateLocus(locus);
    logLik = logLikOld;
    epsilon_pos_reject++;
  }
  else{
    logLik = logLikNew;
    epsilon_pos_accept++;
  }
} // update_epsilon_pos()



void Network::update_k(Node * nodeSave)
{
  double logPropDiff = 0.0;
  int randNode = 0;
  int cntRngShots = 0; // required to prevent infinit loop in case all nodes are founders
  do{
    randNode = (int)gsl_rng_uniform_int(this->rGsl, nodes.size());
      cntRngShots++;
  }while(nodes[randNode]->getAncestors()->empty() & (cntRngShots<=nodes.size()));

  nodeSave->copyNode(nodes[randNode]);

  if(!nodes[randNode]->getAncestors()->empty()){
      logPropDiff = change_k_anc(nodes[randNode]);
  }
  //else if(!nodes[randNode]->getSources()->empty()){
  //  change_k_src(nodes[randNode]);
  //}

  double logLikNew = updateLogLik(nodes[randNode]);

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik) + logPropDiff), 1.0)){
    nodes[randNode]->copyNode(nodeSave);
    change_k_reject++;
  }
  else{
    logLik = logLikNew;
    change_k_accept++;
  }
} // update_k()


double Network::updateLogLik(Node * nodeIn)
{
  double logLik_new = logLik;
  logLik_new -= nodeIn->getLogLik();

  nodeIn->resetLogLik();
  nodeIn->clearCache();
  invalidateLoci();
  nodeIn->calcLogLik(Pr_vec, *sources[LOCAL_SOURCE], validLoci, diffusion_coef, tau_s, tau_l, prop_reported_travel);
  validateLoci();
  logLik_new += nodeIn->getLogLik();

  return logLik_new;
} // updateLogLik()



bool Network::update_make_founder(Node * nodeSave)
{
  bool updated = 0;
  double logPropDiff, logLikNew;
  int randNode = gsl_rng_uniform_int(this->rGsl, nodes.size());

  if(settings["BELIEVE_TRAVELHIST"] > 0.0){
    if((nodes[randNode]->getHistory() == 0)){
      return updated;
    }
    else{
      if(nodes[randNode]->isFounder()){
        return updated;
      }
      else{
        nodes[randNode]->removeAllAncestors();
        nodes[randNode]->assignSource(chooseSource(), 1);
        updated = 1;
      }
    }
  }
  else{
    nodeSave->copyNode(nodes[randNode]);

    if(!nodes[randNode]->getAncestors()->empty()){
      logPropDiff = make_founder(nodes[randNode]);

      if(!isCyclic(nodes[randNode])){
        logLikNew = updateLogLik(nodes[randNode]);

        if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik) + logPropDiff), 1)){
          nodes[randNode]->copyNode(nodeSave);
          make_founder_reject++;
        }
        else{
          logLik = logLikNew;
          make_founder_accept++;
        }
        updated = 1;
      }
      else{
        nodes[randNode]->copyNode(nodeSave);
      }
    }
  }

  return updated;
} // update_make_founder()



bool Network::update_make_orphan(Node * nodeSave)
{
  bool updated = 0;
  double logPropDiff, logLikNew;
  int randNode = gsl_rng_uniform_int(rGsl, nodes.size());

  if(nodes[randNode]->isOrphan()){
    return updated;
  }

  if(nodes[randNode]->isFounder() && settings["BELIEVE_TRAVELHIST"] > 0.0){
    return updated;
  }

  logPropDiff = make_orphan(nodes[randNode]);

  if(!isCyclic(nodes[randNode])){
    logLikNew = updateLogLik(nodes[randNode]);

    if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik) + logPropDiff), 1)){
      nodes[randNode]->copyNode(nodeSave);
      make_orphan_reject++;
    }
    else{
      logLik = logLikNew;
      make_orphan_accept++;
    }
    updated = 1;
  }
  else{
    nodes[randNode]->copyNode(nodeSave);
  }

  return updated;
} // update_make_orphan()


void Network::update_lambda()
{
  double lambda_old = lambda;
  double lambda_new = change_lambda();
  double logLikNew;
  double logPriorDiff;
  double logLikOld;

  copyProbabilities();

  lambda = lambda_new;
  logLikOld = logLik;
  recalculateProbs();
  invalidateLoci();
  logLikNew = calcLogLik();

  logPriorDiff = logPrior_lambda(lambda_new) - logPrior_lambda(lambda_old);

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik + logPriorDiff)), 1.0)){
    lambda = lambda_old;
    invalidateLoci();
    revertProbabilities();
    logLik = logLikOld;
    lambda_reject++;
  }
  else{
    logLik = logLikNew;
    lambda_accept++;
  }
} // update_lambda()


void Network::update_mu()
{
  double mu_old = mu;
  double mu_new = change_mu();
  double logLikNew;
  double logPriorDiff;
  double logLikOld;

  copyProbabilities();

  mu = mu_new;
  logLikOld = logLik;
  recalculateProbs();
  invalidateLoci();
  logLikNew = calcLogLik();

  logPriorDiff = logPrior_mu(mu_new) - logPrior_mu(mu_old);

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik + logPriorDiff)), 1.0)){
    mu = mu_old;
    revertProbabilities();
    invalidateLoci();
    logLik = logLikOld;
    mu_reject++;
  }
  else{
    logLik = logLikNew;
    mu_accept++;
  }
} // update_mu()


void Network::update_diffusion_coef()
{
  double
    diffusion_coef_old = diffusion_coef,
    diffusion_coef_new = change_diffusion_coef(),
    logLikNew,
    logLikOld,
    logPriorDiff,
    logProposalDiff;

  diffusion_coef = diffusion_coef_new;
  logLikOld = logLik;

  if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
  {
    copyProbabilities();
    recalculateProbs();
  }

  logLikNew = calcLogLik();
  logPriorDiff = logPrior_diffusion(diffusion_coef_new) - logPrior_diffusion(diffusion_coef_old);
  logProposalDiff = logProposal_diffusion(diffusion_coef_old, diffusion_coef_new, 2.5) - logProposal_diffusion(diffusion_coef_new, diffusion_coef_old, 2.5);

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik + logPriorDiff) - logProposalDiff), 1.0)){
    diffusion_coef = diffusion_coef_old;
    if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
    {
      revertProbabilities();
    }
    logLik = logLikOld;
    diffusion_coef_reject++;
  }
  else{
    logLik = logLikNew;
    diffusion_coef_accept++;
  }
}


void Network::update_tau_l()
{
  double tau_l_old = tau_l,
  tau_l_new = change_tau_l(),
  logLikNew,
  logLikOld,
  logPriorDiff,
  logProposalDiff;

  logLikOld = logLik;
  tau_l = tau_l_new;

  if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
  {
    copyProbabilities();
    recalculateProbs();
  }

  logLikNew = calcLogLik();
  logPriorDiff = logPrior_tau_l(tau_l_new) - logPrior_tau_l(tau_l_old);
  // logProposalDiff = logProposal_tau_l(tau_l_old, tau_l_new, 0.000001) - logProposal_tau_l(tau_l_new, tau_l_old, 0.000001);
  logProposalDiff = logProposal_tau_l(tau_l_old, tau_l_new, 0.25) - logProposal_tau_l(tau_l_new, tau_l_old, 0.25);

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik + logPriorDiff) - logProposalDiff), 1.0))
  {
    tau_l = tau_l_old;
    if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
    {
      revertProbabilities();
    }
    logLik = logLikOld;
    tau_l_reject++;
  }
  else
  {
    logLik = logLikNew;
    tau_l_accept++;
  }
}// Network::update_tau_l()



void Network::update_tau_s()
{
  double tau_s_old = tau_s,
  tau_s_new = change_tau_s(),
  logLikNew,
  logLikOld,
  logPriorDiff,
  logProposalDiff;

  logLikOld = logLik;
  tau_s = tau_s_new;

  if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
  {
    copyProbabilities();
    recalculateProbs();
  }

  logLikNew = calcLogLik();
  logPriorDiff = logPrior_tau_s(tau_s_new) - logPrior_tau_s(tau_s_old);
  // logProposalDiff = logProposal_tau_s(tau_s_old, tau_s_new, 0.000001) - logProposal_tau_s(tau_s_new, tau_s_old, 0.000001);
  logProposalDiff = logProposal_tau_s(tau_s_old, tau_s_new, 0.25) - logProposal_tau_s(tau_s_new, tau_s_old, 0.25);

  if(gsl_rng_uniform(this->rGsl) > fmin(exp(beta * (logLikNew - logLik + logPriorDiff) - logProposalDiff), 1.0))
  {
    tau_s = tau_s_old;
    if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
    {
      revertProbabilities();
    }
    logLik = logLikOld;
    tau_s_reject++;
  }
  else
  {
    logLik = logLikNew;
    tau_s_accept++;
  }
}// Network::update_tau_s()



bool Network::update_remove_anc(Node * nodeSave)
{
  bool updated = 0;

  int randNode = gsl_rng_uniform_int(this->rGsl, nodes.size());

  nodeSave->copyNode(nodes[randNode]);

  if(!nodes[randNode]->getAncestors()->empty()){
    double logPropDiff = remove_anc(nodes[randNode]);

    if(!isCyclic(nodes[randNode])){
      double logLikNew = updateLogLik(nodes[randNode]);
      double random_float = gsl_rng_uniform(this->rGsl);
      if(random_float > fmin(exp(beta * (logLikNew - logLik) + logPropDiff), 1.0)){
        nodes[randNode]->copyNode(nodeSave);
        remove_anc_reject++;
      }
      else{
        logLik = logLikNew;
        remove_anc_accept++;
      }
      updated = 1;
    }
    else{
      nodes[randNode]->copyNode(nodeSave);
    }
  }

  return updated;
} // update_remove_anc()



bool Network::update_switch_anc_anc(Node * nodeSave)
{
  bool updated = 0;

  int randNode = gsl_rng_uniform_int(this->rGsl, nodes.size());

  nodeSave->copyNode(nodes[randNode]);

  if(!nodes[randNode]->getAncestors()->empty()){
    double logPropDiff = switch_anc_anc(nodes[randNode]);

    if(logPropDiff < 100){
      if(!isCyclic(nodes[randNode])){
        double logLikNew = updateLogLik(nodes[randNode]);
        double random_float = gsl_rng_uniform(this->rGsl);
        if(random_float > fmin(exp(beta * (logLikNew - logLik) + logPropDiff), 1.0)){
          nodes[randNode]->copyNode(nodeSave);
          switch_anc_anc_reject++;
        }
        else{
          logLik = logLikNew;
          switch_anc_anc_accept++;
        }
        updated = 1;
      }
      else{
        nodes[randNode]->copyNode(nodeSave);
      }
    }
    else{
      nodes[randNode]->copyNode(nodeSave);
    }
  }

  return updated;
} // update_switch_anc_anc()



void Network::writeHeader()
{
  {
    ofstream fileScalars, fileNetwork;
    fileScalars.open(scalars_destination);

    fileScalars <<
      "logLik,";

    std::vector<std::string>::iterator locItr;

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++)
      fileScalars << "epsilon_neg_" << *locItr << ",";

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++)
      fileScalars << "epsilon_pos_" << *locItr << ",";

    fileScalars <<
      "lambda," <<
      "mu," <<
      "diffusion_coef," <<
      "tau_s," <<
      "tau_l," <<
      "epsilon_neg_accept," <<
      "epsilon_pos_accept," <<
      "lambda_accept," <<
      "mu_accept," <<
      "diffusion_coef_accept," <<
      "tau_s_accept," <<
      "tau_l_accept," <<
      "change_k_accept," <<
      "switch_anc_accept," <<
      "add_anc_accept," <<
      "make_founder_accept," <<
      "make_orphan_accept" <<
      endl;

    fileScalars.close();
    fileNetwork.open(network_destination);
    fileNetwork.close();
  }
} // writeHeader()



void Network::writeNetwork() {
  ofstream file;
  file.open(network_destination, ios::app);

  if (file.is_open()) {
    vector<Node *>::iterator itrNode;
    vector<pair<Node *,int> >::iterator itrAnc;
    vector<pair<Source *,int> >::iterator itrSrc;

    for(itrNode = nodes.begin(); itrNode != nodes.end(); itrNode++){
      if(!(*itrNode)->getAncestors()->empty()){
        for(itrAnc = (*itrNode)->getAncestors()->begin();
            itrAnc != (*itrNode)->getAncestors()->end();
            itrAnc++)
          file <<
            itrAnc->first->getCode() << "-" <<
            itrAnc->second << "-" <<
            (*itrNode)->getCode() << ";";
      }
      if(!(*itrNode)->getSources()->empty()){
        for(itrSrc = (*itrNode)->getSources()->begin();
            itrSrc != (*itrNode)->getSources()->end();
            itrSrc++)
          file <<
            itrSrc->first->getName() << "-" <<
            itrSrc->second << "-" <<
            (*itrNode)->getCode() << ";";
      }
    }
    file << endl;

    file.close();
  }
} // writeNetwork()



void Network::writeScalars(){
  ofstream file;
  file.open(scalars_destination, ios::app);

  file <<
    logLik << ",";

  std::vector<std::string>::iterator locItr;
  locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++)
    file << epsilon_neg[*locItr] << ",";

  locItr = locus_order.begin();
  for(; locItr != locus_order.end(); locItr++)
    file << epsilon_pos[*locItr] << ",";

  file <<
    lambda << "," <<
    mu << "," <<
    diffusion_coef << "," <<
    tau_s << "," <<
    tau_l << "," <<
    double(epsilon_neg_accept) / double(epsilon_neg_accept + epsilon_neg_reject) << "," <<
    double(epsilon_pos_accept) / double(epsilon_pos_accept + epsilon_pos_reject) << "," <<
    double(lambda_accept) / double(lambda_reject + lambda_accept) << "," <<
    double(mu_accept) / double(mu_accept + mu_reject) << "," <<
    double(diffusion_coef_accept) / double(diffusion_coef_accept + diffusion_coef_reject) << "," <<
    double(tau_s_accept) / double(tau_s_accept + tau_s_reject) << "," <<
    double(tau_l_accept) / double(tau_l_accept + tau_l_reject) << "," <<
    double(change_k_accept) / double(change_k_accept + change_k_reject) << "," <<
    double(switch_anc_anc_accept) / double(switch_anc_anc_accept + switch_anc_anc_reject) << "," <<
    double(add_anc_accept) / double(add_anc_accept + add_anc_reject) << "," <<
    double(make_founder_accept) / double(make_founder_accept + make_founder_reject) << "," <<
    double(make_orphan_accept) / double(make_orphan_accept + make_orphan_reject) <<
    endl;

    epsilon_neg_accept = epsilon_neg_reject = 0;
    epsilon_pos_accept = epsilon_pos_reject = 0;
    lambda_accept = lambda_reject = 0;
    mu_accept = mu_reject = 0;
    diffusion_coef_accept = diffusion_coef_reject = 0;
    tau_s_accept = tau_s_reject = 0;
    tau_l_accept = tau_l_reject = 0;
    change_k_accept = change_k_reject = 0;
    switch_anc_anc_accept = switch_anc_anc_reject = 0;
    add_anc_accept = add_anc_reject = 0;
    remove_anc_accept = remove_anc_reject = 0;
    make_founder_accept = make_founder_reject = 0;
    make_orphan_accept = make_orphan_reject = 0;

  file.close();
} // writeScalars()
