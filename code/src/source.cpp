#include "source.h"



void Source::setLocusFreqs(string locus, vector<double> * freqs){
  frequencies[locus] = freqs;
} // setLocusFreqs()



void Source::copyProbabilities(){
  Pr_vec_copy = Pr_vec;
} // copyProbabilities()



void Source::revertProbabilities() {
  Pr_vec = Pr_vec_copy;
} // revertProbabilities()



void Source::initPrVecs() {
  std::vector<std::string>::iterator locItr;
  vector<double>::iterator alleleItr;
  string locus;
  int s_i = 0;
  Pr_vec.resize((settings["MAX_K"]) * 5 * k_stride, 0);
  Pr_vec_copy.resize((settings["MAX_K"]) * 5 * k_stride, 0);
} // initPrVecs()



void Source::setPrVecs(map<string, vector<double> >* Pr_g_obs_offspring, vector<map<string, vector<double> > >* Pr_g_inheritance){
  std::vector<std::string>::iterator locItr;
  string locusLabel;
  int s_i = 0;
  double freq_true;
  for(int k = 0; k < settings["MAX_K"]; k++){
    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++){

      locusLabel = *locItr;
      for(int a = 0; a < alleleCount[locusLabel]; a++){

        freq_true =
          ((*frequencies[locusLabel])[a] - (*Pr_g_obs_offspring)[locusLabel][2]) /
          (1.0 - (*Pr_g_obs_offspring)[locusLabel][1] - (*Pr_g_obs_offspring)[locusLabel][2]);
        // XXS = unk | unk | src
        // 0 = 0S
        // 0
        Pr_vec[s_i] =
          log((*Pr_g_obs_offspring).at(locusLabel)[1] * freq_true +
          (*Pr_g_obs_offspring).at(locusLabel)[0] * (1.0 - freq_true));
        // 1 = 1S

        // 1
        Pr_vec[s_i + 1] =
          log((*Pr_g_obs_offspring).at(locusLabel)[3] * freq_true +
          (*Pr_g_obs_offspring).at(locusLabel)[2] * (1.0 - freq_true));

        // XOS = unk | 0 | src
        // 0 = 0S

        // 2
        Pr_vec[s_i + 2] =
          log((*Pr_g_obs_offspring).at(locusLabel)[0] * (1.0 - freq_true));
        // 1 = 1S
        // 3
        Pr_vec[s_i + 3] =
          log((*Pr_g_obs_offspring).at(locusLabel)[2] * (1.0 - freq_true));
        // OS = 0 | src
        // 0 = 0S
        // 4
        Pr_vec[s_i + 4] =
          log((1.0 - freq_true));

        s_i += 5;
      }
    }
  }
} // setPrVecs()



void Source::setPrVecs(const map<string, vector<double> > &Pr_g_obs_offspring, const vector<map<string, vector<double> > > &Pr_g_inheritance, const std::string &locus){
  std::vector<std::string>::iterator locItr;
  int s_i = 0;
  double freq_true;

  for(int k = 0; k < settings["MAX_K"]; k++){

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++) {
      if(locus != *locItr) {
        s_i += 5 * alleleCount[*locItr];
      } else {
        for(int a = 0; a < alleleCount[locus]; a++) {

          freq_true =
            ((*frequencies[locus])[a] - Pr_g_obs_offspring.at(locus)[2]) /
            (1.0 - Pr_g_obs_offspring.at(locus)[1] - Pr_g_obs_offspring.at(locus)[2]);
          // XXS = unk | unk | src
          // 0 = 0S
          // 0
          Pr_vec[s_i] =
            log(Pr_g_obs_offspring.at(locus)[1] * freq_true +
            Pr_g_obs_offspring.at(locus)[0] * (1.0 - freq_true));
          // 1 = 1S

          // 1
          Pr_vec[s_i + 1] =
            log(Pr_g_obs_offspring.at(locus)[3] * freq_true +
            Pr_g_obs_offspring.at(locus)[2] * (1.0 - freq_true));

          // XOS = unk | 0 | src
          // 0 = 0S

          // 2
          Pr_vec[s_i + 2] =
            log(Pr_g_obs_offspring.at(locus)[0] * (1.0 - freq_true));
          // 1 = 1S
          // 3
          Pr_vec[s_i + 3] =
            log(Pr_g_obs_offspring.at(locus)[2] * (1.0 - freq_true));
          // OS = 0 | src
          // 0 = 0S
          // 4
          Pr_vec[s_i + 4] =
            log((1.0 - freq_true));
          s_i += 5;
        }
      }
    }
  }
} // setPrVecs()
