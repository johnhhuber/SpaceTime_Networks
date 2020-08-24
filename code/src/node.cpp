#include "node.h"
#include <typeinfo>

void Node::calcEdgeProbs(double beta)
{
  map<Node *, int>::iterator genDistItr = genDists.begin();
  map<Node *, int>::iterator timeDistItr = timeDists.begin();
  map<Node *, std::tuple<double, double> >::iterator spaceDistItr = spaceDists.begin();

  double edgeProb;
  double edgeProbSum = 0.0;
  for(; timeDistItr != timeDists.end(); genDistItr++, timeDistItr++, spaceDistItr++){
    edgeProb = 0.0;
    if(settings["LOGLIKGEN_TOGGLE"] > 0.0)
    {
      edgeProb += gsl_ran_gaussian_pdf(genDistItr->second, settings["SIGMA_gendist"] / pow(beta, settings["PROPOSAL_TEMP_FACTOR"]));
    }
    if(settings["LOGLIKTIME_TOGGLE"] > 0.0)
    {
      // edgeProb += gsl_ran_gaussian_pdf(timeDistItr->second, settings["TIMEDIST_weight"] * settings["SIGMA_time"] / pow(beta, settings["PROPOSAL_TEMP_FACTOR"]));
      edgeProb += 1.0;
    }
    if(settings["LOGLIKSPACE_TOGGLE"] > 0.0)
    {
      // edgeProb += gsl_ran_gaussian_pdf(std::get<0>(spaceDistItr->second), settings["SIGMA_space"] / pow(beta, settings["PROPOSAL_TEMP_FACTOR"]));
      // edgeProb += gsl_ran_gaussian_pdf(std::get<1>(spaceDistItr->second), settings["SIGMA_space"] / pow(beta, settings["PROPOSAL_TEMP_FACTOR"]));
      edgeProb += 1.0;
    }
    if(timeDistItr != timeDists.end()){
      edgeProbs.insert(pair<Node *, double>(timeDistItr->first, edgeProb));
      edgeProbSum += edgeProb;
    }
  }

  map<Node *, double>::iterator probItr = edgeProbs.begin();
  for(; probItr != edgeProbs.end(); probItr++){
    probItr->second = probItr->second / edgeProbSum;
  }
} // calcEdgeProbs()



int Node::calcGenDist(Node * nodeIn)
{
  int genDist = 0;

  map<string,vector<bool> *>::iterator locusItr_i, locusItr_j;
  vector<bool>::iterator alleleItr_i, alleleItr_j;

  for(locusItr_i = loci.begin(), locusItr_j = nodeIn->getLoci()->begin();
      locusItr_i != loci.end() && locusItr_j != nodeIn->getLoci()->end();
      locusItr_i++, locusItr_j++){
    for(alleleItr_i = (*locusItr_i).second->begin(), alleleItr_j = (*locusItr_j).second->begin();
        alleleItr_i != (*locusItr_i).second->end() && alleleItr_j != (*locusItr_j).second->end();
        alleleItr_i++, alleleItr_j++){
      if(*alleleItr_i != *alleleItr_j)
        genDist++;
    }
  }

  return genDist;
} // calcGenDist()



void Node::calcLogLik(const std::vector<double> &Pr_vec, const Source &sourceLocal, const std::map<std::string, bool> &validLoci, double diffusion_coef, double tau_s, double tau_l, double prop_reported_travel)
{
  if(settings["LOGLIKGEN_TOGGLE"] > 0.0){
    calcLogLikGen(Pr_vec, sourceLocal, validLoci);
    assert(!std::isnan(logLikGen));
  }

  if(settings["LOGLIKTIME_TOGGLE"] > 0.0){
    calcLogLikTime();
    assert(!std::isnan(logLikTime));
  }

  if(settings["LOGLIKSPACE_TOGGLE"] > 0.0){
    if(!(settings["BELIEVE_TRAVELHIST"] > 0.0) || ((settings["BELIEVE_TRAVELHIST"] > 0.0) & !this->isFounder()))
    {
      calcLogLikSpace(diffusion_coef);
      assert(!std::isnan(logLikSpace));
    }
  }

  if(settings["LOGLIKTRAVELHIST_TOGGLE"] > 0.0)
  {
    calcLogLikTravelHist(tau_s, tau_l, prop_reported_travel);
    assert(!std::isnan(logLikTravelHist));
  }
} // calcLogLik()



void Node::clearCache() {
  srcLocusCache.clear();
  locusCache.clear();
}// clearCache()



void Node::clearLocusCache() {
  locusCache.clear();
}// clearLocusCache()



void Node::clearSrcLocusCache() {
  srcLocusCache.clear();
}// clearSrcLocusCache()



void Node::calcLogLikGen(const std::vector<double> &Pr_vec, const Source &sourceLocal, const std::map<std::string, bool> &validLoci)
{
  vector<pair<Node *,int> >::iterator itrAnc;
  vector<pair<Source *,int> >::iterator itrSrc;

  std::vector<std::string>::iterator locItr;
  std::vector<bool> locus_alleles;

  // vector<double> src_Pr_vec;

  logLikGen = 0.0;
  double locusLik;
  int s_i = 0;
  int a_i = 0;
  std::tuple<int, int, std::string> key;
  std::tuple<int, std::string, std::string> srcKey;
  int code;

  string locus;
  int binary_offspring_parent;
  int binary_parent;
  int k;
  int allele_present;
  double numerator_E1E2 [4];
  double numerator_E3 [4];
  double denominator [3];
  // src_Pr_vec = *sourceLocal->get_Pr_vec();

  for(itrAnc = ancestors.begin(); itrAnc != ancestors.end(); itrAnc++){
    k = itrAnc->second;
    s_i = 0;
    a_i = 0;
    code = itrAnc->first->getCode();

    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++){
      locus = *locItr;
      key = std::make_tuple(k, code, locus);

      if(lociData[locus]){
        if(validLoci.at(locus)) {
          s_i += 5 * alleleCount[locus];
          a_i += 22 * alleleCount[locus];
          logLikGen += locusCache[key];
        } else {
          locus_alleles = *loci[locus];
          // reset values of terms in numerator and denominator
          for(int i = 0; i < 4; i++){
            numerator_E1E2[i] = 0.0;
            numerator_E3[i] = 0.0;
            if(i < 3) {
              denominator[i] = 0.0;
            }
          }
          if((*(itrAnc->first->getLociData()))[locus]){

            // vector<bool> allele_data = (*itrAnc->first->getLoci()->find(locus)->second);

            // loop through each allele at the present locus
            for(unsigned a = 0; a < alleleCount[locus]; a++){
              // determine similarity of offspring and parent alleles
              binary_offspring_parent = 2 * int(locus_alleles[a]) + (*(*itrAnc->first->getLoci())[locus])[a];
              binary_parent = (*(*itrAnc->first->getLoci())[locus])[a];
              // add probability to each numerator and denominator term
              numerator_E1E2[0] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + binary_offspring_parent);
              numerator_E1E2[1] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 8 + binary_offspring_parent);
              numerator_E1E2[2] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 4 + binary_offspring_parent);
              numerator_E1E2[3] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 12 + binary_offspring_parent);

              numerator_E3[0] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + binary_parent);
              numerator_E3[1] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 8 + binary_parent);
              numerator_E3[2] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 4 + binary_parent);
              numerator_E3[3] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 12 + binary_parent);

              denominator[0] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 16 + binary_parent);
              denominator[1] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 18 + binary_parent);
              denominator[2] += Pr_vec.at(a_i + (k - 1) * 22 * k_stride + 20 + binary_parent);

              s_i += 5;
              a_i += 22;
            }

            // calculate locus-level probability and add to overall gen likelihood
            locusLik =
              log(
                exp(numerator_E1E2[0]) -
                exp(numerator_E1E2[1]) -
                exp(numerator_E1E2[2]) +
                exp(numerator_E1E2[3])) -
              log(1.0 -
                exp(denominator[0]) -
                exp(denominator[1]) +
                exp(denominator[2]) -
                exp(numerator_E3[0]) +
                exp(numerator_E3[1]) +
                exp(numerator_E3[2]) -
                exp(numerator_E3[3]));


          } else {
            // loop through each allele at the present locus
            for(unsigned a = 0; a < alleleCount[locus]; a++){
              allele_present = int(locus_alleles[a]);

              /*
                This is only working right now because k always is 1.  Will need to determine the stride length and multiply accordingly
                by K to calculate correct offset. Or actually maybe not, since the calculations currently dont take K into account.
              */
              numerator_E1E2[0] += sourceLocal.get_Pr_vec().at(s_i + (k - 1) * 5 * k_stride + allele_present);
              numerator_E1E2[1] += sourceLocal.get_Pr_vec().at(s_i + (k - 1) * 5 * k_stride + 2 + allele_present);
              numerator_E3[0] += sourceLocal.get_Pr_vec().at(s_i + (k - 1) * 5 * k_stride);
              numerator_E3[1] += sourceLocal.get_Pr_vec().at(s_i + (k - 1) * 5 * k_stride + 2);
              denominator[0] += sourceLocal.get_Pr_vec().at(s_i + (k - 1) * 5 * k_stride + 4);


              s_i += 5;
              a_i += 22;

            }

            // calculate locus-level probability and add to overall gen likelihood
            locusLik =
              log(
                exp(numerator_E1E2[0]) -
                exp(numerator_E1E2[1])) -
              log(1.0 -
                exp(denominator[0]) -
                exp(numerator_E3[0]) +
                exp(numerator_E3[1]));
          }
          locusCache[key] = locusLik;
          logLikGen += locusLik;
        }
      } else {
        s_i += 5 * alleleCount[locus];
        a_i += 22 * alleleCount[locus];
      }
    }
  }

  for(itrSrc = sources.begin(); itrSrc != sources.end(); itrSrc++){
    s_i = 0;
    k = itrSrc->second;
    locItr = locus_order.begin();
    for(; locItr != locus_order.end(); locItr++){
      locus = *locItr;
      srcKey = std::make_tuple(k, itrSrc->first->getName(), locus);

      if(lociData[locus]){
        if(validLoci.at(locus)) {
          s_i += 5 * alleleCount[locus];
          logLikGen += srcLocusCache[srcKey];
        } else {
          numerator_E1E2[0] = 0.0;
          numerator_E1E2[1] = 0.0;
          numerator_E3[0] = 0.0;
          numerator_E3[1] = 0.0;
          denominator[0] = 0.0;

          // loop through each allele at the present locus
          for(unsigned a = 0; a < alleleCount[locus]; a++){
            allele_present = int((*loci[locus])[a]);
            // add probability to each numerator and denominator term
            numerator_E1E2[0] += itrSrc->first->get_Pr_vec()[s_i + (k - 1) * 5 * k_stride + allele_present];
            numerator_E1E2[1] += itrSrc->first->get_Pr_vec()[s_i + (k - 1) * 5 * k_stride + 2 + allele_present];
            numerator_E3[0] += itrSrc->first->get_Pr_vec()[s_i + (k - 1) * 5 * k_stride];
            numerator_E3[1] += itrSrc->first->get_Pr_vec()[s_i + (k - 1) * 5 * k_stride + 2];
            denominator[0] += itrSrc->first->get_Pr_vec()[s_i + (k - 1) * 5 * k_stride + 4];
            s_i += 5;
          }

          // calculate locus-level probability and add to overall gen likelihood
          locusLik =
            log(
              exp(numerator_E1E2[0]) -
              exp(numerator_E1E2[1])) -
            log(1.0 -
              exp(denominator[0]) -
              exp(numerator_E3[0]) +
              exp(numerator_E3[1]));
          srcLocusCache[srcKey] = locusLik;
          logLikGen += locusLik;
        }
      } else {
        s_i += 5 * alleleCount[locus];
      }
    }
  }
} // calcLogLikGen()



void Node::calcLogLikTime()
{
  int SI, k;
  std::string transmission_scenario;
  vector<pair<Node *,int> >::iterator itrAnc;
  vector<pair<Source *, int> >::iterator itrSrc;

  logLikTime = 0.0;

  if(sources.empty())
  {
    for(itrAnc = ancestors.begin(); itrAnc != ancestors.end(); itrAnc++)
    {
      SI = time - itrAnc->first->getTime();
      k = itrAnc->second;
      transmission_scenario = itrAnc->first->getType() + type;

      if(transmission_scenario == "STST")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_STS_1) && (SI < MAX_time_diff_STS_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STS_1, TIME_p_STS_1, TIME_k_STS_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_STS_2) && (SI < MAX_time_diff_STS_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STS_2, TIME_p_STS_2, TIME_k_STS_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_STS_3) && (SI < MAX_time_diff_STS_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STS_3, TIME_p_STS_3, TIME_k_STS_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "STAT" )
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_STA_1) && (SI < MAX_time_diff_STA_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STA_1, TIME_p_STA_1, TIME_k_STA_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_STA_2) && (SI < MAX_time_diff_STA_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STA_2, TIME_p_STA_2, TIME_k_STA_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_STA_3) && (SI < MAX_time_diff_STA_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STA_3, TIME_p_STA_3, TIME_k_STA_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "STAU")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_STA_1) && (SI < MAX_time_diff_STA_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STA_1, TIME_p_STA_1, TIME_k_STA_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_STA_2) && (SI < MAX_time_diff_STA_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STA_2, TIME_p_STA_2, TIME_k_STA_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_STA_3) && (SI < MAX_time_diff_STA_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_STA_3, TIME_p_STA_3, TIME_k_STA_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "ATAT")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_ATA_1) && (SI < MAX_time_diff_ATA_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATA_1, TIME_p_ATA_1, TIME_k_ATA_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_ATA_2) && (SI < MAX_time_diff_ATA_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATA_2, TIME_p_ATA_2, TIME_k_ATA_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_ATA_3) && (SI < MAX_time_diff_ATA_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATA_3, TIME_p_ATA_3, TIME_k_ATA_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "ATAU")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_ATA_1) && (SI < MAX_time_diff_ATA_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATA_1, TIME_p_ATA_1, TIME_k_ATA_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_ATA_2) && (SI < MAX_time_diff_ATA_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATA_2, TIME_p_ATA_2, TIME_k_ATA_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_ATA_3) && (SI < MAX_time_diff_ATA_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATA_3, TIME_p_ATA_3, TIME_k_ATA_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "AUAT")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_AUA_1) && (SI < MAX_time_diff_AUA_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUA_1, TIME_p_AUA_1, TIME_k_AUA_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_AUA_2) && (SI < MAX_time_diff_AUA_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUA_2, TIME_p_AUA_2, TIME_k_AUA_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_AUA_3) && (SI < MAX_time_diff_AUA_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUA_3, TIME_p_AUA_3, TIME_k_AUA_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "AUAU")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_AUA_1) && (SI < MAX_time_diff_AUA_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUA_1, TIME_p_AUA_1, TIME_k_AUA_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_AUA_2) && (SI < MAX_time_diff_AUA_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUA_2, TIME_p_AUA_2, TIME_k_AUA_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_AUA_3) && (SI < MAX_time_diff_AUA_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUA_3, TIME_p_AUA_3, TIME_k_AUA_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "ATST")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_ATS_1) && (SI < MAX_time_diff_ATS_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATS_1, TIME_p_ATS_1, TIME_k_ATS_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_ATS_2) && (SI < MAX_time_diff_ATS_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATS_2, TIME_p_ATS_2, TIME_k_ATS_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_ATS_3) && (SI < MAX_time_diff_ATS_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_ATS_3, TIME_p_ATS_3, TIME_k_ATS_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
      else if(transmission_scenario == "AUST")
      {
        switch(k)
        {
          case 1:
            if((SI >= -TIME_shift_AUS_1) && (SI < MAX_time_diff_AUS_1))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUS_1, TIME_p_AUS_1, TIME_k_AUS_1));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 2:
            if((SI >= -TIME_shift_AUS_2) && (SI < MAX_time_diff_AUS_2))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUS_2, TIME_p_AUS_2, TIME_k_AUS_2));
            else
              logLikTime += MIN_loglik_time;
            break;
          case 3:
            if((SI >= -TIME_shift_AUS_3) && (SI < MAX_time_diff_AUS_3))
              logLikTime += log(gsl_ran_negative_binomial_pdf(SI + TIME_shift_AUS_3, TIME_p_AUS_3, TIME_k_AUS_3));
            else
              logLikTime += MIN_loglik_time;
            break;
        }
      }
    }
  }
  else
  {
    // want to have an if statement for the travel history
    for(itrSrc = sources.begin(); itrSrc != sources.end(); itrSrc++)
    {
      if(itrSrc->first->getName() == "local" || itrSrc->first->getName() == "s")
      {
        if(type == "ST")
          logLikTime += log(LIK_founder_time_S);
        if(type == "AT")
          logLikTime += log(LIK_founder_time_A);
        if(type == "AU")
          logLikTime += log(LIK_founder_time_A);
      }
    }
  }
}// calcLogLikTime()



void Node::calcLogLikSpace(double diffusion_coef)
{
  double sigma;
  int k;
  double dist_lat, dist_lon;
  std::string caseHistory;
  vector<pair<Node *,int> >::iterator itrAnc;
  vector<pair<Source *,int> >::iterator itrSrc;

  logLikSpace = 0.0;

  if(sources.empty())
  {
    // confirm that offspring does not have missing coordinates
    if(!(fabs(lat - MISSING_VALUE) < std::numeric_limits<double>::epsilon()) & !(fabs(lon - MISSING_VALUE) < std::numeric_limits<double>::epsilon()))
    {
      for(itrAnc = ancestors.begin(); itrAnc != ancestors.end(); itrAnc++)
      {
        k = itrAnc->second;
        //SI = time - itrAnc->first->getTime();
        caseHistory = itrAnc->first->getType() + type[0];
        sigma = calcSpaceSigma(itrAnc->first, diffusion_coef, k);

        // check whether ancestor has missing coordinates
        // if(std::isnan(itrAnc->first->getLat()) || std::isnan(itrAnc->first->getLon()))
        if((fabs(itrAnc->first->getLat() - MISSING_VALUE) < std::numeric_limits<double>::epsilon()) || (fabs(itrAnc->first->getLon() - MISSING_VALUE) < std::numeric_limits<double>::epsilon()))
        {
          if(sigma > 0.0)
          {
            logLikSpace += INTERCEPT_UNKNOWNCOORDS + SLOPE_UNKNOWNCOORDS * log(sigma);
          }
          else
          {
            logLikSpace += MIN_loglik_space;
          }
        }
        else
        {
          dist_lat = std::get<0>(spaceDists.find(itrAnc->first)->second);
          dist_lon = std::get<1>(spaceDists.find(itrAnc->first)->second);

          if(sigma > 0.0)
          {
            logLikSpace += log(gsl_ran_gaussian_pdf(dist_lat, sigma));
            logLikSpace += log(gsl_ran_gaussian_pdf(dist_lon, sigma));
          }
          else
          {
            logLikSpace += MIN_loglik_space;
          }
        }
      }
    }
    else
    {
      for(itrAnc = ancestors.begin(); itrAnc != ancestors.end(); itrAnc++)
      {
        k = itrAnc->second;
        caseHistory = itrAnc->first->getType() + type[0];
        sigma = calcSpaceSigma(itrAnc->first, diffusion_coef, k);

        if(sigma > 0.0)
        {
          logLikSpace += INTERCEPT_UNKNOWNCOORDS + SLOPE_UNKNOWNCOORDS * log(sigma);
        }
        else
        {
          logLikSpace += MIN_loglik_space;
        }
      }
    }
  }
  else
  {
    // want an if statement to toggle the believe travel history
    for(itrSrc = sources.begin(); itrSrc != sources.end(); itrSrc++)
    {
      logLikSpace += INTERCEPT_UNKNOWNLOCAL + SLOPE_UNKNOWNLOCAL * log(diffusion_coef);
    }
  }
} // calcLogLikSpace()



void Node::calcLogLikTravelHist(double tau_s, double tau_l, double prop_reported_travel)
{
  vector<pair<Node *, int> >::iterator itrAnc;
  vector<pair<Source *,int> >::iterator itrSrc;

  logLikTravelHist = 0.0;

  if(!(fabs(history - MISSING_VALUE) < std::numeric_limits<double>::epsilon()))
  {
    if(sources.empty())
    {
      for(itrAnc = ancestors.begin(); itrAnc != ancestors.end(); itrAnc++)
      {
        switch(history)
        {
          case 1: logLikTravelHist += log(tau_l);
            break;
          case 0: logLikTravelHist += log(1 - tau_l);
            break;
        }
      }
    }
    else
    {
      for(itrSrc = sources.begin(); itrSrc != sources.end(); itrSrc++)
      {
        if(itrSrc->first->getName() == "local")
        {
          switch(history)
          {
            case 1: logLikTravelHist += log(tau_l);
              break;
            case 0: logLikTravelHist += log(1 - tau_l);
              break;
          }
        }
        else
        {
          switch(history)
          {
            case 1: logLikTravelHist += log(tau_s);
              break;
            case 0: logLikTravelHist += log(1 - tau_s);
              break;
          }
        }
      }
    }
  }else{
    if(sources.empty())
    {
      for(itrAnc = ancestors.begin(); itrAnc != ancestors.end(); itrAnc++)
      {
        logLikTravelHist += log(prop_reported_travel * tau_l + (1 - prop_reported_travel) * (1 - tau_l));
      }
    }
    else
    {
      for(itrSrc = sources.begin(); itrSrc != sources.end(); itrSrc++)
      {
        logLikTravelHist += log(prop_reported_travel * tau_s + (1 - prop_reported_travel) * (1 - tau_s));
      }
    }
  }
}// calcLogLikTravelHist())



int Node::calcTimeDist(Node * nodeIn)
{
  return nodeIn->getTime() - time;
} // calcTimeDist()



std::tuple<double,double> Node::calcSpaceDist(Node * nodeIn)
{

  if(std::isnan(lon) || std::isnan(lat) || std::isnan(nodeIn->getLon()) || std::isnan(nodeIn->getLat()))
  {
    return std::make_tuple(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
  }
  else
  {
    const double EARTH_RADIUS = 6371.0;
    double dLat = (lat - nodeIn->getLat()) * (M_PI / 180.0);
    double dLon = (lon - nodeIn->getLon()) * (M_PI / 180.0);

    double dx = EARTH_RADIUS * dLat;
    double dy = EARTH_RADIUS * dLon  * cos(M_PI * lat / 180.0);

    return std::make_tuple(dx, dy);
  }
} // calcSpaceDist()



double Node::calcSpaceSigma(Node * nodeIn, double diffusion_coef, int k)
{
  std::string caseHistory = nodeIn->getType() + type[0];
  double variance = 0.0;
  int index;
  if(caseHistory == "STS")
  {
    index = time - nodeIn->getTime() + MAX_time_IDP_S - 1;
    switch(k)
    {
      case 1:
        if(index < 0 || index >= sigma_STS_k_1.size()){
          variance = 0.0;
        }
        else
          variance = diffusion_coef * sigma_STS_k_1[index];
        break;
      case 2:
        if(index < 0 || index >= sigma_STS_k_2.size())
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_STS_k_2[index];
        break;
      case 3:
        if(index < 0 || index >= sigma_STS_k_3.size())
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_STS_k_3[index];
        break;
    }
  }

  else if(caseHistory == "STA")
  {
    index = time - nodeIn->getTime() + MAX_time_IDP_S - 1;
    switch(k)
    {
      case 1:
        if(index < 0 || index >= static_cast<int>(sigma_STA_k_1.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_STA_k_1[index];
        break;
      case 2:
        if(index < 0 || index >= static_cast<int>(sigma_STA_k_2.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_STA_k_2[index];
        break;
      case 3:
        if(index < 0 || index >= static_cast<int>(sigma_STA_k_3.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_STA_k_3[index];
        break;
    }
  }

  else if(caseHistory == "ATS")
  {
    index = time - nodeIn->getTime() + MAX_time_IDP_A - 1;
    switch(k)
    {
      case 1:
        if(index < 0 || index >= static_cast<int>(sigma_ATS_k_1.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_ATS_k_1[index];
        break;
      case 2:
        if(index < 0 || index >= static_cast<int>(sigma_ATS_k_2.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_ATS_k_2[index];
        break;
      case 3:
        if(index < 0 || index >= static_cast<int>(sigma_ATS_k_3.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_ATS_k_3[index];
        break;
    }
  }
  else if(caseHistory == "ATA")
  {
    index = time - nodeIn->getTime() + MAX_time_IDP_A - 1;
    switch(k)
    {
      case 1:
        if(index < 0 || index >= static_cast<int>(sigma_ATA_k_1.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_ATA_k_1[index];
        break;
      case 2:
        if(index < 0 || index >= static_cast<int>(sigma_ATA_k_2.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_ATA_k_2[index];
        break;
      case 3:
        if(index < 0 || index >= static_cast<int>(sigma_ATA_k_3.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_ATA_k_3[index];
        break;
    }
  }
  else if(caseHistory == "AUS")
  {
    index = time - nodeIn->getTime() + MAX_time_IDP_A - 1;
    switch(k)
    {
      case 1:
        if(index < 0 || index >= static_cast<int>(sigma_AUS_k_1.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_AUS_k_1[index];
        break;
      case 2:
        if(index < 0 || index >= static_cast<int>(sigma_AUS_k_2.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_AUS_k_2[index];
        break;
      case 3:
        if(index < 0 || index >= static_cast<int>(sigma_AUS_k_3.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_AUS_k_3[index];
        break;
    }
  }
  else
  {
    index = time - nodeIn->getTime() + MAX_time_IDP_A - 1;
    switch(k)
    {
      case 1:
        if(index < 0 || index >= static_cast<int>(sigma_AUA_k_1.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_AUA_k_1[index];
        break;
      case 2:
        if(index < 0 || index >= static_cast<int>(sigma_AUA_k_2.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_AUA_k_2[index];
        break;
      case 3:
        if(index < 0 || index >= static_cast<int>(sigma_AUA_k_3.size()))
          variance = 0.0;
        else
          variance = diffusion_coef * sigma_AUA_k_3[index];
        break;
    }
  }
  return sqrt(variance);
} // calcSpaceSigma()



void Node::checkForNAs()
{
  bool present;

  map<string,vector<bool> *>::iterator locus_itr;
  vector<bool>::iterator allele_itr;

  for(locus_itr = loci.begin(); locus_itr != loci.end(); locus_itr++){
    present = 0;
    for(allele_itr = locus_itr->second->begin(); allele_itr != locus_itr->second->end(); allele_itr++){
      if(*allele_itr){
        present = 1;
      }
    }
    lociData.insert(pair<string,bool>(locus_itr->first, present));
  }
} // checkForNAs()



void Node::copyNode(Node * nodeIn)
{
  ancestors = *(nodeIn->getAncestors());
  code = nodeIn->getCode();
  edgeProbs = *(nodeIn->getEdgeProbs());
  genDists = *(nodeIn->getGenDists());
  timeDists = *(nodeIn->getTimeDists());
  spaceDists = *(nodeIn->getSpaceDists());
  loci = *(nodeIn->getLoci());
  logLikGen = nodeIn->getLogLikGen();
  logLikTime = nodeIn->getLogLikTime();
  logLikSpace = nodeIn->getLogLikSpace();
  logLikTravelHist = nodeIn->getLogLikTravelHist();
  sources = *(nodeIn->getSources());
  locusCache = nodeIn->getLocusCache();
  srcLocusCache = nodeIn->getSrcLocusCache();
  time = nodeIn->getTime();
  lat = nodeIn->getLat();
  lon = nodeIn->getLon();
  history = nodeIn->getHistory();
} // copyNode()



bool Node::hasAncSrc()
{
  if(ancestors.size() == 0 && sources.size() == 0)
    return 0;
  else
    return 1;
} // hasAncSrc()



bool Node::isFounder()
{
  if(!sources.empty()){
    if(sources.begin()->first->getName() != LOCAL_SOURCE){
      return 1;
    }
  }
  return 0;
} // isFounder()



bool Node::isOrphan()
{
  if(!sources.empty()){
    if(sources.begin()->first->getName() == LOCAL_SOURCE){
      return 1;
    }
  }
  return 0;
} // isOrphan()



bool Node::isAncestor(Node * nodeIn)
{
  if(!ancestors.empty()){
    for(vector<pair<Node *,int> >::iterator itrVec = ancestors.begin();
        itrVec != ancestors.end();
        itrVec++){
      if(nodeIn == (*itrVec).first)
        return 1;
    }
  }
  return 0;
} // isAncestor()



void Node::loadTime(
  int n_node)
{
  string line;

  vector<string> filenames;
  filenames.push_back("../data/times.txt");

  vector<string>::iterator itrLoc;
  for(itrLoc = filenames.begin();
      itrLoc != filenames.end();
      itrLoc++){
    ifstream NewInfile(itrLoc->c_str());

    for(short i = 0; i < n_node; i++)
      getline(NewInfile, line);

    this->time = stoi(line, nullptr);
  }
} // loadTime()



void Node::removeAncestor(int toRemove)
{
  if(!ancestors.empty()){
    vector<pair<Node *,int> >::iterator vecItr = ancestors.begin();
    for(int i = 0; i < toRemove; i++)
      vecItr++;

    ancestors.erase(vecItr);
  }
} // removeAncestor()



void Node::removeAllAncestors()
{
  if(!ancestors.empty()){
    ancestors.clear();
  }
} // removeAllAncestors()



void Node::removeSource()
{
  if(!sources.empty())
    sources.pop_back();
} // removeSource()



void Node::resetAll()
{
  ancestors.clear();
  code = 0;
  edgeProbs.clear();
  genDists.clear();
  timeDists.clear();
  spaceDists.clear();
  loci.clear();
  logLikGen = 0.0;
  logLikTime = 0.0;
  logLikSpace = 0.0;
  logLikTravelHist = 0.0;
  sources.clear();
  time = 0;
  lat = 0.0;
  lon = 0.0;
  history = 0;
} // resetAll()



void Node::setGenDist(Node * nodeIn, int distIn)
{
  genDists.insert(pair<Node *, int>(nodeIn, distIn));
} // setGenDist()



void Node::set_k_anc(Node * ancIn, int kIn)
{
  vector<pair<Node *,int> >::iterator itrVec = ancestors.begin();
  for(; (*itrVec).first != ancIn; itrVec++){;}

  ancestors.erase(itrVec);

  ancestors.push_back(pair<Node *,int>(ancIn, kIn));
} // set_k_anc()



void Node::set_k_src(Source * srcIn, int kIn)
{
  vector<pair<Source *,int> >::iterator itrVec = sources.begin();
  for(; (*itrVec).first != srcIn; itrVec++){;}

  sources.erase(itrVec);

  sources.push_back(pair<Source *,int>(srcIn, kIn));
} // set_k_src()



void Node::setTimeDist(Node * nodeIn, int distIn)
{
  timeDists.insert(pair<Node *, int>(nodeIn, distIn));
} // setTimeDist()


bool Node::operator<(const Node& rhs) {
  return code < rhs.getCode();
};


void Node::setSpaceDist(Node * nodeIn, std::tuple<double, double> spaceDistIn)
{
  spaceDists.insert(pair<Node *, std::tuple<double, double> >(nodeIn, spaceDistIn));
} // setSpaceDist()
