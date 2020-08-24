#include "main.h"
#include "hyperparams.h"
#include "network.h"
#include "node.h"
#include "source.h"
#include <chrono>



double LIK_founder_time_S;
double LIK_founder_time_A;

vector<double> sigma_AUA_k_1;
vector<double> sigma_AUA_k_2;
vector<double> sigma_AUA_k_3;
vector<double> sigma_AUS_k_1;
vector<double> sigma_AUS_k_2;
vector<double> sigma_AUS_k_3;
vector<double> sigma_STA_k_1;
vector<double> sigma_STA_k_2;
vector<double> sigma_STA_k_3;
vector<double> sigma_ATA_k_1;
vector<double> sigma_ATA_k_2;
vector<double> sigma_ATA_k_3;
vector<double> sigma_STS_k_1;
vector<double> sigma_STS_k_2;
vector<double> sigma_STS_k_3;
vector<double> sigma_ATS_k_1;
vector<double> sigma_ATS_k_2;
vector<double> sigma_ATS_k_3;

map<string,double> settings;

int WRITE_STATUS;

void calc_LIK_founder_time();
void loadSettings(std::string, std::map<string,double>&);
void loadSigma(std::vector<std::string> , int , std::vector<double>& , std::vector<double>& , std::vector<double>& );
void loadSigmaFile(std::string , std::vector<double>& , std::vector<double>& , std::vector<double>& );
void writeHeader(string);
void writeSwaps(string, bool, double, double, double, double);

int main(int argc, char* argv[])
{
  gsl_rng *rGsl = gsl_rng_alloc(gsl_rng_taus);
  string node_file_path;
  string sources_file_path;
  string network_file_path;
  string scalars_file_path;
  string parameters_file_path;
  string network_destination;
  string scalars_destination;
  string swaps_destination;
  string sigma_AUA_file_path = "variancesAUA.csv";
  string sigma_AUS_file_path = "variancesAUS.csv";
  string sigma_STA_file_path = "variancesSTA.csv";
  string sigma_ATA_file_path = "variancesATA.csv";
  string sigma_STS_file_path = "variancesSTS.csv";
  string sigma_ATS_file_path = "variancesATS.csv";
  string settings_file_path;

  int rng_seed = 4;
  map<string, float> epsilon_neg;
  map<string, float> epsilon_pos;
  float lambda;
  float mu;
  float diffusion_coef;
  bool hotstart = false;
  string HOTSTART_FLAG = "--hotstart";

  if (argc > 0) {
    if (argc == 8) {
      node_file_path = argv[1];
      sources_file_path = argv[2];
      settings_file_path = argv[3];
      network_destination = argv[4];
      scalars_destination = argv[5];
      swaps_destination = argv[6];
      rng_seed = 292745 + std::stoi(argv[7]);
    } else if (argc == 4) {
      node_file_path = argv[1];
      sources_file_path = argv[2];
      settings_file_path = argv[3];
    } else if (argc == 6) {
      node_file_path = argv[1];
      sources_file_path = argv[2];
      settings_file_path = argv[3];
      network_file_path = argv[4];
      parameters_file_path = argv[5];
     } else if ((argc == 9) && (argv[8] == HOTSTART_FLAG)) {
      node_file_path = argv[1];
      sources_file_path = argv[2];
      settings_file_path = argv[3];
      network_destination = argv[4];
      scalars_destination = argv[5];
      swaps_destination = argv[6];
      rng_seed = 292745 + std::stoi(argv[7]);
      hotstart = true;
    } else {
      cout << endl;
      cout << "Signature of method is as follows: " << endl;
      cout << "mc3 nodes_data_file sources_data_file network_destination_file scalars_destination_file swaps_destination_file seed" << endl;
      cout << endl;
      cout << "Otherwise, files at these destinations will be used:" << endl;
      cout << "\tNodes: " << node_file_path << endl;
      cout << "\tSources: " << sources_file_path << endl;
      cout << "\tNetwork Destination: " << network_destination << endl;
      cout << "\tScalars Destination: " << scalars_destination << endl;
      cout << "\tSwaps Destination: " << swaps_destination << endl;
      cout << endl;
      exit(0);
    }
  }


  double beta;
  calc_LIK_founder_time();

  loadSigmaFile(sigma_AUA_file_path, sigma_AUA_k_1, sigma_AUA_k_2, sigma_AUA_k_3);
  loadSigmaFile(sigma_AUS_file_path, sigma_AUS_k_1, sigma_AUS_k_2, sigma_AUS_k_3);
  loadSigmaFile(sigma_STA_file_path, sigma_STA_k_1, sigma_STA_k_2, sigma_STA_k_3);
  loadSigmaFile(sigma_ATA_file_path, sigma_ATA_k_1, sigma_ATA_k_2, sigma_ATA_k_3);
  loadSigmaFile(sigma_STS_file_path, sigma_STS_k_1, sigma_STS_k_2, sigma_STS_k_3);
  loadSigmaFile(sigma_ATS_file_path, sigma_ATS_k_1, sigma_ATS_k_2, sigma_ATS_k_3);

  if(argc != 6){
    unsigned int threads;
    double temp_beta;

    size_t networkIndex1;
    size_t networkIndex2;

    Network * network;
    Network * rNetwork1;
    Network * rNetwork2;
    vector<Network *> networks;
    vector<Network *>::iterator itrNetwork;

    loadSettings(settings_file_path, settings);

    gsl_rng_set(rGsl, rng_seed);

    if(settings["N_CHAINS"] == 0){
      threads = omp_get_num_procs();
      omp_set_num_threads(threads);
    } else {
      omp_set_num_threads(settings["N_CHAINS"]);
      threads = settings["N_CHAINS"];
    }

    for(unsigned int i = 0; i < threads; i++){
      beta = double(1) / double(1 + settings["DELTA_T"]*(i));
      Network* network = new Network(i, beta, node_file_path, sources_file_path, network_file_path, network_destination, scalars_destination, settings, gsl_rng_get(rGsl));
      if(hotstart) {
        network->hotstartInitialize();
      } else {
        network->initialize();
      }
      networks.push_back(network);
    }
    // cout << "All networks initialized. Random number generator seed set to " << rng_seed << "." << endl;

    double accept_prob, accept_num = 0;
    bool accept;
    if(!hotstart) {
      writeHeader(swaps_destination);
    }


    auto start = std::chrono::system_clock::now();
    auto start_time = std::chrono::system_clock::to_time_t(start);

    cout << "Start Time: " << std::ctime(&start_time) << endl;

    for(int i = 0; i < int(settings["SAMPLE_NUMBER"]) * int(settings["SAMPLE_FREQUENCY"]); i++){
      std::chrono::high_resolution_clock::time_point loop_start = std::chrono::high_resolution_clock::now();
      if(i % int(settings["SAMPLE_FREQUENCY"]) == 0){
        WRITE_STATUS = 1;
        if( (i / int(settings["SAMPLE_FREQUENCY"])) % int(settings["SAMPLE_FREQUENCY"]) == 0){
          cout << "Sampled  " <<  i / int(settings["SAMPLE_FREQUENCY"]) << " of " << int(settings["SAMPLE_NUMBER"]) << endl;
        }
      } else {
        WRITE_STATUS = 0;
      }

      // cout << "Parallel Loop Starting" << endl;
      #pragma omp parallel for private(network)
      for(size_t j = 0; j < networks.size(); j++){
        network = networks[j];
        network->simulate();
      }
      // cout << "Done Simulating" << endl;

      std::chrono::high_resolution_clock::time_point current = std::chrono::high_resolution_clock::now();

      std::chrono::duration<double> elapsed_seconds = current - loop_start;

      float comps_per_sec = (settings["SWAP_FREQUENCY"] * networks.size()) / elapsed_seconds.count();

      // cout << "Computations/s: " << comps_per_sec << "\n" << endl;


      if(networks.size() > 1){
        networkIndex1 = gsl_rng_uniform_int(rGsl, networks.size());
        networkIndex2 = gsl_rng_uniform_int(rGsl, networks.size());

        while(networkIndex2 == networkIndex1){
          networkIndex2 = gsl_rng_uniform_int(rGsl, networks.size());
        }
        // cout << "Selecting Network... " << networkIndex1 << " " << networkIndex2 << endl;
        rNetwork1 = networks[networkIndex1];
        rNetwork2 = networks[networkIndex2];

        accept = 0;
        accept_prob = fmin(1, exp(
          rNetwork1->getBeta() * rNetwork2->getLogLik() + rNetwork2->getBeta() * rNetwork1->getLogLik() -
          rNetwork1->getBeta() * rNetwork1->getLogLik() - rNetwork2->getBeta() * rNetwork2->getLogLik()));
        if(gsl_rng_uniform(rGsl) < accept_prob){
          temp_beta = rNetwork1->getBeta();
          rNetwork1->setBeta(rNetwork2->getBeta());
          rNetwork2->setBeta(temp_beta);
          accept = 1;
          accept_num++;
        }

        if(WRITE_STATUS){
          writeSwaps(swaps_destination, accept, accept_prob, rNetwork2->getBeta(), rNetwork1->getBeta(), accept_num / int(settings["SAMPLE_FREQUENCY"]));
          accept_num = 0;
        }
      }
    }
  } else {
    // cout << "Running Calculator Mode." << endl;
    loadSettings(settings_file_path, settings);
  	Network* network = new Network(1, 1, node_file_path, sources_file_path, network_file_path, parameters_file_path, settings, gsl_rng_get(rGsl));
    cout << network->calcLogLik() << endl;
  }

  return(0);
}// main()



void calc_LIK_founder_time()
{
	LIK_founder_time_S = 0.0;
	for(int t = 0; t < MAX_time_diff_AUS_1; t++)
		LIK_founder_time_S += pow(gsl_ran_negative_binomial_pdf(t, TIME_p_AUS_1, TIME_k_AUS_1), 2);

  LIK_founder_time_A = 0.0;
  for(int t = 0; t < MAX_time_diff_AUA_1; t++)
    LIK_founder_time_A += pow(gsl_ran_negative_binomial_pdf(t, TIME_p_AUA_1, TIME_k_AUA_1), 2);
} // calcLIK_founder_time()



void loadSettings(std::string settings_file_path, std::map<string,double>& settings)
{
  ifstream file(settings_file_path);
  char settings_delim = ',';
  string line;
  vector<string> settings_data;
  string setting_name;
  double setting_value;

  if(!file.is_open()) {
    perror("Error While Opening Settings File");
  }

  while(getline(file, line)) {
    settings_data = util::split(line, settings_delim);
    setting_name = settings_data[0];
    setting_value = stod(settings_data[1], nullptr);
    settings.insert(pair<string,double>(setting_name, setting_value));
  }

  file.close();
} // loadSettings()



void loadSigma(std::vector<std::string> sigma_data, int code, std::vector<double>& sigma_k_1, std::vector<double>& sigma_k_2, std::vector<double>& sigma_k_3)
{
  double sigma;
  switch(code)
  {
    case 1:
      for(auto data = sigma_data.begin(); data != sigma_data.end(); data++)
      {
        sigma = stod(*data, nullptr);
        sigma_k_1.push_back(sigma);
      }
      break;
    case 2:
      for(auto data = sigma_data.begin(); data != sigma_data.end(); data++)
      {
        sigma = stod(*data, nullptr);
        sigma_k_2.push_back(sigma);
      }
      break;
    case 3:
      for(auto data = sigma_data.begin(); data != sigma_data.end(); data++)
      {
        sigma = stod(*data, nullptr);
        sigma_k_3.push_back(sigma);
      }
      break;
  }
}// loadSigma()



void loadSigmaFile(std::string sigma_file_path, std::vector<double>& sigma_k_1, std::vector<double>& sigma_k_2, std::vector<double>& sigma_k_3)
{
  ifstream file(sigma_file_path);
  int code = 1;
  std::string line;
  char delim = ',';
  std::vector<std::string> sigma_data;
  if (!file.is_open())
  {
    perror("Error While Opening File");
  }
  while(getline(file, line))
  {
    sigma_data = util::split(line, delim);
    loadSigma(sigma_data, code, sigma_k_1, sigma_k_2, sigma_k_3);
    code++;
  }
  file.close();
}// loadSigmaFile()



void writeHeader(string swaps_destination)
{
  {
    ofstream fileSwaps;
    fileSwaps.open(swaps_destination);

    fileSwaps <<
      "accepted," <<
      "accept_prob," <<
      "beta1," <<
      "beta2," <<
      "accept_rate" << endl;

    fileSwaps.close();
  }
} // writeHeader()



void writeSwaps(string swaps_destination, bool accepted, double prob, double beta1, double beta2, double rate){
  ofstream file;

  file.open(swaps_destination, ios::app);

  file << accepted << "," <<
    prob << "," <<
    beta1 << "," <<
    beta2 << "," <<
    rate << endl;

  file.close();
} // writeScalars()
