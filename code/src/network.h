#ifndef _NETWORK_H
#define _NETWORK_H



#include "source.h"
#include "node.h"
#include "main.h"
#include "hyperparams.h"
#include "util.h"
#include <assert.h>
#include <string>
#include <limits>


using namespace std;


class Network
{
	unsigned int id;
	gsl_rng *rGsl;

	map<string, double> settings;

	unsigned int change_k_accept;
	unsigned int change_k_reject;

	map<string, double> epsilon_neg;
	unsigned int epsilon_neg_accept;
	unsigned int epsilon_neg_reject;

	map<string, double> epsilon_pos;
	unsigned int epsilon_pos_accept;
	unsigned int epsilon_pos_reject;

	double lambda;
	unsigned int lambda_accept;
	unsigned int lambda_reject;

	double logLik;

	unsigned max_alleles;

	double mu;
	unsigned int mu_accept;
	unsigned int mu_reject;

	double nonLocalSources;
	int numSources;

	double diffusion_coef;
	unsigned int diffusion_coef_accept;
	unsigned int diffusion_coef_reject;

	double tau_s;
	unsigned int tau_s_accept;
	unsigned int tau_s_reject;

	double tau_l;
	unsigned int tau_l_accept;
	unsigned int tau_l_reject;

	double prop_reported_travel;
	int num_reported_travel;
	int num_travel_history;

	std::vector<double> log_proposal_k;

	std::vector<std::string> locus_order;
	std::map<std::string, int> alleleCount;
	std::map<std::string, bool> validLoci;
	int k_stride;

	map<string, double> prevMax;
	map<string, double> prevMin;

	unsigned int switch_anc_anc_accept;
	unsigned int switch_anc_anc_reject;

	unsigned int add_anc_accept;
	unsigned int add_anc_reject;

	unsigned int remove_anc_accept;
	unsigned int remove_anc_reject;

	unsigned int make_founder_accept;
	unsigned int make_founder_reject;

	unsigned int make_orphan_accept;
	unsigned int make_orphan_reject;

	double beta;

	string node_file_path;
	string sources_file_path;
	string network_file_path;
	string scalars_file_path;
	string parameters_file_path;
	string network_destination;
	string scalars_destination;

	vector<Node *> nodes;
	vector<double> proposal_probs;
	map<string,Source *> sources;

	map<string, vector<double> > Pr_g_obs_offspring;
	vector<map<string, vector<double> > > Pr_g_inheritance;
	vector<map<string, vector<vector<double> > *> > Pr_g_true_parent;

	vector<map<string, vector<vector<double> > > > Pr_XXXX;
	vector<map<string, vector<vector<double> > > > Pr_XXOX;
	vector<map<string, vector<vector<double> > > > Pr_XOXX;
	vector<map<string, vector<vector<double> > > > Pr_XOOX;
	vector<map<string, vector<vector<double> > > > Pr_XOX;
	vector<map<string, vector<vector<double> > > > Pr_OXX;
	vector<map<string, vector<vector<double> > > > Pr_OOX;

	vector<double> Pr_vec;

	map<string, vector<double> > Pr_g_obs_offspring_copy;
	vector<map<string, vector<double> > > Pr_g_inheritance_copy;
	vector<map<string, vector<vector<double> > *> > Pr_g_true_parent_copy;

	vector<map<string, vector<vector<double> > > > Pr_XXXX_copy;
	vector<map<string, vector<vector<double> > > > Pr_XXOX_copy;
	vector<map<string, vector<vector<double> > > > Pr_XOXX_copy;
	vector<map<string, vector<vector<double> > > > Pr_XOOX_copy;
	vector<map<string, vector<vector<double> > > > Pr_XOX_copy;
	vector<map<string, vector<vector<double> > > > Pr_OXX_copy;
	vector<map<string, vector<vector<double> > > > Pr_OOX_copy;

  	std::vector<double> Pr_vec_copy;

	double add_anc(Node *);
	void calcGenDists();
	void calcEdgeProbs();
	void calcSourcesMinMax();
	void calcTimeDists();
	void calcSpaceDists();
	map<string, double> change_epsilon_neg(const std::string &);
	map<string, double> change_epsilon_pos(const std::string &);
	double change_lambda();
	double change_mu();
	double change_diffusion_coef();
	double change_tau_s();
	double change_tau_l();
	double change_k_anc(Node *);
	void change_k_src(Node *);
	int getRandom_k();
	Source* chooseSource();
	bool findNode(Node *);
	void hotstartLoadNetwork();
	void hotstartLoadParameters();
	void initialize_epsilon_neg();
	void initialize_epsilon_pos();
  void initialize_k_proposal();
	void validateLoci();
	void invalidateLoci();
	void validateLocus(const std::string &);
	void invalidateLocus(const std::string &);
	void initializeParams();
	bool isAcyclic(Node *);
	bool isCyclic(Node *, Node *);
	bool isCyclic(Node *);
	void loadNodesFile();
	void loadNetworkFile();
	void loadParametersFile();
	void loadScalarsFile();
	void loadNode(vector<string> *, vector<string> *, int);
	void loadProposalProbs();
	void loadSourcesFile();
	void loadSource(string, string, vector<double> *);
	double logPrior_diffusion(double);
	double logPrior_epsilon_neg(double);
	double logPrior_epsilon_pos(double);
	double logPrior_lambda(double);
	double logPrior_mu(double);
	double logPrior_tau_l(double);
	double logPrior_tau_s(double);
	double logProposal_k(int);
	double logProposal_diffusion(double, double, double);
	double logProposal_tau_l(double, double, double);
	double logProposal_tau_s(double, double, double);
	double make_founder(Node *);
	double make_orphan(Node *);
	void randomTree();
	double remove_anc(Node *);
	double runif_interval(double, double);
	double switch_anc_anc(Node *);

	void recalculateProbs();
	void recalculateProbs(const std::string &);
	void copyProbabilities();
	void revertProbabilities();
	void setPrGObsOffspring();
	void setPrGObsOffspring(const std::string &);
	void setPrGInheritance();
	void setPrGInheritance(const std::string &);
	void setPrGTrueParent();
	void setPrGTrueParent(const std::string &);
	void setPrGSourcesAll();
	void setPrVecs();
	void setPrVecs(const std::string &);
	void initPrVecs();
	void initPrGObsOffspring();
	void initPrGTrueParent();

	void update(Node *);
	bool update_add_anc(Node *);
	void update_epsilon_neg();
	void update_epsilon_pos();
	void update_k(Node *);
	void update_lambda();
	bool update_make_founder(Node *);
	bool update_make_orphan(Node *);
	void update_mu();
	void update_diffusion_coef();
	void update_tau_s();
	void update_tau_l();
	bool update_remove_anc(Node *);
	bool update_switch_anc_anc(Node *);
	void writeHeader();
	void writeNetwork();
	void writeScalars();

public:

	// Regular Init
	Network(unsigned int id, double beta, string node_file_path, string sources_file_path, string network_file_path, string network_destination, string scalars_destination, map<string,double> settingsIn, long unsigned int seed)
	{
		cout << "Generated Network " << id << " with seed " << seed << endl;
		this->settings = settingsIn;
		this->id = id;
		this->beta = beta;
		this->node_file_path = node_file_path;
		this->sources_file_path = sources_file_path;
		this->network_destination = network_destination;
		this->scalars_destination = scalars_destination;
		this->nonLocalSources = 0;
		this->numSources = 0;
		this->k_stride = 0;
		this->rGsl = gsl_rng_alloc(gsl_rng_taus);
		this->settings = settingsIn;
		gsl_rng_set(this->rGsl, seed);
	};

	// Calculator Init
	Network(unsigned int id, double beta, string node_file_path, string sources_file_path, string network_file_path, string parameters_file_path, map<string,double> settingsIn, long unsigned int seed)
	{
		cout << "Generated Network " << id << " with seed " << seed << endl;
		this->settings = settingsIn;
		this->id = id;
		this->beta = beta;
		this->node_file_path = node_file_path;
		this->sources_file_path = sources_file_path;
		this->network_file_path = network_file_path;
		this->parameters_file_path = parameters_file_path;
		this->nonLocalSources = 0;
		this->numSources = 0;
		this->k_stride = 0;
		this->loadSourcesFile();
		this->loadNodesFile();
		this->loadNetworkFile();
		this->loadParametersFile();
    this->initialize_k_proposal();
		this->initPrGObsOffspring();
		this->initPrGTrueParent();
		this->initPrVecs();
		this->calcEdgeProbs();
		this->invalidateLoci();
		this->recalculateProbs();
		this->rGsl = gsl_rng_alloc(gsl_rng_taus);
		gsl_rng_set(this->rGsl, seed);
	};

	double calcLogLik();
	int getNumNodes(){return nodes.size();};
	int getID(){return id;};
	double getBeta(){return beta;};
	double getLogLik(){return logLik;};
	void setBeta(double);
	vector<Node *> getNodes(){return nodes;};
	void initialize();
	void hotstartInitialize();
	void simulate();
	double updateLogLik(Node *);
}; // class Network



#endif // _NETWORK_H
