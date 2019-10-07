#include "cstm.hpp"

CSTM::CSTM() {
        _ndim_d = NDIM_D;
        _sigma_u = SIGMA_U;
        _sigma_phi = SIGMA_PHI;
        _sigma_alpha0 = SIGMA_ALPHA;
        _gamma_alpha_a = GAMMA_ALPHA_A;
        _gamma_alpha_b = GAMMA_ALPHA_B;
        _g0 = NULL;
        _word_vectors = NULL;
        _doc_vectors = NULL;
        _n_k = NULL;
        _word_count = NULL;
        _sum_n_k = NULL;
        _Zi = NULL;
        _log_likelihood_first_term = NULL;
        _tmp_vec = NULL;
        _vocabrary_size = 0;
        _num_documents = 0;
        _sum_word_frequency = 0;
        _ignore_word_count = 0;
        _standard_normal_distribution = normal_distribution<double>(0, 1);
}

CSTM::~CSTM() {}
