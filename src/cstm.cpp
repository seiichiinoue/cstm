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
    _vocabulary_size = 0;
    _num_documents = 0;
    _sum_word_frequency = 0;
    _ignore_word_count = 0;
    _standard_normal_distribution = normal_distribution<double>(0, 1);
}

CSTM::~CSTM() {}

void CSTM::add_word(id word_id, int doc_id) {
    int *count = _n_k[doc_id];
    count[word_id] += 1;
    _sum_word_frequency += 1;
}
double CSTM::generate_noise_from_standard_normal_distribution() {
    return _standard_normal_distribution(sampler::minstd);
}
double CSTM::generate_noise_doc() {
    return _noise_doc(sampler::minstd);
}
double CSTM::generate_noise_word() {
    return _noise_word(sampler::minstd);
}
double *CSTM::generate_vector() {
    double *vec = new double[_ndim_d];
    for (int i=0; i<_ndim_d; ++i) {
            vec[i] = generate_noise_from_standard_normal_distribution();
    }
    return vec;
}
double *CSTM::draw_word_vector(double *old_vec) {
    for (int i=0; i<_ndim_d; ++i) {
        _tmp_vec[i] = old_vec[i]+generate_noise_word();
    }
    return _tmp_vec;
}
double *CSTM::draw_doc_vector(double *old_vec) {
    for (int i=0; i<_ndim_d; ++i) {
        _tmp_vec[i] = old_vec[i]+generate_noise_doc();
    }
    return _tmp_vec;
}
double CSTM::draw_alpha0(double old_alpha0) {
    double z = _noise_alpha0(sampler::minstd);
    return old_alpha0 * exp(z);
}
double CSTM::sum_alpha_word_given_doc(int doc_id) {
    double sum = 0;
    for (id word_id=0; word_id<_vocabulary_size; ++word_id) {
        sum += compute_alpha_word_given_doc(word_id, doc_id);
    }
    return sum;
}
double CSTM::compute_alpha_word_given_doc(id word_id, int doc_id) {
    double *word_vec = _word_vectors[word_id];
    double *doc_vec = _doc_vectors[doc_id];
    double g0 = get_g0_of_word(word_id);
    return _compute_alpha_word(word_vec, doc_vec, g0);
}
double CSTM::compute_reduced_log_probability_document(id word_id, int doc_id) {
    double log_pw = 0;
    double Zi = _Zi[doc_id];
    int n_k = get_word_count_in_doc(word_id, doc_id);
    if (n_k == 0) {
        return _compute_reduced_log_probability_document(word_id, doc_id, n_k, Zi, 0);
    }
    double alpha_k = compute_alpha_word_given_doc(word_id, doc_id);
    return _compute_reduced_log_probability_document(word_id, doc_id, n_k, Zi, alpha_k);
}