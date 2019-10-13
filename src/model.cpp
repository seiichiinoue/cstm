#include <iostream>
#include <thread>
#include <string>
#include <set>
#include <unordered_set>
#include <unordered_map> 

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "cstm.hpp"
#include "vocab.hpp"

template<typename T>
struct multiset_comparator {
    bool operator()(const pair<id, T> &a, const pair<id, T> &b) {
        return a.second > b.second;
    }
};

void split_string_by(const string &str, char delim, vector<string> &elems) {
    elems.clear();
    string item;
    for (char ch : str) {
        if (ch == delim) {
            if (!item.empty()) {
                elems.push_back(item);
            }
            item.clear();
        } else {
            item += ch;
        }
    }
    if (!item.empty()) {
        elems.push_back(item);
    }
}
void split_word_by(const wstring &str, wchar_t delim, vector<wstring> &elems) {
    elems.clear();
    wstring item;
    for (wchar_t ch : str) {
        if (ch == delim) {
            if (!item.empty()) {
                elems.push_back(item);
            }
            item.clear();
        } else {
            item += ch;
        }
    }
    if (!item.empty()) {
        elems.push_back(item);
    }
}

class CSTMTrainer {
public:
    CSTM *_cstm;
    Vocab *_vocab;
    vector<vector<vector<id>>> _dataset;
    vector<unordered_set<id>> _word_ids_in_doc;
    vector<int> _sum_word_frequency;    // word frequency per document
    vector<id> _random_word_ids;
    vector<int> _random_doc_ids;
    unordered_map<id, unordered_set<int>> _docs_containing_word;    // word -> [doc_ids]
    unordered_map<id, int> _word_frequency;
    unordered_map<string, int> _doc_filename_to_id;
    unordered_map<int, string> _doc_id_to_filename;
    
    double *_old_vec_copy;
    double *_new_vec_copy;
    double **_old_vec_copy_thread;
    double **_new_vec_copy_thread;
    double *_old_alpha_words;
    double *_Zi_cache;
    
    int _ndim_d;
    int _ignored_vocabulary_size;
    std::thread *_doc_threads;

    // stat
    // times of acceptance 
    int _num_acceptance_doc;    
    int _num_acceptance_word;
    int _num_acceptance_alpha0;
    // times of rejection 
    int _num_rejection_doc;
    int _num_rejection_word;
    int _num_rejection_alpha0;
    // times of sampling
    int _num_word_vec_sampled;
    int _num_doc_vec_sampled;
    
    int _random_sampling_word_index;
    int _random_sampling_doc_index;
    unordered_map<id, int> _num_updates_word;
    unordered_map<int, int> _num_updates_doc;

    CSTMTrainer() {
        setlocale(LC_CTYPE, "ja_JP.UTF-8");
        _cstm = new CSTM();
        _vocab = new Vocab();
        _old_vec_copy = NULL;
        _new_vec_copy = NULL;
        _old_vec_copy_thread = NULL;
        _new_vec_copy_thread = NULL;
        _old_alpha_words = NULL;
        _Zi_cache = NULL;
        _doc_threads = NULL;
        _ndim_d = NULL;
        _ignored_vocabulary_size = 0;
        reset_statistics();
        _random_sampling_doc_index = 0;
        _random_sampling_word_index = 0;
    }
    ~CSTMTrainer() {}

    void reset_statistics() {
        int _num_acceptance_doc = 0;    
        int _num_acceptance_word = 0;
        int _num_acceptance_alpha0 = 0;
        int _num_rejection_doc = 0;
        int _num_rejection_word = 0;
        int _num_rejection_alpha0 = 0;
        int _num_word_vec_sampled = 0;
        int _num_doc_vec_sampled = 0;
    }

    void prepare() {
        int num_docs = _dataset.size();
        int vocabulary_size = _word_frequency.size();
        // CSTM
        _cstm->initialize(_ndim_d, vocabulary_size, num_docs);
        for (int doc_id=0; doc_id<num_docs; ++doc_id) {
            vector<vector<id>> &dataset = _dataset[doc_id];
            for (int data_index=0; data_index<dataset.size(); ++data_index) {
                vector<id> &word_ids = dataset[data_index];
                for (const id word_id : word_ids) {
                    _cstm->add_word(word_id, doc_id);
                }
            }
            _num_updates_doc[doc_id] = 0;
            _random_doc_ids.push_back(doc_id);
        }
        _cstm->prepare();
        // Zi
        for (int doc_id=0; doc_id<num_docs; ++doc_id) {
            _cstm->update_Zi(doc_id);
        }
        // random sampling of words
        for (id word_id=0; word_id<vocabulary_size; ++word_id) {
            _num_updates_word[word_id] = 0;
            int count = _cstm->get_word_count(word_id);
            if (count <= _cstm->get_ignore_word_count()) {
                _ignored_vocabulary_size += 1;
                continue;
            }
            _random_word_ids.push_back(word_id);
        }
        std::shuffle(_random_word_ids.begin(), _random_word_ids.end(), sampler::mt);
        std::shuffle(_random_doc_ids.begin(), _random_doc_ids.end(), sampler::mt);
    }
    int add_document(string filepath) {
        wifstream ifs(filepath.c_str());
        assert(ifs.fail() == false);
        // add document
        int doc_id = _dataset.size();
        _dataset.push_back(vector<vector<id>>());
        _word_ids_in_doc.push_back(unordered_set<id>());
        _sum_word_frequency.push_back(0);
        // read file
        wstring sentence;
        vector<wstring> sentences;
        while (getline(ifs, sentence) && !sentence.empty()) {
            sentences.push_back(sentence);
        }
        for (wstring &sentence : sentences) {
            vector<wstring> words;
            split_word_by(sentence, L' ', words);
            add_sentence_to_doc(words, doc_id);
        }
        vector<string> components;
        split_string_by(filepath, '/', components);
        string filename = components.back();
        _doc_filename_to_id[filename] = doc_id;
        _doc_id_to_filename[doc_id] = filename;
        return doc_id;
    }
    void add_sentence_to_doc(vector<wstring> &words, int doc_id) {
        if (words.size() > 0) {
            vector<vector<id>> &dataset = _dataset[doc_id];
            _sum_word_frequency[doc_id] += words.size();
            vector<id> word_ids;
            for (auto word : words) {
                if (word.size() == 0) {
                    continue;
                }
                id word_id = _vocab->add_string(word);
                word_ids.push_back(word_id);
                unordered_set<int> &docs = _docs_containing_word[word_id];
                docs.insert(doc_id);
                unordered_set<id> &word_set = _word_ids_in_doc[doc_id];
                word_set.insert(word_id);
                _word_frequency[word_id] += 1;
            }
            dataset.push_back(word_ids);
        }
    }
    bool is_doc_contain_word(int doc_id, id word_id) {
        unordered_set<int> &set = _docs_containing_word[word_id];
        auto itr = set.find(doc_id);
        return itr != set.end();
    }
    int get_num_documents() {
        return _dataset.size();
    }
    int get_vocabulary_size() {
        return _word_frequency.size();
    }
    int get_ignored_vocabulary_size() {
        return _ignored_vocabulary_size;
    }
    int get_ndim_d() {
        return _cstm->_ndim_d;
    }
    int get_sum_word_frequency() {
        return std::accumulate(_sum_word_frequency.begin(), _sum_word_frequency.end(), 0);
    }
    int get_num_word_vec_sampled() {
        return _num_word_vec_sampled;
    }
    int get_num_doc_vec_sampled() {
        return _num_doc_vec_sampled;
    }
    double get_alpha0() {
        return _cstm->_alpha0;
    }
    double get_mh_acceptance_rate_for_doc_vector() {
        return _num_acceptance_doc / (double)(_num_acceptance_doc + _num_rejection_doc);
    }
    double get_mh_acceptance_rate_for_word_vector() {
        return _num_acceptance_word / (double)(_num_acceptance_word + _num_rejection_word);
    }
    double get_mh_acceptance_rate_for_alpha0() {
        return _num_acceptance_alpha0 / (double)(_num_acceptance_alpha0 + _num_rejection_alpha0);
    }
    double *get_word_vector(id word_id) {
        double *old_vec = _cstm->get_word_vector(word_id);
        std::memcpy(_old_vec_copy, old_vec, _cstm->_ndim_d * sizeof(double));
        return _old_vec_copy;
    }
    double *get_doc_vector(int doc_id) {
        double *old_vec = _cstm->get_doc_vector(doc_id);
        std::memcpy(_old_vec_copy, old_vec, _cstm->_ndim_d * sizeof(double));
        return _old_vec_copy;
    }
    double *draw_word_vector(double *old_vec) {
        double *new_vec = _cstm->draw_word_vector(old_vec);
        std::memcpy(_new_vec_copy, new_vec, _cstm->_ndim_d * sizeof(double));
        return _new_vec_copy;
    }
    double *draw_doc_vector(double *old_vec) {
        double *new_vec = _cstm->draw_doc_vector(old_vec);
        std::memcpy(_new_vec_copy, new_vec, _cstm->_ndim_d * sizeof(double));
        return _new_vec_copy;
    }
    void set_ignore_word_count(int count){
        _cstm->set_ignore_word_count(count);
    }
    void set_ndim_d(int ndim_d){
        _ndim_d = ndim_d;
    }
    void set_alpha0(double alpha0){
        _cstm->_alpha0 = alpha0;
    }
    void set_sigma_u(double sigma_u){
        _cstm->_sigma_u = sigma_u;
    }
    void set_sigma_phi(double sigma_phi){
        _cstm->_sigma_phi = sigma_phi;
    }
    void set_sigma_alpha0(double sigma_alpha0){
        _cstm->_sigma_alpha0 = sigma_alpha0;
    }
    void set_gamma_alpha_a(double gamma_alpha_a){
        _cstm->_gamma_alpha_a = gamma_alpha_a;
    }
    void set_gamma_alpha_b(double gamma_alpha_b){
        _cstm->_gamma_alpha_b = gamma_alpha_b;
    }
    double compute_log_likelihood_data() {
        double log_pw = 0;
        int n = 0;
        for (int doc_id=0; doc_id<get_num_documents(); ++doc_id) {
            unordered_set<id> &word_ids = _word_ids_in_doc[doc_id];
            log_pw += _cstm->compute_log_probability_document_given_words(doc_id, word_ids);
        }
        return log_pw;
    }
    double compute_perplexity() {
        double log_pw = 0;
        int n = 0;
        for (int doc_id=0; doc_id<get_num_documents(); ++doc_id) {
            unordered_set<id> &word_ids = _word_ids_in_doc[doc_id];
            log_pw += _cstm->compute_log_probability_document_given_words(doc_id, word_ids);
        }
        return exp(-log_pw / get_sum_word_frequency());
    }
    void update_all_Zi() {
        for (int doc_id=0; doc_id<get_num_documents(); ++doc_id) {
            _cstm->update_Zi(doc_id);
        }
    }
    void perform_mh_sampling_word() {
        // choose word vector for update
        int limit = (int)((get_vocabulary_size() - _ignored_vocabulary_size) / (double)get_num_documents());
        if (_random_sampling_word_index + limit >= _random_word_ids.size()) {
            std::shuffle(_random_word_ids.begin(), _random_word_ids.end(), sampler::mt);
            _random_sampling_word_index = 0;
        }
        for (int i=0; i<limit; ++i) {
            id word_id = _random_word_ids[i + _random_sampling_word_index];
            double *old_vec = get_word_vector(word_id);
            double *new_vec = draw_word_vector(old_vec);
            accept_word_vector_if_needed(new_vec, old_vec, word_id);
            _num_word_vec_sampled += 1;
            _num_updates_word[word_id] += 1;
        }
        _random_sampling_word_index += limit;
    }
    bool accept_word_vector_if_needed(double *new_word_vec, double *old_word_vec, id word_id) {
        auto itr = _docs_containing_word.find(word_id);
        assert(itr != _docs_containing_word.end());
        unordered_set<int> &docs = itr->second;
        assert(docs.size() > 0);
        // likelihood of old word vec
        double log_pw_old = 0;
        for (int doc_id=0; doc_id<get_num_documents(); ++doc_id) {
            double old_alpha_word = _cstm->compute_alpha_word_given_doc(word_id, doc_id);
            double old_Zi = _cstm->get_Zi(doc_id);
            int n_k = _cstm->get_word_count_in_doc(word_id, doc_id);
            log_pw_old += _cstm->_compute_reduced_log_probability_document(word_id, doc_id, n_k, old_Zi, old_alpha_word);
            _old_alpha_words[doc_id] = old_alpha_word;
        }
        // likelihood of new word vec
        double g0 = _cstm->get_g0_of_word(word_id);
        double log_pw_new = 0;
        for (int doc_id=0; doc_id<get_num_documents(); ++doc_id) {
            double *doc_vec = _cstm->get_doc_vector(doc_id);
            double new_alpha_word = _cstm->_compute_alpha_word(new_word_vec, doc_vec, g0);
            double old_alpha_word = _old_alpha_words[doc_id];
            // simplify calculation of Zi
            double old_Zi = _cstm->get_Zi(doc_id);
            double new_Zi = old_Zi - old_alpha_word + new_alpha_word;
            int n_k = _cstm->get_word_count_in_doc(word_id, doc_id);
            log_pw_new += _cstm->_compute_reduced_log_probability_document(word_id, doc_id, n_k, new_Zi, new_alpha_word);
            _Zi_cache[doc_id] = new_Zi;
        }
        // prior distribution
        double log_prior_old = _cstm->compute_log_prior_vector(old_word_vec);
        double log_prior_new = _cstm->compute_log_prior_vector(new_word_vec);
        // acceptance rate
        double log_acceptance_rate = log_pw_new + log_prior_new - log_pw_old - log_prior_old;
        double acceptance_ratio = std::min(1.0, exp(log_acceptance_rate));
        double bernoulli = sampler::uniform(0, 1);
        if (bernoulli <= acceptance_ratio) {
            _num_acceptance_word += 1;
            // set new vector
            _cstm->set_word_vector(word_id, new_word_vec);
            for (int doc_id=0; doc_id<get_num_documents(); ++doc_id) {
                _cstm->set_Zi(doc_id, _Zi_cache[doc_id]);
            }
            return true;
        }
        _num_rejection_word += 1;
        return false;
    }
    void perform_mh_sampling_document() {
        // choose doc vector for update
        if (_random_sampling_doc_index >= _random_doc_ids.size()) {
            std::shuffle(_random_doc_ids.begin(), _random_doc_ids.end(), sampler::mt);
            _random_sampling_doc_index = 0;
        }
        int doc_id = _random_doc_ids[_random_sampling_doc_index];
        double *old_vec = get_doc_vector(doc_id);
        double *new_vec = draw_doc_vector(old_vec);
        accept_document_vector_if_needed(new_vec, old_vec, doc_id);
        _num_doc_vec_sampled += 1;
        _num_updates_doc[doc_id] += 1;
        _random_sampling_doc_index += 1;
        return;
    }
    bool accept_document_vector_if_needed(double *new_doc_vec, double *old_doc_vec, int doc_id) {
        double original_Zi = _cstm->get_Zi(doc_id);
        // likelihood of old doc vector
        double log_pw_old = _cstm->compute_log_probability_document(doc_id);
        // likelihood of new doc vector
        _cstm->set_doc_vector(doc_id, new_doc_vec);
        _cstm->update_Zi(doc_id);
        double log_pw_new = _cstm->compute_log_probability_document(doc_id);
        // prior distribution
        double log_prior_old = _cstm->compute_log_prior_vector(old_doc_vec);
        double log_prior_new = _cstm->compute_log_prior_vector(new_doc_vec);
        // acceptance rate
        double log_acceptance_rate = log_pw_new + log_prior_new - log_pw_old - log_prior_old;
        double acceptance_ratio = std::min(1.0, exp(log_acceptance_rate));
        double bernoulli = sampler::uniform(0, 1);
        if (bernoulli <= acceptance_ratio) {
            _num_acceptance_doc += 1;
            return true;
        }
        // undo
        _cstm->set_doc_vector(doc_id, old_doc_vec);
        _cstm->set_Zi(doc_id, original_Zi);
        _num_rejection_doc += 1;
        return false;
    }
    void perform_mh_sampling_alpha0() {
        int doc_id = sampler::uniform_int(0, _cstm->_num_documents - 1);
        double old_alpha0 = _cstm->get_alpha0();
        double new_alpha0 = _cstm->draw_alpha0(old_alpha0);
        accept_alpha0_if_needed(new_alpha0, old_alpha0);
    }
    bool accept_alpha0_if_needed(double new_alpha0, double old_alpha0) {
        int num_docs = _dataset.size();
        // likelihood of old a0
        double log_pw_old = 0;
        for (int doc_id=0; doc_id<num_docs; ++doc_id) {
            log_pw_old += _cstm->compute_log_probability_document(doc_id);
            _Zi_cache[doc_id] = _cstm->get_Zi(doc_id);
        }
        // likelihood of new a0
        _cstm->set_alpha0(new_alpha0);
        update_all_Zi();
        double log_pw_new = 0;
        for (int doc_id=0; doc_id<num_docs; ++doc_id) {
            log_pw_new += _cstm->compute_log_probability_document(doc_id);
        }
        // prior distribution
        double log_prior_old = _cstm->compute_log_prior_alpha0(old_alpha0);
        double log_prior_new = _cstm->compute_log_prior_alpha0(new_alpha0);
        // acceptance rate
        double log_acceptance_rate = log_pw_new + log_prior_new - log_pw_old - log_prior_old;
        double acceptance_ratio = std::min(1.0, exp(log_acceptance_rate));
        double bernoulli = sampler::uniform(0, 1);
        if (bernoulli <= acceptance_ratio) {
            _num_acceptance_alpha0 += 1;
            return true;
        }
        _num_rejection_alpha0 += 1;
        // undo
        _cstm->set_alpha0(old_alpha0);
        for (int doc_id=0; doc_id<num_docs; ++doc_id) {
            _cstm->set_Zi(doc_id, _Zi_cache[doc_id]);
        }
        return false;
    }
};

DEFINE_int32(NDIM, 20, "Number of hidden size");
DEFINE_int32(IGNORE, 10, "Number of ignore word");
DEFINE_int32(EPOCH, 100000, "Num of epoch");
DEFINE_string(DATA, "data/wiki-wakati.txt", "Data file");

int main() {
    CSTMTrainer trainer;
    trainer.set_ndim_d(FLAGS_NDIM);
    trainer.set_sigma_u(0.02);
    trainer.set_sigma_phi(0.04);
    trainer.set_sigma_alpha0(0.2);
    trainer.set_gamma_alpha_a(5);
    trainer.set_gamma_alpha_b(500);
    trainer.set_ignore_word_count(FLAGS_IGNORE);
    
    // read file
    int doc_id = trainer.add_document(FLAGS_DATA);
    // initialize
    trainer.prepare();
    // training
    for (int i=0; i<FLAGS_EPOCH; ++i) {
        trainer.perform_mh_sampling_document();
        trainer.perform_mh_sampling_word();
        trainer.perform_mh_sampling_alpha0();
        trainer.reset_statistics();
    }
    // save model
    return 0;
}