#include <iostream>
#include <thread>
#include <string>
#include <set>
#include <unordered_set>
#include <unordered_map> 
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
    int _num_threads;
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
        _num_threads = 1;
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

};