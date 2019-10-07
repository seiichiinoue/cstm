#include <unordered_set>
#include <cassert>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <numeric>
#include <limits>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <boost/random.hpp>
#include "sampler.hpp"

#define PI 3.14159265358979323846
#define NDIM_D 2                // 文書,単語ベクトルの次元数
#define SIGMA_U 0.01            // 文書ベクトののランダムウォーク幅
#define SIGMA_PHI 0.02          // 単語ベクトルのランダムウォーク幅
#define SIGMA_ALPHA 0.2         // a0のランダムウォーク幅
#define GAMMA_ALPHA_A 5         // a0のガンマ事前分布のハイパーパラメータ
#define GAMMA_ALPHA_B 1         // a0がガンマ事前分のはハイパーパラメータ

class CSTM {
public:
    int **_n_k;                 // 文書ごとの単語の出現頻度
    int *_sum_n_k;              // 文書ごとの単語の出現頻度の総和
    int *_word_count;
    double *_Zi;
    double *_g0;                // 単語のデフォルト確率
    double **_word_vectors;     // 単語ベクトル
    double **_doc_vectors;      // 文書ベクトル
    int _ndim_d;
    int _num_documents;
    int _vocabrary_size;
    int _sum_word_frequency;    // 全単語の出現回の総和
    int _ignore_word_count;
    double _sigma_u;
    double _sigma_phi;
    double _sigma_alpha0;
    double _gamma_alpha_a;
    double _gamma_alpha_b;
    double _alpha0;
    double *_tmp_vec;
    double *_log_likelihood_first_term;
    normal_distribution<double> _standard_normal_distribution;
    normal_distribution<double> _noise_word;
    normal_distribution<double> _noise_doc;
    normal_distribution<double> _noise_alpha0;

    CSTM();
    ~CSTM();

};