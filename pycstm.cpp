#include <boost/python.hpp>
#include "src/model.cpp"

class PyCSTM {
public:
    CSTM *_cstm;
    Vocab *_vocab;
    unordered_map<id, int> _word_frequency;
    vector<unordered_set<id>> _word_ids_in_doc;
    vector<int> _sum_word_frequency;
    unordered_map<id, unordered_set<int>> _docs_containing_word;
    double *_vec_copy;
    unordered_map<string, int> _doc_filename_to_id;
    unordered_map<int, string> _doc_id_to_filename;
    PyCSTM(string filename) {
        int ndim_d = get_ndim_d();
        _vec_copy = new double[ndim_d];
    }
    ~PyCSTM() {
        delete _cstm;
        delete _vocab;
        delete[] _vec_copy;
    }
    bool load(string filename) {
        std::ifstream ifs(filename);
        if (ifs.good()) {
            _vocab = new Vocab();
            _cstm = new CSTM();
            boost::archive::binary_iarchive iarchive(ifs);
            iarchive >> *_vocab;
            iarchive >> *_cstm;
            iarchive >> _word_frequency;
            iarchive >> _word_ids_in_doc;
            iarchive >> _docs_containing_word;
            iarchive >> _sum_word_frequency;
            iarchive >> _doc_filename_to_id;
            iarchive >> _doc_id_to_filename;
            return true;
        }
        return false;
    }
    int get_num_documents() {
        return _cstm->_num_documents;
    }
    int get_vocabulary_size() {
        return _cstm->_vocabulary_size;
    }
    int get_ndim_d() {
        return _cstm->_ndim_d;
    }
    int get_sum_word_frequency() {
        return std::accumulate(_sum_word_frequency.begin(), _sum_word_frequency.end(), 0);
    }
    double get_alpha0() {
        return _cstm->_alpha0;
    }
    int get_doc_id_by_filename(string filename) {
        auto itr = _doc_filename_to_id.find(filename);
        return itr->second;
    }
    string get_doc_filename_by_id(int doc_id) {
        auto itr = _doc_id_to_filename.find(doc_id);
        return itr->second;
    }
    double *get_word_vector(id word_id) {
        double *original_vec = _cstm->get_word_vector(word_id);
        std::memcpy(_vec_copy, original_vec, _cstm->_ndim_d * sizeof(double));
        return _vec_copy;
    }
    double *get_doc_vector(int doc_id) {
        double *original_vec = _cstm->get_doc_vector(doc_id);
        std::memcpy(_vec_copy, original_vec, _cstm->_ndim_d * sizeof(double));
        return _vec_copy;
    }
    python::list _convert_vector_to_list(double *vector) {
        python::list vector_list;
        for (int i=0; i<_cstm->_ndim_d; ++i) {
            vector_list.append(vector[i]);
        }
        return vector_list;
    }
    python::list get_word_vector_by_word(wstring str) {
        id word_id = _vocab->get_word_id(str);
        return get_word_vector_by_id(word_id);
    }
    python::list get_word_vector_by_id(id word_id) {
        double *vector = get_word_vector(word_id);
        return _convert_vector_to_list(vector);
    }
    python::list get_word_vectors() {
        python::list vector_array;
        for (id word_id=0; word_id<get_vocabulary_size(); ++word_id) {
            double *vector = get_word_vector(word_id);
            vector_array.append(_convert_vector_to_list(vector));
        }
        return vector_array;
    }
    python::list get_doc_vectors() {
        python::list vector_array;
        for (int doc_id=0; doc_id<get_num_documents(); ++doc_id) {
            python::list vector_list;
            double *vector = get_doc_vector(doc_id);
            for (int i=0; i<_cstm->_ndim_d; ++i) {
                vector_list.append(vector[i]);
            }
            vector_array.append(vector_list);
        }
        return vector_array;
    }
};

BOOST_PYTHON_MODULE(pycstm) {
    python::class_<PyCSTM>("cstm")
    .def("get_num_documents", &PyCSTM::get_num_documents)
    .def("get_vocabulary_size", &PyCSTM::get_vocabulary_size)
    .def("get_ndim_d", &PyCSTM::get_ndim_d)
    .def("get_sum_word_frequency", &PyCSTM::get_sum_word_frequency)
    .def("get_alpha0", &PyCSTM::get_alpha0)
    .def("get_doc_id_by_filename", &PyCSTM::get_doc_id_by_filename)
    .def("get_doc_filename_by_id", &PyCSTM::get_doc_filename_by_id)
    .def("get_word_vector_by_word", &PyCSTM::get_word_vector_by_word)
    .def("get_word_vectors", &PyCSTM::get_word_vectors)
    .def("get_doc_vectors", &PyCSTM::get_doc_vectors)
    .def("get_word_vector_by_id", &PyCSTM::get_word_vector_by_id)
    .def("load", &PyCSTM::load);
}