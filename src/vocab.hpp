#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <string>
#include <fstream>
#include "common.hpp"
using namespace std;

namespace cstm {
    class Vocab {
    public:
        unordered_map<id, wstring> _string_by_word_id;
        unordered_map<id, id> _hash_to_id;
        hash<wstring> _hash_func;
        Vocab() {}
        id add_string(wstring &str) {
            id hash = hash_string(str);
            auto itr = _hash_to_id.find(hash);
            if (itr==_hash_to_id.end()) {
                id word_id = _hash_to_id.size();
                _string_by_word_id[word_id] = str;
                _hash_to_id[hash] = word_id;
                return word_id;
            }
            return itr->second;
        }
        id get_word_id(wstring &str) {
            id hash = hash_string(str);
            auto itr = _hash_to_id.find(hash);
            return itr->second;
        }
        id hash_string(wstring &str) {
            return (id)_hash_func(str);
        }
        wstring word_id_to_string(id word_id) {
            auto itr = _string_by_word_id.find(word_id);
            return itr->second;
        }
        wstring token_ids_to_sentence(vector<id> &token_ids) {
            wstring sentence = L"";
            for (const auto &word_id : token_ids) {
                wstring word = word_id_to_string(word_id);
                sentence += word;
                sentence += L" ";
            }
            return sentence;
        }
        int word_exists(wstring &str) {
            id hash = hash_string(str);
            auto itr = _hash_to_id.find(hash);
            return itr != _hash_to_id.end();
        }
        int num_words() {
            return _string_by_word_id.size();
        }
        // serialize
        template<class Archive>
        void serialize(Archive &archive, unsigned int version) {
            archive & _string_by_word_id;
            archive & _hash_to_id;
        }
    };
}