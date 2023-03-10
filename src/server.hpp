#pragma once

#include <iostream>
#include <helib/helib.h>
#include <memory>
#include <string>
#include <vector>
#include "globals.hpp"
#include "comparator.hpp"
#include "tools.hpp"

#define MAX_NUMBER_BITS 4
#define NOISE_THRES 2
#define WARN false

using namespace std;

class Server{
public:
    
    //Setup
    Server(const helib::Context &context);
    void GenData(int _num_rows, int _num_cols);
    void SetData(vector<vector<unsigned long>> &db);
    void SetColumnHeaders(vector<string> &headers);
    
    //Querries
    helib::Ctxt CountingQuery(bool conjunctive, vector<pair<int, int>>& query);
    pair<helib::Ctxt, helib::Ctxt> MAFQuery(int snp, bool conjunctive, vector<pair<int, int>> &query);
    vector<helib::Ctxt> DistrubtionQuery(vector<pair<int, int>>& prs_params);
    pair<helib::Ctxt, helib::Ctxt> SimilarityQuery(int target_column, vector<helib::Ctxt>& d, int threshold);

    
    void AddOneMod2(helib::Ctxt& a);
    helib::Ctxt MultiplyMany(vector<helib::Ctxt>& v);
    helib::Ctxt AddMany(vector<helib::Ctxt>& v);
    helib::Ctxt SquashCtxt(helib::Ctxt& ciphertext, int num_data_entries = 10);
    helib::Ctxt SquashCtxtLogTime(helib::Ctxt& ciphertext);
    helib::Ctxt EQTest(unsigned long a, helib::Ctxt& b);
    vector<vector<helib::Ctxt>> filter(vector<pair<int, int>>& query);
    
    //Encrypt / Decrypt Methods
    helib::Ptxt<helib::BGV> DecryptPlaintext(helib::Ctxt ctxt);
    vector<long> Decrypt(helib::Ctxt ctxt);
    helib::Ctxt Encrypt(unsigned long a);
    helib::Ctxt Encrypt(vector<unsigned long> a);
    helib::Ctxt GetAnyElement();
    
    void PrintContext();
    void PrintEncryptedDB(bool with_headers);
    int GetSlotSize();

    int StorageOfOneElement();
    
private:
    const helib::Context* context;
    helib::SecKey secret_key;
    const helib::PubKey& public_key;

    unique_ptr<he_cmp::Comparator> comparator;
    
    bool db_set;
    
    int num_rows;
    int num_cols;
    int num_compressed_rows;
    int num_slots;
    
    vector<vector<helib::Ctxt>> encrypted_db; 
    vector<string> column_headers;
    
    int one_over_two;
    int neg_three_over_two;
    int neg_one_over_two;

    int plaintext_modulus;
};