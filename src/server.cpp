#include "server.hpp"
#include "tools.hpp"

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------

//                                                         HELPER FUNCTIONS

// ------------------------------------------------------------------------------------------------------------------------


template<typename T, typename Allocator>
void print_vector(const vector<T, Allocator>& vect, int num_entries)
{
    cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++){
        cout << ", " << vect[i];
    } 
    cout << endl;
}

Server::Server(const helib::Context &context): secret_key(context), public_key(secret_key){
    this->context = &context;
        
    secret_key.GenSecKey();
    helib::addSome1DMatrices(secret_key);

    const helib::EncryptedArray& ea = context.getEA();
    num_slots = ea.size();
    plaintext_modulus = context.getP();

    one_over_two = get_inverse(1,2,plaintext_modulus);
    neg_three_over_two = get_inverse(-3,2,plaintext_modulus);
    neg_one_over_two = get_inverse(-1,2,plaintext_modulus);

    db_set = false;
    
    comparator = unique_ptr<he_cmp::Comparator>(new he_cmp::Comparator(context, he_cmp::UNI, 1, 1, secret_key, false));
}

void Server::GenData(int _num_rows, int _num_cols){
    num_rows = _num_rows;
    num_cols = _num_cols;
    
    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;
    
    encrypted_db = vector<vector<helib::Ctxt>>();
    for(int i = 0; i < num_cols; i++){
        vector<helib::Ctxt> cipher_vector = vector<helib::Ctxt>();
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ptxt<helib::BGV> ptxt(*context);
            
            helib::Ctxt ctxt(public_key);
            
            public_key.Encrypt(ctxt, ptxt);
            
            cipher_vector.push_back(ctxt);
        }
        encrypted_db.push_back(cipher_vector);
    }
    
    db_set = true;
    
}

void Server::SetData(vector<vector<unsigned long>> &db){
    num_cols = db.size();
    if (num_cols == 0){
        throw invalid_argument("ERROR: DB has zero columns! THIS DOES NOT WORK!");
    }
    
    num_rows = db[0].size();
    
    num_compressed_rows = num_rows % num_slots == 0 ? num_rows / num_slots : (num_rows / num_slots) + 1;
    
    encrypted_db = vector<vector<helib::Ctxt>>();
    for(int i = 0; i < num_cols; i++){
        vector<helib::Ctxt> cipher_vector = vector<helib::Ctxt>();
        for (int j = 0; j < num_compressed_rows; j++){
            
            helib::Ptxt<helib::BGV> ptxt(*context);
            
            int entries_left = min(num_slots, num_rows - (j * num_slots));
            for (int k = 0; k < entries_left; k++){
                ptxt[k] = db[i][j*num_slots + k];
            }
            
            helib::Ctxt ctxt(public_key);
            
            public_key.Encrypt(ctxt, ptxt);
            
            cipher_vector.push_back(ctxt);
        }
        encrypted_db.push_back(cipher_vector);
    }
    

    db_set = true;
}

void Server::SetColumnHeaders(vector<string> &headers){
    column_headers = vector<string>();
    for (int i = 0; i < headers.size(); i++){
        column_headers.push_back(headers[i]);
    }
}


helib::Ctxt Server::CountingQuery(bool conjunctive, vector<pair<int, int>>& query){
    if (!db_set){
        throw invalid_argument("ERROR: DB needs to be set to run query");
    }
    
    vector<vector<helib::Ctxt>> cols = filter(query);

    int num_columns = cols[0].size();

    vector<helib::Ctxt> filter_results;
    if (conjunctive){
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else{
        for (int i = 0; i < num_compressed_rows; i++){
            for (int j = 0; j < num_columns; j++){
                AddOneMod2(cols[i][j]);
            }
        }
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (int j = 0; j < num_compressed_rows; j++){
            AddOneMod2(filter_results[j]);
        }
    }

    if (constants::DEBUG){
        print_vector(Decrypt(filter_results[0]));
    }

    helib::Ctxt result = AddMany(filter_results);
    result = SquashCtxt(result);
    return result;
    
}

pair<helib::Ctxt, helib::Ctxt> Server::MAFQuery(int snp, bool conjunctive, vector<pair<int, int>> &query){
    vector<vector<helib::Ctxt>> cols = filter(query);
    int num_columns = cols[0].size();

    vector<helib::Ctxt> filter_results;
    if (conjunctive){
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
    }
    else{
        for (int i = 0; i < num_compressed_rows; i++){
            for (int j = 0; j < num_columns; j++){
                AddOneMod2(cols[i][j]);
            }
        }
        for (int j = 0; j < num_compressed_rows; j++){
            helib::Ctxt temp = MultiplyMany(cols[j]);
            filter_results.push_back(temp);
        }
        for (int j = 0; j < num_compressed_rows; j++){
            AddOneMod2(filter_results[j]);
        }
    }

    vector<helib::Ctxt> indv_MAF = vector<helib::Ctxt>();

    for (int i = 0; i < num_compressed_rows; i++){
        helib::Ctxt clone = encrypted_db[snp][i];
        clone *= filter_results[i];
        indv_MAF.push_back(clone);
    }

    helib::Ctxt freq = AddMany(indv_MAF);
    helib::Ctxt number_of_patients = AddMany(filter_results);

    freq = SquashCtxt(freq);
    number_of_patients = SquashCtxt(number_of_patients);

    number_of_patients.multByConstant(NTL::ZZX(2));
    return pair(freq, number_of_patients);
}

vector<helib::Ctxt> Server::DistrubtionQuery(vector<pair<int, int>>& prs_params){
    
    vector<helib::Ctxt> scores;

    for(int j = 0; j < num_compressed_rows; j++){
        vector<helib::Ctxt> indvs_scores;
        for(pair<int, int> i : prs_params){
            helib::Ctxt temp = encrypted_db[i.first][j];

            temp.multByConstant(NTL::ZZX(i.second));
            indvs_scores.push_back(temp);
        }
        helib::Ctxt score = AddMany(indvs_scores);
        scores.push_back(score);
    }
    return scores;
}


pair<helib::Ctxt, helib::Ctxt> Server::SimilarityQuery(int target_column, vector<helib::Ctxt>& d, int threshold){
    // Compute Normalized Score

    vector<vector<helib::Ctxt>> normalized_scores = vector<vector<helib::Ctxt>>();

    for (int j = 0; j < num_compressed_rows; j++){
        vector<helib::Ctxt> temp = vector<helib::Ctxt>();
        for (size_t i = 0; i < d.size(); i++){
            helib::Ctxt clone = encrypted_db[i][j];
            clone -= d[i];
            clone.square();
            temp.push_back(clone);
        }
        normalized_scores.push_back(temp);
    }
    vector<helib::Ctxt> scores = vector<helib::Ctxt>();

    for (int j = 0; j < num_compressed_rows; j++){
        scores.push_back(AddMany(normalized_scores[j]));
    } 
    if (constants::DEBUG){
        cout << "After scoring:" << endl;
        for (int j = 0; j < num_compressed_rows; j++){
            print_vector(Decrypt(scores[j]));
        } 
    }

    helib::Ctxt thres = Encrypt((unsigned long)threshold);
    vector<helib::Ctxt> predicate = vector<helib::Ctxt>();
    for (int j = 0; j < num_compressed_rows; j++){
        helib::Ctxt res = scores[j];
        comparator->compare(res, scores[j], thres);
        predicate.push_back(res);
    }

    if (constants::DEBUG){
        cout << "After thres:" << endl;
        for (int j = 0; j < num_compressed_rows; j++){
            print_vector(Decrypt(predicate[j]));
        } 
    }
    
    vector<helib::Ctxt> inverse_predicate = vector<helib::Ctxt>();
    for (int j = 0; j < num_compressed_rows; j++){
        helib::Ctxt inv = inverse_predicate[j];
        AddOneMod2(inv);
        inverse_predicate.push_back(inv);
    }

    for (int j = 0; j < num_compressed_rows; j++){
        predicate[j] *= encrypted_db[target_column][j];
        inverse_predicate[j] *= encrypted_db[target_column][j];
    }

    helib::Ctxt count_with = AddMany(predicate);
    helib::Ctxt count_without = AddMany(inverse_predicate);

    count_with = SquashCtxt(count_with);
    count_without = SquashCtxt(count_without);
    
    return pair(count_with, count_without);
}



void Server::AddOneMod2(helib::Ctxt& a){
    //   0 -> 1
    //   1 -> 0
    // f(x)-> -x+1
    
    a.negate();
    a.addConstant(NTL::ZZX(1));
}

helib::Ctxt Server::MultiplyMany(vector<helib::Ctxt>& v){
    int num_entries = v.size();
    int depth = ceil(log2(num_entries));

    for (int d = 0; d < depth; d++){
        int jump_factor = pow(2, d);
        int skip_factor = 2 * jump_factor;

        for (int i = 0; i < num_entries; i+= skip_factor){            
            v[i].multiplyBy(v[i + jump_factor]);
        }
     }
     return v[0];
}

helib::Ctxt Server::AddMany(vector<helib::Ctxt>& v){
    int num_entries = v.size();
    int depth = ceil(log2(num_entries));

    for (int d = 0; d < depth; d++){
        int jump_factor = pow(2, d);
        int skip_factor = 2 * jump_factor;

        for (int i = 0; i < num_entries; i+= skip_factor){
            v[i] += v[i + jump_factor];
        }
     }
     return v[0];
}

helib::Ctxt Server::SquashCtxt(helib::Ctxt& ciphertext, int num_data_elements){
    const helib::EncryptedArray& ea = context->getEA();

    helib::Ctxt result = ciphertext;
    
    for (int i = 1; i < num_data_elements; i++) {
        ea.rotate(ciphertext, -(1));
        result += ciphertext;
    }
    return result;
}

helib::Ctxt Server::SquashCtxtLogTime(helib::Ctxt& ciphertext){
    const helib::EncryptedArray& ea = context->getEA();

    helib::Ctxt result = ciphertext;

    int depth = ceil(log2(num_slots));

    int shift = num_slots;
    int decrease;
    
    for (int i = 0; i < depth; i++) {
        if ((shift % 2) == 1){
            decrease = (shift - 1) / 2;
            shift = (shift + 1) / 2;
        }
        else{
            decrease = shift / 2;
            shift = shift - decrease;
        }

        if (shift == 0)
            continue;

        ea.rotate(result, -(shift));
        ciphertext += result;
        result = ciphertext;
    }
    return ciphertext;
}

helib::Ctxt Server::EQTest(unsigned long a, helib::Ctxt& b){

    helib::Ctxt clone = b;
    helib::Ctxt result = b;
    
    switch (a){
        case 0:
        {
            //f(x) = x^2 / 2 - 3/2 x + 1
            // 0 -> 1
            // 1 -> 0
            // 2 -> 0
            result.square();

            result.multByConstant(NTL::ZZX(one_over_two));
            clone.multByConstant(NTL::ZZX(neg_three_over_two));
            
            result += clone;

            result.addConstant(NTL::ZZX(1));
            return result;
        }
        case 1:
        {
            //f(x) = -x^2 + 2x 
            // 0 -> 0
            // 1 -> 1
            // 2 -> 0
            clone.square();

            result.multByConstant(NTL::ZZX(2));

            result -= clone;
            return result;
        }
        case 2:{
            //f(x) = x^2 / 2 - x / 2
            // 0 -> 0
            // 1 -> 0
            // 2 -> 1
           
            result.square();

            result.multByConstant(NTL::ZZX(one_over_two)); 
            clone.multByConstant(NTL::ZZX(neg_one_over_two));
            
            result += clone;
            return result;
        }
        default:
            cout << "Can't use a value of a other than 0, 1, or 2" << endl;
            throw invalid_argument("ERROR: invalid value for EQTest");
    }
}

vector<vector<helib::Ctxt>> Server::filter(vector<pair<int, int>>& query){
    vector<vector<helib::Ctxt>> feature_cols;

    for(int j = 0; j < num_compressed_rows; j++){
        vector<helib::Ctxt> indv_vector;
        for(pair<int, int> i : query){
            indv_vector.push_back(EQTest(i.second, encrypted_db[i.first][j]));
            if (constants::DEBUG == 2){
                cout << "checking equality to " << i.second << endl;
                cout << "original:";
                print_vector(Decrypt(encrypted_db[i.first][j]));
                cout << "result  :"; 
                print_vector(Decrypt(EQTest(i.second, encrypted_db[i.first][j])));
            }
        }
        feature_cols.push_back(indv_vector);
    }
    return feature_cols;
}

vector<long> Server::Decrypt(helib::Ctxt ctxt){
    helib::Ptxt<helib::BGV> new_plaintext_result(*context);
    secret_key.Decrypt(new_plaintext_result, ctxt);
    
    vector<helib::PolyMod> poly_mod_result = new_plaintext_result.getSlotRepr();
    
    vector<long> result = vector<long>(num_slots);
    
    for (int i = 0; i < num_slots; i++){
        result[i] = (long)poly_mod_result[i];
    }
    
    return result; 
}

helib::Ptxt<helib::BGV> Server::DecryptPlaintext(helib::Ctxt ctxt){
    
    if (constants::DEBUG && ctxt.capacity() < 2){
        cout << "NOISE BOUNDS EXCEEDED!!!" << endl;
    }

    helib::Ptxt<helib::BGV> new_plaintext_result(*context);
    secret_key.Decrypt(new_plaintext_result, ctxt);
    
    return new_plaintext_result;
}

helib::Ctxt Server::Encrypt(unsigned long a){
    helib::Ptxt<helib::BGV> ptxt(*context);
    
    for (int i = 0; i < num_slots; i++)
        ptxt[i] = a;

    helib::Ctxt ctxt(public_key);
    public_key.Encrypt(ctxt, ptxt);
    
    return ctxt; 
}

helib::Ctxt Server::Encrypt(vector<unsigned long> a){
    if (a.size() > num_slots){
        throw invalid_argument("Trying to encrypt vector with too many elements");
    }
    helib::Ptxt<helib::BGV> ptxt(*context);
    
    for (size_t i = 0; i < a.size(); ++i) {
        ptxt[i] = a[i];
    }

    helib::Ctxt ctxt(public_key);
    public_key.Encrypt(ctxt, ptxt);
    
    return ctxt; 
}

helib::Ctxt Server::GetAnyElement(){
    return encrypted_db[0][0];
}


void Server::PrintContext(){
    context->printout();
    cout << endl;
    cout << "Security: " << context->securityLevel() << endl;  
    cout << "Num slots: " << num_slots << endl;
}

void Server::PrintEncryptedDB(bool with_headers){
	if (with_headers){
        vector<int> string_length_count = vector<int>();

        cout << "|";
        for (int i = 0; i < num_cols; i++){
            cout << column_headers[i] << "|";
            string_length_count.push_back(column_headers[i].length());
        }
        cout << endl;
        cout << "--------------";
        for (int j = 0; j < num_compressed_rows; j++){
            
            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (int i = 0; i < num_cols; i++){
                temp_storage.push_back(Decrypt(encrypted_db[i][j]));
            }
            for (int jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++){

                cout << endl << "|";
                
                for (int i = 0; i < num_cols; i++){
                    if(i>0){
                        for (int space = 0; space < string_length_count[i]; space++){
                            cout << " ";
                        }
                        
                    }
                    
                    cout << temp_storage[i][jj];
                    
                    if(i==num_cols-1){
                        for (int space = 0; space < string_length_count[i]; space++){
                            cout << " ";
                        }
                        cout << "|";
                    }
                }
            }
        }
    }
    else{
        for (int j = 0; j < num_compressed_rows; j++){
            
            vector<vector<long>> temp_storage = vector<vector<long>>();
            for (int i = 0; i < num_cols; i++){
                temp_storage.push_back(Decrypt(encrypted_db[i][j]));
            }
            for (int jj = 0; jj < min(num_slots, num_rows - (j * num_slots)); jj++){

                cout << endl << "|";
                
                for (int i = 0; i < num_cols; i++){
                    if(i>0){
                        cout << " ";
                    }
                    
                    cout << temp_storage[i][jj];
                    
                    if(i==num_cols-1){
                        cout << "|";
                    }
                }
            }
        }
    }
    cout << endl;
}

int Server::GetSlotSize(){
    return num_slots;
}

// IMPORTED FROM HELIB SOURCE CODE
inline long estimateCtxtSize(const helib::Context& context, long offset)
{
  // Return in bytes.

  // We assume that the size of each element in the DCRT is BINIO_64BIT

  // sizeof(BINIO_EYE_CTXT_BEGIN) = 4;
  // BINIO_32BIT = 4
  // sizeof(long) = BINIO_64BIT = 8
  // xdouble = s * sizeof(long) = 2 * BINIO_64BIT = 16

  // We assume that primeSet after encryption is context.ctxtPrimes
  // We assume we have exactly 2 parts after encryption
  // We assume that the DCRT prime set is the same as the ctxt one

  long size = 0;

  // Header metadata
  size += 24;

  // Begin eye-catcher
  size += 4;

  // Begin Ctxt metadata
  // 64 = header_size = ptxtSpace (long) + intFactor (long) + ptxtMag (xdouble)
  //                    + ratFactor (xdouble) + noiseBound (xdouble)
  size += 64;

  // primeSet.write(str);
  // size of set (long) + each prime (long)
  size += 8 + context.getCtxtPrimes().card() * 8;

  // Begin Ctxt content size
  // write_raw_vector(str, parts);
  // Size of the parts vector (long)
  size += 8;

  long part_size = 0;
  // Begin CtxtPart size

  // skHandle.write(str);
  // powerOfS (long) + powerOfX (long) + secretKeyID (long)
  part_size += 24;

  // Begin DCRT size computation

  // this->DoubleCRT::write(str);
  // map.getIndexSet().write(str);
  // size of set (long) + each prime (long)
  part_size += 8 + context.getCtxtPrimes().card() * 8;

  // DCRT data write as write_ntl_vec_long(str, map[i]);
  // For each prime in the ctxt modulus chain
  //    size of DCRT column (long) + size of each element (long) +
  //    size of all the slots (column in DCRT) (PhiM long elements)
  long dcrt_size = (8 + 8 * context.getPhiM()) * context.getCtxtPrimes().card();

  part_size += dcrt_size;

  // End DCRT size
  // End CtxtPart size

  size += 2 * part_size; // 2 * because we assumed 2 parts
  // End Ctxt content size

  // End eye-catcher
  size += 4;

  return size + offset;
}

int Server::StorageOfOneElement(){
    if (!db_set){
        throw invalid_argument("ERROR: DB needs to be set to get storage cost");
    }

    return estimateCtxtSize(*context, 0);
}
