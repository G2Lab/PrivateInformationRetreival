/* Copyright (C) 2019-2021 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

// This is a sample program for education purposes only.
// It attempts to show the various basic mathematical
// operations that can be performed on both ciphertexts
// and plaintexts.

#include <iostream>

#include <helib/helib.h>
#include "client.hpp"
#include "server.hpp"
#include "globals.hpp"


template<typename T, typename Allocator>
void print_vector(const vector<T, Allocator>& vect, int num_entries)
{
    cout << vect[0];
    for (int i = 1; i < min((int)vect.size(), num_entries); i++){
        cout << ", " << vect[i];
    } 
    cout << endl;
}

int main()
{
    /*  Example of BGV scheme  */
    
    
    std::cout << "Initialising context object..." << std::endl;

    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(constants::M)
                               .p(constants::P)
                               .r(constants::R)
                               .bits(constants::BITS)
                               .c(constants::C)
                               .build();

    Server server = Server(context);
    Client client = Client(context);
    
    server.PrintContext();
    
    vector<vector<unsigned long>> fake_db = vector<vector<unsigned long>> {
          {0, 0, 0, 1, 1, 1, 2, 2, 2, 0},
          {0, 1, 2, 0, 1, 2, 0, 1, 2, 1},
          {0, 0, 0, 0, 1, 1, 0, 1, 0, 0}
     };
    vector<string> headers = vector<string>{"snp1", "snp2","ALS"};
    
    server.SetData(fake_db);
    server.SetColumnHeaders(headers);

    cout << "Printing DB" << endl;
    cout << "-----------------------------------------------------" << endl;
    server.PrintEncryptedDB(true);

    cout << "Running sample queries:" << endl;
    cout << "-----------------------------------------------------" << endl;
    
    vector<pair<int, int>> query;

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 1)" << endl;
    query = vector<pair<int, int>>{pair(0,0), pair(1,1)};
    helib::Ctxt result = server.CountingQuery(true, query);
    cout << "Count: " << server.Decrypt(result)[0] << endl;

    cout << "Running Counting query (snp 0 = 0 and snp 1 = 2)" << endl;
    query = vector<pair<int, int>>{pair(0,0), pair(1,2)};
    result = server.CountingQuery(true, query);
    cout << "Count: " << server.Decrypt(result)[0] << endl;

    cout << "Running Counting query (snp 0 = 0 or snp 1 = 1)" << endl;
    query = vector<pair<int, int>>{pair(0,0), pair(1,1)};
    result = server.CountingQuery(false, query);
    cout << "Count: " << server.Decrypt(result)[0] << endl;

    cout << "Running MAF query filter (snp 0 = 0 or snp 1 = 1), target snp = 0" << endl;
    query = vector<pair<int, int>>{pair(0,0), pair(1,1)};
    pair<helib::Ctxt, helib::Ctxt> result_pair = server.MAFQuery(0, false, query);
    cout << "Nom: " << server.Decrypt(result_pair.first)[0] << endl;
    cout << "Dom: " << server.Decrypt(result_pair.second)[0] << endl;
    double AF = (double)(server.Decrypt(result_pair.first)[0]) / (double)(server.Decrypt(result_pair.second)[0]);
    cout << "Computed MAF: " << min(AF, 1 - AF) << endl;


    cout << "Running Distribution query (alpha: snp 0 = 2 and snp 1 = 5)" << endl;
    query = vector<pair<int, int>>{pair(0,2), pair(1,5)};
    vector<helib::Ctxt> results_distribution = server.DistrubtionQuery(query);
    print_vector(server.Decrypt(results_distribution[0]));



    cout << "Running Similarity query (d: snp 0 = 2 and snp 1 = 2, target = 0)" << endl;
    vector<helib::Ctxt> d = vector<helib::Ctxt>();
    d.push_back(server.Encrypt(2));
    d.push_back(server.Encrypt(2));
    pair<helib::Ctxt, helib::Ctxt> result_sim = server.SimilarityQuery(0, d, 2);
    cout << "Count with target:   " << server.Decrypt(result_sim.first)[0] << endl;
    cout << "Count without target:" << server.Decrypt(result_sim.second)[0] << endl;
    
    cout << "Printing DB" << endl;
    cout << "-----------------------------------------------------" << endl;
    server.PrintEncryptedDB(true);

    return 0;
}