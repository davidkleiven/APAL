#include "track_value_logger.hpp"
#include <fstream>
#include <string>
#include <iostream>

using namespace std;
bool read_from_file(){
    string fname("track_values_logger.txt");

    ofstream out(fname);

    if (!out.good()){
        throw runtime_error("Could not open file in read_from_file!");
    }

    out << "# Iter,Item1,Item2,Item3,Item4,Item5\n";
    out << "1,0.5,0.2,0.6,0.8,1.0\n";
    out.close();

    TrackValueLogger logger(fname);
    const vector<string> from_file = logger.get_keys();

    vector<string> expected;
    expected.push_back("Item1");
    expected.push_back("Item2");
    expected.push_back("Item3");
    expected.push_back("Item4");
    expected.push_back("Item5");

    if (expected.size() != from_file.size()){
        cout << "Read " << from_file.size() << " from file. Expected " << expected.size() << endl;
        return false;
    }

    for (unsigned int i=0;i<expected.size();i++){
        if (expected[i] != from_file[i]){
            cout << "Key " << i << ": " << expected[i] << "!=" << from_file[i] << endl;
            return false;
        }
    }
    return true;
}

void init_keys_from_entry(){
    track_value_row_t entry;
    entry["Value1"] = 0.1;
    entry["Value2"] = -1.2;
    entry["Value3"] = 0.4;

    string fname("init_key_from_entry.txt");
    TrackValueLogger logger(fname);
    logger.log(1, entry);
}

void track_values_append(){
    string fname("test_track_values_append.txt");

    ofstream out(fname);

    if (!out.good()){
        throw runtime_error("track_values_append could not open file!");
    }

    out << "# Iter,Value1,Value2\n";
    out << "1,0.2,0.3\n";
    out.close();

    TrackValueLogger logger(fname);

    track_value_row_t entry;
    entry["Value1"] = -0.2;
    entry["Value2"] = 9.9;
    logger.log(2, entry);
}