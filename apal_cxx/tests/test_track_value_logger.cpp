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

    out << "# Item1,Item2,Item3,Item4,Item5\n";
    out << "0.5,0.2,0.6,0.8,1.0\n";
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