#include "track_value_logger.hpp"
#include <iostream>
using namespace std;

TrackValueLogger::TrackValueLogger(const string &fname): fname(fname){
    logfile.open(fname, ios::out | ios::app);

    if (!logfile.good()){
        throw runtime_error("TrackValueLogger could not open logfile!");
    }

    init_keys_from_file();
}

TrackValueLogger::~TrackValueLogger(){
    logfile.close();
}

void TrackValueLogger::init_keys_from_file(){
    if (keys.size() != 0){
        throw runtime_error("init_keys_from_file should never be called if keys are already initialised!");
    }

    ifstream read_stream(fname);

    if (!read_stream.good()){
        throw runtime_error("Could not open logfile for reading!");
    }

    string line;
    getline(read_stream, line);
    read_stream.close();

    if (line.substr(0, 2) != "# "){
        return;
    }

    // Erase first # and the following whitespace
    line.erase(0, 2);

    unsigned int pos = line.find(",");
    unsigned int loop_indx = 0;
    unsigned int max_loop_indx = 100;
    while (loop_indx++ < max_loop_indx){
        int pos = line.find(",");
        if (pos != std::string::npos){
            keys.push_back(line.substr(0, pos));
            line.erase(0, pos+1);
        }
        else{
            keys.push_back(line);
            break;
        }
    }

    if (loop_indx >= max_loop_indx){
        throw runtime_error("Entered potentially infinite loop while reading keys from file!");
    }
}

void TrackValueLogger::init_keys_from_entry(const track_value_row_t &entry){
    if (keys.size() != 0){
        throw runtime_error("init_keys_from_entry should never be called if keys are already initialised!");
    }

    for (auto iter=entry.begin(); iter != entry.end(); ++iter){
        keys.push_back(iter->first);
    }
}