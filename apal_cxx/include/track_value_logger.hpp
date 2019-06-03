#ifndef TRACK_VALUE_LOGGER_H
#define TRACK_VALUE_LOGGER_H

#include <string>
#include <map>
#include <fstream>
#include <vector>

typedef std::map<std::string, double> track_value_row_t;
typedef std::vector<std::string> keys_t;

class TrackValueLogger{
public:
    TrackValueLogger(const std::string &fname);
    ~TrackValueLogger();

    /** Log item to file */
    void log(unsigned int iter, const track_value_row_t &entry);

    /** Get all the keys */
    const keys_t& get_keys() const{return keys;};
private:
    std::string fname;
    std::ofstream logfile;
    keys_t keys;

    bool keys_match();
    void init_keys_from_file();
    void init_keys_from_entry(const track_value_row_t &entry);
    void write_keys();
};
#endif