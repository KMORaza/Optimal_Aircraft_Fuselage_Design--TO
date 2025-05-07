#include "output.h"
#include <sstream>

Output::Output(const std::string& filename) {
    file_.open(filename);
}

void Output::writeHeader(const std::vector<std::string>& headers) {
    if (!file_.is_open()) return;
    for (size_t i = 0; i < headers.size(); ++i) {
        file_ << headers[i];
        if (i < headers.size() - 1) file_ << ",";
    }
    file_ << "\n";
}

void Output::writeData(const std::vector<double>& data) {
    if (!file_.is_open()) return;
    for (size_t i = 0; i < data.size(); ++i) {
        file_ << data[i];
        if (i < data.size() - 1) file_ << ",";
    }
    file_ << "\n";
}
