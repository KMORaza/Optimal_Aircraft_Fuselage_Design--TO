#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include <vector>
#include <fstream>

class Output {
public:
    Output(const std::string& filename);
    void writeHeader(const std::vector<std::string>& headers);
    void writeData(const std::vector<double>& data);

private:
    std::ofstream file_;
};

#endif // OUTPUT_H
