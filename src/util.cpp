
#include <boost/numeric/ublas/vector.hpp>

#include "scsi/util.h"

void numeric_table::read(std::istream &strm)
{
    typedef boost::numeric::ublas::vector<double> vector_t;

    std::string rawline;

    std::vector<double> values;

    std::vector<vector_t> table;

    unsigned line = 0;
    while(std::getline(strm, rawline)) {
        line++;
        values.clear();

        size_t cpos = rawline.find_first_not_of(" \t");
        if(cpos==rawline.npos)
            continue; // skip blank and comment lines

        cpos = rawline.find_last_not_of("\r\n");
        if(cpos!=rawline.npos)
            rawline = rawline.substr(0, cpos+1);

        if(rawline[0]=='%') {
            // headers
            size_t P = rawline.find_first_not_of(" \t", 1);

            colnames_t cols;

            unsigned col=0;
            while(P!=rawline.npos) {
                size_t E = rawline.find_first_of(" \t", P);

                const std::string& ent = rawline.substr(P, E==rawline.npos ? E : E-P);
                cols[ent] = col++;

                P = rawline.find_first_not_of(" \t", E);
            }

            colnames.swap(cols);

        } else {
            // data

            std::istringstream lstrm(rawline);

            while(lstrm.good()) {
                double val;
                lstrm>>val;
                values.push_back(val);
            }

            if(!lstrm.eof() && lstrm.fail())
                throw std::runtime_error(std::string("Error parsing data line '")+rawline+"'");

            if(!table.empty() && table.back().size()!=values.size())
                throw std::runtime_error("Line w/ different # of elements");

            table.push_back(vector_t(values.size()));

            std::copy(values.begin(),
                      values.end(),
                      table.back().begin());
        }
    }

    if(!strm.eof() && strm.fail())
        throw std::runtime_error(std::string("Error parsing line '")+rawline+"'");
    else if(table.empty()) {
        this->table.clear();
    } else {
        value_t result(table.size(), table.front().size());

        for(size_t r=0; r<result.size1(); r++) {
            const vector_t& R=table[r];
            for(size_t c=0; c<result.size2(); c++) {
                std::copy(R.begin(),
                          R.end(),
                          result.find2(2, r,0));
            }
        }

        this->table.swap(result);
    }
}

