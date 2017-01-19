
#include <fstream>

#include <ctime>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread/mutex.hpp>
#define BOOST_FILESYSTEM_DYN_LINK
#include <boost/filesystem.hpp>

#include "flame/util.h"

void numeric_table::readvec(std::vector<double> vec, int numcol)
{
	value_t table(int(vec.size()/numcol),numcol);
	for (unsigned i=0; i < vec.size(); i++)
	{
		int ii,jj;
		jj=i%numcol;
		ii=int((i-jj)/numcol);
		table(ii, jj) = vec[i];
    }
    this->table.swap(table);
}

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
                if(!lstrm.fail())
                    values.push_back(val);
            }

            if(!lstrm.eof() && lstrm.fail())
                throw std::runtime_error(SB()<<"Error parsing data line "<<line<<" '"<<rawline<<"'");

            if(!table.empty() && table.back().size()!=values.size())
                throw std::runtime_error(SB()<<"Line "<<line<<" w/ different # of elements");

            table.push_back(vector_t(values.size()));

            std::copy(values.begin(),
                      values.end(),
                      table.back().begin());
        }
    }

    if(!strm.eof() && strm.fail())
        throw std::runtime_error(SB()<<"Error parsing line "<<line<<" '"<<rawline<<"'");
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

struct numeric_table_cache::Pvt {
    boost::mutex lock;

    struct Value {
        std::time_t lastmod;
        typedef boost::shared_ptr<numeric_table> table_pointer;
        table_pointer table;
    };

    typedef std::map<std::string, Value> cache_t;
    cache_t cache;
};

numeric_table_cache::numeric_table_cache()
    :pvt(new Pvt)
{}

numeric_table_cache::~numeric_table_cache() {}

numeric_table_cache::table_pointer numeric_table_cache::fetch(const std::string& path)
{
    Pvt::Value::table_pointer ret;

    boost::filesystem::path P(path);
    if(!P.is_absolute())
        throw std::logic_error("numeric_table_cache want's absolute paths");

    std::time_t mtime = boost::filesystem::last_write_time(P);

    boost::mutex::scoped_lock L(pvt->lock);

    Pvt::cache_t::const_iterator it = pvt->cache.find(path);
    if(it==pvt->cache.end() || mtime>it->second.lastmod) {
        // miss

        ret.reset(new numeric_table);

        std::ifstream strm(path.c_str());

        ret->read(strm);

        Pvt::Value v;
        v.lastmod = boost::filesystem::last_write_time(P); // fetch again to avoid some, but not all, races
        v.table = ret;

        pvt->cache[path] = v;
    } else {
        ret = it->second.table;
    }

    return ret;
}

void numeric_table_cache::clear()
{
    boost::mutex::scoped_lock L(pvt->lock);
    pvt->cache.clear();
}

 //TODO: worry about global ctor order or calls before main() ???
static numeric_table_cache ntc_single;

numeric_table_cache* numeric_table_cache::get()
{
    return &ntc_single;
}
