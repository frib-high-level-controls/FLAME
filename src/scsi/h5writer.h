#ifndef H5WRITER_H
#define H5WRITER_H

#include <boost/noncopyable.hpp>

#include "base.h"

struct H5StateWriter : public boost::noncopyable
{
    struct Pvt;
    Pvt *pvt;
public:
    H5StateWriter();
    H5StateWriter(const char *);
    H5StateWriter(const std::string&);
    ~H5StateWriter();

    void open(const char *);
    void open(const std::string&);
    void close();

    void clear();

    void prepare(const StateBase *);

    void append(const StateBase *);

    void setAttr(const char *name, const char *val);
    void setAttr(const char *name, const std::string& val) { setAttr(name, val.c_str()); }

    static void dontPrint();
};

#endif // H5WRITER_H
