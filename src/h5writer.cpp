

#include <sstream>

#include <boost/foreach.hpp>

#include <H5Cpp.h>

#include "flame/h5writer.h"
#include "flame/h5loader.h"

// H5::Exception doesn't derive from std::exception
// so translate to some type which does.
// TODO: sub-class mixing H5::Exception and std::exception?
#define CATCH() catch(H5::Exception& he) { \
    std::ostringstream strm; \
    strm<<"H5 Error "<<he.getDetailMsg(); \
    throw std::runtime_error(strm.str()); \
    }

namespace {
struct StateElement {
    unsigned idx;
    StateBase::ArrayInfo info;
    H5::DataSet dset;
    std::vector<hsize_t> shape;
    H5::DataSpace memspace;
};
}

struct H5StateWriter::Pvt {
    H5::H5File file;
    H5::Group group;
    std::vector<StateElement> elements;
};

H5StateWriter::H5StateWriter() :pvt(new Pvt) {}

H5StateWriter::H5StateWriter(const char *spec) :pvt(new Pvt)
{
    open(spec);
}

H5StateWriter::H5StateWriter(const std::string& spec) :pvt(new Pvt)
{
    open(spec);
}

H5StateWriter::~H5StateWriter()
{
    try{
        close();
    } catch(std::runtime_error& e) {
        std::cerr<<"H5StateWriter is ignoring exception in dtor : "<<e.what()<<"\n";
    }
    delete pvt;
}

void H5StateWriter::open(const char *spec)
{
    open(std::string(spec));
}

void H5StateWriter::open(const std::string& spec)
{
    try {
        close();
        /* The provided spec may contain both file path and group(s)
         * seperated by '/' which is ambigious as the file path
         * may contain '/' as well...
         * so do as h5ls does and strip off from the right hand side until
         * and try to open while '/' remain.
         */
        size_t sep = spec.npos;

        while(true) {
            sep = spec.find_last_of('/', sep-1);

            std::string fname(spec.substr(0, sep));

            try {
                pvt->file.openFile(fname, H5F_ACC_RDWR|H5F_ACC_CREAT);
            } catch(H5::FileIException& e) {
                if(sep==spec.npos) {
                    // no more '/' so this is failure
                    throw std::runtime_error("Unable to open file");
                }
                continue; // keep trying
            } CATCH()

            if(sep!=spec.npos) {
                std::string group(spec.substr(sep+1));
                // delete group if it already exists
                try{
                    H5G_stat_t stat;
                    pvt->file.getObjinfo(group, stat);
                    pvt->file.unlink(group);
                } catch(H5::FileIException&) {
                    // ignore non-existant
                }

                pvt->group = pvt->file.createGroup(group);
            } else {
                //TODO: cleanup root?

                pvt->group = pvt->file.openGroup("/");
            }

            return;
        }
    } CATCH()
}

void H5StateWriter::close()
{
    try{
        pvt->group.close();
        pvt->file.close();
    }CATCH()
}

void H5StateWriter::prepare(const StateBase *RS)
{
    try {
        // hack since getArray() is non-const
        // we won't actually modify the state
        StateBase *S = const_cast<StateBase*>(RS);

        assert(pvt->elements.empty());

        for(unsigned idx=0; true; idx++) {
            StateBase::ArrayInfo info;

            if(!S->getArray(idx, info))
                break;

            H5::DataType dtype;
            switch(info.type) {
            case StateBase::ArrayInfo::Double:
                dtype = H5::DataType(H5::PredType::NATIVE_DOUBLE);
                break;
            case StateBase::ArrayInfo::Sizet:
                if(sizeof(size_t)==8)
                    dtype = H5::DataType(H5::PredType::NATIVE_UINT64);
                else if(sizeof(size_t)==4)
                    dtype = H5::DataType(H5::PredType::NATIVE_UINT32);
                else
                    throw std::logic_error("unsupported size_t");
                break;
            default:
                continue; // TODO string
            }

            StateElement elem;
            elem.idx = idx;
            elem.info = info;

            // first dim is simulation "time"
            std::vector<hsize_t> dims   (info.ndim+1),
                                 maxdims(info.ndim+1, H5S_UNLIMITED);
            std::copy(info.dim,
                      info.dim+info.ndim,
                      dims.begin()+1);

            elem.shape = dims; // copy

            // size w/ first dim==0
            dims[0] = 0;

            H5::DataSpace dspace(dims.size(), &dims[0], &maxdims[0]);

            dims[0] = 10; // chunk size in "time" steps
            // other chunk sizes are multiple of initial size
            H5::DSetCreatPropList props;
            props.setChunk(dims.size(), &dims[0]);

            // memspace is simple from origin to [1,shape]
            dims[0] = 1;
            elem.memspace = H5::DataSpace(dims.size(), &dims[0]);

            elem.dset = pvt->group.createDataSet(info.name, dtype, dspace, props);

            pvt->elements.push_back(elem);
        }

        if(pvt->elements.empty()) {
            throw std::logic_error("state type has not elements to store?");
        }
    } CATCH()
}

void H5StateWriter::append(const StateBase *RS)
{
    try {
        StateBase *S = const_cast<StateBase*>(RS);

        if(pvt->elements.empty())
            prepare(RS);

        BOOST_FOREACH(StateElement& elem, pvt->elements)
        {
            StateBase::ArrayInfo info;

            if(!S->getArray(elem.idx, info))
                throw std::logic_error("can't re-fetch state parameter?");

            assert((elem.info.ndim==info.ndim) && (elem.info.type==info.type));

            size_t index = elem.shape[0]; // we will write data[index,...]
            elem.shape[0]++;

            std::copy(info.dim,
                      info.dim+info.ndim,
                      elem.shape.begin()+1);

            elem.dset.extend(&elem.shape[0]); // resize

            // filespace is hyper from [index,0...] to [index,shape]
            std::vector<hsize_t> start(elem.shape.size(), 0);
            start[0] = index;

            elem.shape[0] = 1; // reuse as count

            H5::DataSpace filespace(elem.dset.getSpace());
            filespace.selectHyperslab(H5S_SELECT_SET, &elem.shape[0], &start[0]);

            H5::DataType dtype(elem.dset.getDataType());

            elem.dset.write(info.ptr, dtype, elem.memspace, filespace);

            elem.shape[0] = index+1;
        }
    } CATCH()
}

void H5StateWriter::setAttr(const char *name, const char *val)
{
    try {
        if(H5Aexists(pvt->group.getId(), name)>0)
        //if(pvt->group.attrExists(name)) // H5Aexists was added in 1.8.0, c++ wrapper wasn't added until later...
        {
            pvt->group.removeAttr(name);
        }

        {
            H5::DataSpace scalar;
            H5::StrType dtype(H5::PredType::C_S1, std::max((size_t)16u, strlen(val)));
            H5::Attribute attr(pvt->group.createAttribute(name, dtype, scalar));
            attr.write(dtype, val);
        }
    } CATCH()
}

void H5StateWriter::dontPrint()
{
    try {
        H5::Exception::dontPrint();
    }CATCH()
}
