
#include <sstream>

#include <H5Cpp.h>

#include "scsi/h5loader.h"

// H5::Exception doesn't derive from std::exception
// so translate to some type which does.
// TODO: sub-class mixing H5::Exception and std::exception?
#define CATCH() catch(H5::Exception& he) { \
    std::ostringstream strm; \
    strm<<"H5 Error "<<he.getDetailMsg(); \
    throw std::runtime_error(strm.str()); \
    }

struct H5Loader::Pvt {
    H5::H5File file;
    H5::Group group;
};


H5Loader::H5Loader() :pvt(new Pvt) {}

H5Loader::H5Loader(const char *spec) :pvt(new Pvt)
{
    open(spec);
}

H5Loader::H5Loader(const std::string& spec) :pvt(new Pvt)
{
    open(spec);
}

H5Loader::~H5Loader()
{
    try{
        close();
    } catch(...) {
        delete pvt;
        throw;
    }
    delete pvt;
}

void H5Loader::open(const char *spec)
{
    open(std::string(spec));
}

void H5Loader::open(const std::string& spec)
{
    close();
    size_t sep = spec.find_first_of(':');
    if(sep==0)
        throw std::runtime_error("Spec. missing file name");

    std::string fname(spec.substr(0,sep));
    std::string path("/");

    if(sep!=spec.npos) {
        path = spec.substr(sep+1);
    }

    try {
        pvt->file.openFile(fname, H5F_ACC_RDONLY);
    } catch(H5::FileIException& e) {
        throw std::runtime_error("Unable to open file");
    } CATCH()

    try {
        pvt->group = pvt->file.openGroup(path);
    } catch(H5::FileIException& e) {
        throw std::runtime_error("Unable to open group");
    } CATCH()
}

void H5Loader::close()
{
    try{
        pvt->group.close();
        pvt->file.close();
    }CATCH()
}

H5Loader::matrix_t
H5Loader::load(const char * setname)
{
    H5::DataSet dset;
    try{
        dset = pvt->group.openDataSet(setname);
    }catch(H5::GroupIException& e) {
        throw std::runtime_error("H5 group does not exist");
    }CATCH()

    try{
        H5::DataSpace fspace(dset.getSpace());

        size_t N = fspace.getSimpleExtentNdims();
        if(N>2)
            throw std::runtime_error("Can't load > 2d as matrix");

        hsize_t fsize[2] = {1,1};
        fspace.getSimpleExtentDims(fsize);

        H5::DataSpace mspace(2, fsize);

        matrix_t ret(fsize[0], fsize[1]);

        fspace.selectAll();
        mspace.selectAll();

        matrix_t::array_type& storage(ret.data());
        dset.read(&storage[0], H5::PredType::NATIVE_DOUBLE, mspace, fspace);

        return ret;
    }CATCH()
}

H5Loader::matrix_t
H5Loader::load(const std::string& set)
{
    return load(set.c_str());
}

void H5Loader::dontPrint()
{
    try {
        H5::Exception::dontPrint();
    }CATCH()
}
