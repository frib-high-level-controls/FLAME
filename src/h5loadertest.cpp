
#include <typeinfo>
#include <iostream>

#include "scsi/h5loader.h"

int main(int argc, char *argv[])
{
    if(argc<2) {
        std::cerr<<"Usage: "<<argv[0]<<" <h5file>[/<h5/path]/<dsetname>\n";
        return 1;
    }

    try{
        H5Loader::dontPrint();

        std::string spec(argv[1]);
        size_t sep = spec.find_last_of('/');
        if(sep==spec.npos) {
            std::cerr<<"Missing dataset name (eg. 'file.h5/dataset')\n";
            return 1;
        }

        H5Loader loader(spec.substr(0,sep));

        H5Loader::matrix_t data(loader.load(spec.substr(sep+1)));

        std::cout<<"[\n";
        for(size_t i=0; i<data.size1(); i++) {
            std::cout<<"[ ";
            for(size_t j=0; j<data.size2(); j++) {
                std::cout<<data(i,j)<<", ";
            }
            std::cout<<"]\n";
        }
        std::cout<<"]\n";
    } catch(std::exception& e){
        std::cerr<<"Error : "<<e.what()<<"\n";
        return 1;
    }
    return 0;
}
