
#include <typeinfo>
#include <iostream>

#include "scsi/h5loader.h"

int main(int argc, char *argv[])
{
    if(argc<3) {
        std::cerr<<"Usage: "<<argv[0]<<" <h5file>[:<h5/path] <dsetname>\n";
        return 1;
    }

    try{
        H5Loader::dontPrint();

        H5Loader loader(argv[1]);

        H5Loader::matrix_t data(loader.load(argv[2]));

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
