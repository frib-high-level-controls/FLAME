
#include <stdio.h>

#include "scsi/config.h"

extern int glps_debug;

int main(int argc, char *argv[])
{
    FILE *in = stdin;
    if(argc>1) {
        in = fopen(argv[1], "r");
        if(!in) {
            fprintf(stderr, "Failed to open %s\n", argv[1]);
            return 2;
        }
    }

    glps_debug = 1;

    std::auto_ptr<Config> conf;

    try{
        GLPSParser P;
        conf.reset(P.parse(in));
        fprintf(stderr, "Parsing succeeds\n");
    }catch(std::exception& e) {
        fprintf(stderr, "Parse error: %s\n", e.what());
        fclose(in);
        return 1;
    }

    std::cerr<<"Generic AST:\n";
    std::cerr << *conf;
    std::cerr<<"GLPS:\n";
    GLPSPrint(std::cout, *conf);

    fprintf(stderr, "Done\n");
    fclose(in);
    return 0;
}
