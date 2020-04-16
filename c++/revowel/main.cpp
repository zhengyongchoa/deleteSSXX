#include <iostream>
#include <fstream>
#include <vector>
#include "delsx.hpp"
using namespace std ;
const char* filename = "/Users/momo/Desktop/cpp/revowel/audio/g1.pcm";
const char*  outfilename = "/Users/momo/Desktop/out/0324/0408/c_out_g1.pcm";


int writeToFile(char *buf, int size) ;

int main() {

    long size;

    std::ifstream ifstream1;
    ifstream1.open(filename, ios::in||ios::binary);

    cout<<ifstream1.is_open()<<endl;
    size = ifstream1.tellg();
    cout<<size<<endl;
    ifstream1.seekg(0, ios::beg);

    std::vector<char> data(size, 0);
    std::vector<short> outputData(size/2, 0);
    ifstream1.read(data.data(), size);
    delsxlInit();
    delsxprocess((short*)data.data(), (int)size/2 ,outputData);
    delsxTerminate();
    writeToFile((char*)outputData.data(), size);
    ifstream1.close();

    return 0;
}

int writeToFile(char *buf, int size) {
    FILE * file = fopen(outfilename, "a");
    if(!file){
        cout<<"error"<<errno<<endl;
    }
    fwrite(buf, size, 1, file);
    fclose(file);
    return 0;
}