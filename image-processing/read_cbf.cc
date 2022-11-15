#include <iostream>
#include <fstream>
using namespace std;

static int length;
static int length1 ;
static char *buffer;
static int allocated = 0;

size_t allocate_read_cbf(int w, int h)
{
   length =w*h;
   length1 = length-1;
   buffer = new char [length*4+5000];
   allocated = 1;
}
size_t filesize(char *filename)
{
    ifstream fi(filename, ios::ate);
    size_t size = fi.tellg();
    return size;
}

int read_file(char filename[], int *datai)
{
    if(! allocated) cerr <<" you have to allocate read_cbf"<<endl;
    ifstream fi(filename);
    
    int size = filesize(filename);
    fi.read(buffer, size);
    int pos;   //?
    for(int i = 0; i < size-4; i++){
        if(buffer[i] =='\x0c'){
            if(buffer[i+1]=='\x1a'){
                if(buffer[i+2] =='\x04'){
                    if(buffer[i+3] =='\xd5'){
                        pos = i+4;
                        i = size-4;
                    }
                }
            }
        }
    }
    int value = 0;
    int count = 0;
    char c1 = 0x80;
    short s1 = 0x8000;
    while(pos < size){
        if(buffer[pos]!=c1){
            value +=  buffer[pos];
            datai[count++] = value;
            if(count == length1) return 0;
            pos++;
        }else{
            pos++;
            short *temp = reinterpret_cast<short *>(buffer + pos);
            if((*temp) != s1){
                value += *temp;
                datai[count++] = value;
                if(count == length1) return 0;
                pos += 2;
            }else{
                pos+=2;
                int *ltemp = reinterpret_cast<int *>(buffer + pos);
                value += *ltemp;
                datai[count++] = value;
                if(count == length1) return 0;
                pos += 4;
            }
        }
    }
    cout << count << endl;
}
#ifdef DEBUG_MAIN
int main(int argc, char *argv[])
{
    allocate_read_cbf(619, 487);
    buffer = new char [length*4];
    int *datai = new int[length];
    read_file(argv[1], datai);  //input the cbf fila
    ofstream fo("temp.dat");
    fo.write(reinterpret_cast<char *> (datai), sizeof(int)*length);
    fo.close();
}
#endif
