#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

int main(){
    char c;
    ifstream readFile;
    readFile.open("Crancknicolson_1D.cpp");
    while(readFile.peek()!=EOF){
        if (readFile.peek()=='\n') cout << "\\n";
        if(readFile.peek()==' ') cout <<" ";
        readFile >> c;
        cout << c;
        
    }
    readFile.close();
    return 0;
}
