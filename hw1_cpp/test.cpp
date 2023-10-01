// program named mainreturn.cpp
#include <iostream>
using namespace std;

// defining main with arguments
int main(int argc, char **argv)
{
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
        cout << argv[i] << "\n";

    return 0;
}
