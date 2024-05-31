/*
 *  main.cpp
 *
 *  Created by Patrick on 17.04.19.
 *  Modified by Yijun Wang and Jonas Luther on 19.01.24.
 *
 */
#include "driver.h"
#include <iostream>

int main(int argc, char** argv)
{
    // read inputs
    if (argc != 2) {
        cerr << "ERROR: Missing input file: control.txt" << endl;
    } else{
        string dir = argv[1] + string("/");
        Driver d;
        d.run2D_acti(dir);
    }

    return 0;
}
