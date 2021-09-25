//
// Created by 159-egi on 9/15/17.
//

#ifndef AMR_CONSOLE_H
#define AMR_CONSOLE_H

#include <allstd.h>

namespace ConsoleLog {
    string prd(const double x, const int decDigits, const int width);
    string center(const string s, const size_t w);
    string right(const string s, const size_t w);
    string left(const string s, const size_t w);
}

#endif //AMR_CONSOLE_H
