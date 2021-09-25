//
// Created by 159-egi on 9/15/17.
//

#include "io/console.h"

#include <iostream>
#include <sstream>

namespace ConsoleLog {

/* Convert double to string with specified number of places after the decimal
   and left padding. */
string prd(const double x, const int decDigits, const int width) {
    std::stringstream ss;
    ss << std::fixed << std::right;
    ss.fill(' ');              // fill space around displayed #
    ss.width(width);          // set  width around displayed #
    ss.precision(decDigits); // set # places after decimal
    ss << x;
    return ss.str();
}

/*! Center-aligns string within a field of width w. Pads with blank spaces
    to enforce alignment. */
string center(const string s, size_t w) {
    stringstream ss, spaces;

    if(w < s.size()) {
        w = s.size();
    }

    size_t padding = w - s.size();                // count excess room to pad

    for(size_t i=0; i < padding/2; ++i)
        spaces << " ";
    ss << spaces.str() << s << spaces.str();    // format with padding
    if(padding>0 && padding%2!=0)               // if odd #, add 1 space
        ss << " ";
    return ss.str();
}

/* Right-aligns string within a field of width w. Pads with blank spaces
   to enforce alignment. */
string right(const string s, size_t w) {
    stringstream ss, spaces;

    if(w < s.size()) {
        w = s.size();
    }

    size_t padding = w - s.size();   // count excess room to pad

    for(size_t i=0; i < padding; ++i)
        spaces << " ";
    ss << spaces.str() << s;       // format with padding
    return ss.str();
}

/*! Left-aligns string within a field of width w. Pads with blank spaces
    to enforce alignment. */
string left(const string s, size_t w) {
    stringstream ss, spaces;

    if(w < s.size()) {
        w = s.size();
    }

    size_t padding = w - s.size(); // count excess room to pad
    for(size_t i=0; i < padding; ++i)
        spaces << " ";
    ss << s << spaces.str();     // format with padding
    return ss.str();
}

}