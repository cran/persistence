//
//  eccezioni.h
//  communityDetection
//
//  Created by Alessandro Avellone on 19/02/24.
//

#ifndef eccezioni_h
#define eccezioni_h


#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#define throw_line(arg) throw eccezioni(arg, __FILE__, __LINE__);

class eccezioni : public std::runtime_error {
    std::string msg;
public:
    eccezioni(const std::string &arg, const char *file, int line) :
    std::runtime_error(arg) {
        std::ostringstream o;
        o << file << ":" << line << ": " << arg;
        msg = o.str();
    }
    ~eccezioni() throw() {}
    const char *what() const throw() {
        return msg.c_str();
    }
};


#endif /* eccezioni_h */
