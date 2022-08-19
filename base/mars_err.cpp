#include "mars_err.hpp"
#include <stdio.h>
#include <iostream>
// Handling err.h errors on windows.
#ifdef _WIN32
namespace mars {

    void errorx(int select_err,
                const char* format,
                const char* file,
                int line,
                const char* function,
                const char* message) {
        char buffer[1024];
        int nchars = snprintf(buffer, 1024, format, file, line, function, message);
        if (nchars <= 1024) {
        }
        std::cerr << "Windows does not support err.h now.";
        // errx(select_err, format, file, line, function, message);
    }

    void errorx(int select_err, const char* format) { std::cerr << "Windows does not support err.h now."; }

    void warningx(const char* format, const char* message) {
        std::cerr << "Windows does not support err.h now.";
        // warnx(format, message);
    }
}  // namespace mars
#else
#include <err.h>
namespace mars {

    void errorx(int select_err,
                const char* format,
                const char* file,
                int line,
                const char* function,
                const char* message) {
        errx(select_err, format, file, line, function, message);
    }

    void errorx(int select_err, const char* format) { errx(select_err, format); }

    void warningx(const char* format, const char* message) { warnx(format, message); }
}  // namespace mars
#endif