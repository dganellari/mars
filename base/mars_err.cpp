#include "mars_err.hpp"
#ifdef _WIN32
namespace mars {

    void errorx(int select_err,
                const char* format,
                const char* file,
                int line,
                const char* function,
                const char* message) {
        // errx(select_err, format, file, line, function, message);
    }

    void warningx(const char* format, const char* message) {
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

    void warningx(const char* format, const char* message) { warnx(format, message); }
}  // namespace mars
#endif