#ifndef MARS_ERR_HPP
#define MARS_ERR_HPP

namespace mars {
    void errorx(int select_err,
                const char* format,
                const char* file,
                int line,
                const char* function,
                const char* message);

    void warningx(const char* format, const char* message);
}  // namespace mars
#endif