#ifndef GENERATION_MARS_MEMORY_HPP_
#define GENERATION_MARS_MEMORY_HPP_

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

namespace mars {
    namespace memory {

        // values are in KB.

        inline int read_line(char* line) {
            int i = strlen(line);
            const char* p = line;
            while (*p < '0' || *p > '9') p++;
            line[i - 3] = '\0';
            i = atoi(p);
            return i;
        }

        inline int get_virtual_memory() {
            FILE* file = fopen("/proc/self/status", "r");
            int result = -1;
            char line[128];

            while (fgets(line, 128, file) != NULL) {
                if (strncmp(line, "VmSize:", 7) == 0) {
                    result = read_line(line);
                    break;
                }
            }
            fclose(file);
            return result;
        }

        inline int get_physical_memory() {
            FILE* file = fopen("/proc/self/status", "r");
            int result = -1;
            char line[128];

            while (fgets(line, 128, file) != NULL) {
                if (strncmp(line, "VmRSS:", 6) == 0) {
                    result = read_line(line);
                    break;
                }
            }
            fclose(file);
            return result;
        }

    }  // namespace memory
}  // namespace mars

#endif /* GENERATION_MARS_MEMORY_HPP_ */
