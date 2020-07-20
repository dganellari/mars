#ifndef MARS_INSTANCE_HPP
#define MARS_INSTANCE_HPP

namespace mars {

    class MARS {
    public:
        static void init(int argc, char* argv[]);
        static int finalize();
    };
}  // namespace mars

#endif  // MARS_INSTANCE_HPP
