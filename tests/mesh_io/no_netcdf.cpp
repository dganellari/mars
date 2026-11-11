#include "mars_read_exodus_mesh.hpp"
int main()
{
    int rejected = 0;
    try
    {
        mars::readExodusMeshWithElementPartitioning<8>("unused.e", 0, 1);
    }
    catch (const std::runtime_error&)
    {
        ++rejected;
    }
    try
    {
        mars::visitExodusSideSets<4>("unused.e", [](const auto&, const auto&) {});
    }
    catch (const std::runtime_error&)
    {
        ++rejected;
    }
    try
    {
        mars::readExodusSideSetsHex8("unused.e", 0);
    }
    catch (const std::runtime_error&)
    {
        ++rejected;
    }
    return rejected == 3 ? 0 : 1;
}
