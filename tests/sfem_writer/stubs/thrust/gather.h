#pragma once
namespace thrust
{
inline constexpr int device = 0;
template<class T>
T* device_pointer_cast(T* p)
{
    return p;
}
template<class Map, class Input, class Output>
void gather(int, Map first, Map last, Input input, Output output)
{
    while (first != last)
        *output++ = input[*first++];
}
} // namespace thrust
