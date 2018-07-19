#ifndef MARS_STREAM_HPP
#define MARS_STREAM_HPP

#include <ostream>
#include <istream>

namespace mars {

	template<typename T>
	inline void read(T &val, std::istream &is)
	{
		static_assert(std::is_trivial<T>::value, "implement read for your type");
		is.read((char *)&val, sizeof(T));
	}

	template<typename T>
	inline void write(const T &val, std::ostream &os)
	{
		static_assert(std::is_trivial<T>::value, "implement write for your type");
		os.write((const char *)&val, sizeof(T));
	}


	template<typename T>
	inline void read(T *val, const std::size_t n, std::istream &is)
	{
		static_assert(std::is_trivial<T>::value, "implement read for your type");
		is.read((char *)val, n*sizeof(T));
	}

	template<typename T>
	inline void write(const T *val, const std::size_t n, std::ostream &os)
	{
		static_assert(std::is_trivial<T>::value, "implement write for your type");
		os.write((const char *)val, n*sizeof(T));
	}
}

#endif //MARS_STREAM_HPP
