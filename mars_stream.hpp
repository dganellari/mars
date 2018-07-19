#ifndef MARS_STREAM_HPP
#define MARS_STREAM_HPP

#include <ostream>
#include <istream>
#include <set>

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


	template<typename T>
	inline void read(std::set<T> &val, std::istream &is)
	{
		std::size_t n;
		read(n, is);

		for(std::size_t i = 0; i < n; ++i) {
			T v;
			read(v, is);
			val.insert(v);
		}
	}

	template<typename T>
	inline void write(const std::set<T> &val, std::ostream &os)
	{
		std::size_t n = val.size();
		write(n, os);

		for(const auto &v : val) {
			write(v, os);
		}
	}


	template<typename T>
	inline void read(std::vector<T> &val, std::istream &is)
	{
		std::size_t n;
		read(n, is);
		val.resize(n);

		for(std::size_t i = 0; i < n; ++i) {
			read(val[i], is);
		}
	}

	template<typename T>
	inline void write(const std::vector<T> &val, std::ostream &os)
	{
		std::size_t n = val.size();
		write(n, os);

		for(const auto &v : val) {
			write(v, os);
		}
	}

}

#endif //MARS_STREAM_HPP
