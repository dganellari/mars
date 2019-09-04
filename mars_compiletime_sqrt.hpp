#ifndef MARS_COMPILETIME_Sqrt__HPP
#define MARS_COMPILETIME_Sqrt__HPP



namespace mars{


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Compile-time Square-Root based on ratios Sqrt_{X}=Sqrt_{X/1}
///// For implementation details: 
///// https://stackoverflow.com/questions/36321295/rational-approximation-of-square-root-of-stdratio-at-compile-time
/////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;

using ZeroR = ratio<0>;
using OneR = ratio<1>;
template <typename R> using Square = ratio_multiply<R, R>;

// Find the largest integer N such that Predicate<N>::value is true.
template <template <intmax_t N> class Predicate, typename Enabled = void>
struct BinarySearch {
  template <intmax_t N>
  struct SafeDouble_ {
    const intmax_t static value = 2 * N;
    static_assert(value > 0, "Overflows when computing 2 * N");
  };

  template <intmax_t Lower, intmax_t Upper, typename Enabled1 = void>
  struct DoubleSidedSearch_ : DoubleSidedSearch_<Lower, Lower+(Upper-Lower)/2> {};

  template <intmax_t Lower, intmax_t Upper>
  struct DoubleSidedSearch_<Lower, Upper, typename enable_if<Upper-Lower==1>::type> : integral_constant<intmax_t, Lower> {};

  template <intmax_t Lower, intmax_t Upper>
  struct DoubleSidedSearch_<Lower, Upper, typename enable_if<(Upper-Lower>1 && Predicate<Lower+(Upper-Lower)/2>::value)>::type>
      : DoubleSidedSearch_<Lower+(Upper-Lower)/2, Upper> {};

  template <intmax_t Lower, typename Enabled1 = void>
  struct SingleSidedSearch_ : DoubleSidedSearch_<Lower, SafeDouble_<Lower>::value> {};

  template <intmax_t Lower>
  struct SingleSidedSearch_<Lower, typename enable_if<Predicate<SafeDouble_<Lower>::value>::value>::type>
      : SingleSidedSearch_<SafeDouble_<Lower>::value> {};

  const static intmax_t value = SingleSidedSearch_<1>::value;
};

template <template <intmax_t N> class Predicate>
struct BinarySearch<Predicate, typename enable_if<!Predicate<1>::value>::type> : integral_constant<intmax_t, 0> {};

// Find largest integer N such that N<=Sqrt_(R)
template <typename R>
struct IntegerR {
  template <intmax_t N> using Predicate_ = ratio_less_equal<ratio<N>, ratio_divide<R, ratio<N>>>;
  const static intmax_t value = BinarySearch<Predicate_>::value;
};

template <typename R>
struct IsPerfectSquare {
  const static intmax_t DenSqrt__ = IntegerR<ratio<R::den>>::value;
  const static intmax_t NumSqrt__ = IntegerR<ratio<R::num>>::value;
  const static bool value = DenSqrt__ * DenSqrt__ == R::den && NumSqrt__ * NumSqrt__ == R::num;
  using Sqrt_ = ratio<NumSqrt__, DenSqrt__>;
};

// Represents Sqrt_(P)-Q.
template <typename Tp, typename Tq>
struct Remainder {
  using P = Tp;
  using Q = Tq;
};

// Represents 1/R = I + Rem where R is a Remainder.
template <typename R>
struct Reciprocal {
  using P_ = typename R::P;
  using Q_ = typename R::Q;
  using Den_ = ratio_subtract<P_, Square<Q_>>;
  using A_ = ratio_divide<Q_, Den_>;
  using B_ = ratio_divide<P_, Square<Den_>>;
  const static intmax_t I_ = (A_::num + IntegerR<ratio_multiply<B_, Square<ratio<A_::den>>>>::value) / A_::den;
  using I = ratio<I_>;
  using Rem = Remainder<B_, ratio_subtract<I, A_>>;
};

// Expands Sqrt_(R) to continued fraction:
// f(x)=C1+1/(C2+1/(C3+1/(...+1/(Cn+x)))) = (U*x+V)/(W*x+1) and Sqrt_(R)=f(Rem).
// The error |f(Rem)-V| = |(U-W*V)x/(W*x+1)| <= |U-W*V|*Rem <= |U-W*V|/I' where
// I' is the integer part of reciprocal of Rem.
template <typename R, intmax_t N>
struct ContinuedFraction {
  template <typename T>
  using Abs_ = typename conditional<ratio_less<T, ZeroR>::value, ratio_subtract<ZeroR, T>, T>::type;

  using Last_ = ContinuedFraction<R, N-1>;
  using Reciprocal_ = Reciprocal<typename Last_::Rem>;
  using Rem = typename Reciprocal_::Rem;
  using I_ = typename Reciprocal_::I;
  using Den_ = ratio_add<typename Last_::W, I_>;
  using U = ratio_divide<typename Last_::V, Den_>;
  using V = ratio_divide<ratio_add<typename Last_::U, ratio_multiply<typename Last_::V, I_>>, Den_>;
  using W = ratio_divide<OneR, Den_>;
  using Error = Abs_<ratio_divide<ratio_subtract<U, ratio_multiply<V, W>>, typename Reciprocal<Rem>::I>>;
};

template <typename R>
struct ContinuedFraction<R, 1> {
  using U = OneR;
  using V = ratio<IntegerR<R>::value>;
  using W = ZeroR;
  using Rem = Remainder<R, V>;
  using Error = ratio_divide<OneR, typename Reciprocal<Rem>::I>;
};

template <typename R, typename Eps, intmax_t N=1, typename Enabled = void>
struct Sqrt__ : Sqrt__<R, Eps, N+1> {};

template <typename R, typename Eps, intmax_t N>
struct Sqrt__<R, Eps, N, typename enable_if<ratio_less_equal<typename ContinuedFraction<R, N>::Error, Eps>::value>::type> {
  using type = typename ContinuedFraction<R, N>::V;
};

template <typename R, typename Eps, typename Enabled = void>
struct Sqrt_ {
  static_assert(ratio_greater_equal<R, ZeroR>::value, "R can't be negative");
};

template <typename R, typename Eps>
struct Sqrt_<R, Eps, typename enable_if<ratio_greater_equal<R, ZeroR>::value && IsPerfectSquare<R>::value>::type> {
  using type = typename IsPerfectSquare<R>::Sqrt_;
};

template <typename R, typename Eps>
struct Sqrt_<R, Eps, typename enable_if<(ratio_greater_equal<R, ZeroR>::value && !IsPerfectSquare<R>::value)>::type> : Sqrt__<R, Eps> {};

// Test finding Sqrt_(N/D) with error 1/Eps
template <intmax_t N, intmax_t D, intmax_t Eps>
void test() {
  using T = typename Sqrt_<ratio<N, D>, ratio<1, Eps>>::type;
  constexpr double M=static_cast<Real>(T::num) / static_cast<Real>(T::den);
  cout << "Sqrt_(" << N << "/" << D << ") ~ " << T::num << "/" << T::den << "=="<<M << ", "
       << "error=" << abs(sqrt(N/(double)D) - T::num/(double)T::den) << ", "
       << "eps=" << 1/(double)Eps << endl;
}

// template <intmax_t N>
// static constexpr auto Sqrt=static_cast<Real>( Sqrt_<ratio<N, 1>, ratio<1, 10000000000000000>>::type::num) 
//                          / static_cast<Real>( Sqrt_<ratio<N, 1>, ratio<1, 10000000000000000>>::type::den);


template <intmax_t N,intmax_t D>
static constexpr auto Sqrt=static_cast<Real>( Sqrt_<ratio<N, D>, ratio<1, 10000000000000000>>::type::num) 
                         / static_cast<Real>( Sqrt_<ratio<N, D>, ratio<1, 10000000000000000>>::type::den);


template<typename T>
constexpr Real Sqrt_Real(const T& D)
{
  return Sqrt<static_cast<intmax_t>(D*10000000000000000),10000000000000000,10000000000000000>;
}

template <intmax_t N >
void test2() {
  // cout << "Sqrt_(" << N/D<< ") ~ " << Sqrt<N,D> << ", "<<endl;
  // cout     << "error sqrt=" << abs(sqrt(N) - Sqrt<N>) << ", "<<endl;
  // cout     << "error ^2=" << N - Sqrt<N>*Sqrt<N> << ", "<<endl;
}

template <intmax_t N, intmax_t D>
void testAll() {
  test<N, D, 10000>();
  test<N, D, 10000000000>();
  test<N, D, 10000000000000>();
  test<N, D, 10000000000000000>();
}

template <intmax_t N>
void testAll2() {
  test2<N>();
}






}


#endif