#ifndef MARS_EXPRESSIONS_HPP
#define MARS_EXPRESSIONS_HPP
#include "mars_base.hpp"


namespace mars {

///////////////////////////////////////////////////////////////////////////////////////
/////////////////////    EXPRESSIONS FOR TEMPLATE EXPRESSIONS      ////////////////////
///////////////////////////////////////////////////////////////////////////////////////
template<typename Derived>
class Expression
{
public:
        inline constexpr Derived &derived() { return static_cast<Derived &>(*this); }
        inline constexpr const Derived &derived() const { return static_cast<const Derived &>(*this); }
        inline constexpr operator Derived &() {return derived();}
        inline constexpr  operator const Derived &() const {return derived();}

         ///@return a string with the name of the class
        virtual std::string getClass() const {
            return "Expression of ";
        }
};

}
#endif //MARS_EXPRESSIONS_HPP