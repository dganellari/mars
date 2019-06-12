#include "mars_base.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following abstract base class for entities                                                                         ////////
////// so that std::vector<*Entity<Integer Dim, Integer ManifoldDim, Integer EntityDim>                                                 ////////
////// with components pointing to different entities                                                                                   ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace mars{
    
class BaseEntity{
public:
        BaseEntity(){};
        virtual ~BaseEntity()=default;
};













}