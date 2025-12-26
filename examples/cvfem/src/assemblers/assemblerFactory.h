#include "asmCvFem.h"
#include "asmCvFemInlineSF.h"
#include "asmCvFemTeamScratch.h"
#include "asmCvFemTeamBalance.h"
// #include "asmEdge.h"

#include "FEhex.h"

enum class AsmType
{
    NONE,
    CVFEM,
    CVFEM_TEAM_SCRATCH,
    CVFEM_TEAM_BALANCE,
    // EDGE,
    CVFEM_INLINE_SF
};

AsmType strToAsm(const std::string& s)
{
    std::map<std::string, AsmType> map = {
        {"cvfem", AsmType::CVFEM},
        {"cvfemTeamScratch", AsmType::CVFEM_TEAM_SCRATCH},
        {"cvfemTeamBalance", AsmType::CVFEM_TEAM_BALANCE},
        // {"edge", AsmType::EDGE},
        {"cvfemInlineSF", AsmType::CVFEM_INLINE_SF}};

    const auto res = map.find(s);

    if (res != map.end())
    {
        return res->second;
    }
    else
    {
        std::stringstream errMsg;
        errMsg << "Assembly type " << s
               << " not found! Please use one of the following options: \n";
        for (auto i = map.begin(); i != map.end(); ++i)
        {
            errMsg << i->first << "\n";
        }
        std::cerr << errMsg.str() << std::endl;
        return AsmType::NONE;
    }
}

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv = true,
          bool isShifted = true>
std::unique_ptr<
    AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>>
createAssemblerInstance(const std::string& type,
                        const stk::mesh::BulkData& bulk)
{
    switch (strToAsm(type))
    {
        case AsmType::CVFEM:
            return std::make_unique<AsmCvFem<scalar,
                                             label,
                                             BLOCKSIZE,
                                             SPATIAL_DIM,
                                             includeAdv,
                                             isShifted>>(bulk);
        case AsmType::CVFEM_TEAM_SCRATCH:
            return std::make_unique<AsmCvFemTeamScratch<scalar,
                                                        label,
                                                        BLOCKSIZE,
                                                        SPATIAL_DIM,
                                                        includeAdv,
                                                        isShifted>>(bulk);
        // case AsmType::EDGE:
        //     return std::make_unique<AsmEdge<scalar,
        //                                     label,
        //                                     BLOCKSIZE,
        //                                     SPATIAL_DIM,
        //                                     includeAdv,
        //                                     isShifted>>(bulk);
        case AsmType::CVFEM_INLINE_SF:
            return std::make_unique<AsmCvFemInlineSF<scalar,
                                                     label,
                                                     BLOCKSIZE,
                                                     SPATIAL_DIM,
                                                     FEhex<scalar>,
                                                     includeAdv,
                                                     isShifted>>(bulk);

        case AsmType::CVFEM_TEAM_BALANCE:
            return std::make_unique<AsmCvFemTeamBalance<scalar,
                                                        label,
                                                        BLOCKSIZE,
                                                        SPATIAL_DIM,
                                                        FEhex<scalar>,
                                                        includeAdv,
                                                        isShifted>>(bulk);
        default:
            std::cerr << "Assembly type: " << type << " is not supported!"
                      << std::endl;
            return nullptr;
    }
}