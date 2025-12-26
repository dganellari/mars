class eigenFEM
{
private:
    Eigen::Matrix3d invJac_;

public:
    eigenFEM(/* args */);
    ~eigenFEM();
};

eigenFEM::eigenFEM(/* args */)
{
}

eigenFEM::~eigenFEM()
{
}

template <typename CoordViewType>
KOKKOS_INLINE_FUNCTION Eigen::Matrix3d
calcInvJac(const ConstHexDerivView& referenceGradWeights,
           const CoordViewType& coords,
           const unsigned ip) const
{
    constexpr int dim = SPATIAL_DIM;
    // TODO: HEX ONLY!!
    // TODO: HEX ONLY!!
    // TODO: HEX ONLY!!
    constexpr int nodesPerElem = 8;

    Eigen::Matrix<T, nodesPerElem, dim> myRefGrad;
    for (int n = 0; n < nodesPerElem; ++n)
    {
        for (int d = 0; d < dim; ++d)
        {
            myRefGrad(n, d) = referenceGradWeights(ip, n, d);
        }
    }

    Eigen::Matrix<T, nodesPerElem, dim> myCoords;
    for (int n = 0; n < nodesPerElem; ++n)
    {
        for (int d = 0; d < dim; ++d)
        {
            myCoords(n, d) = coords(n, d)._data.get();
        }
    }

    return (myRefGrad.transpose() * myCoords).inverse();

    // Eigen::Matrix3d myJac = myRefGrad.transpose() * myCoords;
    // Eigen::Matrix3d myInvJac = myJac.inverse();
    // return myJac.inverse();
    // Eigen::Vector3d oneWeights =
    //     myInvJac * myRefGrad(1, Eigen::placeholders::all).transpose();
}