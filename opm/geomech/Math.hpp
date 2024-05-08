#pragma once
// copy from dumux
namespace Opm
{
template <class Scalar>
Dune::FieldVector<Scalar, 3>
crossProduct(const Dune::FieldVector<Scalar, 3>& vec1, const Dune::FieldVector<Scalar, 3>& vec2)
{
    return {vec1[1] * vec2[2] - vec1[2] * vec2[1],
            vec1[2] * vec2[0] - vec1[0] * vec2[2],
            vec1[0] * vec2[1] - vec1[1] * vec2[0]};
}
} // namespace Opm