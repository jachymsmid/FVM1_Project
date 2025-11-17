#include <array>

template
<
  std::size_t Size,
  class RealNumber
>
using Vector = std::array< RealNumber, Size >;

template
<
  std::size_t Size,
  class RealNumber
>
struct EulersEquations {
    RealNumber gamma;

    Vector< Size, RealNumber > operator()(const Vector< Size, RealNumber > &u) const {
        Vector< Size, RealNumber > flux;

        RealNumber rho = u[0];
        RealNumber mom = u[1];
        RealNumber E   = u[2];

        RealNumber v = mom / rho;
        RealNumber p = (gamma - 1.0) * (E - 0.5 * rho * v * v);

        flux[0] = mom;
        flux[1] = mom * v + p;
        flux[2] = v * (E + p);

        return flux;
    }
};
