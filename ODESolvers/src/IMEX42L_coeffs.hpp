#ifndef CARPETX_ODESOLVERS_IMEX42L_COEFFS_HPP
#define CARPETX_ODESOLVERS_IMEX42L_COEFFS_HPP

namespace ODESolvers {

// Pareschi and Russo, Journal of Scientific Computing 25 (2005) 129-155,
// Table 6 (arXiv:1009.2757). IMEX-SSP3(4,3,3): the explicit part is Shu's
// SSPRK3, padded with a leading stage whose flux weight is zero, and the
// implicit part is an L-stable DIRK. Every implicit diagonal is
// alpha = 0.24169426078821, so each of the four stages solves a step of
// size alpha*dt. The first stage is not an explicit evaluation of the
// stiff source.
//
// Explicit nodes are (0, 0, 1, 1/2). Implicit nodes are
// (alpha, 0, 1, 1/2). The two abscissae differ, so the flux and the
// source of one stage are evaluated at different times. Weights are
// (0, 1/6, 1/6, 2/3). The printed decimals are the ones in Table 6.
// beta and eta are those decimals. The last implicit row's third entry
// is the expression 1/2 - beta - eta - alpha from that table.
//
// The implicit matrix has one negative entry, -alpha, on the second
// stage (i = 1). A later stage predictor can leave the range of the
// data. It stays a combination of solved stages, and it does not
// contain an explicit stiff rate.
//
// The solvers load imex42l in imex_tableaux.hxx. These functions are the
// standalone copy of the same decimals. A static assert next to imex42l
// requires the two copies to be identical.

template <typename T> static constexpr auto imex42l_alpha() -> T {
  return T(0.24169426078821);
}

template <typename T> static constexpr auto imex42l_beta() -> T {
  return T(0.06042356519705);
}

template <typename T> static constexpr auto imex42l_eta() -> T {
  return T(0.12915286960590);
}

template <typename T> static constexpr auto imex42l_b(const int i) -> T {
  if (i == 0)
    return T(0);
  if (i == 3)
    return T(2) / T(3);
  return T(1) / T(6);
}

template <typename T> static constexpr auto imex42l_c_exp(const int i) -> T {
  if (i == 2)
    return T(1);
  if (i == 3)
    return T(1) / T(2);
  return T(0);
}

template <typename T> static constexpr auto imex42l_c_imp(const int i) -> T {
  if (i == 0)
    return imex42l_alpha<T>();
  if (i == 2)
    return T(1);
  if (i == 3)
    return T(1) / T(2);
  return T(0);
}

template <typename T>
static constexpr auto imex42l_a_exp(const int i, const int j) -> T {
  if (i == 2 && j == 1)
    return T(1);
  if (i == 3 && (j == 1 || j == 2))
    return T(1) / T(4);
  return T(0);
}

template <typename T>
static constexpr auto imex42l_a_imp(const int i, const int j) -> T {
  const T alpha = imex42l_alpha<T>();
  if (i == j)
    return alpha;
  if (i == 1 && j == 0)
    return -alpha;
  if (i == 2 && j == 1)
    return T(1) - alpha;
  if (i == 3 && j == 0)
    return imex42l_beta<T>();
  if (i == 3 && j == 1)
    return imex42l_eta<T>();
  if (i == 3 && j == 2)
    return T(1) / T(2) - imex42l_beta<T>() - imex42l_eta<T>() - alpha;
  return T(0);
}

} // namespace ODESolvers

#endif // CARPETX_ODESOLVERS_IMEX42L_COEFFS_HPP
