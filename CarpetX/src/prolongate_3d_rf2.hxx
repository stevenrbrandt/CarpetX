#ifndef PROLONGATE_3D_RF2_HXX
#define PROLONGATE_3D_RF2_HXX

#include "driver.hxx"

#include <AMReX_Interpolater.H>

#include <cassert>
#include <iostream>

namespace CarpetX {

enum class centering_t : int { vertex = 0, cell = 1 };
std::ostream &operator<<(std::ostream &os, centering_t cent);
constexpr auto VC = centering_t::vertex;
constexpr auto CC = centering_t::cell;

enum class interpolation_t { poly, cons, hermite };
std::ostream &operator<<(std::ostream &os, interpolation_t cent);
constexpr auto POLY = interpolation_t::poly;
constexpr auto CONS = interpolation_t::cons;
constexpr auto HERMITE = interpolation_t::hermite;

template <centering_t CENTI, centering_t CENTJ, centering_t CENTK,
          interpolation_t INTPI, interpolation_t INTPJ, interpolation_t INTPK,
          int ORDERI, int ORDERJ, int ORDERK>
class prolongate_3d_rf2 final : public amrex::Interpolater {

  // Centering must be vertex (0) or cell (1)
  static_assert(CENTI == VC || CENTI == CC, "");
  static_assert(CENTJ == VC || CENTJ == CC, "");
  static_assert(CENTK == VC || CENTK == CC, "");

  // Conservative must be one of the possible choices
  static_assert(INTPI == POLY || INTPI == CONS || INTPI == HERMITE, "");
  static_assert(INTPJ == POLY || INTPJ == CONS || INTPJ == HERMITE, "");
  static_assert(INTPK == POLY || INTPK == CONS || INTPK == HERMITE, "");

  // Order must be nonnegative
  static_assert(ORDERI >= 0, "");
  static_assert(ORDERJ >= 0, "");
  static_assert(ORDERK >= 0, "");

  static constexpr std::array<centering_t, dim> indextype() {
    return {CENTI, CENTJ, CENTK};
  }
  static constexpr std::array<interpolation_t, dim> conservative() {
    return {INTPI, INTPJ, INTPK};
  }
  static constexpr std::array<int, dim> order() {
    return {ORDERI, ORDERJ, ORDERK};
  }

public:
  virtual ~prolongate_3d_rf2() override;

  virtual amrex::Box CoarseBox(const amrex::Box &fine, int ratio) override;
  virtual amrex::Box CoarseBox(const amrex::Box &fine,
                               const amrex::IntVect &ratio) override;

  virtual void interp(const amrex::FArrayBox &crse, int crse_comp,
                      amrex::FArrayBox &fine, int fine_comp, int ncomp,
                      const amrex::Box &fine_region,
                      const amrex::IntVect &ratio,
                      const amrex::Geometry &crse_geom,
                      const amrex::Geometry &fine_geom,
                      amrex::Vector<amrex::BCRec> const &bcr, int actual_comp,
                      int actual_state, amrex::RunOn gpu_or_cpu) override;

  virtual void interp_face(const amrex::FArrayBox &crse, int crse_comp,
                           amrex::FArrayBox &fine, int fine_comp, int ncomp,
                           const amrex::Box &fine_region,
                           const amrex::IntVect &ratio,
                           const amrex::IArrayBox &solve_mask,
                           const amrex::Geometry &crse_geom,
                           const amrex::Geometry &fine_geom,
                           amrex::Vector<amrex::BCRec> const &bcr,
                           int actual_bccomp, amrex::RunOn gpu_or_cpu) override;
};

////////////////////////////////////////////////////////////////////////////////

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c000_o0;
extern prolongate_3d_rf2<VC, VC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c001_o0;
extern prolongate_3d_rf2<VC, CC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c010_o0;
extern prolongate_3d_rf2<VC, CC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c011_o0;
extern prolongate_3d_rf2<CC, VC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c100_o0;
extern prolongate_3d_rf2<CC, VC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c101_o0;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c110_o0;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c111_o0;

extern prolongate_3d_rf2<VC, VC, VC, CONS, CONS, CONS, 1, 1, 1>
    prolongate_cons_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, CONS, CONS, CONS, 1, 1, 2>
    prolongate_cons_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, CONS, CONS, CONS, 1, 2, 1>
    prolongate_cons_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, CONS, CONS, CONS, 1, 2, 2>
    prolongate_cons_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, CONS, CONS, CONS, 2, 1, 1>
    prolongate_cons_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, CONS, CONS, CONS, 2, 1, 2>
    prolongate_cons_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, CONS, 2, 2, 1>
    prolongate_cons_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 2, 2, 2>
    prolongate_cons_3d_rf2_c111_o1;

// Prolongation operators for discrete differential forms:
// interpolating (non-conservative) for vertex centred directions,
// conservative for cell centred directions.
extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_ddf_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 0>
    prolongate_ddf_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 0, 1>
    prolongate_ddf_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 0, 0>
    prolongate_ddf_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 0, 1, 1>
    prolongate_ddf_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 0, 1, 0>
    prolongate_ddf_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 0, 0, 1>
    prolongate_ddf_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_ddf_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_ddf_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 3, 3, 2>
    prolongate_ddf_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 3, 2, 3>
    prolongate_ddf_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 3, 2, 2>
    prolongate_ddf_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 2, 3, 3>
    prolongate_ddf_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 2, 3, 2>
    prolongate_ddf_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 2, 2, 3>
    prolongate_ddf_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 2, 2, 2>
    prolongate_ddf_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5>
    prolongate_ddf_3d_rf2_c000_o5;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 5, 5, 4>
    prolongate_ddf_3d_rf2_c001_o5;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 5, 4, 5>
    prolongate_ddf_3d_rf2_c010_o5;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 5, 4, 4>
    prolongate_ddf_3d_rf2_c011_o5;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 4, 5, 5>
    prolongate_ddf_3d_rf2_c100_o5;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 4, 5, 4>
    prolongate_ddf_3d_rf2_c101_o5;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 4, 4, 5>
    prolongate_ddf_3d_rf2_c110_o5;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 4, 4, 4>
    prolongate_ddf_3d_rf2_c111_o5;

} // namespace CarpetX

#endif // #ifndef PROLONGATE_3D_RF2_HXX
