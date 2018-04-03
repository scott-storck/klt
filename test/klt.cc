// Karhunen-Loève Transform
// S.P.Storck 07-20-2017

#include <algorithm>
#include <complex>
#include <cstring>
#include <stdexcept>
#include <string>
#include <primitive.h> // XM
#include "klt.hh"

//==============================================================================
// MAIN
//==============================================================================
void mainroutine()
{
  // Args
  int arg = 1;
  const std::string in_fname = m_apick(arg++);
  const std::string kltb_fname = m_apick(arg++);
  const std::string kltc_fname = m_apick(arg++);
  const std::string eval_fname = m_apick(arg++);
  const int in_len = std::max(m_lpick(arg++), 2);
  const double in_olap_factor = std::max(std::min(m_dpick(arg++), 0.999999), 0.0);
  const int acm_order = std::max(std::min(m_lpick(arg++), in_len), 2);
  const int num_eig = std::max(std::min(m_lpick(arg++), acm_order), 1);
  const double out_olap_factor = std::max(std::min(m_dpick(arg++), 0.999999), 0.0);

  // Switches
#if KLT_SUPPORT_WIN
  const int window = m_get_switch_def("WIN", 0);
#endif
#if KLT_SUPPORT_EVALN
  const int eval_normalized = m_get_switch_def("EVALN", 1);
#endif

  // Compute xfer/cons lens
  const int in_clen = in_len * (1.0 - in_olap_factor);
  const int in_olen = in_len - in_clen;
  const int out_len = acm_order * (1.0 - out_olap_factor);;

  // Input file
  CPHEADER in_hcb;
  m_init(in_hcb, in_fname, "1000", "CF", 0);
  m_open(in_hcb, HCBF_INPUT);

  // KLT basis functions output file
  CPHEADER kltb_hcb;
  m_init(kltb_hcb, kltb_fname, "2000", "CF", 0);
  kltb_hcb.xstart = in_hcb.xstart;
  kltb_hcb.xdelta = in_hcb.xdelta;
  kltb_hcb.xunits = in_hcb.xunits;
  kltb_hcb.subsize = out_len;
  kltb_hcb.ystart = in_hcb.xstart;
  kltb_hcb.ydelta = (in_hcb.xdelta * in_clen) / num_eig;
  kltb_hcb.yunits = in_hcb.xunits;
  m_open(kltb_hcb, HCBF_OUTPUT + HCBF_OPTIONAL);

  // KLT coeffs output file
  CPHEADER kltc_hcb;
  m_init(kltc_hcb, kltc_fname, "2000", "CF", 0);
  kltc_hcb.xstart = 0;
  kltc_hcb.xdelta = 1;
  kltc_hcb.xunits = 0;
  kltc_hcb.subsize = num_eig;
  kltc_hcb.ystart = in_hcb.xstart;
  kltc_hcb.ydelta = in_hcb.xdelta * in_clen;
  kltc_hcb.yunits = 1;
  m_open(kltc_hcb, HCBF_OUTPUT + HCBF_OPTIONAL);

  // Eigenvalues output file
  CPHEADER eval_hcb;
  m_init(eval_hcb, eval_fname, "2000", "SF", 0);
  eval_hcb.xstart = 0;
  eval_hcb.xdelta = 1;
  eval_hcb.xunits = 0;
  eval_hcb.subsize = num_eig;
  eval_hcb.ystart = in_hcb.xstart;
  eval_hcb.ydelta = in_hcb.xdelta * in_clen;
  eval_hcb.yunits = 1;
  m_open(eval_hcb, HCBF_OUTPUT + HCBF_OPTIONAL);

  try {
    // Create KLT object
    KLT klt(in_len,
#if KLT_SUPPORT_WIN
            window,
#endif
#if KLT_SUPPORT_EVALN
            eval_normalized
#endif
            acm_order, num_eig);

    // Begin pipe section
    m_sync();

    // Main loop...
    while (m_do(in_len, in_hcb.xfer_len)) {
      // Read input file...
      in_hcb.cons_len = in_clen;
      int ngot = 0;
      m_grabx(in_hcb, klt.in_buf, ngot);
      if (ngot < 1) {
        break;
      } else if (ngot < in_len) {
        memset(&(klt.in_buf[ngot]), 0, (in_len - ngot) * in_hcb.bpa);
      }

      // KLT
      try {
        klt.transform();
      } catch (std::runtime_error err) {
        m_warning(err.what());
      }

      // Write output files...
      if (eval_hcb.open) {
        m_filad(eval_hcb, klt.eval_buf, 1);
      }
      if (kltc_hcb.open)
        m_filad(kltc_hcb, klt.kltc_buf, 1);
      if (kltb_hcb.open)
        m_filad(kltb_hcb, klt.kltb_buf, num_eig);
    } // end while (main loop)

    // Done
    m_close(in_hcb);
    if (kltb_hcb.open)
      m_close(kltb_hcb);
    if (kltc_hcb.open)
      m_close(kltc_hcb);
    if (eval_hcb.open)
      m_close(eval_hcb);
  }
  catch (std::runtime_error err) {
    m_close(in_hcb);
    if (kltb_hcb.open)
      m_close(kltb_hcb);
    if (kltc_hcb.open)
      m_close(kltc_hcb);
    if (eval_hcb.open)
      m_close(eval_hcb);
    m_error(err.what());
  }

} // end mainroutine()
