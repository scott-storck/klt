// Karhunen-Loève Transform Library
// S.P.Storck 07-20-2017

#include "klt.hh"

#if KLT_SUPPORT_WIN
#include <cmath>
#endif
#include <complex>
#include <cstdlib>
#include <cstring>
#if (KLT_DEBUG & KLT_DEBUG_VERBOSE) || (KLT_DEBUG & KLT_DEBUG_FINE)
#include <iostream>
#endif
#include <sstream>
#include <stdexcept>
#include <mkl_lapacke.h> // MKL
#if KLT_SUPPORT_WIN && (KLT_DEBUG & KLT_DEBUG_WIN)
#include <primitive.h>
#endif

static const size_t ALIGN = 128;

//---------------------------------------------------------------------------
// Constructor
//---------------------------------------------------------------------------
KLT::KLT(int in_len,
#if KLT_SUPPORT_WIN
         int window,
#endif
#if KLT_SUPPORT_EVALN
         int eval_normalized,
#endif
         int acm_order,
         int num_eig) :
  in_len(in_len),
#if KLT_SUPPORT_WIN
  window(window),
#endif
#if KLT_SUPPORT_EVALN
  eval_normalized(eval_normalized),
#endif
  acm_order(acm_order),
  num_eig(num_eig),
  in_buf(NULL),
#if KLT_SUPPORT_WIN
  win_buf(NULL),
#endif
  ac_buf(NULL),
  eval_buf(NULL),
  kltc_buf(NULL),
  kltb_buf(NULL),
  d_buf(NULL),
  e_buf(NULL),
  tau_buf(NULL),
  ib_buf(NULL),
  is_buf(NULL),
  if_buf(NULL)
{
#if KLT_DEBUG & KLT_DEBUG_FINE
  std::cout<<"in_len="<<in_len<<
    "  acm_order="<<acm_order<<
    "  num_eig="<<num_eig<<std::endl;
#endif
  // Allocate buffers...
  // Input
  if (posix_memalign(reinterpret_cast<void**>(&in_buf),
                     ALIGN, in_len*sizeof(std::complex<float>))) {
    std::ostringstream oss;
    oss << "Failed to allocate in_buf (size " << in_len << ")";
    throw std::runtime_error(oss.str());
  }
#if KLT_SUPPORT_WIN
  // Window
  if (posix_memalign(reinterpret_cast<void**>(&win_buf),
                     ALIGN, in_len*sizeof(float))) {
    std::ostringstream oss;
    oss << "Failed to allocate win_buf (size " << in_len << ")";
    throw std::runtime_error(oss.str());
  }
#endif
  // Auto-correlation toeplitz matrix (lower triangular packed, col major order)
  const size_t ac_size = ((acm_order + 1) * acm_order) / 2;
  if (posix_memalign(reinterpret_cast<void**>(&ac_buf),
                     ALIGN, ac_size*sizeof(std::complex<float>))) {
    std::ostringstream oss;
    oss << "Failed to allocate ac_buf (size " << ac_size << ")";
    throw std::runtime_error(oss.str());
  }
  // Eigenvalues
  // This vector must be allocated to size acm_order rather than
  // just num_eig because LAPACKE_cstein() requires an oversized
  // buffer.
  if (posix_memalign(reinterpret_cast<void**>(&eval_buf),
                     ALIGN, acm_order*sizeof(float))) {
    std::ostringstream oss;
    oss << "Failed to allocate eval_buf (size " << acm_order << ")";
    throw std::runtime_error(oss.str());
  }
#if KLT_DEBUG & KLT_DEBUG_FINE
  std::cout<<"eval_buf="<<eval_buf<<"["<<acm_order<<"]"<<std::endl;
#endif
  // KLT coeffs
  if (posix_memalign(reinterpret_cast<void**>(&kltc_buf),
                     ALIGN, num_eig*sizeof(std::complex<float>))) {
    std::ostringstream oss;
    oss << "Failed to allocate kltc_buf (size " << num_eig << ")";
    throw std::runtime_error(oss.str());
  }
  // KLT basis funcs (eigenvectors)
  const size_t kltb_size = acm_order * num_eig;
  if (posix_memalign(reinterpret_cast<void**>(&kltb_buf),
                     ALIGN, kltb_size*sizeof(std::complex<float>))) {
    std::ostringstream oss;
    oss << "Failed to allocate kltb_buf (size " << kltb_size << ")";
    throw std::runtime_error(oss.str());
  }
  // Eigendecomp temp buffers
  if (posix_memalign(reinterpret_cast<void**>(&d_buf),
                     ALIGN, acm_order*sizeof(float))) {
    std::ostringstream oss;
    oss << "Failed to allocate d_buf (size " << acm_order << ")";
    throw std::runtime_error(oss.str());
  }
  const size_t e_size = acm_order - 1;
  if (posix_memalign(reinterpret_cast<void**>(&e_buf),
                     ALIGN, e_size*sizeof(float))) {
    std::ostringstream oss;
    oss << "Failed to allocate e_buf (size " << e_size << ")";
    throw std::runtime_error(oss.str());
  }
  const size_t tau_size = acm_order - 1;
  if (posix_memalign(reinterpret_cast<void**>(&tau_buf),
                     ALIGN, tau_size*sizeof(std::complex<float>))) {
    std::ostringstream oss;
    oss << "Failed to allocate tau_buf (size " << tau_size << ")";
    throw std::runtime_error(oss.str());
  }
  if (posix_memalign(reinterpret_cast<void**>(&ib_buf),
                     ALIGN, acm_order*sizeof(int))) {
    std::ostringstream oss;
    oss << "Failed to allocate ib_buf (size " << acm_order << ")";
    throw std::runtime_error(oss.str());
  }
  if (posix_memalign(reinterpret_cast<void**>(&is_buf),
                     ALIGN, acm_order*sizeof(int))) {
    std::ostringstream oss;
    oss << "Failed to allocate is_buf (size " << acm_order << ")";
    throw std::runtime_error(oss.str());
  }
  if (posix_memalign(reinterpret_cast<void**>(&if_buf),
                     ALIGN, num_eig*sizeof(int))) {
    std::ostringstream oss;
    oss << "Failed to allocate if_buf (size " << num_eig << ")";
    throw std::runtime_error(oss.str());
  }
#if KLT_SUPPORT_WIN
  // Create window
  if (window) {
    init_window();
#if KLT_DEBUG & KLT_DEBUG_WIN
    CPHEADER win_hcb;
    m_init(win_hcb, "klt_win", "1000", "SF", 0);
    win_hcb.xstart = 0;
    win_hcb.xdelta = 1;
    win_hcb.xunits = 1;
    m_open(win_hcb, HCBF_OUTPUT);
    m_filad(win_hcb, win_buf, in_len);
    m_close(win_hcb);
#endif
  }
#endif
}


//---------------------------------------------------------------------------
// Destructor
//---------------------------------------------------------------------------
KLT::~KLT()
{
  if (in_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free in_buf"<<std::endl;
#endif
    free(in_buf);
  }
#if KLT_SUPPORT_WIN
  if (win_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free win_buf"<<std::endl;
#endif
    free(win_buf);
  }
#endif
  if (ac_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free ac_buf"<<std::endl;
#endif
    free(ac_buf);
  }
  if (eval_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free eval_buf"<<std::endl;
#endif
    free(eval_buf);
  }
  if (kltc_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free kltc_buf"<<std::endl;
#endif
    free(kltc_buf);
  }
  if (kltb_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free kltb_buf"<<std::endl;
#endif
    free(kltb_buf);
  }
  if (d_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free d_buf"<<std::endl;
#endif
    free(d_buf);
  }
  if (e_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free e_buf"<<std::endl;
#endif
    free(e_buf);
  }
  if (tau_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free tau_buf"<<std::endl;
#endif
    free(tau_buf);
  }
  if (ib_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free ib_buf"<<std::endl;
#endif
    free(ib_buf);
  }
  if (is_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free is_buf"<<std::endl;
#endif
    free(is_buf);
  }
  if (if_buf != NULL) {
#if KLT_DEBUG & KLT_DEBUG_FINE
    std::cout<<"Free if_buf"<<std::endl;
#endif
    free(if_buf);
  }
}


//---------------------------------------------------------------------------
// Transform in_buf
//   If an error occurrs, throws std::runtime_error and sets all of eval_buf,
//   kltc_buf, and kltb_buf to 0.0f.
//---------------------------------------------------------------------------
void KLT::transform()
{
#if KLT_SUPPORT_WIN
  // Apply window
  if (window) {
    apply_window();
  }
#endif
  // Auto-corr matrix
  acorr_matrix();
  // Eigendecomp
  eigendecomp();
#if KLT_SUPPORT_EVALN
  // Normalize eigenvalues?
  if (eval_normalized) {
    float max_eval = 0.0f;
    for (int eidx=0; eidx<num_eig; eidx++) {
      max_eval = std::max(eval_buf[eidx], max_eval);
    }
    const float inv_max_eval = 1.0f / max_eval;
    for (int eidx=0; eidx<num_eig; eidx++) {
      eval_buf[eidx] *= inv_max_eval;
    }
  }
#endif
  // Compute KLT coeffs
  for (int cidx=0; cidx<num_eig; cidx++) {
    kltc_buf[cidx] = std::complex<float>(0.0f, 0.0f);
    for (int tidx=0; tidx<acm_order; tidx++) {
      kltc_buf[cidx] += in_buf[tidx] * std::conj(kltb_buf[cidx*acm_order+tidx]);
    }
  }
  // Apply coeffs to KLT basis funcions
  for (int cidx=0; cidx<num_eig; cidx++) {
    for (int tidx=0; tidx<acm_order; tidx++) {
      kltb_buf[cidx*acm_order+tidx] *= kltc_buf[cidx];
    }
  }
}


//-----------------------------------------------------------------------------
// Auto-correlation Matrix
//-----------------------------------------------------------------------------
void KLT::acorr_matrix()
{
  for (int aidx=0; aidx < acm_order; aidx++) {
    ac_buf[aidx] = std::complex<float>(0.0f, 0.0f);
    for (int w2idx=0; w2idx < in_len - aidx; w2idx++) {
      const int w1idx = w2idx + aidx;
      ac_buf[aidx] += in_buf[w1idx] * std::conj(in_buf[w2idx]);
    }
  }
#if KLT_DEBUG & KLT_DEBUG_VERBOSE
  std::cout<<"ac: ";
  for (int aidx=0; aidx < acm_order; aidx++) {
    std::cout<<ac_buf[aidx]<<"  ";
  }
  std::cout<<std::endl;
#endif
  // Toeplitz (lower triangular packed)
  int aoidx = acm_order;
  for (int col_len = acm_order-1; col_len > 0; col_len--) {
    for (int aiidx = 0; aiidx < col_len; aiidx++) {
      ac_buf[aoidx++] = ac_buf[aiidx];
    }
  }
}


//-----------------------------------------------------------------------------
// Compute eigenvalues (eval_buf) & eigenvectors (kltb_buf) for ac_buf
//   WARNING: ac_buf is modified in the process.
//   If an error occurrs, throws std::runtime_error and sets all of eval_buf,
//   kltc_buf, and kltb_buf to 0.0f.
//-----------------------------------------------------------------------------
void KLT::eigendecomp()
{
  static const int major_order = LAPACK_COL_MAJOR;
  static const char uplo = 'L';
  int info;
  info = LAPACKE_chptrd(major_order, uplo, acm_order,
                        reinterpret_cast<lapack_complex_float*>(ac_buf),
                        d_buf, e_buf,
                        reinterpret_cast<lapack_complex_float*>(tau_buf));
  if (info) {
    memset(eval_buf, 0, num_eig*sizeof(float));
    memset(kltb_buf, 0, num_eig*acm_order*sizeof(std::complex<float>));
    std::ostringstream oss;
    oss << "LAPACKE_chptrd() failed, info=" << info;
    throw std::runtime_error(oss.str());
  } else {
    static const char range = 'I';
    static const char order = 'B';
    static const float vl = 0.0f;
    static const float vu = 0.0f;
    static const float abstol = 0.0f;
    const int il = acm_order - num_eig + 1;
    const int iu = acm_order;
    int num_eig_found;
    int nsplit;
    info = LAPACKE_sstebz(range, order, acm_order, vl, vu, il, iu, abstol,
                          d_buf, e_buf, &num_eig_found, &nsplit, eval_buf,
                          ib_buf, is_buf);
#if KLT_DEBUG & KLT_DEBUG_FINE
    if (nsplit != 1) {
      std::cout<<" WARNING: nsplit="<<nsplit<<std::endl;
    }
    if (num_eig_found != num_eig) {
      std::cout<<" WARNING: num_eig_found="<<num_eig_found<<std::endl;
    }
#endif
    if (info) {
      memset(eval_buf, 0, num_eig*sizeof(float));
      memset(kltb_buf, 0, num_eig*acm_order*sizeof(std::complex<float>));
      std::ostringstream oss;
      oss << "LAPACKE_sstebz() failed, info=" << info;
      throw std::runtime_error(oss.str());
    } else {
      info = LAPACKE_cstein(major_order, acm_order, d_buf, e_buf, num_eig_found, eval_buf,
                            ib_buf, is_buf,
                            reinterpret_cast<lapack_complex_float*>(kltb_buf),
                            acm_order, if_buf);
      if (info) {
        memset(eval_buf, 0, num_eig*sizeof(float));
        memset(kltb_buf, 0, num_eig*acm_order*sizeof(std::complex<float>));
        std::ostringstream oss;
        oss << "LAPACKE_cstein() failed, info=" << info;
        throw std::runtime_error(oss.str());
      } else {
        static const char side = 'L';
        static const char trans = 'N';
        info = LAPACKE_cupmtr(major_order, side, uplo, trans, acm_order, num_eig_found,
                              reinterpret_cast<lapack_complex_float*>(ac_buf),
                              reinterpret_cast<lapack_complex_float*>(tau_buf),
                              reinterpret_cast<lapack_complex_float*>(kltb_buf),
                              acm_order);
        if (info) {
          memset(eval_buf, 0, num_eig*sizeof(float));
          memset(kltb_buf, 0, num_eig*acm_order*sizeof(std::complex<float>));
          std::ostringstream oss;
          oss << "LAPACKE_cupmtr() failed, info=" << info;
          throw std::runtime_error(oss.str());
        }
      }
    }
  }
#if KLT_DEBUG & KLT_DEBUG_VERBOSE
  std::cout<<"eval: ";
  for (int eidx=0; eidx<num_eig; eidx++) {
    std::cout<<eval_buf[eidx]<<",  ";
  }
  std::cout<<std::endl;
  for (int vidx=0; vidx<num_eig; vidx++) {
    std::cout<<"evec["<<vidx<<"]: ";
    for (int tidx=0; tidx<acm_order; tidx++) {
      std::cout<<kltb_buf[vidx*acm_order+tidx]<<",  ";
    }
    std::cout<<std::endl;
  }
#endif
}


#if KLT_SUPPORT_WIN
//-----------------------------------------------------------------------------
// Create window
//-----------------------------------------------------------------------------
void KLT::init_window()
{
  static const bool PERIODIC = false;
  static const double TWO_PI = 2.0 * M_PI;
  const int num_out = in_len;
  const int win_len = in_len + static_cast<int>(PERIODIC);
  const double z1 = TWO_PI / win_len;
  for (int oidx=0; oidx<num_out; oidx++) {
    const double z = z1 * oidx;
    // HFT90D Flat-top Window
    win_buf[oidx] =
      1.0 -
      1.942604 * cos(z) +
      1.340318 * cos(2.0 * z) -
      0.440811 * cos(3.0 * z) +
      0.043097 * cos(4.0 * z);
  }
}


//-----------------------------------------------------------------------------
// Apply window to in_buf in-place
//-----------------------------------------------------------------------------
void KLT::apply_window()
{
#if KLT_DEBUG & KLT_DEBUG_VERBOSE
  std::cout<<"in: ";
  for (int widx=0; widx < in_len; widx++) {
    std::cout<<in_buf[widx]<<"  ";
  }
  std::cout<<std::endl;
#endif
  for (int widx=0; widx < in_len; widx++) {
    in_buf[widx] *= win_buf[widx];
  }
#if KLT_DEBUG & KLT_DEBUG_VERBOSE
  std::cout<<"iw: ";
  for (int widx=0; widx < in_len; widx++) {
    std::cout<<in_buf[widx]<<"  ";
  }
  std::cout<<std::endl;
#endif
}
#endif // KLT_SUPPORT_WIN
