// Karhunen-Loève Transform Library
// S.P.Storck 07-20-2017

#ifndef __KLT_HH__

// DEBUG bitmask values
#define KLT_DEBUG_NONE 0
#define KLT_DEBUG_WIN 1
#define KLT_DEBUG_VERBOSE 2
#define KLT_DEBUG_FINE 4

// DEBUG bitmask
#define KLT_DEBUG (KLT_DEBUG_NONE)
// Compile KLT window code?
#define KLT_SUPPORT_WIN 0
// Support normalized eigenvalue output
#define KLT_SUPPORT_EVALN 0

#include <complex>

class KLT
{
public:
  //---------------------------------------------------------------------------
  // Constructor
  //   in_len: in_buf size.
  //   window:  apply window? (support may be compiled out)
  //   eval_normalized: normalize eigenvalues? (support may be compiled out)
  //   acm_order: auto-correlation matrix order.
  //   num_eig: number (from largest eigenvalue down) of eigenvalues/vectors to
  //           solve (must be <= acm_order).
  //   If an error occurrs, throws std::runtime_error.
  //---------------------------------------------------------------------------
  KLT(int in_len,
#if KLT_SUPPORT_WIN
      int window,
#endif
#if KLT_SUPPORT_EVALN
      int eval_normalized,
#endif
      int acm_order,
      int num_eig);

  //---------------------------------------------------------------------------
  // Destructor
  //---------------------------------------------------------------------------
  ~KLT();

  //---------------------------------------------------------------------------
  // Transform contents of in_buf, update eval_buf, kltc_buf, and kltb_buf.
  //   If an error occurrs, throws std::runtime_error and sets all of eval_buf,
  //   kltc_buf, and kltb_buf to 0.0f.
  //---------------------------------------------------------------------------
  void transform();

  //---------------------------------------------------------------------------
  // Input/Output Buffers
  //   in_buf: input buffer (size in_len).
  //   eval_buf: output eigenvalues row vector (size 1 x num_eig).
  //   kltc_buf: output KLT coeffs row vector (size 1 x num_eig).
  //   kltb_buf: output KLT basis functions matrix (size acm_order x num_eig)
  //             (already weighted by the KLT coeffs).  The basis functions are
  //             output as column vectors in column-major order (equivalent to
  //             row vectors in row-major order).
  //---------------------------------------------------------------------------
  std::complex<float>* in_buf;
  float* eval_buf;
  std::complex<float>* kltc_buf;
  std::complex<float>* kltb_buf;


private:
#if KLT_SUPPORT_WIN
  void init_window();
  void apply_window();
#endif
  void acorr_matrix();
  void eigendecomp();

  //---------------------------------------------------------------------------
  // Config
  //---------------------------------------------------------------------------
  const int in_len;
#if KLT_SUPPORT_WIN
  const int window;
#endif
#if KLT_SUPPORT_EVALN
  const int eval_normalized;
#endif
  const int acm_order;
  const int num_eig;

  //---------------------------------------------------------------------------
  // Internal/temp buffers
  //   win_buf: window buffer (size in_len).
  //   ac_buf: temp buffer (size ((acm_order+1)*acm_order)/2).
  //   d_buf: temp buffer (size acm_order).
  //   e_buf: temp buffer (size (acm_order-1)).
  //   tau_buf: temp buffer (size (acm_order-1)).
  //   ib_buf: temp buffer (size acm_order).
  //   is_buf: temp buffer (size acm_order).
  //   if_buf: temp buffer (size num_eig).
  //---------------------------------------------------------------------------
#if KLT_SUPPORT_WIN
  float* win_buf;
#endif
  std::complex<float>* ac_buf;
  float* d_buf;
  float* e_buf;
  std::complex<float>* tau_buf;
  int* ib_buf;
  int* is_buf;
  int* if_buf;
};

#endif // __KLT_HH__
