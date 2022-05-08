/*************************************************************************
ALGLIB 3.18.0 (source code generated 2021-10-25)
Copyright (c) Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
#ifndef _solvers_pkg_h
#define _solvers_pkg_h
#include "ap.h"
#include "alglibinternal.h"
#include "alglibmisc.h"
#include "linalg.h"

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (DATATYPES)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
#if defined(AE_COMPILE_POLYNOMIALSOLVER) || !defined(AE_PARTIAL_BUILD)
typedef struct
{
    double maxerr;
} polynomialsolverreport;
#endif
#if defined(AE_COMPILE_DIRECTDENSESOLVERS) || !defined(AE_PARTIAL_BUILD)
typedef struct
{
    double r1;
    double rinf;
} densesolverreport;
typedef struct
{
    double r2;
    ae_matrix cx;
    ae_int_t n;
    ae_int_t k;
} densesolverlsreport;
#endif
#if defined(AE_COMPILE_DIRECTSPARSESOLVERS) || !defined(AE_PARTIAL_BUILD)
typedef struct
{
    ae_int_t terminationtype;
    ae_int_t nmv;
    ae_int_t iterationscount;
    double r2;
} sparsesolverreport;
#endif
#if defined(AE_COMPILE_ITERATIVESPARSE) || !defined(AE_PARTIAL_BUILD)
typedef struct
{
    ae_int_t n;
    ae_vector x0;
    double epsf;
    ae_int_t maxits;
    ae_int_t algotype;
    ae_int_t gmresk;
    ae_bool xrep;
    ae_bool running;
    ae_bool userterminationneeded;
    ae_vector b;
    ae_vector xf;
    ae_int_t repiterationscount;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    double repr2;
    ae_int_t requesttype;
    ae_vector x;
    ae_vector ax;
    double reply1;
    ae_vector wrkb;
    sparsematrix convbuf;
    fblsgmresstate gmressolver;
    rcommstate rstate;
} sparsesolverstate;
#endif
#if defined(AE_COMPILE_LINCG) || !defined(AE_PARTIAL_BUILD)
typedef struct
{
    ae_vector rx;
    ae_vector b;
    ae_int_t n;
    ae_int_t prectype;
    ae_vector cx;
    ae_vector cr;
    ae_vector cz;
    ae_vector p;
    ae_vector r;
    ae_vector z;
    double alpha;
    double beta;
    double r2;
    double meritfunction;
    ae_vector x;
    ae_vector mv;
    ae_vector pv;
    double vmv;
    ae_vector startx;
    double epsf;
    ae_int_t maxits;
    ae_int_t itsbeforerestart;
    ae_int_t itsbeforerupdate;
    ae_bool xrep;
    ae_bool xupdated;
    ae_bool needmv;
    ae_bool needmtv;
    ae_bool needmv2;
    ae_bool needvmv;
    ae_bool needprec;
    ae_int_t repiterationscount;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    ae_bool running;
    ae_vector tmpd;
    rcommstate rstate;
} lincgstate;
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nmv;
    ae_int_t terminationtype;
    double r2;
} lincgreport;
#endif
#if defined(AE_COMPILE_LINLSQR) || !defined(AE_PARTIAL_BUILD)
typedef struct
{
    normestimatorstate nes;
    ae_vector rx;
    ae_vector b;
    ae_int_t n;
    ae_int_t m;
    ae_int_t prectype;
    ae_vector ui;
    ae_vector uip1;
    ae_vector vi;
    ae_vector vip1;
    ae_vector omegai;
    ae_vector omegaip1;
    double alphai;
    double alphaip1;
    double betai;
    double betaip1;
    double phibari;
    double phibarip1;
    double phii;
    double rhobari;
    double rhobarip1;
    double rhoi;
    double ci;
    double si;
    double theta;
    double lambdai;
    ae_vector d;
    double anorm;
    double bnorm2;
    double dnorm;
    double r2;
    ae_vector x;
    ae_vector mv;
    ae_vector mtv;
    double epsa;
    double epsb;
    double epsc;
    ae_int_t maxits;
    ae_bool xrep;
    ae_bool xupdated;
    ae_bool needmv;
    ae_bool needmtv;
    ae_bool needmv2;
    ae_bool needvmv;
    ae_bool needprec;
    ae_int_t repiterationscount;
    ae_int_t repnmv;
    ae_int_t repterminationtype;
    ae_bool running;
    ae_bool userterminationneeded;
    ae_vector tmpd;
    ae_vector tmpx;
    rcommstate rstate;
} linlsqrstate;
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nmv;
    ae_int_t terminationtype;
} linlsqrreport;
#endif
#if defined(AE_COMPILE_NLEQ) || !defined(AE_PARTIAL_BUILD)
typedef struct
{
    ae_int_t n;
    ae_int_t m;
    double epsf;
    ae_int_t maxits;
    ae_bool xrep;
    double stpmax;
    ae_vector x;
    double f;
    ae_vector fi;
    ae_matrix j;
    ae_bool needf;
    ae_bool needfij;
    ae_bool xupdated;
    rcommstate rstate;
    ae_int_t repiterationscount;
    ae_int_t repnfunc;
    ae_int_t repnjac;
    ae_int_t repterminationtype;
    ae_vector xbase;
    double fbase;
    double fprev;
    ae_vector candstep;
    ae_vector rightpart;
    ae_vector cgbuf;
} nleqstate;
typedef struct
{
    ae_int_t iterationscount;
    ae_int_t nfunc;
    ae_int_t njac;
    ae_int_t terminationtype;
} nleqreport;
#endif

}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{

#if defined(AE_COMPILE_POLYNOMIALSOLVER) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************

*************************************************************************/
class _polynomialsolverreport_owner
{
public:
    _polynomialsolverreport_owner();
    _polynomialsolverreport_owner(const _polynomialsolverreport_owner &rhs);
    _polynomialsolverreport_owner& operator=(const _polynomialsolverreport_owner &rhs);
    virtual ~_polynomialsolverreport_owner();
    alglib_impl::polynomialsolverreport* c_ptr();
    alglib_impl::polynomialsolverreport* c_ptr() const;
protected:
    alglib_impl::polynomialsolverreport *p_struct;
};
class polynomialsolverreport : public _polynomialsolverreport_owner
{
public:
    polynomialsolverreport();
    polynomialsolverreport(const polynomialsolverreport &rhs);
    polynomialsolverreport& operator=(const polynomialsolverreport &rhs);
    virtual ~polynomialsolverreport();
    double &maxerr;

};
#endif

#if defined(AE_COMPILE_DIRECTDENSESOLVERS) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************

*************************************************************************/
class _densesolverreport_owner
{
public:
    _densesolverreport_owner();
    _densesolverreport_owner(const _densesolverreport_owner &rhs);
    _densesolverreport_owner& operator=(const _densesolverreport_owner &rhs);
    virtual ~_densesolverreport_owner();
    alglib_impl::densesolverreport* c_ptr();
    alglib_impl::densesolverreport* c_ptr() const;
protected:
    alglib_impl::densesolverreport *p_struct;
};
class densesolverreport : public _densesolverreport_owner
{
public:
    densesolverreport();
    densesolverreport(const densesolverreport &rhs);
    densesolverreport& operator=(const densesolverreport &rhs);
    virtual ~densesolverreport();
    double &r1;
    double &rinf;

};


/*************************************************************************

*************************************************************************/
class _densesolverlsreport_owner
{
public:
    _densesolverlsreport_owner();
    _densesolverlsreport_owner(const _densesolverlsreport_owner &rhs);
    _densesolverlsreport_owner& operator=(const _densesolverlsreport_owner &rhs);
    virtual ~_densesolverlsreport_owner();
    alglib_impl::densesolverlsreport* c_ptr();
    alglib_impl::densesolverlsreport* c_ptr() const;
protected:
    alglib_impl::densesolverlsreport *p_struct;
};
class densesolverlsreport : public _densesolverlsreport_owner
{
public:
    densesolverlsreport();
    densesolverlsreport(const densesolverlsreport &rhs);
    densesolverlsreport& operator=(const densesolverlsreport &rhs);
    virtual ~densesolverlsreport();
    double &r2;
    real_2d_array cx;
    ae_int_t &n;
    ae_int_t &k;

};
#endif

#if defined(AE_COMPILE_DIRECTSPARSESOLVERS) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
This structure is a sparse solver report (both direct and iterative solvers
use this structure).

Following fields can be accessed by users:
* TerminationType (specific error codes depend on the solver  being  used,
  with positive values ALWAYS signaling  that something useful is returned
  in X, and negative values ALWAYS meaning critical failures.
* NMV - number of matrix-vector products performed (0 for direct solvers)
* IterationsCount - inner iterations count (0 for direct solvers)
* R2 - squared residual
*************************************************************************/
class _sparsesolverreport_owner
{
public:
    _sparsesolverreport_owner();
    _sparsesolverreport_owner(const _sparsesolverreport_owner &rhs);
    _sparsesolverreport_owner& operator=(const _sparsesolverreport_owner &rhs);
    virtual ~_sparsesolverreport_owner();
    alglib_impl::sparsesolverreport* c_ptr();
    alglib_impl::sparsesolverreport* c_ptr() const;
protected:
    alglib_impl::sparsesolverreport *p_struct;
};
class sparsesolverreport : public _sparsesolverreport_owner
{
public:
    sparsesolverreport();
    sparsesolverreport(const sparsesolverreport &rhs);
    sparsesolverreport& operator=(const sparsesolverreport &rhs);
    virtual ~sparsesolverreport();
    ae_int_t &terminationtype;
    ae_int_t &nmv;
    ae_int_t &iterationscount;
    double &r2;

};
#endif

#if defined(AE_COMPILE_ITERATIVESPARSE) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
This object stores state of the sparse linear solver object.

You should use ALGLIB functions to work with this object.
Never try to access its fields directly!
*************************************************************************/
class _sparsesolverstate_owner
{
public:
    _sparsesolverstate_owner();
    _sparsesolverstate_owner(const _sparsesolverstate_owner &rhs);
    _sparsesolverstate_owner& operator=(const _sparsesolverstate_owner &rhs);
    virtual ~_sparsesolverstate_owner();
    alglib_impl::sparsesolverstate* c_ptr();
    alglib_impl::sparsesolverstate* c_ptr() const;
protected:
    alglib_impl::sparsesolverstate *p_struct;
};
class sparsesolverstate : public _sparsesolverstate_owner
{
public:
    sparsesolverstate();
    sparsesolverstate(const sparsesolverstate &rhs);
    sparsesolverstate& operator=(const sparsesolverstate &rhs);
    virtual ~sparsesolverstate();

};
#endif

#if defined(AE_COMPILE_LINCG) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
This object stores state of the linear CG method.

You should use ALGLIB functions to work with this object.
Never try to access its fields directly!
*************************************************************************/
class _lincgstate_owner
{
public:
    _lincgstate_owner();
    _lincgstate_owner(const _lincgstate_owner &rhs);
    _lincgstate_owner& operator=(const _lincgstate_owner &rhs);
    virtual ~_lincgstate_owner();
    alglib_impl::lincgstate* c_ptr();
    alglib_impl::lincgstate* c_ptr() const;
protected:
    alglib_impl::lincgstate *p_struct;
};
class lincgstate : public _lincgstate_owner
{
public:
    lincgstate();
    lincgstate(const lincgstate &rhs);
    lincgstate& operator=(const lincgstate &rhs);
    virtual ~lincgstate();

};


/*************************************************************************

*************************************************************************/
class _lincgreport_owner
{
public:
    _lincgreport_owner();
    _lincgreport_owner(const _lincgreport_owner &rhs);
    _lincgreport_owner& operator=(const _lincgreport_owner &rhs);
    virtual ~_lincgreport_owner();
    alglib_impl::lincgreport* c_ptr();
    alglib_impl::lincgreport* c_ptr() const;
protected:
    alglib_impl::lincgreport *p_struct;
};
class lincgreport : public _lincgreport_owner
{
public:
    lincgreport();
    lincgreport(const lincgreport &rhs);
    lincgreport& operator=(const lincgreport &rhs);
    virtual ~lincgreport();
    ae_int_t &iterationscount;
    ae_int_t &nmv;
    ae_int_t &terminationtype;
    double &r2;

};
#endif

#if defined(AE_COMPILE_LINLSQR) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
This object stores state of the LinLSQR method.

You should use ALGLIB functions to work with this object.
*************************************************************************/
class _linlsqrstate_owner
{
public:
    _linlsqrstate_owner();
    _linlsqrstate_owner(const _linlsqrstate_owner &rhs);
    _linlsqrstate_owner& operator=(const _linlsqrstate_owner &rhs);
    virtual ~_linlsqrstate_owner();
    alglib_impl::linlsqrstate* c_ptr();
    alglib_impl::linlsqrstate* c_ptr() const;
protected:
    alglib_impl::linlsqrstate *p_struct;
};
class linlsqrstate : public _linlsqrstate_owner
{
public:
    linlsqrstate();
    linlsqrstate(const linlsqrstate &rhs);
    linlsqrstate& operator=(const linlsqrstate &rhs);
    virtual ~linlsqrstate();

};


/*************************************************************************

*************************************************************************/
class _linlsqrreport_owner
{
public:
    _linlsqrreport_owner();
    _linlsqrreport_owner(const _linlsqrreport_owner &rhs);
    _linlsqrreport_owner& operator=(const _linlsqrreport_owner &rhs);
    virtual ~_linlsqrreport_owner();
    alglib_impl::linlsqrreport* c_ptr();
    alglib_impl::linlsqrreport* c_ptr() const;
protected:
    alglib_impl::linlsqrreport *p_struct;
};
class linlsqrreport : public _linlsqrreport_owner
{
public:
    linlsqrreport();
    linlsqrreport(const linlsqrreport &rhs);
    linlsqrreport& operator=(const linlsqrreport &rhs);
    virtual ~linlsqrreport();
    ae_int_t &iterationscount;
    ae_int_t &nmv;
    ae_int_t &terminationtype;

};
#endif

#if defined(AE_COMPILE_NLEQ) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************

*************************************************************************/
class _nleqstate_owner
{
public:
    _nleqstate_owner();
    _nleqstate_owner(const _nleqstate_owner &rhs);
    _nleqstate_owner& operator=(const _nleqstate_owner &rhs);
    virtual ~_nleqstate_owner();
    alglib_impl::nleqstate* c_ptr();
    alglib_impl::nleqstate* c_ptr() const;
protected:
    alglib_impl::nleqstate *p_struct;
};
class nleqstate : public _nleqstate_owner
{
public:
    nleqstate();
    nleqstate(const nleqstate &rhs);
    nleqstate& operator=(const nleqstate &rhs);
    virtual ~nleqstate();
    ae_bool &needf;
    ae_bool &needfij;
    ae_bool &xupdated;
    double &f;
    real_1d_array fi;
    real_2d_array j;
    real_1d_array x;

};


/*************************************************************************

*************************************************************************/
class _nleqreport_owner
{
public:
    _nleqreport_owner();
    _nleqreport_owner(const _nleqreport_owner &rhs);
    _nleqreport_owner& operator=(const _nleqreport_owner &rhs);
    virtual ~_nleqreport_owner();
    alglib_impl::nleqreport* c_ptr();
    alglib_impl::nleqreport* c_ptr() const;
protected:
    alglib_impl::nleqreport *p_struct;
};
class nleqreport : public _nleqreport_owner
{
public:
    nleqreport();
    nleqreport(const nleqreport &rhs);
    nleqreport& operator=(const nleqreport &rhs);
    virtual ~nleqreport();
    ae_int_t &iterationscount;
    ae_int_t &nfunc;
    ae_int_t &njac;
    ae_int_t &terminationtype;

};
#endif

#if defined(AE_COMPILE_POLYNOMIALSOLVER) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
Polynomial root finding.

This function returns all roots of the polynomial
    P(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
Both real and complex roots are returned (see below).

INPUT PARAMETERS:
    A       -   array[N+1], polynomial coefficients:
                * A[0] is constant term
                * A[N] is a coefficient of X^N
    N       -   polynomial degree

OUTPUT PARAMETERS:
    X       -   array of complex roots:
                * for isolated real root, X[I] is strictly real: IMAGE(X[I])=0
                * complex roots are always returned in pairs - roots occupy
                  positions I and I+1, with:
                  * X[I+1]=Conj(X[I])
                  * IMAGE(X[I]) > 0
                  * IMAGE(X[I+1]) = -IMAGE(X[I]) < 0
                * multiple real roots may have non-zero imaginary part due
                  to roundoff errors. There is no reliable way to distinguish
                  real root of multiplicity 2 from two  complex  roots  in
                  the presence of roundoff errors.
    Rep     -   report, additional information, following fields are set:
                * Rep.MaxErr - max( |P(xi)| )  for  i=0..N-1.  This  field
                  allows to quickly estimate "quality" of the roots  being
                  returned.

NOTE:   this function uses companion matrix method to find roots. In  case
        internal EVD  solver  fails  do  find  eigenvalues,  exception  is
        generated.

NOTE:   roots are not "polished" and  no  matrix  balancing  is  performed
        for them.

  -- ALGLIB --
     Copyright 24.02.2014 by Bochkanov Sergey
*************************************************************************/
void polynomialsolve(const real_1d_array &a, const ae_int_t n, complex_1d_array &x, polynomialsolverreport &rep, const xparams _xparams = alglib::xdefault);
#endif

#if defined(AE_COMPILE_DIRECTDENSESOLVERS) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
Dense solver for A*x=b with N*N real matrix A and N*1 real vectorx  x  and
b. This is "slow-but-feature rich" version of the  linear  solver.  Faster
version is RMatrixSolveFast() function.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations. It is also very significant on small matrices.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, RMatrixSolveFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolve(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix, x  and  b  are  vectors.  This is a "fast" version of  linear
solver which does NOT provide  any  additional  functions  like  condition
number estimation or iterative refinement.

Algorithm features:
* efficient algorithm O(N^3) complexity
* no performance overhead from additional functionality

If you need condition number estimation or iterative refinement, use  more
feature-rich version - RMatrixSolve().

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 16.03.2015 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const real_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

Similar to RMatrixSolve() but solves task with multiple right parts (where
b and x are NxM matrices). This is  "slow-but-robust"  version  of  linear
solver with additional functionality  like  condition  number  estimation.
There also exists faster version - RMatrixSolveMFast().

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* optional iterative refinement
* O(N^3+M*N^2) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations. It also very significant on small matrices.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, RMatrixSolveMFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is ill conditioned or singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvem(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

Similar to RMatrixSolve() but solves task with multiple right parts (where
b and x are NxM matrices). This is "fast" version of linear  solver  which
does NOT offer additional functions like condition  number  estimation  or
iterative refinement.

Algorithm features:
* O(N^3+M*N^2) complexity
* no additional functionality, highest performance

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix given by its LU decomposition, x and b are real vectors.  This
is "slow-but-robust" version of the linear LU-based solver. Faster version
is RMatrixLUSolveFast() function.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! RMatrixLUSolveFast() function.

INPUT PARAMETERS
    LUA     -   array[N,N], LU decomposition, RMatrixLU result
    P       -   array[N], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[N], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolve(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix given by its LU decomposition, x and b are real vectors.  This
is "fast-without-any-checks" version of the linear LU-based solver. Slower
but more robust version is RMatrixLUSolve() function.

Algorithm features:
* O(N^2) complexity
* fast algorithm without ANY additional checks, just triangular solver

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvefast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

Similar to RMatrixLUSolve() but solves  task  with  multiple  right  parts
(where b and x are NxM matrices). This  is  "robust-but-slow"  version  of
LU-based solver which performs additional  checks  for  non-degeneracy  of
inputs (condition number estimation). If you need  best  performance,  use
"fast-without-any-checks" version, RMatrixLUSolveMFast().

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver.
           !
           ! This performance penalty is especially apparent when you  use
           ! ALGLIB parallel capabilities (condition number estimation  is
           ! inherently  sequential).  It   also   becomes significant for
           ! small-scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! RMatrixLUSolveMFast() function.

INPUT PARAMETERS
    LUA     -   array[N,N], LU decomposition, RMatrixLU result
    P       -   array[N], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvem(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

Similar to RMatrixLUSolve() but solves  task  with  multiple  right parts,
where b and x are NxM matrices.  This is "fast-without-any-checks" version
of LU-based solver. It does not estimate  condition number  of  a  system,
so it is extremely fast. If you need better detection  of  near-degenerate
cases, use RMatrixLUSolveM() function.

Algorithm features:
* O(M*N^2) complexity
* fast algorithm without ANY additional checks, just triangular solver

INPUT PARAMETERS:
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N,M]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvemfast(const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where BOTH ORIGINAL A AND ITS
LU DECOMPOSITION ARE KNOWN. You can use it if for some  reasons  you  have
both A and its LU decomposition.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolve(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

Similar to RMatrixMixedSolve() but  solves task with multiple right  parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolvem(const real_2d_array &a, const real_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Complex dense solver for A*X=B with N*N  complex  matrix  A,  N*M  complex
matrices  X  and  B.  "Slow-but-feature-rich"   version   which   provides
additional functions, at the cost of slower  performance.  Faster  version
may be invoked with CMatrixSolveMFast() function.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3+M*N^2) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, CMatrixSolveMFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, const bool rfs, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Complex dense solver for A*X=B with N*N  complex  matrix  A,  N*M  complex
matrices  X  and  B.  "Fast-but-lightweight" version which  provides  just
triangular solver - and no additional functions like iterative  refinement
or condition number estimation.

Algorithm features:
* O(N^3+M*N^2) complexity
* no additional time consuming functions

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N,M]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 16.03.2015 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Complex dense solver for A*x=B with N*N complex matrix A and  N*1  complex
vectors x and b. "Slow-but-feature-rich" version of the solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear  system
           ! and  performs  iterative   refinement,   which   results   in
           ! significant performance penalty  when  compared  with  "fast"
           ! version  which  just  performs  LU  decomposition  and  calls
           ! triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, CMatrixSolveFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolve(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Complex dense solver for A*x=B with N*N complex matrix A and  N*1  complex
vectors x and b. "Fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^3) complexity
* no additional time consuming features, just triangular solver

INPUT PARAMETERS:
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N complex A given by its  LU  decomposition,
and N*M matrices X and B (multiple right sides).   "Slow-but-feature-rich"
version of the solver.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver.
           !
           ! This performance penalty is especially apparent when you  use
           ! ALGLIB parallel capabilities (condition number estimation  is
           ! inherently  sequential).  It   also   becomes significant for
           ! small-scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! CMatrixLUSolveMFast() function.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvem(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N complex A given by its  LU  decomposition,
and N*M matrices X and B (multiple  right  sides).  "Fast-but-lightweight"
version of the solver.

Algorithm features:
* O(M*N^2) complexity
* no additional time-consuming features

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N,M]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvemfast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Complex dense linear solver for A*x=b with complex N*N A  given  by its LU
decomposition and N*1 vectors x and b. This is  "slow-but-robust"  version
of  the  complex  linear  solver  with  additional  features   which   add
significant performance overhead. Faster version  is  CMatrixLUSolveFast()
function.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems.
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! CMatrixLUSolveFast() function.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolve(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Complex dense linear solver for A*x=b with N*N complex A given by  its  LU
decomposition and N*1 vectors x and b. This is  fast  lightweight  version
of solver, which is significantly faster than CMatrixLUSolve(),  but  does
not provide additional information (like condition numbers).

Algorithm features:
* O(N^2) complexity
* no additional time-consuming features, just triangular solver

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is exactly singular (ill conditioned matrices
                        are not recognized).
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * info>0    =>  overwritten by solution
                * info=-3   =>  filled by zeros

NOTE: unlike  CMatrixLUSolve(),  this   function   does   NOT   check  for
      near-degeneracy of input matrix. It  checks  for  EXACT  degeneracy,
      because this check is easy to do. However,  very  badly  conditioned
      matrices may went unnoticed.


  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvefast(const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver. Same as RMatrixMixedSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolvem(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver. Same as RMatrixMixedSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or exactly singular.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolve(const complex_2d_array &a, const complex_2d_array &lua, const integer_1d_array &p, const ae_int_t n, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A,  and
N*M vectors X and B. It is "slow-but-feature-rich" version of the solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant   performance   penalty  when
           ! compared with "fast" version  which  just  performs  Cholesky
           ! decomposition and calls triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, SPDMatrixSolveMFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or non-SPD.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvem(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A,  and
N*M vectors X and B. It is "fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming features

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N,M], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvemfast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense linear solver for A*x=b with N*N real  symmetric  positive  definite
matrix A,  N*1 vectors x and b.  "Slow-but-feature-rich"  version  of  the
solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant   performance   penalty  when
           ! compared with "fast" version  which  just  performs  Cholesky
           ! decomposition and calls triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, SPDMatrixSolveFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    matrix is very badly conditioned or non-SPD.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolve(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense linear solver for A*x=b with N*N real  symmetric  positive  definite
matrix A,  N*1 vectors x and  b.  "Fast-but-lightweight"  version  of  the
solver.

Algorithm features:
* O(N^3) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming features like condition number estimation

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or non-SPD
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N], it contains:
                * info>0    =>  solution
                * info=-3   =>  filled by zeros

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvefast(const real_2d_array &a, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*M vectors X and B. It  is  "slow-but-
feature-rich" version of the solver which estimates  condition  number  of
the system.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver. Amount of  overhead  introduced  depends  on  M  (the
           ! larger - the more efficient).
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! SPDMatrixCholeskySolveMFast() function.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or badly conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0 contains solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvem(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, real_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*M vectors X and B. It  is  "fast-but-
lightweight" version of  the  solver  which  just  solves  linear  system,
without any additional functions.

Algorithm features:
* O(M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional functionality

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[N,M], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or badly conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N]:
                * for info>0 overwritten by solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvemfast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*x=b with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*1 real vectors x and b. This is "slow-
but-feature-rich"  version  of  the  solver  which,  in  addition  to  the
solution, performs condition number estimation.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! SPDMatrixCholeskySolveFast() function.

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[N], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0  - solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolve(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, densesolverreport &rep, real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*x=b with N*N symmetric positive definite matrix A given
by its Cholesky decomposition, and N*1 real vectors x and b. This is "fast-
but-lightweight" version of the solver.

Algorithm features:
* O(N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional features

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[N], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * for info>0  - overwritten by solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvefast(const real_2d_array &cha, const ae_int_t n, const bool isupper, const real_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B, with N*N Hermitian positive definite matrix A  and
N*M  complex  matrices  X  and  B.  "Slow-but-feature-rich" version of the
solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver.
           !
           ! This performance penalty is especially apparent when you  use
           ! ALGLIB parallel capabilities (condition number estimation  is
           ! inherently  sequential).  It   also   becomes significant for
           ! small-scale problems (N<100).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! HPDMatrixSolveMFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvem(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B, with N*N Hermitian positive definite matrix A  and
N*M complex matrices X and B. "Fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming features like condition number estimation

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly  singular or is not positive definite.
                        B is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[0..N-1]:
                * overwritten by solution
                * zeros, if problem was not solved

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvemfast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*x=b, with N*N Hermitian positive definite matrix A, and
N*1 complex vectors  x  and  b.  "Slow-but-feature-rich"  version  of  the
solver.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant   performance   penalty  when
           ! compared with "fast" version  which  just  performs  Cholesky
           ! decomposition and calls triangular solver.
           !
           ! This  performance  penalty  is  especially  visible  in   the
           ! multithreaded mode, because both condition number  estimation
           ! and   iterative    refinement   are   inherently   sequential
           ! calculations.
           !
           ! Thus, if you need high performance and if you are pretty sure
           ! that your system is well conditioned, we  strongly  recommend
           ! you to use faster solver, HPDMatrixSolveFast() function.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolve(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*x=b, with N*N Hermitian positive definite matrix A, and
N*1 complex vectors  x  and  b.  "Fast-but-lightweight"  version  of   the
solver without additional functions.

Algorithm features:
* O(N^3) complexity
* matrix is represented by its upper or lower triangle
* no additional time consuming functions

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or not positive definite
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[0..N-1]:
                * overwritten by solution
                * zeros, if A is exactly singular (diagonal of its LU
                  decomposition has exact zeros).

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 17.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvefast(const complex_2d_array &a, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition and N*M complex matrices X  and  B.  This is
"slow-but-feature-rich" version of the solver which, in  addition  to  the
solution, estimates condition number of the system.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which  results  in  significant  performance   penalty   when
           ! compared with "fast"  version  which  just  calls  triangular
           ! solver. Amount of  overhead  introduced  depends  on  M  (the
           ! larger - the more efficient).
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large Cholesky decomposition.  However,  if  you call
           ! this  function  many  times  for  the same  left  side,  this
           ! overhead BECOMES significant. It  also   becomes  significant
           ! for small-scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! HPDMatrixCholeskySolveMFast() function.


INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                HPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[N,M], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0 contains solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvem(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, densesolverreport &rep, complex_2d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*X=B with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition and N*M complex matrices X  and  B.  This is
"fast-but-lightweight" version of the solver.

Algorithm features:
* O(M*N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time-consuming features

INPUT PARAMETERS
    CHA     -   array[N,N], Cholesky decomposition,
                HPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[N,M], right part
    M       -   right part size

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task was solved
    B       -   array[N]:
                * for info>0 overwritten by solution
                * for info=-3 filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvemfast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_2d_array &b, const ae_int_t m, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*x=b with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition, and N*1 complex vectors x and  b.  This  is
"slow-but-feature-rich" version of the solver  which  estimates  condition
number of the system.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

IMPORTANT: ! this function is NOT the most efficient linear solver provided
           ! by ALGLIB. It estimates condition  number  of  linear system,
           ! which results in 10-15x  performance  penalty  when  compared
           ! with "fast" version which just calls triangular solver.
           !
           ! This performance penalty is insignificant  when compared with
           ! cost of large LU decomposition.  However,  if you  call  this
           ! function many times for the same  left  side,  this  overhead
           ! BECOMES significant. It  also  becomes significant for small-
           ! scale problems (N<50).
           !
           ! In such cases we strongly recommend you to use faster solver,
           ! HPDMatrixCholeskySolveFast() function.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    Rep     -   additional report, following fields are set:
                * rep.r1    condition number in 1-norm
                * rep.rinf  condition number in inf-norm
    X       -   array[N]:
                * for info>0  - solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolve(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, densesolverreport &rep, complex_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver for A*x=b with N*N Hermitian positive definite matrix A given
by its Cholesky decomposition, and N*1 complex vectors x and  b.  This  is
"fast-but-lightweight" version of the solver.

Algorithm features:
* O(N^2) complexity
* matrix is represented by its upper or lower triangle
* no additional time-consuming features

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is is exactly singular or ill conditioned
                        B is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved
    B       -   array[N]:
                * for info>0  - overwritten by solution
                * for info=-3 - filled by zeros

  -- ALGLIB --
     Copyright 18.03.2015 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvefast(const complex_2d_array &cha, const ae_int_t n, const bool isupper, const complex_1d_array &b, ae_int_t &info, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Dense solver.

This subroutine finds solution of the linear system A*X=B with non-square,
possibly degenerate A.  System  is  solved in the least squares sense, and
general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
returned. If A is non-degenerate, solution in the usual sense is returned.

Algorithm features:
* automatic detection (and correct handling!) of degenerate cases
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..NRows-1,0..NCols-1], system matrix
    NRows   -   vertical size of A
    NCols   -   horizontal size of A
    B       -   array[0..NCols-1], right part
    Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
                considered  zero.  Set  it to 0.0, if you don't understand
                what it means, so the solver will choose good value on its
                own.

OUTPUT PARAMETERS
    Info    -   return code:
                * -4    SVD subroutine failed
                * -1    if NRows<=0 or NCols<=0 or Threshold<0 was passed
                *  1    if task is solved
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1,0..M-1], it contains:
                * solution of A*X=B (even for singular A)
                * zeros, if SVD subroutine failed

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R2        reciprocal of condition number: 1/cond(A), 2-norm.
* N         = NCols
* K         dim(Null(A))
* CX        array[0..N-1,0..K-1], kernel of A.
            Columns of CX store such vectors that A*CX[i]=0.

  ! FREE EDITION OF ALGLIB:
  !
  ! Free Edition of ALGLIB supports following important features for  this
  ! function:
  ! * C++ version: x64 SIMD support using C++ intrinsics
  ! * C#  version: x64 SIMD support using NET5/NetCore hardware intrinsics
  !
  ! We  recommend  you  to  read  'Compiling ALGLIB' section of the ALGLIB
  ! Reference Manual in order  to  find  out  how to activate SIMD support
  ! in ALGLIB.

  ! COMMERCIAL EDITION OF ALGLIB:
  !
  ! Commercial Edition of ALGLIB includes following important improvements
  ! of this function:
  ! * high-performance native backend with same C# interface (C# version)
  ! * multithreading support (C++ and C# versions)
  ! * hardware vendor (Intel) implementations of linear algebra primitives
  !   (C++ and C# versions, x86/x64 platform)
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvels(const real_2d_array &a, const ae_int_t nrows, const ae_int_t ncols, const real_1d_array &b, const double threshold, ae_int_t &info, densesolverlsreport &rep, real_1d_array &x, const xparams _xparams = alglib::xdefault);
#endif

#if defined(AE_COMPILE_DIRECTSPARSESOLVERS) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
Sparse linear solver for A*x=b with N*N  sparse  real  symmetric  positive
definite matrix A, N*1 vectors x and b.

This solver  converts  input  matrix  to  SKS  format,  performs  Cholesky
factorization using  SKS  Cholesky  subroutine  (works  well  for  limited
bandwidth matrices) and uses sparse triangular solvers to get solution  of
the original system.

INPUT PARAMETERS
    A       -   sparse matrix, must be NxN exactly
    IsUpper -   which half of A is provided (another half is ignored)
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    X       -   array[N], it contains:
                * rep.terminationtype>0    =>  solution
                * rep.terminationtype=-3   =>  filled by zeros
    Rep     -   solver report, following fields are set:
                * rep.terminationtype - solver status; >0 for success,
                  set to -3 on failure (degenerate or non-SPD system).

  -- ALGLIB --
     Copyright 26.12.2017 by Bochkanov Sergey
*************************************************************************/
void sparsespdsolvesks(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Sparse linear solver for A*x=b with N*N  sparse  real  symmetric  positive
definite matrix A, N*1 vectors x and b.

This solver  converts  input  matrix  to  CRS  format,  performs  Cholesky
factorization using supernodal Cholesky  decomposition  with  permutation-
reducing ordering and uses sparse triangular solver to get solution of the
original system.

INPUT PARAMETERS
    A       -   sparse matrix, must be NxN exactly
    IsUpper -   which half of A is provided (another half is ignored)
    B       -   array[N], right part

OUTPUT PARAMETERS
    X       -   array[N], it contains:
                * rep.terminationtype>0    =>  solution
                * rep.terminationtype=-3   =>  filled by zeros
    Rep     -   solver report, following fields are set:
                * rep.terminationtype - solver status; >0 for success,
                  set to -3 on failure (degenerate or non-SPD system).

  -- ALGLIB --
     Copyright 26.12.2017 by Bochkanov Sergey
*************************************************************************/
void sparsespdsolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Sparse linear solver for A*x=b with N*N real  symmetric  positive definite
matrix A given by its Cholesky decomposition, and N*1 vectors x and b.

IMPORTANT: this solver requires input matrix to be in  the  SKS  (Skyline)
           or CRS (compressed row storage) format. An  exception  will  be
           generated if you pass matrix in some other format.

INPUT PARAMETERS
    A       -   sparse NxN matrix stored in CRs or SKS format, must be NxN
                exactly
    IsUpper -   which half of A is provided (another half is ignored)
    B       -   array[N], right part

OUTPUT PARAMETERS
    X       -   array[N], it contains:
                * rep.terminationtype>0    =>  solution
                * rep.terminationtype=-3   =>  filled by zeros
    Rep     -   solver report, following fields are set:
                * rep.terminationtype - solver status; >0 for success,
                  set to -3 on failure (degenerate or non-SPD system).

  -- ALGLIB --
     Copyright 26.12.2017 by Bochkanov Sergey
*************************************************************************/
void sparsespdcholeskysolve(const sparsematrix &a, const bool isupper, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Sparse linear solver for A*x=b with general (nonsymmetric) N*N sparse real
matrix A, N*1 vectors x and b.

This solver converts input matrix to CRS format, performs LU factorization
and uses sparse triangular solvers to get solution of the original system.

INPUT PARAMETERS
    A       -   sparse matrix, must be NxN exactly, any storage format
    N       -   size of A, N>0
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    X       -   array[N], it contains:
                * rep.terminationtype>0    =>  solution
                * rep.terminationtype=-3   =>  filled by zeros
    Rep     -   solver report, following fields are set:
                * rep.terminationtype - solver status; >0 for success,
                  set to -3 on failure (degenerate system).

  -- ALGLIB --
     Copyright 26.12.2017 by Bochkanov Sergey
*************************************************************************/
void sparsesolve(const sparsematrix &a, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Sparse linear solver for A*x=b with general (nonsymmetric) N*N sparse real
matrix A given by its LU factorization, N*1 vectors x and b.

IMPORTANT: this solver requires input matrix  to  be  in  the  CRS  sparse
           storage format. An exception will  be  generated  if  you  pass
           matrix in some other format (HASH or SKS).

INPUT PARAMETERS
    A       -   LU factorization of the sparse matrix, must be NxN exactly
                in CRS storage format
    P, Q    -   pivot indexes from LU factorization
    N       -   size of A, N>0
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    X       -   array[N], it contains:
                * rep.terminationtype>0    =>  solution
                * rep.terminationtype=-3   =>  filled by zeros
    Rep     -   solver report, following fields are set:
                * rep.terminationtype - solver status; >0 for success,
                  set to -3 on failure (degenerate system).

  -- ALGLIB --
     Copyright 26.12.2017 by Bochkanov Sergey
*************************************************************************/
void sparselusolve(const sparsematrix &a, const integer_1d_array &p, const integer_1d_array &q, const real_1d_array &b, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);
#endif

#if defined(AE_COMPILE_ITERATIVESPARSE) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
Solving sparse symmetric linear system A*x=b using GMRES(k) method. Sparse
symmetric A is given by its lower or upper triangle.

NOTE: use SparseSolveGMRES() to solve system with nonsymmetric A.

This function provides convenience API for an 'expert' interface  provided
by SparseSolverState class. Use SparseSolver  API  if  you  need  advanced
functions like providing initial point, using out-of-core API and so on.

INPUT PARAMETERS:
    A       -   sparse symmetric NxN matrix in any sparse storage  format.
                Using CRS format is recommended because it avoids internal
                conversion.
                An exception will be generated if  A  is  not  NxN  matrix
                (where  N  is  a  size   specified  during  solver  object
                creation).
    IsUpper -   whether upper or lower triangle of A is used:
                * IsUpper=True  => only upper triangle is used and lower
                                   triangle is not referenced at all
                * IsUpper=False => only lower triangle is used and upper
                                   triangle is not referenced at all
    B       -   right part, array[N]
    K       -   k parameter for  GMRES(k), k>=0.  Zero  value  means  that
                algorithm will choose it automatically.
    EpsF    -   stopping condition, EpsF>=0. The algorithm will stop  when
                residual will decrease below EpsF*|B|. Having EpsF=0 means
                that this stopping condition is ignored.
    MaxIts  -   stopping condition, MaxIts>=0.  The  algorithm  will  stop
                after performing MaxIts iterations. Zero  value  means  no
                limit.

NOTE: having both EpsF=0 and MaxIts=0 means that stopping criteria will be
      chosen automatically.

OUTPUT PARAMETERS:
    X       -   array[N], the solution
    Rep     -   solution report:
                * Rep.TerminationType completion code:
                    * -5    CG method was used for a matrix which  is  not
                            positive definite
                    * -4    overflow/underflow during solution
                            (ill conditioned problem)
                    *  1    ||residual||<=EpsF*||b||
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            best point found is returned
                    *  8    the  algorithm  was  terminated   early  with
                            SparseSolverRequestTermination() being called
                            from other thread.
                * Rep.IterationsCount contains iterations count
                * Rep.NMV contains number of matrix-vector calculations
                * Rep.R2 contains squared residual

  -- ALGLIB --
     Copyright 25.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolvesymmetricgmres(const sparsematrix &a, const bool isupper, const real_1d_array &b, const ae_int_t k, const double epsf, const ae_int_t maxits, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Solving sparse linear system A*x=b using GMRES(k) method.

This function provides convenience API for an 'expert' interface  provided
by SparseSolverState class. Use SparseSolver  API  if  you  need  advanced
functions like providing initial point, using out-of-core API and so on.

INPUT PARAMETERS:
    A       -   sparse NxN matrix in any sparse storage format. Using  CRS
                format   is   recommended   because   it  avoids  internal
                conversion.
                An exception will be generated if  A  is  not  NxN  matrix
                (where  N  is  a  size   specified  during  solver  object
                creation).
    B       -   right part, array[N]
    K       -   k parameter for  GMRES(k), k>=0.  Zero  value  means  that
                algorithm will choose it automatically.
    EpsF    -   stopping condition, EpsF>=0. The algorithm will stop  when
                residual will decrease below EpsF*|B|. Having EpsF=0 means
                that this stopping condition is ignored.
    MaxIts  -   stopping condition, MaxIts>=0.  The  algorithm  will  stop
                after performing MaxIts iterations. Zero  value  means  no
                limit.

NOTE: having both EpsF=0 and MaxIts=0 means that stopping criteria will be
      chosen automatically.

OUTPUT PARAMETERS:
    X       -   array[N], the solution
    Rep     -   solution report:
                * Rep.TerminationType completion code:
                    * -5    CG method was used for a matrix which  is  not
                            positive definite
                    * -4    overflow/underflow during solution
                            (ill conditioned problem)
                    *  1    ||residual||<=EpsF*||b||
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            best point found is returned
                    *  8    the  algorithm  was  terminated   early  with
                            SparseSolverRequestTermination() being called
                            from other thread.
                * Rep.IterationsCount contains iterations count
                * Rep.NMV contains number of matrix-vector calculations
                * Rep.R2 contains squared residual

  -- ALGLIB --
     Copyright 25.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolvegmres(const sparsematrix &a, const real_1d_array &b, const ae_int_t k, const double epsf, const ae_int_t maxits, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function initializes sparse linear iterative solver object.

This solver can be used  to  solve  nonsymmetric  and  symmetric  positive
definite NxN (square) linear systems.

The solver provides  'expert'  API  which  allows  advanced  control  over
algorithms being used, including ability to get progress report, terminate
long-running solver from other thread, out-of-core solution and so on.

NOTE: there are also convenience  functions  that  allows  quick  one-line
      access to the solvers:
      * SparseSolveCG() to solve SPD linear systems
      * SparseSolveGMRES() to solve unsymmetric linear systems.

NOTE: if you want to solve MxN (rectangular) linear problem  you  may  use
      LinLSQR solver provided by ALGLIB.

USAGE (A is given by the SparseMatrix structure):

    1. User initializes algorithm state with SparseSolverCreate() call
    2. User  selects   algorithm  with one of the SparseSolverSetAlgo???()
       functions. By default, GMRES(k) is used with automatically chosen k
    3. Optionally, user tunes solver parameters, sets starting point, etc.
    4. Depending on whether system is symmetric or not, user calls:
       * SparseSolverSolveSymmetric() for a  symmetric system given by its
         lower or upper triangle
       * SparseSolverSolve() for a nonsymmetric system or a symmetric  one
         given by the full matrix
    5. User calls SparseSolverResults() to get the solution

    It is possible to call SparseSolverSolve???() again to  solve  another
    task with same dimensionality but different matrix and/or  right  part
    without reinitializing SparseSolverState structure.

USAGE (out-of-core mode):

    1. User initializes algorithm state with SparseSolverCreate() call
    2. User  selects   algorithm  with one of the SparseSolverSetAlgo???()
       functions. By default, GMRES(k) is used with automatically chosen k
    3. Optionally, user tunes solver parameters, sets starting point, etc.
    4. After that user should work with out-of-core interface  in  a  loop
       like one given below:

        > alglib.sparsesolveroocstart(state)
        > while alglib.sparsesolverooccontinue(state) do
        >     alglib.sparsesolveroocgetrequestinfo(state, out RequestType)
        >     alglib.sparsesolveroocgetrequestdata(state, out X)
        >     if RequestType=0 then
        >         [calculate  Y=A*X, with X=R^N]
        >     alglib.sparsesolveroocsendresult(state, in Y)
        > alglib.sparsesolveroocstop(state, out X, out Report)

INPUT PARAMETERS:
    N       -   problem dimensionality (fixed at start-up)

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolvercreate(const ae_int_t n, sparsesolverstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets the solver algorithm to GMRES(k).

NOTE: if you do not need advanced functionality of the  SparseSolver  API,
      you   may   use   convenience   functions   SparseSolveGMRES()   and
      SparseSolveSymmetricGMRES().

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    K       -   GMRES parameter, K>=0:
                * recommended values are in 10..100 range
                * larger values up to N are possible but have little sense
                  - the algorithm will be slower than any dense solver.
                * values above N are truncated down to N
                * zero value means that  default  value  is  chosen.  This
                  value is 50 in the current version, but  it  may  change
                  in future ALGLIB releases.

  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolversetalgogmres(const sparsesolverstate &state, const ae_int_t k, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets starting point.
By default, zero starting point is used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    X       -   starting point, array[N]

OUTPUT PARAMETERS:
    State   -   new starting point was set

  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolversetstartingpoint(const sparsesolverstate &state, const real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets stopping criteria.

INPUT PARAMETERS:
    EpsF    -   algorithm will be stopped if norm of residual is less than
                EpsF*||b||.
    MaxIts  -   algorithm will be stopped if number of iterations is  more
                than MaxIts.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
If  both  EpsF  and  MaxIts  are  zero then small EpsF will be set to small
value.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void sparsesolversetcond(const sparsesolverstate &state, const double epsf, const ae_int_t maxits, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Procedure for  the  solution of A*x=b with sparse symmetric A given by its
lower or upper triangle.

This function will work with any solver algorithm  being   used,  SPD  one
(like CG) or not (like GMRES). Using unsymmetric solvers (like  GMRES)  on
SPD problems is suboptimal, but still possible.

NOTE: the  solver  behavior is ill-defined  for  a  situation  when a  SPD
      solver is used on indefinite matrix. It  may solve the problem up to
      desired precision (sometimes, rarely)  or  return  with  error  code
      signalling violation of underlying assumptions.

INPUT PARAMETERS:
    State   -   algorithm state
    A       -   sparse symmetric NxN matrix in any sparse storage  format.
                Using CRS format is recommended because it avoids internal
                conversion.
                An exception will be generated if  A  is  not  NxN  matrix
                (where  N  is  a  size   specified  during  solver  object
                creation).
    IsUpper -   whether upper or lower triangle of A is used:
                * IsUpper=True  => only upper triangle is used and lower
                                   triangle is not referenced at all
                * IsUpper=False => only lower triangle is used and upper
                                   triangle is not referenced at all
    B       -   right part, array[N]

RESULT:
    This function returns no result.
    You can get the solution by calling SparseSolverResults()

  -- ALGLIB --
     Copyright 25.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolversolvesymmetric(const sparsesolverstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Procedure for the solution of A*x=b with sparse nonsymmetric A

IMPORTANT: this function will work with any solver algorithm  being  used,
           symmetric solver like CG,  or  not.  However,  using  symmetric
           solvers on nonsymmetric problems is  dangerous.  It  may  solve
           the problem up  to  desired  precision  (sometimes,  rarely) or
           terminate with error code signalling  violation  of  underlying
           assumptions.

INPUT PARAMETERS:
    State   -   algorithm state
    A       -   sparse NxN matrix in any sparse storage  format.
                Using CRS format is recommended because it avoids internal
                conversion.
                An exception will be generated if  A  is  not  NxN  matrix
                (where  N  is  a  size   specified  during  solver  object
                creation).
    B       -   right part, array[N]

RESULT:
    This function returns no result.
    You can get the solution by calling SparseSolverResults()

  -- ALGLIB --
     Copyright 25.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolversolve(const sparsesolverstate &state, const sparsematrix &a, const real_1d_array &b, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Sparse solver results.

This function must be called after calling one of the SparseSolverSolve()
functions.

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[N], solution
    Rep     -   solution report:
                * Rep.TerminationType completion code:
                    * -5    CG method was used for a matrix which  is  not
                            positive definite
                    * -4    overflow/underflow during solution
                            (ill conditioned problem)
                    *  1    ||residual||<=EpsF*||b||
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            best point found is returned
                    *  8    the  algorithm  was  terminated   early  with
                            SparseSolverRequestTermination() being called
                            from other thread.
                * Rep.IterationsCount contains iterations count
                * Rep.NMV contains number of matrix-vector calculations
                * Rep.R2 contains squared residual
s
  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void sparsesolverresults(const sparsesolverstate &state, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function turns on/off reporting during out-of-core processing.

When the solver works in the out-of-core mode, it  can  be  configured  to
report its progress by returning current location. These location  reports
are implemented as a special kind of the out-of-core request:
* SparseSolverOOCGetRequestInfo() returns -1
* SparseSolverOOCGetRequestData() returns current location
* SparseSolverOOCGetRequestData1() returns squared norm of the residual
* SparseSolverOOCSendResult() shall NOT be called

This function has no effect when SparseSolverSolve() is used because  this
function has no method of reporting its progress.

NOTE: when used with GMRES(k), this function reports progress  every  k-th
      iteration.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

  -- ALGLIB --
     Copyright 01.10.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolversetxrep(const sparsesolverstate &state, const bool needxrep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function initiates out-of-core mode of the sparse solver.  It  should
be used in conjunction with other out-of-core-related  functions  of  this
subspackage in a loop like one given below:

> alglib.sparsesolveroocstart(state)
> while alglib.sparsesolverooccontinue(state) do
>     alglib.sparsesolveroocgetrequestinfo(state, out RequestType)
>     alglib.sparsesolveroocgetrequestdata(state, out X)
>     if RequestType=0 then
>         [calculate  Y=A*X, with X=R^N]
>     alglib.sparsesolveroocsendresult(state, in Y)
> alglib.sparsesolveroocstop(state, out X, out Report)

INPUT PARAMETERS:
    State       -   solver object

  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolveroocstart(const sparsesolverstate &state, const real_1d_array &b, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function performs iterative solution of  the  linear  system  in  the
out-of-core mode. It should be used in conjunction with other out-of-core-
related functions of this subspackage in a loop like one given below:

> alglib.sparsesolveroocstart(state)
> while alglib.sparsesolverooccontinue(state) do
>     alglib.sparsesolveroocgetrequestinfo(state, out RequestType)
>     alglib.sparsesolveroocgetrequestdata(state, out X)
>     if RequestType=0 then
>         [calculate  Y=A*X, with X=R^N]
>     alglib.sparsesolveroocsendresult(state, in Y)
> alglib.sparsesolveroocstop(state, out X, out Report)

  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
bool sparsesolverooccontinue(const sparsesolverstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function is used to retrieve information  about  out-of-core  request
sent by the solver:
* RequestType=0  means that matrix-vector products A*x is requested
* RequestType=-1 means that solver reports its progress; this  request  is
  returned only when reports are activated wit SparseSolverSetXRep().

This function returns just request type; in order  to  get contents of the
trial vector, use sparsesolveroocgetrequestdata().

It should be used in conjunction with other out-of-core-related  functions
of this subspackage in a loop like one given below:

> alglib.sparsesolveroocstart(state)
> while alglib.sparsesolverooccontinue(state) do
>     alglib.sparsesolveroocgetrequestinfo(state, out RequestType)
>     alglib.sparsesolveroocgetrequestdata(state, out X)
>     if RequestType=0 then
>         [calculate  Y=A*X, with X=R^N]
>     alglib.sparsesolveroocsendresult(state, in Y)
> alglib.sparsesolveroocstop(state, out X, out Report)

INPUT PARAMETERS:
    State           -   solver running in out-of-core mode

OUTPUT PARAMETERS:
    RequestType     -   type of the request to process:
                        * 0   for matrix-vector product A*x, with A  being
                          NxN system matrix  and X being N-dimensional
                          vector
                        *-1   for location and residual report


  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolveroocgetrequestinfo(const sparsesolverstate &state, ae_int_t &requesttype, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function is used  to  retrieve  vector  associated  with  out-of-core
request sent by the solver to user code. Depending  on  the  request  type
(returned by the SparseSolverOOCGetRequestInfo()) this  vector  should  be
multiplied by A or subjected to another processing.

It should be used in conjunction with other out-of-core-related  functions
of this subspackage in a loop like one given below:

> alglib.sparsesolveroocstart(state)
> while alglib.sparsesolverooccontinue(state) do
>     alglib.sparsesolveroocgetrequestinfo(state, out RequestType)
>     alglib.sparsesolveroocgetrequestdata(state, out X)
>     if RequestType=0 then
>         [calculate  Y=A*X, with X=R^N]
>     alglib.sparsesolveroocsendresult(state, in Y)
> alglib.sparsesolveroocstop(state, out X, out Report)

INPUT PARAMETERS:
    State           -   solver running in out-of-core mode
    X               -   possibly  preallocated   storage;  reallocated  if
                        needed, left unchanged, if large enough  to  store
                        request data.

OUTPUT PARAMETERS:
    X               -   array[N] or larger, leading N elements are  filled
                        with vector X.


  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolveroocgetrequestdata(const sparsesolverstate &state, real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function is used to retrieve scalar value associated with out-of-core
request sent by the solver to user code. In  the  current  ALGLIB  version
this function is used to retrieve squared residual  norm  during  progress
reports.

INPUT PARAMETERS:
    State           -   solver running in out-of-core mode

OUTPUT PARAMETERS:
    V               -   scalar value associated with the current request


  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolveroocgetrequestdata1(const sparsesolverstate &state, double &v, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function is used to send user reply to out-of-core  request  sent  by
the solver. Usually it is product A*x for vector X returned by the solver.

It should be used in conjunction with other out-of-core-related  functions
of this subspackage in a loop like one given below:

> alglib.sparsesolveroocstart(state)
> while alglib.sparsesolverooccontinue(state) do
>     alglib.sparsesolveroocgetrequestinfo(state, out RequestType)
>     alglib.sparsesolveroocgetrequestdata(state, out X)
>     if RequestType=0 then
>         [calculate  Y=A*X, with X=R^N]
>     alglib.sparsesolveroocsendresult(state, in Y)
> alglib.sparsesolveroocstop(state, out X, out Report)

INPUT PARAMETERS:
    State           -   solver running in out-of-core mode
    AX              -   array[N] or larger, leading N elements contain A*x


  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolveroocsendresult(const sparsesolverstate &state, const real_1d_array &ax, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This  function  finalizes out-of-core mode of the linear solver. It should
be used in conjunction with other out-of-core-related  functions  of  this
subspackage in a loop like one given below:

> alglib.sparsesolveroocstart(state)
> while alglib.sparsesolverooccontinue(state) do
>     alglib.sparsesolveroocgetrequestinfo(state, out RequestType)
>     alglib.sparsesolveroocgetrequestdata(state, out X)
>     if RequestType=0 then
>         [calculate  Y=A*X, with X=R^N]
>     alglib.sparsesolveroocsendresult(state, in Y)
> alglib.sparsesolveroocstop(state, out X, out Report)

INPUT PARAMETERS:
    State       -   solver state

OUTPUT PARAMETERS:
    X       -   array[N], the solution.
                Zero-filled on the failure (Rep.TerminationType<0).
    Rep     -   report with additional info:
                * Rep.TerminationType completion code:
                    * -5    CG method was used for a matrix which  is  not
                            positive definite
                    * -4    overflow/underflow during solution
                            (ill conditioned problem)
                    *  1    ||residual||<=EpsF*||b||
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            best point found is returned
                    *  8    the  algorithm  was  terminated   early  with
                            SparseSolverRequestTermination() being called
                            from other thread.
                * Rep.IterationsCount contains iterations count
                * Rep.NMV contains number of matrix-vector calculations
                * Rep.R2 contains squared residual

  -- ALGLIB --
     Copyright 24.09.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolveroocstop(const sparsesolverstate &state, real_1d_array &x, sparsesolverreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This subroutine submits request for termination of the running solver.  It
can be called from some other thread which wants the   solver to terminate
or when processing an out-of-core request.

As result, solver  stops  at  point  which  was  "current  accepted"  when
the termination request was submitted and returns error code 8 (successful
termination).  Such   termination   is  a smooth  process  which  properly
deallocates all temporaries.

INPUT PARAMETERS:
    State   -   solver structure

NOTE: calling this function on solver which is NOT running  will  have  no
      effect.

NOTE: multiple calls to this function are possible. First call is counted,
      subsequent calls are silently ignored.

NOTE: solver clears termination flag on its start, it means that  if  some
      other thread will request termination too soon, its request will went
      unnoticed.

  -- ALGLIB --
     Copyright 01.10.2021 by Bochkanov Sergey
*************************************************************************/
void sparsesolverrequesttermination(const sparsesolverstate &state, const xparams _xparams = alglib::xdefault);
#endif

#if defined(AE_COMPILE_LINCG) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
This function initializes linear CG Solver. This solver is used  to  solve
symmetric positive definite problems. If you want  to  solve  nonsymmetric
(or non-positive definite) problem you may use LinLSQR solver provided  by
ALGLIB.

USAGE:
1. User initializes algorithm state with LinCGCreate() call
2. User tunes solver parameters with  LinCGSetCond() and other functions
3. Optionally, user sets starting point with LinCGSetStartingPoint()
4. User  calls LinCGSolveSparse() function which takes algorithm state and
   SparseMatrix object.
5. User calls LinCGResults() to get solution
6. Optionally, user may call LinCGSolveSparse()  again  to  solve  another
   problem  with different matrix and/or right part without reinitializing
   LinCGState structure.

INPUT PARAMETERS:
    N       -   problem dimension, N>0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgcreate(const ae_int_t n, lincgstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets starting point.
By default, zero starting point is used.

INPUT PARAMETERS:
    X       -   starting point, array[N]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetstartingpoint(const lincgstate &state, const real_1d_array &x, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
function. By default, SolveSparse() uses diagonal preconditioner,  but  if
you want to use solver without preconditioning, you can call this function
which forces solver to use unit matrix for preconditioning.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void lincgsetprecunit(const lincgstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
preconditioner. This preconditioning mode is active by default.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void lincgsetprecdiag(const lincgstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets stopping criteria.

INPUT PARAMETERS:
    EpsF    -   algorithm will be stopped if norm of residual is less than
                EpsF*||b||.
    MaxIts  -   algorithm will be stopped if number of iterations is  more
                than MaxIts.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
If  both  EpsF  and  MaxIts  are  zero then small EpsF will be set to small
value.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetcond(const lincgstate &state, const double epsf, const ae_int_t maxits, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Procedure for solution of A*x=b with sparse A.

INPUT PARAMETERS:
    State   -   algorithm state
    A       -   sparse matrix in the CRS format (you MUST contvert  it  to
                CRS format by calling SparseConvertToCRS() function).
    IsUpper -   whether upper or lower triangle of A is used:
                * IsUpper=True  => only upper triangle is used and lower
                                   triangle is not referenced at all
                * IsUpper=False => only lower triangle is used and upper
                                   triangle is not referenced at all
    B       -   right part, array[N]

RESULT:
    This function returns no result.
    You can get solution by calling LinCGResults()

NOTE: this function uses lightweight preconditioning -  multiplication  by
      inverse of diag(A). If you want, you can turn preconditioning off by
      calling LinCGSetPrecUnit(). However, preconditioning cost is low and
      preconditioner  is  very  important  for  solution  of  badly scaled
      problems.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsolvesparse(const lincgstate &state, const sparsematrix &a, const bool isupper, const real_1d_array &b, const xparams _xparams = alglib::xdefault);


/*************************************************************************
CG-solver: results.

This function must be called after LinCGSolve

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[N], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -5    input matrix is either not positive definite,
                            too large or too small
                    * -4    overflow/underflow during solution
                            (ill conditioned problem)
                    *  1    ||residual||<=EpsF*||b||
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            best point found is returned
                * Rep.IterationsCount contains iterations count
                * NMV countains number of matrix-vector calculations

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgresults(const lincgstate &state, real_1d_array &x, lincgreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets restart frequency. By default, algorithm  is  restarted
after N subsequent iterations.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetrestartfreq(const lincgstate &state, const ae_int_t srf, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets frequency of residual recalculations.

Algorithm updates residual r_k using iterative formula,  but  recalculates
it from scratch after each 10 iterations. It is done to avoid accumulation
of numerical errors and to stop algorithm when r_k starts to grow.

Such low update frequence (1/10) gives very  little  overhead,  but  makes
algorithm a bit more robust against numerical errors. However, you may
change it

INPUT PARAMETERS:
    Freq    -   desired update frequency, Freq>=0.
                Zero value means that no updates will be done.

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetrupdatefreq(const lincgstate &state, const ae_int_t freq, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinCGOptimize().

  -- ALGLIB --
     Copyright 14.11.2011 by Bochkanov Sergey
*************************************************************************/
void lincgsetxrep(const lincgstate &state, const bool needxrep, const xparams _xparams = alglib::xdefault);
#endif

#if defined(AE_COMPILE_LINLSQR) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
This function initializes linear LSQR Solver. This solver is used to solve
non-symmetric (and, possibly, non-square) problems. Least squares solution
is returned for non-compatible systems.

USAGE:
1. User initializes algorithm state with LinLSQRCreate() call
2. User tunes solver parameters with  LinLSQRSetCond() and other functions
3. User  calls  LinLSQRSolveSparse()  function which takes algorithm state
   and SparseMatrix object.
4. User calls LinLSQRResults() to get solution
5. Optionally, user may call LinLSQRSolveSparse() again to  solve  another
   problem  with different matrix and/or right part without reinitializing
   LinLSQRState structure.

INPUT PARAMETERS:
    M       -   number of rows in A
    N       -   number of variables, N>0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE: see also linlsqrcreatebuf()  for  version  which  reuses  previously
      allocated place as much as possible.

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrcreate(const ae_int_t m, const ae_int_t n, linlsqrstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function initializes linear LSQR Solver.  It  provides  exactly  same
functionality as linlsqrcreate(), but reuses  previously  allocated  space
as much as possible.

INPUT PARAMETERS:
    M       -   number of rows in A
    N       -   number of variables, N>0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 14.11.2018 by Bochkanov Sergey
*************************************************************************/
void linlsqrcreatebuf(const ae_int_t m, const ae_int_t n, const linlsqrstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This  function  changes  preconditioning  settings of LinLSQQSolveSparse()
function. By default, SolveSparse() uses diagonal preconditioner,  but  if
you want to use solver without preconditioning, you can call this function
which forces solver to use unit matrix for preconditioning.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetprecunit(const linlsqrstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This  function  changes  preconditioning  settings  of  LinCGSolveSparse()
function.  LinCGSolveSparse() will use diagonal of the  system  matrix  as
preconditioner. This preconditioning mode is active by default.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 19.11.2012 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetprecdiag(const linlsqrstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets optional Tikhonov regularization coefficient.
It is zero by default.

INPUT PARAMETERS:
    LambdaI -   regularization factor, LambdaI>=0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetlambdai(const linlsqrstate &state, const double lambdai, const xparams _xparams = alglib::xdefault);


/*************************************************************************
Procedure for solution of A*x=b with sparse A.

INPUT PARAMETERS:
    State   -   algorithm state
    A       -   sparse M*N matrix in the CRS format (you MUST contvert  it
                to CRS format  by  calling  SparseConvertToCRS()  function
                BEFORE you pass it to this function).
    B       -   right part, array[M]

RESULT:
    This function returns no result.
    You can get solution by calling LinCGResults()

NOTE: this function uses lightweight preconditioning -  multiplication  by
      inverse of diag(A). If you want, you can turn preconditioning off by
      calling LinLSQRSetPrecUnit(). However, preconditioning cost is   low
      and preconditioner is very important for solution  of  badly  scaled
      problems.

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsolvesparse(const linlsqrstate &state, const sparsematrix &a, const real_1d_array &b, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets stopping criteria.

INPUT PARAMETERS:
    EpsA    -   algorithm will be stopped if ||A^T*Rk||/(||A||*||Rk||)<=EpsA.
    EpsB    -   algorithm will be stopped if ||Rk||<=EpsB*||B||
    MaxIts  -   algorithm will be stopped if number of iterations
                more than MaxIts.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE: if EpsA,EpsB,EpsC and MaxIts are zero then these variables will
be setted as default values.

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetcond(const linlsqrstate &state, const double epsa, const double epsb, const ae_int_t maxits, const xparams _xparams = alglib::xdefault);


/*************************************************************************
LSQR solver: results.

This function must be called after LinLSQRSolve

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[N], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    *  1    ||Rk||<=EpsB*||B||
                    *  4    ||A^T*Rk||/(||A||*||Rk||)<=EpsA
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            X contains best point found so far.
                            (sometimes returned on singular systems)
                    *  8    user requested termination via calling
                            linlsqrrequesttermination()
                * Rep.IterationsCount contains iterations count
                * NMV countains number of matrix-vector calculations

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrresults(const linlsqrstate &state, real_1d_array &x, linlsqrreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinCGOptimize().

  -- ALGLIB --
     Copyright 30.11.2011 by Bochkanov Sergey
*************************************************************************/
void linlsqrsetxrep(const linlsqrstate &state, const bool needxrep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function is used to peek into LSQR solver and get  current  iteration
counter. You can safely "peek" into the solver from another thread.

INPUT PARAMETERS:
    S           -   solver object

RESULT:
    iteration counter, in [0,INF)

  -- ALGLIB --
     Copyright 21.05.2018 by Bochkanov Sergey
*************************************************************************/
ae_int_t linlsqrpeekiterationscount(const linlsqrstate &s, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This subroutine submits request for termination of the running solver.  It
can be called from some other thread which wants LSQR solver to  terminate
(obviously, the  thread  running  LSQR  solver can not request termination
because it is already busy working on LSQR).

As result, solver  stops  at  point  which  was  "current  accepted"  when
termination  request  was  submitted  and returns error code 8 (successful
termination).  Such   termination   is  a smooth  process  which  properly
deallocates all temporaries.

INPUT PARAMETERS:
    State   -   solver structure

NOTE: calling this function on solver which is NOT running  will  have  no
      effect.

NOTE: multiple calls to this function are possible. First call is counted,
      subsequent calls are silently ignored.

NOTE: solver clears termination flag on its start, it means that  if  some
      other thread will request termination too soon, its request will went
      unnoticed.

  -- ALGLIB --
     Copyright 08.10.2014 by Bochkanov Sergey
*************************************************************************/
void linlsqrrequesttermination(const linlsqrstate &state, const xparams _xparams = alglib::xdefault);
#endif

#if defined(AE_COMPILE_NLEQ) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************
                LEVENBERG-MARQUARDT-LIKE NONLINEAR SOLVER

DESCRIPTION:
This algorithm solves system of nonlinear equations
    F[0](x[0], ..., x[n-1])   = 0
    F[1](x[0], ..., x[n-1])   = 0
    ...
    F[M-1](x[0], ..., x[n-1]) = 0
with M/N do not necessarily coincide.  Algorithm  converges  quadratically
under following conditions:
    * the solution set XS is nonempty
    * for some xs in XS there exist such neighbourhood N(xs) that:
      * vector function F(x) and its Jacobian J(x) are continuously
        differentiable on N
      * ||F(x)|| provides local error bound on N, i.e. there  exists  such
        c1, that ||F(x)||>c1*distance(x,XS)
Note that these conditions are much more weaker than usual non-singularity
conditions. For example, algorithm will converge for any  affine  function
F (whether its Jacobian singular or not).


REQUIREMENTS:
Algorithm will request following information during its operation:
* function vector F[] and Jacobian matrix at given point X
* value of merit function f(x)=F[0]^2(x)+...+F[M-1]^2(x) at given point X


USAGE:
1. User initializes algorithm state with NLEQCreateLM() call
2. User tunes solver parameters with  NLEQSetCond(),  NLEQSetStpMax()  and
   other functions
3. User  calls  NLEQSolve()  function  which  takes  algorithm  state  and
   pointers (delegates, etc.) to callback functions which calculate  merit
   function value and Jacobian.
4. User calls NLEQResults() to get solution
5. Optionally, user may call NLEQRestartFrom() to  solve  another  problem
   with same parameters (N/M) but another starting  point  and/or  another
   function vector. NLEQRestartFrom() allows to reuse already  initialized
   structure.


INPUT PARAMETERS:
    N       -   space dimension, N>1:
                * if provided, only leading N elements of X are used
                * if not provided, determined automatically from size of X
    M       -   system size
    X       -   starting point


OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state


NOTES:
1. you may tune stopping conditions with NLEQSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use NLEQSetStpMax() function to bound algorithm's steps.
3. this  algorithm  is  a  slightly  modified implementation of the method
   described  in  'Levenberg-Marquardt  method  for constrained  nonlinear
   equations with strong local convergence properties' by Christian Kanzow
   Nobuo Yamashita and Masao Fukushima and further  developed  in  'On the
   convergence of a New Levenberg-Marquardt Method'  by  Jin-yan  Fan  and
   Ya-Xiang Yuan.


  -- ALGLIB --
     Copyright 20.08.2009 by Bochkanov Sergey
*************************************************************************/
void nleqcreatelm(const ae_int_t n, const ae_int_t m, const real_1d_array &x, nleqstate &state, const xparams _xparams = alglib::xdefault);
void nleqcreatelm(const ae_int_t m, const real_1d_array &x, nleqstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets stopping conditions for the nonlinear solver

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsF    -   >=0
                The subroutine finishes  its work if on k+1-th iteration
                the condition ||F||<=EpsF is satisfied
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsF=0 and MaxIts=0 simultaneously will lead to  automatic
stopping criterion selection (small EpsF).

NOTES:

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void nleqsetcond(const nleqstate &state, const double epsf, const ae_int_t maxits, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to NLEQSolve().

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void nleqsetxrep(const nleqstate &state, const bool needxrep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when target function  contains  exp()  or  other  fast
growing functions, and algorithm makes  too  large  steps  which  lead  to
overflow. This function allows us to reject steps that are too large  (and
therefore expose us to the possible overflow) without actually calculating
function value at the x+stp*d.

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void nleqsetstpmax(const nleqstate &state, const double stpmax, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool nleqiteration(const nleqstate &state, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This family of functions is used to launcn iterations of nonlinear solver

These functions accept following parameters:
    state   -   algorithm state
    func    -   callback which calculates function (or merit function)
                value func at given point x
    jac     -   callback which calculates function vector fi[]
                and Jacobian jac at given point x
    rep     -   optional callback which is called after each iteration
                can be NULL
    ptr     -   optional pointer which is passed to func/grad/hess/jac/rep
                can be NULL


  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey

*************************************************************************/
void nleqsolve(nleqstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr) = NULL,
    void *ptr = NULL,
    const xparams _xparams = alglib::xdefault);


/*************************************************************************
NLEQ solver results

INPUT PARAMETERS:
    State   -   algorithm state.

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -4    ERROR:  algorithm   has   converged   to   the
                            stationary point Xf which is local minimum  of
                            f=F[0]^2+...+F[m-1]^2, but is not solution  of
                            nonlinear system.
                    *  1    sqrt(f)<=EpsF.
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations
                * ActiveConstraints contains number of active constraints

  -- ALGLIB --
     Copyright 20.08.2009 by Bochkanov Sergey
*************************************************************************/
void nleqresults(const nleqstate &state, real_1d_array &x, nleqreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
NLEQ solver results

Buffered implementation of NLEQResults(), which uses pre-allocated  buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 20.08.2009 by Bochkanov Sergey
*************************************************************************/
void nleqresultsbuf(const nleqstate &state, real_1d_array &x, nleqreport &rep, const xparams _xparams = alglib::xdefault);


/*************************************************************************
This  subroutine  restarts  CG  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used for reverse communication previously
                allocated with MinCGCreate call.
    X       -   new starting point.
    BndL    -   new lower bounds
    BndU    -   new upper bounds

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void nleqrestartfrom(const nleqstate &state, const real_1d_array &x, const xparams _xparams = alglib::xdefault);
#endif
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (FUNCTIONS)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
#if defined(AE_COMPILE_POLYNOMIALSOLVER) || !defined(AE_PARTIAL_BUILD)
void polynomialsolve(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Complex */ ae_vector* x,
     polynomialsolverreport* rep,
     ae_state *_state);
void _polynomialsolverreport_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _polynomialsolverreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _polynomialsolverreport_clear(void* _p);
void _polynomialsolverreport_destroy(void* _p);
#endif
#if defined(AE_COMPILE_DIRECTDENSESOLVERS) || !defined(AE_PARTIAL_BUILD)
void rmatrixsolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void rmatrixsolvefast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void rmatrixsolvem(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_bool rfs,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void rmatrixsolvemfast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void rmatrixlusolve(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void rmatrixlusolvefast(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void rmatrixlusolvem(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void rmatrixlusolvemfast(/* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void rmatrixmixedsolve(/* Real    */ ae_matrix* a,
     /* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void rmatrixmixedsolvem(/* Real    */ ae_matrix* a,
     /* Real    */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void cmatrixsolvem(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_bool rfs,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void cmatrixsolvemfast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void cmatrixsolve(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void cmatrixsolvefast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void cmatrixlusolvem(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void cmatrixlusolvemfast(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void cmatrixlusolve(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void cmatrixlusolvefast(/* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void cmatrixmixedsolvem(/* Complex */ ae_matrix* a,
     /* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void cmatrixmixedsolve(/* Complex */ ae_matrix* a,
     /* Complex */ ae_matrix* lua,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void spdmatrixsolvem(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void spdmatrixsolvemfast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void spdmatrixsolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void spdmatrixsolvefast(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void spdmatrixcholeskysolvem(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_matrix* x,
     ae_state *_state);
void spdmatrixcholeskysolvemfast(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void spdmatrixcholeskysolve(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void spdmatrixcholeskysolvefast(/* Real    */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void hpdmatrixsolvem(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void hpdmatrixsolvemfast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void hpdmatrixsolve(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void hpdmatrixsolvefast(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void hpdmatrixcholeskysolvem(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_matrix* x,
     ae_state *_state);
void hpdmatrixcholeskysolvemfast(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_matrix* b,
     ae_int_t m,
     ae_int_t* info,
     ae_state *_state);
void hpdmatrixcholeskysolve(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     densesolverreport* rep,
     /* Complex */ ae_vector* x,
     ae_state *_state);
void hpdmatrixcholeskysolvefast(/* Complex */ ae_matrix* cha,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* b,
     ae_int_t* info,
     ae_state *_state);
void rmatrixsolvels(/* Real    */ ae_matrix* a,
     ae_int_t nrows,
     ae_int_t ncols,
     /* Real    */ ae_vector* b,
     double threshold,
     ae_int_t* info,
     densesolverlsreport* rep,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void _densesolverreport_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _densesolverreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _densesolverreport_clear(void* _p);
void _densesolverreport_destroy(void* _p);
void _densesolverlsreport_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _densesolverlsreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _densesolverlsreport_clear(void* _p);
void _densesolverlsreport_destroy(void* _p);
#endif
#if defined(AE_COMPILE_DIRECTSPARSESOLVERS) || !defined(AE_PARTIAL_BUILD)
void sparsespdsolvesks(sparsematrix* a,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparsespdsolve(sparsematrix* a,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparsespdcholeskysolve(sparsematrix* a,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparsesolve(sparsematrix* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparselusolve(sparsematrix* a,
     /* Integer */ ae_vector* p,
     /* Integer */ ae_vector* q,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void initsparsesolverreport(sparsesolverreport* rep, ae_state *_state);
void _sparsesolverreport_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _sparsesolverreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _sparsesolverreport_clear(void* _p);
void _sparsesolverreport_destroy(void* _p);
#endif
#if defined(AE_COMPILE_ITERATIVESPARSE) || !defined(AE_PARTIAL_BUILD)
void sparsesolvesymmetricgmres(sparsematrix* a,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_int_t k,
     double epsf,
     ae_int_t maxits,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparsesolvegmres(sparsematrix* a,
     /* Real    */ ae_vector* b,
     ae_int_t k,
     double epsf,
     ae_int_t maxits,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparsesolvercreate(ae_int_t n,
     sparsesolverstate* state,
     ae_state *_state);
void sparsesolversetalgogmres(sparsesolverstate* state,
     ae_int_t k,
     ae_state *_state);
void sparsesolversetstartingpoint(sparsesolverstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void sparsesolversetcond(sparsesolverstate* state,
     double epsf,
     ae_int_t maxits,
     ae_state *_state);
void sparsesolversolvesymmetric(sparsesolverstate* state,
     sparsematrix* a,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void sparsesolversolve(sparsesolverstate* state,
     sparsematrix* a,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void sparsesolverresults(sparsesolverstate* state,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparsesolversetxrep(sparsesolverstate* state,
     ae_bool needxrep,
     ae_state *_state);
void sparsesolveroocstart(sparsesolverstate* state,
     /* Real    */ ae_vector* b,
     ae_state *_state);
ae_bool sparsesolverooccontinue(sparsesolverstate* state,
     ae_state *_state);
void sparsesolveroocgetrequestinfo(sparsesolverstate* state,
     ae_int_t* requesttype,
     ae_state *_state);
void sparsesolveroocgetrequestdata(sparsesolverstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void sparsesolveroocgetrequestdata1(sparsesolverstate* state,
     double* v,
     ae_state *_state);
void sparsesolveroocsendresult(sparsesolverstate* state,
     /* Real    */ ae_vector* ax,
     ae_state *_state);
void sparsesolveroocstop(sparsesolverstate* state,
     /* Real    */ ae_vector* x,
     sparsesolverreport* rep,
     ae_state *_state);
void sparsesolverrequesttermination(sparsesolverstate* state,
     ae_state *_state);
void _sparsesolverstate_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _sparsesolverstate_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _sparsesolverstate_clear(void* _p);
void _sparsesolverstate_destroy(void* _p);
#endif
#if defined(AE_COMPILE_LINCG) || !defined(AE_PARTIAL_BUILD)
void lincgcreate(ae_int_t n, lincgstate* state, ae_state *_state);
void lincgsetstartingpoint(lincgstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void lincgsetb(lincgstate* state,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void lincgsetprecunit(lincgstate* state, ae_state *_state);
void lincgsetprecdiag(lincgstate* state, ae_state *_state);
void lincgsetcond(lincgstate* state,
     double epsf,
     ae_int_t maxits,
     ae_state *_state);
ae_bool lincgiteration(lincgstate* state, ae_state *_state);
void lincgsolvesparse(lincgstate* state,
     sparsematrix* a,
     ae_bool isupper,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void lincgresults(lincgstate* state,
     /* Real    */ ae_vector* x,
     lincgreport* rep,
     ae_state *_state);
void lincgsetrestartfreq(lincgstate* state,
     ae_int_t srf,
     ae_state *_state);
void lincgsetrupdatefreq(lincgstate* state,
     ae_int_t freq,
     ae_state *_state);
void lincgsetxrep(lincgstate* state, ae_bool needxrep, ae_state *_state);
void lincgrestart(lincgstate* state, ae_state *_state);
void _lincgstate_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _lincgstate_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _lincgstate_clear(void* _p);
void _lincgstate_destroy(void* _p);
void _lincgreport_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _lincgreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _lincgreport_clear(void* _p);
void _lincgreport_destroy(void* _p);
#endif
#if defined(AE_COMPILE_LINLSQR) || !defined(AE_PARTIAL_BUILD)
void linlsqrcreate(ae_int_t m,
     ae_int_t n,
     linlsqrstate* state,
     ae_state *_state);
void linlsqrcreatebuf(ae_int_t m,
     ae_int_t n,
     linlsqrstate* state,
     ae_state *_state);
void linlsqrsetb(linlsqrstate* state,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void linlsqrsetprecunit(linlsqrstate* state, ae_state *_state);
void linlsqrsetprecdiag(linlsqrstate* state, ae_state *_state);
void linlsqrsetlambdai(linlsqrstate* state,
     double lambdai,
     ae_state *_state);
ae_bool linlsqriteration(linlsqrstate* state, ae_state *_state);
void linlsqrsolvesparse(linlsqrstate* state,
     sparsematrix* a,
     /* Real    */ ae_vector* b,
     ae_state *_state);
void linlsqrsetcond(linlsqrstate* state,
     double epsa,
     double epsb,
     ae_int_t maxits,
     ae_state *_state);
void linlsqrresults(linlsqrstate* state,
     /* Real    */ ae_vector* x,
     linlsqrreport* rep,
     ae_state *_state);
void linlsqrsetxrep(linlsqrstate* state,
     ae_bool needxrep,
     ae_state *_state);
void linlsqrrestart(linlsqrstate* state, ae_state *_state);
ae_int_t linlsqrpeekiterationscount(linlsqrstate* s, ae_state *_state);
void linlsqrrequesttermination(linlsqrstate* state, ae_state *_state);
void _linlsqrstate_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _linlsqrstate_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _linlsqrstate_clear(void* _p);
void _linlsqrstate_destroy(void* _p);
void _linlsqrreport_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _linlsqrreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _linlsqrreport_clear(void* _p);
void _linlsqrreport_destroy(void* _p);
#endif
#if defined(AE_COMPILE_NLEQ) || !defined(AE_PARTIAL_BUILD)
void nleqcreatelm(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     nleqstate* state,
     ae_state *_state);
void nleqsetcond(nleqstate* state,
     double epsf,
     ae_int_t maxits,
     ae_state *_state);
void nleqsetxrep(nleqstate* state, ae_bool needxrep, ae_state *_state);
void nleqsetstpmax(nleqstate* state, double stpmax, ae_state *_state);
ae_bool nleqiteration(nleqstate* state, ae_state *_state);
void nleqresults(nleqstate* state,
     /* Real    */ ae_vector* x,
     nleqreport* rep,
     ae_state *_state);
void nleqresultsbuf(nleqstate* state,
     /* Real    */ ae_vector* x,
     nleqreport* rep,
     ae_state *_state);
void nleqrestartfrom(nleqstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void _nleqstate_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _nleqstate_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _nleqstate_clear(void* _p);
void _nleqstate_destroy(void* _p);
void _nleqreport_init(void* _p, ae_state *_state, ae_bool make_automatic);
void _nleqreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic);
void _nleqreport_clear(void* _p);
void _nleqreport_destroy(void* _p);
#endif

}
#endif

