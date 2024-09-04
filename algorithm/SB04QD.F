#include "fintrf.h"
C======================================================================
#if 0
C     
C     SB04QD.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     SB04QD.f
C
C     ZTGSYL solves Sylvester equations, using Slicot's SB04QD 
C     routine
     
C     This is a MEX-file for MATLAB.
C     
C======================================================================

C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
 
C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
 
C     Function declarations:
      mwPointer mxCreateDoubleMatrix
#if MX_HAS_INTERLEAVED_COMPLEX
	mwPointer mxGetDoubles
      mwPointer mxGetInt32s
      integer*4 mxClassIDFromClassName
      mwPointer mxCreateNumericMatrix
#else
	  mwPointer mxGetPr, mxGetPi
#endif
      integer mxIsNumeric
      integer*4 mxIsDouble
      mwPointer mxGetM, mxGetN
      
C     Pointers to input/output mxArrays:
      mwPointer N_pr, M_pr
      mwPointer A_pr, B_pr, C_pr
      mwPointer Cout_pr, Zout_pr, infoout_pr 

C     Arguments for computational routine:
C     Maximum size = numel
      integer*4 one, numel, nmax, mmax
      parameter(one = 1, numel = 10000, nmax = 50, mmax = 180)
      mwSize n, m
      mwSize lda, ldb, ldc, ldz
C      real*8 A(numel), B(numel), C(numel), Z(numel)
C      real*8 A(20,20), B(20,20), C(20,20), Z(20,20)
      real*8, allocatable :: A(:,:), B(:,:), C(:,:), Z(:,:)
C      real*8 A_coeff(3,3), B_coeff(3,3), C_coeff(3,3)
      real*8, allocatable :: A_coeff(:,:), B_coeff(:,:),
     +                       C_coeff(:,:)
C      integer*4   Iwork(numel)
      integer*4, allocatable :: Iwork(:)
C      integer*4   Iwork(80)
C      real*8      dwork(numel)      
C      real*8      dwork(980)      
      real*8, allocatable :: Dwork(:)      
      integer*4   ldwork, info
      real*8      tmp_real_to_int


C     Assign pointers to input
      N_pr = prhs(1)
      M_pr = prhs(2)
      A_pr = prhs(3)
      B_pr = prhs(4)
      C_pr = prhs(5)
 
C-----------------------------------------------------------------------
C --------- Checks -----------------------------------------------------

C     Check for proper number of arguments. 
      if (nrhs .ne. 5) then
         call mexErrMsgIdAndTxt ('MATLAB:SB04QD:nInput',
     +                           'Five inputs required.')
      elseif (nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:SB04QD:nOutput',
     +                           'Too many output arguments.')
      endif

C     N must be numeric
      if (mxIsNumeric(N_pr) .ne. 1) then
            call mexErrMsgIdAndTxt ('MATLAB:SB04QD:NonNumeric',
     +               'N must be a numeric variable.')
      endif
C     M must be numeric
      if (mxIsNumeric(M_pr) .ne. 1) then
            call mexErrMsgIdAndTxt ('MATLAB:SB04QD:NonNumeric',
     +               'M must be a numeric variable.')
      endif

C     Column * row must be smaller than numel
      if(mxGetM(A_pr)*mxGetN(A_pr) .gt. numel) then
         call mexErrMsgIdAndTxt ('MATLAB:SB04QD:mSize',
     +               'Size of A must be <= 1000.')
      endif
C     Column * row must be smaller than numel
      if(mxGetM(B_pr)*mxGetN(B_pr) .gt. numel) then
      call mexErrMsgIdAndTxt ('MATLAB:SB04QD:mSize',
     +               'Size of B must be <= 1000.')
      endif
C     Column * row must be smaller than numel
      if(mxGetM(C_pr)*mxGetN(C_pr) .gt. numel) then
      call mexErrMsgIdAndTxt ('MATLAB:SB04QD:mSize',
     +               'Size of C must be <= 1000.')
      endif

C     Check if matrix inputs are doubles.
      if ((mxIsDouble(A_pr) .ne. 1) .or.           
     +    (mxIsDouble(B_pr) .ne. 1) .or.
     +    (mxIsDouble(C_pr) .ne. 1)
     +    ) then
         call mexErrMsgIdAndTxt ('MATLAB:SB04QD:NonDouble',
     +                           'Matrix inputs must be double.')
      endif

C-----------------------------------------------------------------------
C---------- Assignment of arguments of computational routine -----------

#if MX_HAS_INTERLEAVED_COMPLEX
      N_pr = mxGetDoubles(N_pr)
      M_pr = mxGetDoubles(M_pr)
#else
      N_pr = mxGetPr(N_pr)
      M_pr = mxGetPr(M_pr)
#endif

C     Typecast double to integer
      call mxCopyPtrToReal8(N_pr,tmp_real_to_int,1)
      n = tmp_real_to_int
      call mxCopyPtrToReal8(M_pr,tmp_real_to_int,1)
      m = tmp_real_to_int
      
      lda = nmax
      ldb = mmax
      ldc = nmax
      ldz = mmax
      ldwork = max( 1, 2*nmax*nmax+9*nmax, 5*mmax, nmax+mmax )

C-----------------------------------------------------------------------
C---------- Load the data into Fortran arrays -----------

      allocate(A(lda,n))
      allocate(B(ldb,m))
      allocate(C(ldc,m))
      allocate(Z(ldz,m))
      allocate(Iwork(4*n))
      allocate(Dwork(ldwork))

      allocate(A_coeff(n,n))
      allocate(B_coeff(m,m))
      allocate(C_coeff(n,m))

      A(:,:) = 0
      B(:,:) = 0
      C(:,:) = 0
      Z(:,:) = 0
      DWork(:) = 0
      IWork(:) = 0
      info = 0

#if MX_HAS_INTERLEAVED_COMPLEX    
      A_pr = mxGetDoubles(A_pr)
      call mxCopyPtrToReal8(A_pr,A_coeff,n*n)
      B_pr = mxGetDoubles(B_pr)
      call mxCopyPtrToReal8(B_pr,B_coeff,m*m)
      C_pr = mxGetDoubles(C_pr)
      call mxCopyPtrToReal8(C_pr,C_coeff,n*m)
      
#else      
      call mxCopyPtrToReal8(mxGetPr(A_pr),
     +                          A_coeff,lda*n)
      call mxCopyPtrToReal8(mxGetPr(B_pr),
     +                          B_coeff,ldb*m)
      call mxCopyPtrToReal8(mxGetPr(C_pr),
     +                          C_coeff,ldc*m)
#endif

      A(1:n,1:n) = A_coeff
      B(1:m,1:m) = B_coeff
      C(1:n,1:m) = C_coeff
      
C-----------------------------------------------------------------------
C---------- Call of computational routine -----------
      call SB04QD( n, m, 
     +   A, lda,
     +   B, ldb, 
     +   C, ldc, 
     +   Z, ldz, 
     +   Iwork,
     +   Dwork, ldwork, info )

C-----------------------------------------------------------------------
C---------- Set up output pointers ---------------------------------------------

      Cout_pr = mxCreateDoubleMatrix(ldc, m, 0)
      Zout_pr = mxCreateDoubleMatrix(ldz, m, 0)
      infoout_pr = mxCreateNumericMatrix(1, 1,
     +                mxClassIDFromClassName('int32'), 0)

      plhs(1) = Cout_pr
      plhs(2) = Zout_pr
      plhs(3) = infoout_pr

C-----------------------------------------------------------------------
C---------- Load the output into MATLAB arrays -----------
#if MX_HAS_INTERLEAVED_COMPLEX
C     X / C_out
      call mxCopyReal8ToPtr(C,mxGetDoubles(Cout_pr),ldc*m)
C     Z
      call mxCopyReal8ToPtr(Z,mxGetDoubles(Zout_pr),ldz*m)
C     info
      call mxCopyInteger4ToPtr(info,mxGetInt32s(infoout_pr),
     +                         one)
#else
C     X / C_out
      call mxCopyReal8ToPtr(C,mxGetPr(Cout_pr),
     +                          ldc*m)
C     Z
      call mxCopyReal8ToPtr(Z,mxGetPr(Zout_pr),
     +                          ldz*m)
C     info
      call mxCopyInteger4ToPtr(info,mxGetPr(infoout_pr),
     +                         one)
#endif

      return
      end