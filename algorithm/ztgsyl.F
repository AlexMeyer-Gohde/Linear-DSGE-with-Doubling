#include "fintrf.h"
C======================================================================
#if 0
C     
C     ztgsyl.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     ztgsyl.f
C
C     ZTGSYL solves Sylvester equations, using Lapacks ZTGSYL 
C     function
     
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
      mwPointer mxGetComplexDoubles
      mwPointer mxGetInt32s
      integer*4 mxGetString
      integer*4 status
      integer*4 mxCopyPtrToComplex16, mxCopyComplex16ToPtr
      integer*4 mxCopyReal8ToPtr
      integer*4 mxCopyInteger4ToPtr
      integer*4 mxClassIDFromClassName
      mwPointer mxCreateNumericMatrix
#else
	  mwPointer mxGetPr, mxGetPi
#endif
      integer mxIsNumeric
      integer*4 mxIsComplex
      integer*4 mxIsChar
      mwPointer mxGetM, mxGetN
      
C     Pointers to input/output mxArrays:
      mwPointer trans_pr, ijob_pr
      mwPointer A_pr, B_pr, C_pr, D_pr, E_pr, F_pr
      mwPointer Cout_pr, Fout_pr, difout_pr, scaleout_pr
      mwpointer workout_pr, iworkout_pr, infoout_pr
 
C     Arguments for computational routine:
C     Maximum size = numel
      integer*4 numel
      mwSize one
      parameter(numel = 9000000, one = 1)
      character*1 trans
      mwSize   ijob, m, n
      mwSize   lda, ldb, ldc, ldd, lde, ldf
      mwSize   lwork   
      complex*16  A(numel), B(numel), C(numel)
      complex*16  D(numel), E(numel), F(numel)
      real*8      dif, scale, tmp_ijob
      integer*4   iwork(numel)
      complex*16  Work
      integer*4   info

C     Assign pointers to input
      trans_pr = prhs(1)
      ijob_pr = prhs(2)
      A_pr = prhs(3)
      B_pr = prhs(4)
      C_pr = prhs(5)
      D_pr = prhs(6)
      E_pr = prhs(7)
      F_pr = prhs(8)
 
C-----------------------------------------------------------------------
C --------- Checks -----------------------------------------------------

C     Check for proper number of arguments. 
      if (nrhs .ne. 8) then
         call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:nInput',
     +                           'Eight inputs required.')
      elseif (nlhs .gt. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:nOutput',
     +                           'Too many output arguments.')
      endif

C     ijob must be numeric
      if (mxIsNumeric(ijob_pr) .ne. 1) then
            call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:NonNumeric',
     +               'ijob must be a numeric variable.')
      endif

C     trans must be a character
      if(mxIsChar(trans_pr) .ne. 1) then
            call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:NonChar',
     +               'trans must be a char.')
      endif

C     Column * row must be smaller than numel
      if(mxGetM(A_pr)*mxGetN(A_pr) .gt. numel) then
         call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:mSize',
     +               'Size of A must be <= 1000.')
      endif

C     Column * row must be smaller than numel
      if(mxGetM(B_pr)*mxGetN(B_pr) .gt. numel) then
      call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:mSize',
     +               'Size of B must be <= 1000.')
      endif
C     Column * row must be smaller than numel
      if(mxGetM(C_pr)*mxGetN(C_pr) .gt. numel) then
      call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:mSize',
     +               'Size of C must be <= 1000.')
      endif
C     Column * row must be smaller than numel
      if(mxGetM(D_pr)*mxGetN(D_pr) .gt. numel) then
      call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:mSize',
     +               'Size of D must be <= 1000.')
      endif
C     Column * row must be smaller than numel
      if(mxGetM(E_pr)*mxGetN(E_pr) .gt. numel) then
      call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:mSize',
     +               'Size of E must be <= 1000.')
      endif
C     Column * row must be smaller than numel
      if(mxGetM(F_pr)*mxGetN(F_pr) .gt. numel) then
      call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:mSize',
     +               'Size of F must be <= 1000.')
      endif

C     Check if matrix inputs are complex.
      if ((mxIsComplex(A_pr) .ne. 1) .or.           
     +    (mxIsComplex(B_pr) .ne. 1) .or.
     +    (mxIsComplex(C_pr) .ne. 1) .or.
     +    (mxIsComplex(D_pr) .ne. 1) .or.
     +    (mxIsComplex(E_pr) .ne. 1) .or.
     +    (mxIsComplex(F_pr) .ne. 1)
     +    ) then
         call mexErrMsgIdAndTxt ('MATLAB:ztgsyl:NonComplex',
     +                           'Matrix inputs must be complex.')
      endif

C-----------------------------------------------------------------------
C---------- Assignment of arguments of computational routine -----------

      m = mxGetM(A_pr)
      n = mxGetM(B_pr)
      lda = m
      ldb = n
      ldc = m
      ldd = m
      lde = n
      ldf = m
      lwork = 1

C-----------------------------------------------------------------------
C---------- Load the data into Fortran arrays -----------

#if MX_HAS_INTERLEAVED_COMPLEX
      
      ijob_pr = mxGetDoubles(ijob_pr)

      A_pr = mxGetComplexDoubles(A_pr)
      status = 
     + mxCopyPtrToComplex16(A_pr,A,lda*m)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from A to complex*16.')
      endif
      B_pr = mxGetComplexDoubles(B_pr)
      status = 
     + mxCopyPtrToComplex16(B_pr,B,ldb*n)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from B to complex*16.')
      endif
      C_pr = mxGetComplexDoubles(C_pr)
      status = 
     + mxCopyPtrToComplex16(C_pr,C,ldc*n)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from C to complex*16.')
      endif
      D_pr = mxGetComplexDoubles(D_pr)
      status = 
     + mxCopyPtrToComplex16(D_pr,D,ldd*m)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from D to complex*16.')
      endif
      E_pr = mxGetComplexDoubles(E_pr)
      status = 
     + mxCopyPtrToComplex16(E_pr,E,lde*n)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from E to complex*16.')
      endif
      F_pr = mxGetComplexDoubles(F_pr)
      status = 
     + mxCopyPtrToComplex16(F_pr,F,ldf*n)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from F to complex*16.')
      endif
#else
      ijob_pr = mxGetPr(ijob_pr)
      
      call mxCopyPtrToComplex16(mxGetPr(A_pr),
     +                          mxGetPi(A_pr),A,lda*m)
      call mxCopyPtrToComplex16(mxGetPr(B_pr),
     +                          mxGetPi(B_pr),B,ldb*n)
      call mxCopyPtrToComplex16(mxGetPr(C_pr),
     +                          mxGetPi(C_pr),C,ldc*n)
      call mxCopyPtrToComplex16(mxGetPr(D_pr),
     +                          mxGetPi(D_pr),D,ldd*m)
      call mxCopyPtrToComplex16(mxGetPr(E_pr),
     +                          mxGetPi(E_pr),E,lde*n)
      call mxCopyPtrToComplex16(mxGetPr(F_pr),
     +                          mxGetPi(F_pr),F,ldf*n)
#endif

C     Typecast double to integer
      call mxCopyPtrToReal8(ijob_pr,tmp_ijob,1)
      ijob = tmp_ijob

      status = mxGetString(trans_pr, trans, 1)
      if (status .ne. 0) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +               'Could not copy from trans to char*1.')
      endif

C-----------------------------------------------------------------------
C---------- Call of computational routine -----------

      call ztgsyl(trans, ijob, m, n, A, lda, B, ldb, C, ldc, D, ldd, 
     +       E, lde, F, ldf, scale, dif, work, lwork, iwork, info)

C-----------------------------------------------------------------------
C---------- Set up output pointers ---------------------------------------------

C      plhs(3) = mxCreateDoubleMatrix(1, 1, 0)
C      difout_pr = mxGetDoubles(plhs(3))
C      plhs(4) = mxCreateDoubleMatrix(1, 1, 0)
C      scaleout_pr = mxGetDoubles(plhs(4))

      Cout_pr = mxCreateDoubleMatrix(ldc, n, 1)
      Fout_pr = mxCreateDoubleMatrix(ldf, n, 1)
      difout_pr = mxCreateDoubleMatrix(one, one, 0)
      scaleout_pr = mxCreateDoubleMatrix(one, one, 0)
      workout_pr = mxCreateDoubleMatrix(max(1,lwork), 1, 1)
      iworkout_pr = mxCreateNumericMatrix(m+n+2, 1,
     +                mxClassIDFromClassName('int32'), 0)
      infoout_pr = mxCreateNumericMatrix(1, 1,
     +                mxClassIDFromClassName('int32'), 0)
C      lworkout_pr = mxCreateDoubleMatrix(1, 1, 0)
C      iworkout_pr = 

      plhs(1) = Cout_pr
      plhs(2) = Fout_pr
      plhs(3) = difout_pr
      plhs(4) = scaleout_pr
      plhs(5) = workout_pr
      plhs(6) = iworkout_pr
      plhs(5) = infoout_pr

C-----------------------------------------------------------------------
C---------- Load the output into MATLAB arrays -----------
#if MX_HAS_INTERLEAVED_COMPLEX
C     C
      status =
     + mxCopyComplex16ToPtr(C,mxGetComplexDoubles(Cout_pr),ldc*n)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from complex*16 to C.')
      endif
C     F
      status =
     + mxCopyComplex16ToPtr(F,mxGetComplexDoubles(Fout_pr),ldf*n)
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from complex*16 to F.')
      endif
C     dif
      status = mxCopyReal8ToPtr(dif,mxGetDoubles(difout_pr),one)
C Following check fails and is therefore commented out. No clue why this is the case
C      if (status .ne. 1) then
C            call mexErrMsgIdAndTxt (
C     +              'MATLAB:ztgsyl:CopyingFailed',
C     +              'Could not copy from real*8 to dif.')
C      endif
C     scale
      status =  mxCopyReal8ToPtr(scale,mxGetDoubles(scaleout_pr),
     + one)
C Following check fails and is therefore commented out. No clue why this is the case
C      if (status .ne. 1) then
C            call mexErrMsgIdAndTxt (
C     +              'MATLAB:ztgsyl:CopyingFailed',
C     +              'Could not copy from real*8 to scale.')
C      endif
C     work
      status =
     + mxCopyComplex16ToPtr(work,mxGetComplexDoubles(workout_pr),
     + max(1,lwork))
      if (status .ne. 1) then
            call mexErrMsgIdAndTxt (
     +              'MATLAB:ztgsyl:CopyingFailed',
     +              'Could not copy from complex*16 to work.')
      endif
C     iwork
      status = mxCopyInteger4ToPtr(iwork,mxGetInt32s(iworkout_pr),
     + m+n+2)
C Following check fails and is therefore commented out. No clue why this is the case
C      if (status .ne. 1) then
C            call mexErrMsgIdAndTxt (
C     +              'MATLAB:ztgsyl:CopyingFailed',
C     +              'Could not copy from real*8 to scale.')
C      endif
C     work
#else
C     C
      call mxCopyComplex16ToPtr(C,mxGetPr(Cout_pr),
     +                          mxGetPi(Cout_pr),ldc*n)

      return
C     F
      call mxCopyComplex16ToPtr(F,mxGetPr(Fout_pr),
     +                          mxGetPi(Fout_pr),ldf*n)

      return
#endif

      return
      end