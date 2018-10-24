#include "fintrf.h"

c
c     MATLAB interface for calling a QR with cost routine
c
c     Usage
c
c     [AREF,IPIV,TAU] = xqrmc_m(A,c,k)
c     
c     LAPACK routines are available under BSD license (below)
c     Authors: Univ. of Tennessee, Univ. of California Berkeley,
c     Univ. of Colorado Denver, NAG Ltd.
c
c
c Copyright (c) 1992-2013 The University of Tennessee and The University
c                         of Tennessee Research Foundation.  All rights
c                         reserved.
c Copyright (c) 2000-2013 The University of California Berkeley. All
c                         rights reserved.
c Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
c                         reserved.
c
c     Modifications are also available under BSD license (below)
c     Author: Travis Askham
c
c Copyright (c) 2018 Travis Askham
c 
c $COPYRIGHT$
c 
c Additional copyrights may follow
c 
c $HEADER$
c 
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are
c met:
c 
c - Redistributions of source code must retain the above copyright
c   notice, this list of conditions and the following disclaimer.
c 
c - Redistributions in binary form must reproduce the above copyright
c   notice, this list of conditions and the following disclaimer listed
c   in this license in the documentation and/or other materials
c   provided with the distribution.
c 
c - Neither the name of the copyright holders nor the names of its
c   contributors may be used to endorse or promote products derived from
c   this software without specific prior written permission.
c 
c The copyright holders provide no reassurances that the source code
c provided does not infringe any patent, copyright, or any other
c intellectual property rights of third parties.  The copyright holders
c disclaim any liability to any recipient for claims brought against
c recipient by any third party for infringement of that parties
c intellectual property rights.
c 
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
c OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
c LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
c THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
c (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
c OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.c     
c
c 
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
      mwPointer mxcreatedoublematrix
      mwPointer mxgetpr
      mwPointer A, TAU, JPVT, AREF, c, kptr
      mwSize mxgetm, mxgetn
      mwSignedIndex m, n, n2, minmn, mone, k
      mwSignedIndex lwork, info, lda
      double precision work1(10)
      mwSize numel
      double precision one, zero, ar, br, dk
      double precision, allocatable :: work(:), c2(:)
      character ch1, ch2
      character *100 msg

      mwPointer mxCreateNumericMatrix
      integer*4 mexPrintf
      integer*4 classid, ComplexFlag, iout
      mwPointer ptemp(1)
      
      ch1 = 'N'
      ch2 = 'N'
      one = 1.0
      zero = 0.0

c     get inputs

      A = mxgetpr(prhs(1))
      c = mxgetpr(prhs(2))
      m = mxgetm(prhs(1))
      n = mxgetn(prhs(1))
      n2 = mxgetn(prhs(2))

      lda = m

c     allocate outputs

      mone = 1
      minmn = m
      if (n .lt. m) minmn = n

      plhs(1) = mxcreatedoublematrix(m, n, 0.0)
      plhs(2) = mxcreatenumericmatrix(n, mone, 
     +     mxclassidfromclassname('int64'), 0)
      plhs(3) = mxcreatedoublematrix(minmn, mone, 0.0)


c     check if only partial pivoting desired (i.e. first k)
      if ( nrhs .gt. 2 ) then
         kptr = mxgetpr(prhs(3))
         call qrmc_dcopy_in(dk,%val(kptr),mone,mone)
         k = (dk+0.1d0)
      else
         k = minmn
      end if
      
      
c     initialize outputs

      AREF = mxgetpr(plhs(1))
      call qrmc_dcopy_in(%val(AREF),%val(A),m,n)
      JPVT = mxgetpr(plhs(2))
      call qrmc_iset_to_zero(%val(JPVT),n)
      TAU = mxgetpr(plhs(3))

c     make a copy of c

      allocate(c2(n))
      call qrmc_dcopy_in_vec(c2,%val(c),n)

c     query routine for size of work array

      lwork = -1
      call dqrmc(m,n,%val(AREF),lda,%val(JPVT),%val(TAU),work1,
     +     lwork,%val(c),k,info)
      lwork = work1(1)+1
c      write(msg,'(i20,a)') lwork, '\n'
c      iout = mexPrintf(msg)

c     allocate work array and do actual call

      allocate(work(lwork))
      call dqrmc(m,n,%val(AREF),lda,%val(JPVT),%val(TAU),work,
     +     lwork,c2,k,info)

      deallocate(work,c2)

      return
      end

      subroutine qrmc_dcopy_in(a,b,m,n)
c
c     Utility, copy contents of b into a
c     where a and b are double arrays of 
c     size m x n
c
      mwSignedIndex m,n
      double precision a(m,n), b(m,n)
      mwSignedIndex i, j
      
      do j = 1,n
         do i = 1,m
            a(i,j) = b(i,j)
         enddo
      enddo

      return
      end
      
      subroutine qrmc_dcopy_in_vec(v1,v2,n)
c
c     Utility, copy contents of v2 into v1
c     where v1 and v2 are double arrays of 
c     size n
c
      mwSignedIndex n
      double precision v1(n), v2(n)
      mwSignedIndex i
      
      do i = 1,n
         v1(i) = v2(i)
      enddo

      return
      end
      

      subroutine qrmc_iset_to_zero(jpvt,n)
c
c     utility, set entries of jpvt to zero
c
      mwSignedIndex jpvt(*)
      mwSignedIndex n, i

      do i = 1,n
         jpvt(i) = 0
      enddo
      
      return
      end

      SUBROUTINE dqrmc( M, N, A, LDA, JPVT, TAU, WORK, LWORK, C, K,INFO)
c
c     THIS IS A MODIFICATION OF DGEQP3, ACTING AS A WRAPPER FOR 
c     BLOCKED QRs WITH COST FUNCTIONS
c
c  -- LAPACK computational routine (version 3.7.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     December 2016
c
c     .. Scalar Arguments ..
      mwSignedIndex            INFO, LDA, LWORK, M, N, K
c     ..
c     .. Array Arguments ..
      mwSignedIndex            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ), C(*)
c     ..
c
c  =====================================================================
c
c     .. Parameters ..
      mwSignedIndex            INB, INBMIN, IXOVER
      PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LQUERY
      mwSignedIndex       FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB,
     $     NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN, IONE, IZERO,
     $     iout
c     ..
c     .. External Subroutines ..
      EXTERNAL           DGEQRF, DLAQP2_MC, DLAQPS_MC, DORMQR, DSWAP
      EXTERNAL           XERBLA
c     ..
c     .. External Functions ..
      mwSignedIndex            ILAENV
      DOUBLE PRECISION   DNRM2
      EXTERNAL           ILAENV, DNRM2
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
      character *100 msg      
c     ..
c     .. Executable Statements ..
c
c     Test input arguments
c  ====================
c

      IONE = 1
      IZERO = 0
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )

      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
c
      IF( INFO.EQ.0 ) THEN
         MINMN = MIN( M, N )
         IF( MINMN.EQ.0 ) THEN
            IWS = 1
            LWKOPT = 1
         ELSE
            IWS = 3*N + M + N
            LWKOPT = IWS
         END IF
         WORK( 1 ) = LWKOPT
c
         IF( ( LWORK.LT.IWS ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
c
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQP3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

*
*     Initialize pivots
*
      DO 10 J = 1, N
         JPVT( J ) = J
 10   CONTINUE

c     
c     Initialize partial column norms. The first N elements of work
c     store the exact column norms and pivoted costs.
c     
      DO 20 J = 1, N
         WORK( J ) = DNRM2( M, A( 1, J ), IONE )
         WORK( N+J ) = WORK( J )
 20   CONTINUE

c
c     Use unblocked code to factor.
c     
c     
      CALL DLAQP2_MC( M, N, IZERO, A( 1, 1 ), LDA, JPVT( 1 ),
     $     TAU( 1 ), WORK( 1 ), WORK( N+1 ),
     $     WORK( 2*N+1 ), WORK(3*N+1), C(1), K )

c
      WORK( 1 ) = IWS
      RETURN
c
c     End of DGRMC3
c
      END

      SUBROUTINE DLAQP2_MC( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,
     $                   VN3, WORK, C, K )
*
*     MODIFIED: includes a cost per column now
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      mwSignedIndex            LDA, M, N, OFFSET, IONE, K
*     ..
*     .. Array Arguments ..
      mwSignedIndex            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), VN1( * ), VN2( * ),
     $                   WORK( * ), VN3(*), C(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      mwSignedIndex            I, ITEMP, J, MN, MNK, OFFPI, PVT, II
      DOUBLE PRECISION   AII, TEMP, TEMP2, TOL3Z, CTEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      mwSignedIndex            IDAMAX
      DOUBLE PRECISION   DLAMCH, DNRM2
      EXTERNAL           IDAMAX, DLAMCH, DNRM2
*     ..
*     .. Executable Statements ..
*

      IONE = 1

      MN = MIN( M-OFFSET, N )
      MNK = MIN( MN, K)
      TOL3Z = SQRT(DLAMCH('Epsilon'))
*
*     Compute factorization.
*
      DO 20 I = 1, MNK
*
         OFFPI = OFFSET + I
*
*        Determine ith pivot column and swap if necessary.
*

         DO II = I,N
            VN3(II) = VN1(II) - C(II)
         ENDDO

         PVT = I
         VN3MAX = VN3(I)
         DO II = I+1,N
            IF (VN3(II) .gt. VN3MAX) THEN
               VN3MAX = VN3(II)
               PVT=II
            ENDIF
         ENDDO
c         PVT = ( I-1 ) + IDAMAX( N-I+1, VN3( I ), IONE )
*
         IF( PVT.NE.I ) THEN
            CALL DSWAP( M, A( 1, PVT ), IONE, A( 1, I ), IONE )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( I )
            JPVT( I ) = ITEMP
            CTEMP = C( PVT )
            C( PVT ) = C( I )
            C( I ) = CTEMP
            VN1( PVT ) = VN1( I )
            VN2( PVT ) = VN2( I )
         END IF
*
*        Generate elementary reflector H(i).
*
         IF( OFFPI.LT.M ) THEN
            CALL DLARFG( M-OFFPI+1, A( OFFPI, I ), A( OFFPI+1, I ),IONE,
     $                   TAU( I ) )
         ELSE
            CALL DLARFG( IONE, A( M, I ), A( M, I ), IONE, TAU( I ) )
         END IF
*
         IF( I.LT.N ) THEN
*
*           Apply H(i)**T to A(offset+i:m,i+1:n) from the left.
*
            AII = A( OFFPI, I )
            A( OFFPI, I ) = ONE
            CALL DLARF( 'Left', M-OFFPI+1, N-I, A( OFFPI, I ),IONE,
     $                  TAU( I ), A( OFFPI, I+1 ), LDA, WORK( 1 ) )
            A( OFFPI, I ) = AII
         END IF
*
*        Update partial column norms.
*
         DO 10 J = I + 1, N
            IF( VN1( J ).NE.ZERO ) THEN
*
*              NOTE: The following 4 lines follow from the analysis in
*              Lapack Working Note 176.
*
               TEMP = ONE - ( ABS( A( OFFPI, J ) ) / VN1( J ) )**2
               TEMP = MAX( TEMP, ZERO )
               TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
               IF( TEMP2 .LE. TOL3Z ) THEN
                  IF( OFFPI.LT.M ) THEN
                     VN1( J ) = DNRM2( M-OFFPI, A( OFFPI+1, J ), IONE)
                     VN2( J ) = VN1( J )
                  ELSE
                     VN1( J ) = ZERO
                     VN2( J ) = ZERO
                  END IF
               ELSE
                  VN1( J ) = VN1( J )*SQRT( TEMP )
               END IF
            END IF
   10    CONTINUE
*
   20 CONTINUE
*
      RETURN
*
*     End of DLAQP2_MC
*
      END
