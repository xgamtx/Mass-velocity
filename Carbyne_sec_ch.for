	PROGRAM CARBYNE
	IMPLICIT NONE
	EXTERNAL BESS_NUL, BESS, F02HAF
	INTEGER, PARAMETER :: MAXSUM=50
C	?enei aaeaiee a naoea ii R
	INTEGER, Parameter:: NGRID = 300

C	Iaeneiaeuiia cia?aiea eaaioiaiai ?enea l

C	INTEGER, Parameter:: NLMAX = 12
      DOUBLE PRECISION, Parameter:: evs = 13.6058D0


	CHARACTER, PARAMETER :: JOB*1 = 'V'
C	Constraint: JOB = 'N' or 'V'.
	CHARACTER, PARAMETER :: UPLO*1 = 'U'

	DOUBLE PRECISION, Parameter:: 
	+LightSpeed = 274.0746D0
	DOUBLE PRECISION, Parameter:: 
	+PI = 3.14159265358979323846264338328D0

	DOUBLE COMPLEX, ALLOCATABLE :: AHA(:,:), ASMALL(:,:,:,:,:), 
	+BSMALL(:,:,:,:,:), HamAHA(:,:),WORK(:),
	+H(:,:), BigH(:,:,:)

	DOUBLE COMPLEX :: ResPoM, SIMP_I12,
	+ResPoL, ResPoX, ResPoA, SUM_VAR
	DOUBLE PRECISION :: RES(MAXSUM), A, B, C, EPSI, BESSLJ,
	+BESSLJ1, K, Kstart, Kstop
	DOUBLE PRECISION, ALLOCATABLE :: KAPPA(:), RAD(:,:), BigD(:,:), 
	+RO(:,:), PHI(:,:), RSZC(:,:), PSIALL(:,:,:),PSIEALL(:,:,:),D(:), 
     +RMT(:),
	+PSIRALL(:,:,:),PSIERALL(:,:,:), V(:,:), RWORK(:), DE(:),EBig(:,:)
	INTEGER, ALLOCATABLE :: M(:), N(:), P(:), JRIS(:), NATOMS(:)
	INTEGER :: COUNTX, I, I_, X, L, MSMALL, J,
	+LENH, LEND, IT, NINTERVALS, NPNTS, NLMAX,
	+ II, JJ, INFO, LWORK, I1, I1_, NT

	OPEN(1, FILE = 'nt.txt', FORM = 'FORMATTED')
	READ(1,*) NT
      CLOSE(1)     

C************************************************************************************************************
C	 N?eouaaiea COUNTX, M, N, P
C************************************************************************************************************
      OPEN(1, FILE = 'countx.txt', FORM = 'FORMATTED')
	READ(1,*) COUNTX
	ALLOCATE (KAPPA(COUNTX), M(COUNTX), N(COUNTX), P(COUNTX), 
	+AHA(COUNTX,COUNTX),HamAHA(CountX, CountX),D(COUNTX))
	READ(1,*) A, B, C, EPSI, NPNTS, NLMAX
	DO I = 1, COUNTX
		READ(1,*) M(I), N(I), P(I)
	END DO
	close(1)
C	Ia?aiaiiua aey aeaaiiaeecaoee
	INFO=0
	LWORK = 64 * CountX * 2
	ALLOCATE(RWORK(7*CountX*2), WORK(LWORK))

C************************************************************************************************************
C	 ?AN?AO KAPPA
C************************************************************************************************************

      DO I=1,COUNTX	   
		CALL BESS_NUL(N(I),M(I),RES,MAXSUM)
		KAPPA(I)=RES(N(I))/A
c		WRITE(*,*) I, KAPPA(I)
	END DO

	OPEN(1, FILE = 'kstart_kstop.txt', FORM = 'FORMATTED')
	READ(1, *) kStart, kStop
	CLOSE(1)
C************************************************************************************************************
C	 N?eouaaiea JRIS-iiia? rmt a rad, RMT, RAD
C************************************************************************************************************
	ALLOCATE(RAD(NGRID,NT), V(NGRID,NT), JRIS(NT), RMT(NT),NATOMS(NT)) 
	OPEN(1, FILE = 'JRIS_RAD.txt', FORM = 'FORMATTED')
	DO I=1, NT
		READ(1,*) JRIS(I), RMT(I)
		READ(1,*) RAD(:,I)
	END DO
	CLOSE(1)

C************************************************************************************************************
C	 N?eouaaiea EPSI, NATOMS
C************************************************************************************************************

	OPEN(13, FILE = 'NATOMS.TXT', FORM = 'FORMATTED')
	DO IT=1, NT 
		READ(13, *) NATOMS(IT)
	END DO
	CLOSE(13)

C	IINEA N?EOUAAIE? ?ENEA AOIIIA A YEAIAIOA?IIE ??AEEA NATOMS II?II II?AAAEEOU IANNEAU OEEEIA?E?ANEEO EII?AEIAO

	ALLOCATE(RO(MAXVAL(NATOMS),NT), PHI(MAXVAL(NATOMS),NT),
	+RSZC(MAXVAL(NATOMS),NT))

C************************************************************************************************************
C	 N?EOUAAIEA OEEEIA?E?ANEEO EII?AEIAO, RO - RHO, PHI, RSZC - Z; 
C************************************************************************************************************
	
   	OPEN(24, FILE = 'COORD.txt', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO I=1,NATOMS(IT)
			READ(24,'(3f10.4)')RO(I,IT),PHI(I,IT),RSZC(I,IT)
		END DO
	END DO

	OPEN(15, FILE = 'LENH_LEND.txt')
	READ(15, *) LENH, LEND
C	NLMAX - YOI IAENEIAEUIIA CIA?AIEA L IAEAIUEIA(EAAIOIAIA ?ENEI)
	CLOSE(15)
	NLMAX=NLMAX-1
	NLMAX=1
	ALLOCATE(ASMALL(NLMAX,2*NLMAX+1,COUNTX,NPNTS,NT),
	+BSMALL(NLMAX,2*NLMAX+1,COUNTX,NPNTS,NT))
	ALLOCATE(BIGH(NPNTS,COUNTX,COUNTX), BIGD(NPNTS,COUNTX),
	+H(COUNTX,COUNTX),DE(CountX),
     +EBig(NPNTS,COUNTX))
	NINTERVALS=NPNTS-1
C************************************************************************************************************
C	N?eouaaiea PSIALL(NT, NLMAX, EPNTS), PSIEALL(NT, NLMAX, EPNTS), PSIRALL, PSIERALL
C************************************************************************************************************
	ALLOCATE(PSIALL(NT,NLMAX+1,MAXVAL(JRIS)),
	+PSIEALL(NT,NLMAX+1,MAXVAL(JRIS)),
	+PSIRALL(NT,NLMAX+1,MAXVAL(JRIS)),
     +PSIERALL(NT,NLMAX+1,MAXVAL(JRIS)))
 	OPEN(1, FILE = 'PSIALL.txt', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS(IT)
				READ(1, *) PSIALL(IT, L, X)

				PSIALL(IT,L,X) = PSIALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)

 	OPEN(1, FILE = 'PSIEALL.txt', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS(IT)
				READ(1, *) PSIEALL(IT, L, X)

				PSIEALL(IT, L, X)=PSIEALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)
	OPEN(1, FILE = 'PSIRALL.TXT', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS(IT)
				READ(1, *) PSIRALL(IT, L, X)

				PSIRALL(IT, L, X)=PSIRALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)
	OPEN(12, FILE='SO_EnergiesAll.txt')
	write(12,*)''
	close(12)
	OPEN(1, FILE = 'PSIERALL.TXT', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS(IT)
				READ(1, *) PSIERALL(IT, L, X)

				PSIERALL(IT, L, X)=PSIERALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)


	OPEN(1, FILE = 'WPOT.txt', FORM = 'FORMATTED')
	DO IT = 1, NT
		READ(1,*) V(:,IT)
	END DO
	CLOSE(1)

	OPEN(22,FILE='H_AMNP.TXT', ACCESS = 'DIRECT',
	+RECL = LENH+1000)
	DO J = 1, NPNTS
	 READ(22, REC=J)((H(II,JJ),JJ=1,COUNTX),II=1,COUNTX)
	 DO JJ = 1, COUNTX
	 DO II = 1, COUNTX
	  BIGH(J,II,JJ) = H(II,JJ)
	  END DO
	 END DO
	END	DO
	CLOSE(22)

	OPEN(23,FILE='D_ENERGIES.TXT',
	+ACCESS = 'DIRECT', RECL = LEND+1000)
	DO J=1, NPNTS
	 READ(23,REC=J) D 
	 DO II = 1, COUNTX
	  BIGD(J,II) = D(II)
	 END DO
	END DO
	CLOSE(23)
C    ***************************************************?an?ao A e B iaeaiueeo************************************
	DO I=1, COUNTX
		DO L=1,NLMAX
			DO MSMALL=-L, L
				DO J=1, NPNTS
					DO IT=1, NT
	K=Kstart + (kStop - kStart) * (DBLE(J-1)/DBLE(NINTERVALS))
	ASMALL(L,MSMALL+L+1,I,J,IT)=PSIEALL(IT,L+1,JRIS(IT))*
	+SIMP_I12(EPSI, RMT(IT),K+(2*PI/C)*DBLE(P(I)),KAPPA(I), MSMALL,L,2)
     +-PSIERALL(IT,L+1,JRIS(IT))*
     +SIMP_I12(EPSI,RMT(IT),K+(2*PI/C)*DBLE(P(I)), KAPPA(I), MSMALL,L,1)
	BSMALL(L,MSMALL+L+1,I,J,IT)=PSIRALL(IT,L+1,JRIS(IT))*
	+SIMP_I12(EPSI, RMT(IT),K+(2*PI/C)*DBLE(P(I)),KAPPA(I), MSMALL,L,1)
     +-PSIALL(IT,L+1,JRIS(IT))*SIMP_I12(EPSI, RMT(IT),
     +K+(2*PI/C)*DBLE(P(I)), KAPPA(I), MSMALL,L,2)
					END DO
				END DO
			END DO
		END DO
	END DO
C    ********************************************** ?an?aoa AHA e BHB ********************************************
      DO J=NPNTS,1,-1
		DO I1=1,COUNTX
 			DO I1_=I1,COUNTX
			  AHA(I1,I1_)=DCMPLX(0.0D0,0.0D0)
			  DO I=1,COUNTX
			    DO I_=1,COUNTX
				  SUM_VAR=DCMPLX(0.0d0,0.0d0)

				IF (M(I).EQ.M(I_)) THEN
					SUM_VAR=1/(LightSpeed**2*C*A**2)
					CALL BESS(M(I),KAPPA(I)*A,BESSLJ,BESSLJ1)
					SUM_VAR=SUM_VAR/BESSLJ1
					CALL BESS(M(I_),KAPPA(I_)*A,BESSLJ,BESSLJ1)
					SUM_VAR=SUM_VAR/BESSLJ1
					MSMALL=M(I)
					ResPoA = dcmplx(0.0d0, 0.0d0)
					DO IT = 1, NT
						ResPoL=0
						DO L=ABS(M(I)),NLMAX 
							ResPoM=0
							DO X=1, NATOMS(IT)
								ResPoX=
	+mu1(V(:,IT),PSIALL(IT,L+1,:),RAD(:,IT),JRIS(IT),BIGD(J,I1),
     +BIGD(J,I1_))*
     +DCONJG(ASMALL(L,MSMALL+L+1,I,J,IT))*ASMALL(L,MSMALL+L+1,I_,J,IT)+
	+mu3(V(:,IT),PSIEALL(IT,L+1,:),RAD(:,IT),JRIS(IT),BIGD(J,I1),
     +BIGD(J,I1_))*
     +DCONJG(BSMALL(L,MSMALL+L+1,I,J,IT))*BSMALL(L,MSMALL+L+1,I_,J,IT)+
     +mu2(V(:,IT), PSIEALL(IT,L+1,:), PSIALL(IT,L+1,:),RAD(:,IT),
     +JRIS(IT),BIGD(J,I1),BIGD(J,I1_))*
	+(DCONJG(ASMALL(L,MSMALL+L+1,I,J,IT))*BSMALL(L,MSMALL+L+1,I_,J,IT)+
	+DCONJG(BSMALL(L,MSMALL+L+1,I,J,IT))*ASMALL(L,MSMALL+L+1,I_,J,IT))
		ResPoM=ResPoM+ResPoX*DCMPLX(DCOS((2*PI/C)*DBLE(P(I)-P(I_))*
	+	RSZC(X,IT)),DSIN((2*PI/C)*DBLE(P(I)-P(I_))*RSZC(X,IT)))
							END DO
						ResPoL=ResPoL+ResPoM*(2*L+1)
	+					*FCT(L-ABS(MSMALL))/FCT(L+ABS(MSMALL))
						END DO
					ResPoA = ResPoA + ResPoL * RMT(IT)**4
					END DO
				SUM_VAR = SUM_VAR * ResPoA
				END IF					  
				AHA(I1,I1_)=AHA(I1,I1_)-SUM_VAR*DCONJG(BIGH(J,I,I1))*
	+			BIGH(J,I_,I1_)
			END DO
		END DO
		IF (I1.EQ.I1_) THEN
			AHA(I1, I1_)=AHA(I1, I1_)-1/(LightSpeed**2)*BIGD(J,I1)
	+		*BIGD(J,I1_)+ BIGD(J,I1)
		END IF
			END DO
		END DO
		write(*,*)J,NPNTS

		info=0
		CALL F02HAF(JOB,UPLO,COUNTX,AHA,COUNTX,DE,
	+RWORK,WORK,LWORK,INFO)
		OPEN(100,FILE = 'SO_EnergiesAll.txt', STATUS = 'OLD',
	+FORM = 'FORMATTED', POSITION = 'APPEND')
		EBig(J,:) = DE
		DO I = 1, CountX
			EBig(J,I) = EBig(J,I)*evs
			K= kStart + (kStop - kStart) *(DBLE(J-1)/DBLE(NINTERVALS))
			WRITE(100,*) K, EBig(J,I)
	 
		END DO
		close(100)
	END DO
	CONTAINS
	
      INTEGER FUNCTION FCT(X)
	INTEGER :: X, I, RES
	RES=1
	DO I=2,X
		RES=RES*I
	END DO
	FCT=RES
	END FUNCTION

	DOUBLE PRECISION FUNCTION GetF(E1,E2,V)
	DOUBLE PRECISION :: E1,E2,V
	GetF=(E1-V)*(E2-V)-E1*E2
	END FUNCTION GetF

c     ************************************Eioaa?ae acaoa1 IA?AEI************************************************************

	DOUBLE PRECISION FUNCTION pod1(E1,E2,V,U, rad)
	DOUBLE PRECISION :: E1,E2,V,U, rad
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod1 = GetF(E1,E2,V)*U**2*rad**2

	END FUNCTION pod1

	DOUBLE PRECISION FUNCTION mu1(V, u, rad, JRIS, E1, E2)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), u(:), rad(:), E1, E2, integ_1

	integ_1 = 0.d0

	DO i=1, JRIS-1
		integ_1 = integ_1 + (pod1(E1,E2,V(i),u(i),rad(i))+ 
	+pod1(E1,E2,V(i+1),u(i+1),rad(i+1)))/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	mu1=integ_1

	END FUNCTION mu1

c     ************************************Eioaa?ae acaoa2 *********************************************************
C************************************************************************************************************

	DOUBLE PRECISION FUNCTION pod2(V, UE, U, rad, E1, E2)
	DOUBLE PRECISION :: V, UE, U, rad, E1, E2

C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vdurddu(L,IT,i)

	pod2=GetF(E1,E2,V)*UE*U*rad**2
c	pod2=UE*U*rad
c	Vdurddu = V(i,IT)*PSIEALL(IT, L, i)*PSIERALL(IT, L, i)*RAD(i,IT)

	END FUNCTION pod2

C************************************************************************************************************

	DOUBLE PRECISION FUNCTION mu2(V, UE, U, rad, JRIS, E1, E2)
	INTEGER :: i, JRIS
	DOUBLE PRECISION ::integ_2, V(:), UE(:), U(:), rad(:), E1, E2 

C	I?ioaao?a au?eneyao cia?aiea eioaa?aea VdurdduIntValue(l, IT) iaoiaii i?yiioaieuieeia ii n?aaiaio

	integ_2 = 0.d0

	DO i=1, JRIS-1
		integ_2 = integ_2 + (pod2(V(i), UE(i), U(i), rad(i), E1, E2) 
	+	+ pod2(V(i+1), UE(i+1), U(i+1), rad(i+1), E1, E2))
     +    /2.D0 * (RAD(i+1)-RAD(i))
      ENDDO
	mu2=integ_2
	END FUNCTION mu2	
C************************************************************************************************************

	DOUBLE PRECISION FUNCTION pod3(V, UE, rad, E1, E2)
	DOUBLE PRECISION :: V, ue, rad, E1, E2
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey UnderInt3(L,IT,i)	aey caoa 3 aac eniieuciaaiey i?iecaiaiie ii V
	   
	pod3 = GetF(E1,E2,V) * UE**2 * rad**2

	END FUNCTION pod3

C************************************************************************************************************

	DOUBLE PRECISION FUNCTION mu3(V, UE, rad, JRIS, E1, E2)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), UE(:), rad(:), E1, E2, integ_3

C	I?ioaao?a au?eneyao cia?aiea eioaa?aea SubInt3(l, IT) iaoiaii i?yiioaieuieeia ii n?aaiaio,
C	 eioi?ue aoiaeo a eioaa?ae caoa3

	integ_3 = 0.d0

	DO i=1, JRIS-1
		integ_3 = integ_3 + (pod3(V(i), UE(i), rad(i), E1, E2) 
	+	+ pod3(V(i+1), ue(i+1), rad(i+1), E1, E2))/2.D0 
     +	* (rad(i+1)-rad(i))
      ENDDO
	mu3=integ_3
	END FUNCTION mu3

C************************************************************************************************************

	END PROGRAM