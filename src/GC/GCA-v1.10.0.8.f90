!
!	Group Contribution with Association Equation of State (GCA-EoS)
!	Basada en la GC-EoS de Steen Skjold-J�rgensen con el agregado de una contribuci�n
!   asociativa derivada de la ecuaci�n de Chapman (1989,1990), pero formulada a con-
!   tribuci�n grupal.
!
!
SUBROUTINE GCEOS (NC, NG, NST, iDer, iTemp, T, P, XN, phi, dLPhi, dLPhiT, DLPHIP, CZ, IGz, MTyp, IC,fallo)
!-----------------------------------------------------------------------
!	SUBRUTINA GCEOS calcula el logaritmo de los coeficientes de fugaci-
!	dad y sus derivados con respecto a la temperatura, presi�n y compo-
!	sici�n una mezcla de componentes NC.
!	XN especifica el n�mero de moles de cada componente.
!	IOPT especifica si el c�lculo debe ser hecho para la fase gaseosa
!	(MTYP =- 1) o la fase l��quida (MTYP = 1).
!
!	iDer = 1: derivadas 1� respecto del n�mero de moles
!	     = 2:     "     2�
!
!	iTemp = 1: derivadas respecto de T
!
!	iGz: inicializaci�n del c�lculo. Detalles en subr. ZMAX
!
!
!	PHI(I)    = ln del coeficiente de fugacidad del comp i
!     DLPHIT(I) = derivada del ln(phi i) respecto de T
!     DLPHIP(I) = derivada del ln(phi i) respecto de P
!     DLPHI(I,K)= derivada del ln(phi i) respecto de nk
!-----------------------------------------------------------------------
      IMPLICIT real*8(A-H,O-Z)
!
	parameter(ncm = 15, ngm = 15, ngam = 12, NSM = 24)
	character*20 :: versGCA = '1.10.0.8 (2.0RC0)'
	integer sigma
	common/versSUB/versGCA
!
	dimension x(NCM),xms(NCM),Ps(NCM,NGM),xN(NC),phi(NCM),            &
     &          help3(NCM,NGM),help7(NCM,NGM),dLam1(NCM),dLam2(NCM),     &
     &          dLam3(NCM),Z(2),dYdn(NCM),dYvdn(NCM)
!
	DIMENSION DLPHI(NCM,NCM),HELP8(NCM,NGM),HELP9(NCM,NGM),            &
     &          HELP10(NCM,NGM),DFDNN(NCM,NCM),DPDN(NCM),DFDN(NCM)
!
	DIMENSION DH3DT(NCM,NGM),DH7DT(NCM,NGM),DH2DT(NGM),DH4DT(NGM),     &
     &          DH5DT(NGM),DH6DT(NGM),DLAMT1(NCM),DLAMT2(NCM),           &
     &          DLAMT3(NCM),DYDTN(NCM)
!
	real*8 LAMBDA2,LAMBDA,invXs_m(NSM,NSM)
	dimension aux(NSM),dXs_dn(NSM,NCM),dXs_ds(NSM,NSM),                &
     &          d2Fassoc_dsdT(NSM)


	DIMENSION DLPHIT(NCM),DLPHIP(NCM),DFDNT(NCM)
!
!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),rNyT(NGM),QT,    &
     &           TETA(NGM),HELP6(2,NGM),HELP11(2,NGM),HELP12(2,NGM),     &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),XLAM1,XLAM2,XLAM3
!
	COMMON/COORD/ZZ
!
!	Variables espec��ficas del t�rmino dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)
!
!	Versi�n asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Versi�n asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!	Variables ya calculadas en la subrutina ZMAX
	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
	
	common/ScndDer/dPV, dPDT, dPDn
	logical::fallo
    igz = 0
	PI=3.1415926536D0
	RT=T*R
	dPdV = 0
	dPdT = 0
	dPdn = 0
	dXs_ds(:NST,:NST) = 0.D0
!
!	Calcula el n�mero total de moles, fracciones molares, moles
!	de superficie total, y las fracciones superficie de area
	xNtot = sum(xN(:NC))
!
!	Define la fase y su factor de compresibilidad
	IPHAS=1+(MTYP+1)/2
	if (IGZ /= 0 .AND. MTYP /= 0)Z(IPHAS) = CZ
	CALL ZMAX (Z, MTYP, T, P, XN, IGZ, NC, NG ,NST, IC,fallo)
	IPHAS = 1 + (MTYP*IC + 1)/2
	if (MTYP == 0) IPHAS=1+(IC+1)/2
	CZ=Z(IPHAS)
!
!	C�lculo de coeficientes de fugacidad.
	rNyT(:NG) = rNyT(:NG)*xNtot

	QT = QT*XNTOT               	!�rea total de los grupos atractivos
	V = Z(IPHAS)*RT/P*XNTOT     	!Volumen de la fase
	DPV = -DPDV(IPHAS)*XNTOT/V/V	!(dP/dV)_(T,n): viene del COMMON "ZZZ"
!
	DO 22 I=1,NC

	  XSUM=0.D0
	  DO 2 K=1,NG

	    PS(I,K)=0.D0             	!PS(i,k) = ny(k,i)*q(k)
!	    if (NYN == 0) GOTO 2
	    if (NY(I,K) /= 0) then

!	      PSIK=DFLOAT(NYN)*Q(K)		!Necesita convertir Ny en n�mero decimal, sino FORTRAN es capaz de convertirle Ps(i,k) en un entero...
	      PS(I,K)=DFLOAT(NY(I,K))*Q(K)	!(PS)ik	!Necesita convertir Ny en n�mero decimal, sino FORTRAN es capaz de convertirle Ps(i,k) en un entero...
	      XSUM=XSUM+PS(I,K)

	    endif

2	   enddo
	   XMS(I)=XSUM                   	!MS(i) = SUM(ny(k,i)*q(k), k=1,NG)

22	enddo
!
!	Variable auxiliar para tau_ij:
      DNOM=QT/V/RT
!
	dPdTat = 0.D0
	if (NST  >  0) then
!
!	  Esta subrutina trabaja con volumen total en vez de molar como
!	  la sub. ZMAX
	  sm(:NST) = xnTot*sm(:NST)
	  dXs_dV(iPhas,:NST) = dXs_dV(iPhas,:NST)/xnTot
	  sm_Xs(iPhas,:NST) = sm_Xs(iPhas,:NST)*xnTot

	endif
!
!     Variables auxiliares para el calculo del ln(coef.fug)i y sus derivadas
!     con respecto a T,P y nro. de moles del comp. j
!
!     Contribuci�n atractiva
	DO 6 I=1,NC

	  XMSI=XMS(I)
	  DO 61 K=1,NG

	    HEL3K=0.D0
	    HEL7K=0.D0
	    HEL9=0.D0
	    HEL10=0.D0
	    DO 23 J=1,NG

	      if (ny(i,j) > 0) then

	        TAU=E(IPHAS,J,K)
	        AJK=A(J,K)*DNOM
	        PSIJ=PS(I,J)*TAU
	        HEL7K=HEL7K+PSIJ         !H7ik
	        PTIJ=PSIJ*DLOG(TAU)
	        PSIJ=PSIJ*AJK
	        HEL3K=HEL3K+PSIJ         !H3ik
	        HEL10=HEL10+PTIJ
	        HEL9=HEL9+PTIJ*AJK

	      endif
!
23	    enddo
	    H4=HELP4(IPHAS,K)
	    HEL3K=HEL3K/H4
	    HEL7K=HEL7K/H4
	    HEL9=HEL9/H4
	    HEL10=HEL10/H4
	    HELP3(I,K)=HEL3K        !H3ik/H4k
	    HELP7(I,K)=HEL7K        !H7ik/H4k
	    HELP9(I,K)=HEL9
	    HELP10(I,K)=HEL10
	    HELP8(I,K)=HEL7K-XMSI*(1.D0-HELP6(IPHAS,K))

61	  enddo

6	enddo
!
	if (iTemp /= 0) then
!
!  Derivadas del t�rmino atractivo con respecto de T
	  DH3DT(:NC,:NG)=0.D0
	  DH7DT(:NC,:NG)=0.D0
	  DHEL=0.D0
	  DO 300 J=1,NG

	    H30=0.D0
	    H31=0.D0
	    H32=0.D0
	    H33=0.D0
	    H34=0.D0
	    H35=0.D0
	    H36=0.D0
	    AA=A(J,J)
	    DAT=DADT(J,J)
	    DO 301 K=1,NG

	      TAU=E(IPHAS,K,J)
	      TET=TETA(K)
	      DAKJ=DADT(K,J)
	      ARG=DLOG(TAU)
	      AKJ=A(K,J)
	      DDEL=DAKJ-DAT
	      DARG=0.D0
	      DCRT=DABS(AKJ-AA)
	      if (DCRT  >  1.D-2) DARG=ARG*(DDEL/(AKJ-AA)-1.D0/T)
	      ALF=ALFA(K,J)
	      TETAU=TAU*TET
	      HELPI=TAU*DARG
	      HELP=TETAU*DARG
	      H30=H30+HELP
	      HELPI1=HELPI*AKJ*DNOM
	      HELP1=HELP*AKJ*DNOM
	      H31=H31+HELP1
	      H32=H32+HELP*ARG
	      H33=H33+HELP1*ARG
	      HEL2=TETAU*DAKJ*DNOM
	      HELPI2=TAU*DAKJ*DNOM
	      H34=H34+HEL2
	      H35=H35+HEL2*ARG
	      H36=H36+TETAU*DDEL*ALF*DNOM
	      DO 3010 I=1,NC

	        if (ny(i,k) /= 0) then

	          PSIK=PS(I,K)
	          DH3DT(I,J)=DH3DT(I,J)+HELPI1*PSIK+HELPI2*PSIK  !d(H3ij/H4j)/dT
	          DH7DT(I,J)=DH7DT(I,J)+HELPI*PSIK               !d(H7ij/H4j)/dT

	        endif

3010	      enddo

301	    enddo
	    H2=HELP2(IPHAS,J)
	    H4=HELP4(IPHAS,J)
	    H5=HELP5(IPHAS,J)
	    DH2T=(H31+H34)/H4
	    DH4T=H30/H4
	    DH5T=(H33+H35+H31)/H4
	    H6=HELP6(IPHAS,J)
	    DH6T=-H6/T+(H32+H36)/H4
	    DH2DT(J)=DH2T             !d(H2j/H4j)/dT
	    DH4DT(J)=DH4T             !d(H4j)/dT
	    DH5DT(J)=DH5T             !d(H5j/H4j)/dT
	    DH6DT(J)=DH6T             !d(H6j/H4j)/dT
!
!---------(dP/dT)att = -RT*(d2(F)att/dVdT) - R*(d(F)att/dV)
!
	    HEL=DH5T+DH2T-DH2T*H6-H2*DH6T-(H5+H2-2.D0*H2*H6 )*DH4T
	    DHEL = DHEL + HEL*rNYT(J)*Q(J)    !- 2*V/(ZZ*RT)*(dP/dT)att

300	  enddo
!
	  DPDT=XNTOT*R/V-ZZ/2.D0*RT/V*DHEL  !n*R/V + (dP/dT)att (falta contribucion asociativa y repulsiva)
	  dPatdT = -ZZ/2.D0*RT/V*DHEL
!
!	  Ac� termina el if (iTemp == 0)goto 306
!--------------------------------------------------------------------------

	endif
	dXs_dn(:NST,:NC) = 0D0
	dPdTas = 0D0
306	if (NST > 0) then

	  if (iDer == 2 .OR. iTemp == 1) then
!
!  Calculation of the number of associating site moles derivatives of the non-bonded fraction. This is required for the
!  calculation of d2[ln(phi_k)]/[dni dnj], and for d[ln(phi_k)]/dT and dP
!
	    if (NST == 1) then
	    
	      dXs_ds(1,1) = -dXs_dV(iPhas,1)*V/sm(1)
	
	    elseif (NST == 2 .AND. Delta(1,1) <= 0 .AND. Delta(2,2) <= 0) then
	    
	      dXs_ds(major,major) = -Xs(iPhas,major)**2 * Delta(1,2)/V/2D0*( (2D0 - b_aux(iPhas))/root(iPhas) - 1D0)
	      dXs_ds(major,3-major) = -Xs(iPhas,major)**2 *Delta(1,2)/V/2D0*(1D0 + b_aux(iPhas)/root(iPhas))
	      if (sm(major) == sm(3-major)) then
	  
	        dXs_ds(3-major,3-major) = dXs_ds(major,major)
	        dXs_ds(3-major,major) = dXs_ds(major,3-major)
! 	        dXs_ds(major,3-major) = dXs_ds(major,major)
	  
	      else
	  
	        dXs_ds(3-major,major) = -Xs(iPhas,major)**2 *Delta(1,2)/V*(Xs(iPhas,major) + sm(major)*dXs_ds(major,major))
	        dXs_ds(3-major,3-major) = -Xs(iPhas,major)**2 *Delta(1,2)/V*sm(2)*dXs_ds(major,3-major)
		  
	      endif
	
	    else
	
	      do i = 1, NST
	
	        if (sm(i) >= 1D-16) dXs_ds(i,:NST) = -sm(i)*Xs(iPhas,:NST)*Delta(i,:NST)/V/xnTot
! 	        dXs_ds(i,:NST) = -sm(i)*Xs(iPhas,:NST)*Delta(i,:NST)/V/xnTot
	        dXs_ds(i,i) = dXs_ds(i,i) + (1D0/Xs(iPhas,i) - 1D0 - dot_product( sm_Xs(iPhas,:NST), Delta(i,:NST) )/V)/xnTot
	      
	      enddo
!
!	      Conviene hallar uno u otro primero, dependiendo qu� sea m�s
!	      grande: NC o NST
! 	      dXs_ds(:NST,:NST) = dXs_ds(:NST,:NST)/xnTot
! 	      if (NST < NC) then

	      do i = 1,NST

	        dXs_ds(:NST,i) = -dXs_ds(:NST,i)
	        call LUBksb (H(iPhas,:NST,:NST), NST, NST, indx(iPhas,:NST), dXs_ds(:NST,i))

	      enddo
!
!	      dXs/dn = dXs/ds*ds/dn, ds/dn = sigma

! 	    else
! !
! !	      matriz dg(i)/dn(l):
! 	      dXs_dn(:NST,:NC) = matmul(-dXs_ds(:NST,:NST),dfloat(sigma(:NST,:NC))) !/XNTOT
! 	      do i = 1,NC
! 
! 	        call LUBksb (H(iPhas,:NST,:NST), NST, NST, indx(iPhas,:NST), dXs_dn(:NST,i))
! 
! 	      enddo
! 
	    endif
	    dXs_dn(:NST,:NC) = matmul(dXs_ds(:NST,:NST),dfloat(sigma(:NST,:NC)))

	    if (iTemp == 1) then
!
!  Calculation of the temperature derivative of the Helmholtz free energy:
!  (dP/dT)Vn = -(d2A/dTdV)n
!  (d2A/dTdV)n = RT (d2(A/RT)/dTdV)n + R (d(A/RT)/dV)Tn
!  dFVasoc = d(Aasoc/(nRT))/dv = d(Aasoc/RT)/dV
	      dPdTas = - R*dFVas(iPhas)
	      do i = 1, NST  
! 
	        dPdTas = dPdTas + RT*sm_Xs(iPhas,i)*SUM( sm_Xs(iPhas,:NST)*dDeldT(i,:NST)*(dXs_dV(iPhas,:NST)/Xs(iPhas,:NST) - 5D-1/V) )/V
	     
	      enddo
	      do k = 1, NST

	        aux(:NST) = sm(:NST)*dXs_ds(:NST,k)
	        aux(:NST) = matmul( dDeldT(:NST,:NST), aux(:NST) )
	        d2Fassoc_dsdT(k) = dot_product( sm_Xs(iPhas,:NST), aux(:NST) )
	        d2Fassoc_dsdT(k) = -(d2Fassoc_dsdT(k) + Xs(iPhas,k)*dot_product( sm_Xs(iPhas,:NST), dDeldT(:NST,k) ))/V
	        
	        
! 	        d2Fassoc_dsdT(k) = 0.D0
! 	        do i = 1, NST
! 	          if (sm(i) >= 1E-16) then
! ! 	          if (sm(i) >= 1E-16) then
! 	    
! 	            d2Fassoc_dsdT(k) = d2Fassoc_dsdT(k) - sm_Xs(iPhas,i)*Xs(iPhas,k)*dDeldT(k,i)
! 	            do j = 1, NST
! 	              if (sm(j) >= 1E-16) then
! 	    
! 	                d2Fassoc_dsdT(k) = d2Fassoc_dsdT(k) - sm_Xs(iPhas,i)*sm(j)*dXs_ds(j,k)*dDeldT(i,j)
! 	                
! 	              endif
! 	            enddo
! 	            
! 	          endif
! 	        enddo	  
! 	        d2Fassoc_dsdT(k) = d2Fassoc_dsdT(k)/V

	      enddo	      
!
	    endif

	  endif

	endif
!
!  Contribuci�n repulsiva
	XLAM1=XLAM1*XNTOT    !lambda1
	XLAM2=XLAM2*XNTOT    !lamdda2
	XLAM3=XLAM3*XNTOT    !lamdda3
	TLAM1=0.D0
	TLAM2=0.D0
	TLAM3=0.D0
	DO 9 I=1,NC

	  DIA=D(I)
!	  if (ITEMP == 0) GOTO 311
	  if (iTemp /= 0) then
	    
!  d(dc)/dT y var. aux. para d(lambda_k)/dT
	    DIAT=DT(I)
	    DELDT=DIAT
	    XDT=XN(I)*DIAT
	    DLAMT1(I)=DELDT
	    TLAM1=TLAM1+XDT
	    DELDT=2.D0*DELDT*DIA
	    DLAMT2(I)=DELDT
	    XDT=XDT*DIA
	    TLAM2=TLAM2+XDT
	    DLAMT3(I)=1.5D0*DELDT*DIA
	    TLAM3=TLAM3+XDT*DIA
!
	  endif
311	  continue
	  DLAM1(I)=DIA/XLAM1         !d(lambda1)/dni / lambda1
	  DIAV=DIA*DIA
	  DLAM2(I)=DIAV/XLAM2        !d(lambda2)/dni / lambda2
	  DLAM3(I)=DIAV*DIA

9	enddo
	PI6 = PI/6.D0/V
	Y = 1.D0/(1.D0-PI6*XLAM3)
	Y2 = Y*Y
	DYDV = -Y2*PI6*XLAM3/V
!	if (ITEMP == 0) GOTO 312
	if (iTemp /= 0) then
!
	  TLAM2=TLAM2*2.D0
	  TLAM3=TLAM3*3.0D0
	  DYDT=Y2*PI6*TLAM3        !dY/dT
	  TLAM3=TLAM3/XLAM3        !d(lambda1)/dt / lambda1
	  TLAM2=TLAM2/XLAM2        !d(lambda2)/dt / lambda2
	  TLAM1=TLAM1/XLAM1        !d(lamdda3)/dt / lambda3
	  DYDVT=DYDV*TLAM3+2.D0/Y*DYDV*DYDT !d2Y/dTdV
!
	endif
312	continue
	DO 98 I = 1, NC

	  DYDN(I)=Y2*PI6*DLAM3(I)     	!dY/dni
	  if (iTemp /= 0) then
!
	    DYDTN(I) = 2.D0/Y*DYDN(I)*DYDT + DLAMT3(I)/DLAM3(I)*DYDN(I) !d2Y/dTdni
!
	  endif
298	  DYVDN(I) = (2.D0/Y*DYDV - 1.D0/V)*DYDN(I)	!d2Y/dVdni
	  DLAM3(I) = DLAM3(I)/XLAM3       		!d(lambda3)/dni / lambda3

98	enddo
	R0=XLAM2/XLAM3
	R1=XLAM1*R0
	R3=R0*R0*XLAM2
	R4=Y2-Y-DLOG(Y)
	R6=2.D0*Y-1.D0-1.D0/Y
	R33=2.D0+1.D0/Y2
	R34=Y-1.D0
	if (iTemp /= 0) then

	  R30=TLAM1+TLAM2-TLAM3
	  R31=3.D0*TLAM2-2.D0*TLAM3
!
!	  (dP/dT)rep = -R*d(F)rep/dV - R*T*d2(F)rep/dVdT
	  DPDTR=PREP(IPHAS)/T+DYDVT/DYDV*PREP(IPHAS)-RT*DYDV*(3.D0*R1*R30+R3*R31*R6+R3*R33*DYDT-XNTOT/Y2*DYDT)   !(dP/dT)rep
	  DPDT = DPDT + DPDTR + dPdTas  !n*R/V+(dP/dT)att+(dP/dT)rep+(dP/dT)asoc

	endif
314	continue
!-------------------------------------------------------------------------------
!
!
!  C�lculo del ln(coef.fug)i y sus derivadas con respecto a T,P y nro. de moles del comp. j
	DO 8 I=1,NC

	  dFdnas = 0D0
! 	  if (maxval(sigma(:NST,i)) > 0) dFdnas=1.D0	!lo hallo por productoria, no sumatoria
	  dPdnas=0.D0
	  dFdnTas=0.D0
	  DFDTA=0.D0
	  DFDNI=0.D0
	  DPDNI=0.D0
	  DYDNI=DYDN(I)
	  XMSI=XMS(I)
	  NC1=I
	  DO 54 K = NC1,NC

	    DFDNN(I,K) = 0.D0

54	  enddo
	  if (iTemp /= 0) then
!
	    TLAM1I=DLAMT1(I)/XLAM1
	    TLAM2I=DLAMT2(I)/XLAM2
	    TLAM3I=DLAMT3(I)/XLAM3
	    DYDTNI=DYDTN(I)
!
	  endif
299	  CONTINUE
!
!  Contribucion atractiva
	  DO 50 J = 1,NG

	    TET=TETA(J)
	    PSIJ=PS(I,J)
	    H3IJ=HELP3(I,J)
	    H7IJ=HELP7(I,J)
	    H2=HELP2(IPHAS,J)
	    H4=HELP4(IPHAS,J)
	    H5=HELP5(IPHAS,J)
	    H6=HELP6(IPHAS,J)
	    H8IJ=HELP8(I,J)
	    if (iDer > 1 .OR. iTemp /= 0) then
!
	      H12=HELP12(IPHAS,J)
	      H10IJ=HELP10(I,J)
	      H9IJ=HELP9(I,J)
	      H11=HELP11(IPHAS,J)
	      H21=RNYT(J)*Q(J)
	      H20=H5+H2-H2*H6

	    endif
51	    HDER=-H2*(H7IJ-XMSI*(1.D0-H6))+H3IJ+XMSI*H5
	    DELI=-(HDER*TET+PSIJ*H2)
	    DFDNI=DFDNI+DELI                 !2./zz*d(F)att/dni
	    if (iTemp > 0) then
!
	      DH2T=DH2DT(J)
	      DH4T=DH4DT(J)
	      DH5T=DH5DT(J)
	      DH6T=DH6DT(J)
	      HELPT=PSIJ*(DH2T-H2*DH4T)+TET*(DH3DT(I,J)/H4+XMSI*DH5T-(H3IJ+XMSI*H5)*DH4T)
	      HELPT=HELPT-TET*H8IJ*(DH2T-2.D0*H2*DH4T)-TET*H2*(DH7DT(I,J)/H4-DH4T*XMSI+DH6T*XMSI)
	      DFDTA=DFDTA-DELI/T-HELPT         !2./zz*d2(F)att/dnidT
!
	    endif
316	    CONTINUE
	    if (iDer > 1 .OR. iTemp > 0) then
!
	      DP = -PSIJ*H20
	      DP=DP-H21/QT*(XMSI*(3.D0*H5+H2+H11)+H9IJ+H3IJ)
	      DP=DP+H21/QT*(H6*(H3IJ+2.D0*XMSI*H5)+H2*(H10IJ+XMSI*H12)+H7IJ*(H5+H2)+3.D0*XMSI*H2*H6)
	      DP=DP-H21/QT*2.D0*(H7IJ*H2*H6+XMSI*H2*H6*H6)
	      DPDNI=DPDNI+ZZ/2.D0*DP*RT/V      !-R*T*d2(F)att/dVdni
!
	      if (iDer > 1) then

	        H3H5=(H3IJ+XMSI*H5)/QT
	        H2H8=H2*H8IJ/QT
	        H7H6=TET*H2*(H7IJ+H6*XMSI)
	        H8=TET*H8IJ
	        H2H4=TET*H2
	        DO 53 K=NC1,NC

	          XMSK=XMS(K)
	          PSKJ=PS(K,J)
	          H9KJ=HELP9(K,J)
	          H3KJ=HELP3(K,J)
	          H7KJ=HELP7(K,J)
	          H10KJ=HELP10(K,J)
	          D2=(H3KJ+XMSK*H5)/QT
	          D3=XMSK/QT*(H9IJ+H3IJ)
	          D4=(H7KJ-XMSK*(1.D0-H6))/QT
	          D5=(H9KJ+XMSK*(H11+H5))/QT
	          D6=(H10KJ+XMSK*H12)/QT
	          D7=XMSK*H10IJ/QT
	          DELIK=-(PSIJ*D2+(PSKJ-TET*XMSK)*H3H5)
	          DELIK=DELIK+(PSKJ-TET*XMSK)*H2H8-D4*H7H6+D2*H8-DELI*D4
	          DELIK=DELIK+H2H4*(D7+XMSI*D6)-TET*(D3+XMSI*D5)
	          DFDNN(I,K)=DFDNN(I,K)+DELIK      !2./zz*d2(F)att/dnidnj

53	        enddo

	      endif

	    endif

50	  enddo
	  if (iDer == 2) then

	    do k = NC1, NC

!  Antes lo multiplicaba al sumar contribuciones. Desde la aparici�n de asociaci�n, multiplica 
!  ahora y luego suma... cosas heredadas...
	      dFdnn(i,k) = dFdnn(i,k)*zz/2.D0      !d2(F)att/dnidnj.

	    end do

	  end if

!  Contribuci�n asociativa
	  if (maxval(sigma(:NST,i)) > 0) then

!  Coeficiente de fugacidad:
	    do j = 1, NST
!
!  dF/dn(i)  = SUM(sigma(j,i)*ln(Xs(j), j=1..NST) = ln [PROD( Xs(j)^sigma(j,i), j=1..NST)]
	      dFdnas = dFdnas + dlog(Xs(iPhas,j))*dfloat(sigma(j,i))
! 	      dFdnas = dFdnas*Xs(iPhas,j)**dfloat(sigma(j,i))

	    enddo
! 	    dFdnas = dlog(dFdnas)
	    if (iTemp == 1 .OR. iDer == 2) then
!
!  dPasoc/dn(i) = -RT d2F/(dn(i)dV)
	      dPdnas = -RT*sum(dfloat(sigma(:NST,i))*dXs_dV(iPhas,:NST)/Xs(iPhas,:NST))


	    endif
!
!  d2(F)ass/dnidT
	    if (iTemp == 1) then

! 	      dFdnTas = sum(dfloat(sigma(:NST,i))*dXs_dT(:NST)/Xs(iPhas,:NST))
	      dFdnTas = dot_product( dfloat(sigma(:NST,i)), d2Fassoc_dsdT(:NST) )
! 	      dFdnTas = 0.D0
! 	      do k = 1, NST
! 	
! 	        d2Fassoc_dsdT_k = 0.D0
! 	        if (sigma(k,i) > 0) then
! 	          do j = 1, NST
! 	            if (sm(j) >= 1E-16 .AND. Delta(k,j) >= 1E-16) then
! 	     
! 	              d2Fassoc_dsdT_k = d2Fassoc_dsdT_k + sm_Xs(iPhas,j)*Xs(iPhas,k)*dDeldT(k,j)
! 	              do l = 1, NST
! 	                if (sm(l) >= 1E-16 .AND. Delta(j,l) >= 1E-16) then
! 	      
! 	                  d2Fassoc_dsdT_k = d2Fassoc_dsdT_k + sm_Xs(iPhas,j)*sm(l)*dXs_ds(l,k)*dDeldT(j,l)
! 	                  
! 	                endif
! 	              enddo
! 	              
! 	            endif
! 	          enddo	  
! 	          dFdnTas = dFdnTas - d2Fassoc_dsdT_k*dfloat(sigma(k,i))/V
! 	        endif
! 	      enddo
	      
	    endif
	    if (iDer == 2) then
	      do k = NC1, NC

	        aux_d2fn = sum(dfloat(sigma(:NST,i))*dXs_dn(:NST,k)/Xs(iPhas,:NST))
!
!  d2(F)att/dnidnj+d2(F)asoc/dnidnj
	        dFdnn(i,k) = dFdnn(i,k) + aux_d2fn
!
	      enddo
	    endif

	  endif
!
!  Contribucion repulsiva
	  DL1=DLAM1(I)
	  DL2=DLAM2(I)
	  DL3=DLAM3(I)
	  R2=DL1+DL2-DL3
	  R20=R2+DYDNI/R34
	  R20=R20*R1*3.D0
	  R5=3.D0*DL2-2.D0*DL3
	  DFDNR= R34*R20+R3*R4*R5+R3*R6*DYDNI+DLOG(Y)+XNTOT/Y*DYDNI                                      !d(F)rep/dni
!	  if (ITEMP == 0) GOTO 317
	  if (iTemp > 0) then

	    DL1TN=TLAM1I-TLAM1*DL1
	    DL2TN=TLAM2I-TLAM2*DL2
	    DL3TN=TLAM3I-TLAM3*DL3
	    DFDTR=DYDT*R20+R34*R30*R20+3.D0*R34*R1*(DL1TN+DL2TN-DL3TN+(DYDTNI-DYDT*DYDNI/R34)/R34)
	    DHEL=R31*R4*R5+DYDT*R6*R5+R4*(3.D0*DL2TN-2.D0*DL3TN)+R6*R31*DYDNI + R33*DYDT*DYDNI+R6*DYDTNI
	    DFDTR=DFDTR+R3*DHEL+DYDT/Y-XNTOT/Y*(DYDNI*DYDT/Y-DYDTNI) !d2(F)rep/dnidT
!
!  d2(A^R/RT)/dn(i)dT:
	    DFDNT(I)=DFDTA*ZZ/2.D0 + DFDTR + dFdnTas

	  endif
317	  CONTINUE
!	  IF (IDER == 1 .AND. ITEMP == 0) GOTO 8
	  if (iDer > 1 .OR. iTemp > 0) then
!
!  -R*T*d2(F)rep/dnidV
	    DPREP = PREP(IPHAS)*DYVDN(I)/DYDV - DYDV*RT*(3.D0*R1*R2 + R3*R6*R5 + R3*R33*DYDNI + 1.D0/Y - XNTOT/Y2*DYDNI)
!  dP/dni
	    DPDN(I) = DPDNI + DPREP + RT/V + dPdnas
!	    if (IDER == 1) GOTO 8
	    if (iDer > 1) then
!
	      DO 56 K = NC1, NC

	        DYDNK=DYDN(K)
	        DLAM2K=DLAM2(K)
	        DLAM1K=DLAM1(K)
	        DLAM3K=DLAM3(K)
	        R21=DLAM1K+DLAM2K-DLAM3K
	        R25=3.D0*DLAM2K-2.D0*DLAM3K
	        DREP=R20*(DYDNK+R34*R21)
	        DL1IK=DLAM1K*DL1
	        DL2IK=DLAM2K*DL2
	        DL3IK=DLAM3K*DL3
	        DYIK=DYDNI*DYDNK
	        DREP=DREP + 3.D0*R34*R1*( -DL1IK - DL2IK + DL3IK) - 3.D0*R1*(DYIK/R34 - DYIK*2.D0/Y)
	        DREP=DREP+R4*R5*R25*R3
	        DREP=DREP+R3*R5*R6*DYDNK+R3*R4*(2.D0*DLAM3K*DL3 - 3.D0*DLAM2K*DL2)
	        DREP=DREP+R3*R25*R6*DYDNI+R3*(R33+2.D0*R6/Y)*DYIK
	        DREP=DREP+(DYDNI+DYDNK)/Y+XNTOT*DYIK/Y2          !d2(F)rep/dnidnj
!
!  d2(A^r/RT)/dn(i)dn(k):
	        DFDNN(I,K)=DFDNN(I,K)+DREP   !(dFdnn_at + dFdnn_as) + dFdnn_rep

56 	      enddo
!
	    endif

	  endif
!  d(A^r/RT)/dn(i):
	  DFDN(I)=DFDNI*ZZ/2.D0+DFDNR +dFdnas

8	enddo
	DO 60 I=1,NC

	  PHI(I)=-DLOG(Z(IPHAS))+DFDN(I)          !ln(coef.de fugac)i
!	  if (ITEMP == 0) GOTO 62
	  if (iTemp > 0) then

	    DLPHIT(I)=DFDNT(I)+DPDN(I)*DPDT/DPV/RT+1.D0/T  !d(ln(coef.fug.)i)/dT
	    DLPHIP(I)=-DPDN(I)/DPV/RT - 1.D0/P               !d(ln(coef.fug.)i)/dP

	  endif
!62	  if (IDER == 1) GOTO 60
62	  continue
	  if (iDer > 1) then
	    DO 6100 K=I,NC

	      DLPHI(I,K)=DFDNN(I,K)+DPDN(I)*DPDN(K)/DPV/RT
	      DLPHI(I,K)=1.D0+XNTOT*DLPHI(I,K)         !n*d(ln(coef.fug.)i)/dnj
	      DLPHI(K,I)=DLPHI(I,K)                    !n*d(ln(coef.fug.)j)/dni

6100	    enddo
	  endif
	!por maxwell son iguales
!
!	      ln(coef fug i) = d(ln coef fug mezcla)/dni
!	      ln(coef fug j) = d(ln coef fug mezcla)/dnj
!	      Entonces la derivada de ln(coef fug)i respecto de j o viceversa [derivada de
!	      ln(coef fug)j respecto de i] termina dando la derivada cruzada de la misma funcion
!           LN(COEF FUG MEZCLA) y por Maxwell deben ser iguales.(Ec dif. exacta)
!
60	enddo
99	RETURN

ENDSUBROUTINE
!-----------------------------------------------------------------------
!
!	NTLPY Calcula la entalp��a residual a P y T de una mezcla de NC com-
!	ponentes.
!	Est� definida como
!
!	      r     conf          conf,GI
!	     H   = H   (T,P,n) - H       (T,P,n)
!
!	donde Hconf es la entalp��a configuracional.
!	xN es el vector de n�mero de moles de la mezcla (como siempre).
!	mTyp especifica la fases: -1 gas
!	                           1 l��quido
!	                           0 fase con menos energ��a libre de Gibbs
!
!	Las unidades, como en toda la GCA, est�n en [atm cm^3], ya que
!	R = 82.05 atm cm^3/(mol K).
!
SUBROUTINE NTLPY (NC, NG, NST, T, P, XN, HRES, CZ, IGZ, MTYP, IC)
!
	IMPLICIT real*8(A-H,O-Z)
	parameter(NCM = 15, NGM = 15, NSM = 24)
	integer sigma
!
	DIMENSION XN(NC),Z(2)
!
	real*8 LAMBDA2,LAMBDA,Id(NST,NST)
	dimension assoc_aux(NSM)
!
!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),rNyT(NGM),QT,    &
     &           TETA(NGM),HELP6(2,NGM),HELP11(2,NGM),HELP12(2,NGM),     &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),XLAM1,XLAM2,XLAM3
!
	COMMON/COORD/ZZ
!
!	Variables espec��ficas del t�rmino dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
!	Versi�n asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Versi�n asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!	Variables ya calculadas en la subrutina ZMAX
	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	logical::fallo
!
!
	PI=3.1415926536D0
	RT=T*R
!
!	Nro. total de moles
	xnTot = sum(xN(:NC))
!
!	C�lculo del factor de compresibilidad
	IPHAS=1+(MTYP+1)/2
	if (IGZ /= 0 .AND. MTYP /= 0) Z(IPHAS)=CZ
	CALL ZMAX(Z,MTYP,T,P,XN,IGZ,NC,NG,NST,IC,fallo)
!
!	Variable auxiliar para identificar tipo de fase:
	IPHAS=1+(MTYP*IC+1)/2
	if (MTYP == 0) then

	  IPHAS=1+(IC+1)/2
	  IPHBIS=IPHAS

	END IF
!	Factor compresibilidad de la fase elegida
	CZ=Z(IPHAS)
!
!	Esta subrutina trabaja con volumen total en vez de molar como
!	la sub. ZMAX
!	Nro. de moles del grupo att. k
	rNyT(:NG) = xNtot*rNyT(:NG)
!	�rea grupal total
	QT=QT*XNTOT
	if (NST > 0) then

	  sm(:NST) = xnTot*sm(:NST)
	  sm_Xs(iPhas,:NST) = sm_Xs(iPhas,:NST)*xnTot	  

	endif
!	Volumen de la fase
	V=Z(IPHAS)*RT/P*XNTOT
!-----------------------------------------------------------------------
!	C�lculo de la ent.residual:
!
!	                r       2  dF
!	               H  = -R T  ----  + P V - n R T
!	                           dT
!
!	Contribuci�n atractiva
	DNOM=QT/V/RT
	DADTA=0.D0
	ARES=0.D0
	DO 300 J = 1, NG

	  H30=0.D0
	  H31=0.D0
	  H34=0.D0
	  AA=A(J,J)
	  DAT=DADT(J,J)
	  DO 301 K=1,NG

	    TAU=E(IPHAS,K,J)
	    TET=TETA(K)
	    DAKJ=DADT(K,J)
	    ARG=DLOG(TAU)
	    AKJ=A(K,J)
	    DARG=0.D0
	    DDEL=DAKJ-DAT
	    DCRT=DABS(AKJ-AA)
	    if (DCRT > 1.D-2) DARG=ARG*(DDEL/(AKJ-AA)-1.D0/T)
	    TETAU=TAU*TET
	    HELP=TETAU*DARG
	    H30=H30+HELP
	    HELP1=HELP*AKJ*DNOM
	    H31=H31+HELP1
	    HEL2=TETAU*DAKJ*DNOM
	    H34=H34+HEL2

301	  enddo
	  DH2T=H31+H34
	  DH4T=H30
	  ARES=ARES+RNYT(J)*Q(J)*HELP2(IPHAS,J)               !Fatt=(Aatt)/R/T
	  DADTA=DADTA+RNYT(J)*Q(J)*(DH2T-DH4T*HELP2(IPHAS,J))/HELP4(IPHAS,J)

300	enddo
!	d(Fatt)/dT
	dAdTa=-ZZ/2.D0*(dAdTa-Ares/T)
!	write (*, *) dAdTa
!
!	Contribuci�n asociativa
	dFasdT = 0.D0
	if (NST > 0) then
!
	  assoc_aux(:NST) = matmul( dDeldT(:NST,:NST), sm_Xs(iPhas,:NST) )
	  dFasdT = -dot_product( assoc_aux(:NST),sm_Xs(iPhas,:NST) )/2D0/V

	endif
!	write (*, *) dFasdT
!
!	Contribuci�n repulsiva
	XLAM1=XLAM1*XNTOT        !lambda1
	XLAM2=XLAM2*XNTOT        !lambda2
	XLAM3=XLAM3*XNTOT        !lambda3
	TLAM1=0.D0
	TLAM2=0.D0
	TLAM3=0.D0
	DO 9 I=1,NC

	  DIA=D(I)
	  DIAT=DT(I)
	  XDT=XN(I)*DIAT
	  TLAM1=TLAM1+XDT
	  XDT=XDT*DIA
	  TLAM2=TLAM2+XDT
	  TLAM3=TLAM3+XDT*DIA

9	ENDDO
	PI6=PI/6.D0/V
	Y=1.D0/(1.D0-PI6*XLAM3)
	TLAM2=TLAM2*2.D0
	TLAM3=TLAM3*3.0D0
	DYDT=Y*Y*PI6*TLAM3
	R0=XLAM2/XLAM3
	R1=XLAM1*R0
	R3=R0*R0*XLAM2
	R4=Y*Y-Y-DLOG(Y)
	R6=2.D0*Y-1.D0-1.D0/Y
	R34=Y-1.D0
	R30=TLAM1/XLAM1+TLAM2/XLAM2-TLAM3/XLAM3
	R31=3.D0*TLAM2/XLAM2-2.D0*TLAM3/XLAM3
	!write (*, *) 3.D0*R1*R30*R34+3.D0*R1*DYDT+R3*R6*DYDT+R3*R31*R4+XNTOT/Y*DYDT
!
	DARDT=DADTA + dFasdT + 3.D0*R1*R30*R34+3.D0*R1*DYDT+R3*R6*DYDT+R3*R31*R4+XNTOT/Y*DYDT               !dF/dT
!
	HRES=(-DARDT*T+CZ*XNTOT-XNTOT)*RT    !Hres
!
	RETURN
	ENDsubroutine
!-----------------------------------------------------------------------
!
!
!
subroutine NTLPY_GI (NC, T, Hgi)
!
!	Esta subrutina calcula la entalp��a del gas ideal a T, para dar entalp��-
!	as reales complementando con las entalp��as residuales de NTLPY.
!
!	La subrutina almacena las entalp��as de GI molares dentro del vector Hgi
!
!	La Hgi se halla por correlaci�n. Los par�metros de la misma se leen en
!	en la subrutina PARMOL.
!
!	Francisco, 14/05/11. Durante la versi�n 1.3.3.
!
!	           14/07/11. Limpieza, v-1.3.4
!
	implicit real*8 (a-h, o-z)
	parameter (NCM = 15)
!	character*10 aname(ncm)
	character*10 aname
	dimension Hgi(NC)
	common/NAME/AName(NCM)
	common/GCPROM/PMM(ncm), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM),     &
&	              HHG(NCM)
!
!-----De momento planteo ingresar los par�metros como est�n en la correlaci�n de Passut-Danner: BTU/(molLB R)
	TRan = T*1.8D0	!K -> �R (
	Hgi = 0.D0
	Hgi(:NC) = HHA(:NC) + HHB(:NC)*TRan + HHC(:NC)*TRan**2 + HHD(:NC)*TRan**3 + HHE(:NC)*TRan**4 + HHF(:NC)*TRan**5
!
!	Traspaso de unidades de la GCA: atm*cm3/mol-g
	Hgi(:NC) = Hgi(:NC)*PMM(:NC)*1.0412238D4/454.D0	!(BTU/lb)*(10412 atm cm3/BTU)*(PM lb/mol-lb)*(454 mol-g/mol-lb)
!
	return

endsubroutine
!-----------------------------------------------------------------------
!
!
!
subroutine NTRPY_GI (NC, T, Sgi)
!
!	Esta subrutina calcula la entrop��a del gas ideal a T y 1 atm.
!
!	La subrutina almacena las entalp��as de GI molares dentro del vector Sgi
!
!	La Sgi se halla por correlaci�n. Los par�metros de la misma se leen en
!	en la subrutina PARMOL.
!
!
!	Francisco, 10/03/14. Durante la versi�n 1.9.26
!
	implicit real(8) (a-h,o-z)
	parameter(NCM = 15)
	character(10) aname
	real(8) Sgi(NC)
	common/NAME/AName(NCM)
	common/GCPROM/PMM(NCM), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM),     &
&	              HHG(NCM)
!
!	Correlaci�n de Passut-Danner: BTU/(mol-lb �R)
	TRan = T*1.8D0	!K -> �R
	Sgi = 0.D0
	do i = 1, NC
	  
	  Sgi(i) = HHB(i)*dlog(TRan) + 2.D0*HHC(i)*TRan + 1.5D0*HHD(i)*TRan**2 + 4.D0/3.D0*HHE(i)*TRan**3 &
&	           + 1.25D0*HHF(i)*TRan**4 + HHG(i)
!
!	  Traspaso de unidades de la GCA: atm*cm3/(mol-g K)
	  Sgi(i) = Sgi(i)*PMM(i)*1.0412238D4/454.D0*1.8D0
	  !(Sgi BTU/lb �R)*(10412 atm cm3/BTU)*(PM lb/mol-lb)*(454 mol-g/mol-lb)*(1.8 �R/K)
	 
	enddo
!
	return

endsubroutine

!-----------------------------------------------------------------------
!
!	La subrutina PARAGC eval�a los par�metros dependientes de la tem-
!	peratura
!
!
!	d = di�metros de esfera "blanda"
!     g(i,i) y k(i,j) = par�metros energ�ticos y de interacci�n
!	Delta(i,j) = fuerza de asociaci�n
!
SUBROUTINE PARAGC (T, NC, NG, NST, NTEMP)
!
	IMPLICIT real*8(A-H,O-Z)
	parameter(NCM=15,NGM=15,NGAM=12,NSM=24)
	real*8 kappa
	integer sigma
!
!	Variables espec��ficas del t�rmino dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
!
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
!
	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!
8	DO 1 I=1,NC

	  ARG=-2.D0*TC(I)/T/3.D0
	  DARG=DEXP(ARG)
!	  if (NTEMP == 0)GOTO 1
	  if (nTemp > 0) then

	    DT(I)=DC(I)*DARG*.12D0*ARG/T*1.065655D0

	  endif
	  D(I)=DC(I)*(1.D0-.12D0*DARG)*1.065655D0

1	enddo
	DO 2 I=1,NG

!	  if (T > TSPL(I))GOTO 4
	  if (T < Tspl(i)) then

	    TR=T/TSTR(I)
	    if (nTemp > 0) then
!	    if (NTEMP == 0) GOTO 5

	      DADT(I,I)=GSTR(I)*(G1(I)/TSTR(I)+G2(I)/T)

	    endif
5	    A(I,I)=GSTR(I)*(1.D0+G1(I)*(TR-1.D0)+G2(I)*DLOG(TR))
	    GOTO 2

	  else

	    A(I,I)=GSTR(I)/4.D0/(T/TSPL(I))**EPX(I)
!4	  A(I,I)=GSTR(I)/4.D0/(T/TSPL(I))**EPX(I)
!	  if (NTEMP == 0)GOTO 2
	    if (nTemp > 0) then

	      DADT(I,I)=-EPX(I)*A(I,I)/T

	    endif

	  endif

2	enddo
!2	continue

!	Si hay 1 s�lo grupo, no necesita k(i,j)
	if (NG >= 2) then
	  DO 3 I = 1, NG - 1

	    NG1=I+1
	    DO 30 J=NG1,NG

	      TR=T/(TSTR(I)+TSTR(J))*2.D0
	      XKIJ1=XKIJ(I,J)*(1.D0+AKIJ(I,J)*DLOG(TR))
	      HELP=DSQRT(A(I,I)*A(J,J))
	      A(I,J)=HELP*XKIJ1
!	      if (NTEMP == 0) GOTO 3
	      if (nTemp > 0) then

	        DKIJ=XKIJ(I,J)*AKIJ(I,J)/T
	        DADT(I,J)=DKIJ*HELP+XKIJ1/2.D0/HELP*(A(I,I)*DADT(J,J)+A(J,J)*DADT(I,I))
	        DADT(J,I)=DADT(I,J)

	      endif
!3	  A(J,I)=A(I,J)
	      A(J,I)=A(I,J)

30	    enddo

3	  enddo
	endif
!
	Delta(:NST,:NST) = 0.D0
	if (NST > 0) then
	  do i = 1,NST
	    do j = i,NST

	      Delta(i,j) = kappa(i,j)*(dexp(eps_R(i,j)/T) - 1.D0)
	      Delta(j,i) = Delta(i,j)
	      if (nTemp == 1) then

	        dDeldT(i,j) = -(Delta(i,j) + kappa(i,j))*eps_R(i,j)/T**2
	        dDeldT(j,i) = dDeldT(i,j)

	      endif

	    enddo
	  enddo
	endif
!
      RETURN
ENDsubroutine
!-----------------------------------------------------------------------
!
!	La subrutinas PARGR lee los par�metros grupales dispersivos y aso-
!	ciativos
!
!
!
SUBROUTINE PARGR (NC, NG, NGA, NST, inputFile, outputFile)
!
	IMPLICIT real*8(A-H,O-Z)
	parameter(NCM = 15, NGM = 15, NGAM = 12, NSM = 24)

	DIMENSION X(35),X1(35),X2(35),NNY(NCM,NGM+1),IGROUP(NGM),INT(2*(NGM+1)), &
     &          IDG(NCM,NGM+1),ID(NGM)
!
!	Energ��a y volumen de asociaci�n con 4 sub��ndices: s�lo est�n en
!	esta subrutina.
	dimension enAs(NSM,NGAM,NSM,NGAM),volAs(NSM,NGAM,NSM,NGAM),        &
     &          nyAss(NGAM,NCM),MAssoc(NGA)
	integer sigma, outputFile
	real*8 kappa

	COMMON/COORD/ZZ
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)

	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)

	common/NAME/AName(NCM)
!
!	Comprobaci�n de errores:
	if (NC > NCM) then

	  write (*, *) 'Error: n�mero de compuestos m�ximo sobrepasado'
	  stop

	elseif (NG > NGM) then

	  write (*, *) 'Error: n�mero m�ximo de grupos atractivos sobrepasado'
	  stop

	elseif (NGA > NGAM) then

	  write (*, *) 'Error: n�mero m�ximo de grupos asociativos sobrepasado'
	  stop

	endif
!
!  N�mero de coordinaci�n:
	ZZ=10.D0
!
!  Constante de los gases ideales
	R=82.05D0	!atm cm3/(mol K)
	NGROUP=35
!
      read (inputFile, *) NOWNII, NOWNIJ
!
	do I = 1,NCM

3	  IGROUP(I) = NGROUP + 1

	enddo
!
!  Lectura configuracion grupal atractiva por componente
      DO 1 I=1,NC

!	  JJ=20
	  JJ = NGM*2	!----> deber��a ser algo as��
! 	  write (*, *)  i
	  read (inputFile, *) (INT(J),J=1,JJ)
!	  read (inputFile, *)(int(j), j = 1, 2*NG)
	  NG1=NG+1

	  DO 2 J=1,NG1

	    JJ2=2*J
	    JJ1=JJ2-1
	    IDG(I,J)=INT(JJ1)  !identificacion grupo J perteneciente compuesto I
	    NNY(I,J)=INT(JJ2)  !cantidad grupo J perteneciente compuesto I

2	  enddo

	  II=1
	  IJ=1
!  Ordena en forma creciente el numero de identificacion de los grupos
!7	  if (IDG(I,IJ) == 0) GOTO 1
7	  if (idG(i,ij) /= 0) then

!	    do while (II <= NGM)

!	    if (IDG(I,IJ) > IGROUP(II)) GOTO 4
	      if (idG(i,ij) <= iGroup(ii)) then

!	      if (IDG(I,IJ) == IGROUP(II)) GOTO 5
	        if (idG(i,ij) /= iGroup(ii)) then

	          iTop = NGM - II !creo que deber��a ir esta l��nea en vez de la siguiente.
!	          ITOP=15-II
!	          if (ITOP == 0) GOTO 8
	          if (iTop /= 0) then

	            DO 6 J=1,ITOP

	              J1 = NGM - J     !creo que deber��a ir esta l��nea en vez de la siguiente.
!	              J1 = 15 - J
	              J2=J1+1
	              IGROUP(J2)=IGROUP(J1)

6	            enddo

	          endif
8	          IGROUP(II)=IDG(I,IJ)     !numero de identificacion del grupo

	        endif
5	        IJ=1+IJ

	      endif
4	      II=1+II
!	    if (II <= 15) GOTO 7
	    if (II <= NGM)goto 7 !creo que deber��a ir esta l��nea en lugar de la anterior
!	    enddo

	  endif

1	enddo
!
!	Se genera la configuracion grupal para cada componente teniendo en
!	cuenta los NG grupos
	DO 10 I=1,NC
	  DO 31 J=1,NG

	    NY(I,J)=0
	    JJ=1
!9	    if (IDG(I,JJ) == 0) GOTO 31
9	    if (idG(i,jj) /= 0) then

	      do while(jj <= NG)

!	      if (IDG(I,JJ) /= IGROUP(J)) GOTO 32
	        if (idG(i,jj) == iGroup(j)) then

!	          Cantidad del grupo J en el compuesto I
	          NY(I,J)=NNY(I,JJ)
!	          GOTO 31
	          exit

	        else

32	          JJ=JJ+1
!	          if (JJ <= NG) GOTO 9

	        endif

	      enddo

	    endif

31	  enddo
10 	enddo
!-----------------------------------------------------------------------
!
!	Lectura parametros grupales atractivos puros del banco de datos
	I=0
	J=1
	OPEN (UNIT = 10, FILE = 'FOR010.DAT', STATUS = 'OLD')
!	do while(iGroup(j) > i)

!	  I=I+1
!	  do while(j <= NG)

!11	    READ(10,101) (X(K),K=1,5)
11	    READ (10, *) (X(K),K=1,5)
	    i = i + 1
	  if (IGROUP(J) /= I) GOTO 11

	    ID(J)=I
	    Q(J)=X(1)
	    TSTR(J)=X(2)
	    GSTR(J)=X(3)
	    G1(J)=X(4)
	    G2(J)=X(5)
	    J=J+1
	    if (J <= NG) GOTO 11

!	  enddo

!	enddo
	CLOSE (UNIT = 10)
!-----------------------------------------------------------------------
!
!	Lectura parametros grupales atractivos puros ingresados por el usuario
	if (NOWNII > 0) then	!EQ.0) GOTO 12

	  DO 13 I = 1,NOWNII

!	    write (*, *)  NOWNII, i
	    read (inputFile, *) IOWN, (X(J), J = 1,5)
	    DO 1300 J=1,NG
	      if (IGROUP(J) == iOwn) then	!NE.IOWN) GOTO 13

	        Q(J)=X(1)		! area q del grupo
	        TSTR(J)=X(2)		! temperatura t* del grupo
	        GSTR(J)=X(3)		! parametro g* del grupo
	        G1(J)=X(4)		! parametro g' del grupo
	        G2(J)=X(5)		! parametro g'' del grupo

	      endif
1300	    enddo

13	  enddo

	endif
!-----------------------------------------------------------------------
!
!	Halla temperatura de spline con Newton
12	write (outputFile, 105) (J, J = 1,NC)
	write (outputFile, '(X,75("-"),10("-----"))')
	do I = 1,NG

	  if (G1(I) < -1.D-4) then	!GT.-1.D-4) GOTO 18

	    TR=1.D0-1.D0/G1(I)
	    IT=0
	    abs_AX = 1.
	    do while(abs_AX > 1.D-15)

!  Funci�n a ser igualada a 0 = AX(Tspl)
16	      AX=1.D0+G1(I)*(TR-1.D0)+G2(I)*DLOG(TR)-1.D0/4.D0
!  Derivada de dicha funci�n: d(AX)/d(Tspl)
	      DAX=G1(I)+G2(I)/TR
!  Cociente de Newton:
	      DTR=-AX/DAX
	      TR=TR+DTR
	      IT=IT+1
!	    if (TR < 0.D0 .OR. IT > 10) GOTO 18
	      abs_AX=DABS(AX)
!	      if (AX > 1.D-15) GOTO 16
	      if (TR >= 0D0 .AND. it <= 10) then

	        TSPL(I)=TR*TSTR(I)
	        EPX(I)=-4.D0*(G1(I)*TR+G2(I))
!	        GOTO 15

	      else

! Si Ts est� dando < 0 o no converge, Ts = 1000.
	        TSPL(I)=1.D3
	        epx(i) = 0
	        exit

	      endif

	    enddo

	  else
!
!  Si g' > 0, no puede haber una Ts > 0 (o dif��cilmente la haya)
18	   TSPL(I)=1.D3
	   epx(i) = 0

	  endif
15	  write (outputFile, 106)ID(I), TSTR(I), Q(I), GSTR(I), G1(I), G2(I), Tspl(i), epx(i), (NY(J,I), J = 1, NC)

	enddo
14	continue
!-----------------------------------------------------------------------
105	FORMAT(///, X, 'Pure group atractive parameters:', /, 78X, 'Group configuration per component', /, 3X, &
&	       'GN', 3X, 'T*(K)', 5X, 'q', 3X, 'g*(atm cm6/mol2)', 5X, 'g`', 8X, 'g"', 6X, "Ts(K)", 4X,        &
&	       "exp", X, <NC>I5)
106	FORMAT(1X, I4, 1X, F7.2, F8.4, 4X, F10.1, 2X, 2F10.5, 3X, F7.2, 2X, F5.2, 10I5)

!  Lectura parametros de interaccion binaria grupales atractivos del banco de datos
      I=0
      J=0
!
!          -------        -------       -------
!	do while(j <= NG)

29	    J=J+1
!
	if (J > NG) GOTO 28

!	  do while(i <= iGroup(j))	! > i)

20	    I=I+1
	    II=NGROUP+1-I
	    OPEN (UNIT=11, FILE='FOR011.DAT', STATUS='OLD')
	    OPEN (UNIT=12, FILE='FOR012.DAT', STATUS='OLD')
	    OPEN (UNIT=13, FILE='FOR013.DAT', STATUS='OLD')
	    READ(11,111) (X(K),K=1,II)
	    READ(12,111) (X1(K),K=1,II)
	    READ(13,111) (X2(K),K=1,NGROUP)

	  if (IGROUP(J) /= I) GOTO 20
!	  enddo
	  XKIJ(J,J)=X(1)
	  AKIJ(J,J)=0.D0
	  L=0
	  K=0

27	    L=L+1

	  if (L > NG) GOTO 29

25	      K=K+1
	    if (IGROUP(L) /= K) GOTO 25

	    CABS=DABS(88.888D0-X2(K))
	    if (CABS < 1.D-3) X2(K)=0.D0
	    ALFA(J,L)=X2(K)
!
!          -------        -------       -------
!
	  if (K <= I) GOTO 27

	  KK=K-I+1
	  CABS=DABS(88.888D0-X(KK))
	  IF (CABS > 1.D-3) GOTO 400
	  XKIJ(J,L)=1.D0
	  AKIJ(J,L)=0.D0
	  write (outputFile, 130) IGROUP(J),IGROUP(L)
	  GOTO 401
400	  XKIJ(J,L)=X(KK)
	  AKIJ(J,L)=X1(KK)
401	  AKIJ(L,J)=AKIJ(J,L)
	  XKIJ(L,J)=XKIJ(J,L)
	  GOTO 27

28	CONTINUE
	CLOSE (UNIT=11)
	CLOSE (UNIT=12)
	CLOSE (UNIT=13)
!-----------------------------------------------------------------------
!
!	Lectura parametros de interaccion binaria grupales atractivos
!     ingresados por el usuario
      if (NOWNIJ > 0) then	!EQ.0) GOTO 33
	  DO 21 I=1,NOWNIJ

	    read (inputFile, *) IDI,IDJ,XK,AK,ALF1,ALF2
	    DO 22 J=1,NG

	      if (IGROUP(J) == IDI) INI=J
	      if (IGROUP(J) == IDJ) INJ=J

22	    enddo
	    XKIJ(INI,INJ)=XK			! (k*)ij
	    XKIJ(INJ,INI)=XK			! (k*)ji
	    AKIJ(INI,INJ)=AK			! (k')ij
	    AKIJ(INJ,INI)=AK			! (k')ji
	    ALFA(INI,INJ)=ALF1		! (alfa)ij
	    ALFA(INJ,INI)=ALF2		! (alfa)ji

21	  enddo
	endif
!-----------------------------------------------------------------------
!
!	Impresi�n de par�metros binarios
   33 write (outputFile, 115)
	do i = 1,NG

	  write (outputFile, 116) (XKIJ(I,J),J=1,NG)

	enddo
	write (outputFile, 117)
	do i = 1,NG

	  write (outputFile, 116) (AKIJ(I,J),J=1,NG)

	enddo
	write (outputFile, 118)
	DO 23 I=1,NG

	  write (outputFile, 116) (ALFA(I,J),J=1,NG)

23	enddo
!-----------------------------------------------------------------------
!
!	Lectura de par�metros asociativos ingresados por el usuario.
!	(Te doy los "asociativos del banco de datos, toditito el mont�n
!	a fin de a��o")
!
!	...mentira (30/12/2011)
	NST = 0
	if (NGA > 0) then
!
!	  Vector M = sitios de asociaci�n
	  read (inputFile, *)(MAssoc(I),I=1,NGA)
	  do i = 1,NGA - 1
	    if (MAssoc(i) > MAssoc(i+1)) then

	      write(*,'("Error: los grupos asociativos deben ingresarse en orden creciente")')
	      write (*, *) Massoc(:NGA)
!	      stop

	    endif
	  enddo
!
!	  Lectura de energ��a y volumen de asociaci�n en una misma l��nea, aprovechando las simetr��as.
	  do i = 1,NGA
	    do j = i,NGA
	      do k = 1,Massoc(i)
	        if (Massoc(i) == Massoc(j)) then
	          do l = k,Massoc(j)
!
	            read (inputFile, *)enAs(k,i,l,j),volAs(k,i,l,j)
!	            Simetr��as
	            enAs(l,i,k,j) = enAs(k,i,l,j)
	            enAs(l,j,k,i) = enAs(k,i,l,j)
	            enAs(k,j,l,i) = enAs(k,i,l,j)
	            volAs(l,i,k,j) = volAs(k,i,l,j)
	            volAs(l,j,k,i) = volAs(k,i,l,j)
	            volAs(k,j,l,i) = volAs(k,i,l,j)
!
	          enddo
	        else
	          do l = 1,Massoc(j)
!
	            read (inputFile, *)enAs(k,i,l,j),volAs(k,i,l,j)
!	            Simetr��as
	            enAs(l,i,k,j) = enAs(k,i,l,j)
	            enAs(l,j,k,i) = enAs(k,i,l,j)
	            enAs(k,j,l,i) = enAs(k,i,l,j)
	            volAs(l,i,k,j) = volAs(k,i,l,j)
	            volAs(l,j,k,i) = volAs(k,i,l,j)
	            volAs(k,j,l,i) = volAs(k,i,l,j)
!
	          enddo
	        endif
	      enddo
	    enddo
	  enddo
!
!  El siguiente mamarracho surge de que no me las ingeni� para leer
!  la configuraci�n grupal asociativa y pasarla a notaci�n de sitios
!  corrida (1..NST) en un s�lo lazo.
!		Francisco, 16/07/2011 durante v-1.8.1 (alfa 1)

!  Conversion de variables a notacion por sitio "m" (y no, por sitio
!  "k" en el grupo "i")
	  m1 = 0
	  m2 = 0
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      m2 = m1
	      m1 = m1 + 1      	!Defino actual sitio "k" como sitio "m1"
	      do j = i,NGA      !esta definicion solo es valida gracias a que hasta ahora, solo ha habido
	        if (i == j) then  !un tipo de sitio en cada grupo.
	          do l = k,Massoc(j)

	            m2 = m2 + 1 !defino actual sitio "l" como sitio "m2"
	            eps_R(m1,m2) = enAs(k,i,l,j)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = volAs(k,i,l,j)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        else
	          do l = 1,Massoc(j)

	            m2 = m2 + 1  !defino actual sitio "l" como sitio "m2"
	            eps_R(m1,m2) = enAs(k,i,l,j)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = volAs(k,i,l,j)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        endif
	      enddo

	    enddo
	  enddo
	  NST = sum(Massoc(:NGA))
	  if (NST /= m1) then

	    write (outputFile, '(/,"Error: NST <> m1, verificar",/)')

	  endif
!
        m1 = 0
	  DO I=1,NGA

	    read (inputFile, *) (NYASS(I,J),J=1,NC)

	  enddo
!  Reasignaci�n en notaci�n 1..NST:
	  do i = 1,NC

	    m1 = 0
	    do j = 1,NGA
	      do k = 1, MAssoc(j)

	        m1 = m1 + 1
	        sigma(m1,i) = nyAss(j,i)

	      enddo
	    enddo

	  enddo
!
!  Impresion parametros asociativos
	  write (outputFile, 217) NGA
217	  FORMAT(//,X,"Number of associating groups =",I3)
!
	  write (outputFile, 218)
218	  FORMAT(/,X,'Associating group configuration ',/,26X,'Compound No')

	  write (outputFile, '(" Group    No of sites ",<NC>(I4,X))')((I), I=1,NC)
!
	  write (outputFile, '(X,25("-"),<NC>(5("-")))')
	  DO I=1,NGA

	    write (outputFile, 219) I,MAssoc(I), (nyass(I,J), J=1,NC)
219	    FORMAT(2X,I3,8X,I3,6X,<NC>(I4,X))

	  END DO
!
	  write (outputFile, 220)
 220	  FORMAT(/,X,'Associating energy matrix, [epsilon/R] (K)',/)

!  Escritura de par�metros energ�ticos de asociaci�n en forma de matriz
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write (outputFile, 2210)((enAs(k,i,l,j),l = 1,MAssoc(j)),j=1,NGA)

	    enddo
	  enddo
2210	  format(3X,'|',<NST>(X,F7.1),' |')
!  Escritura del volumen de asociaci�n
	  write (outputFile, '(/," Associating volume matrix, [kappa] (cm3/mol)",/)')
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write (outputFile, 2220)((volAs(k,i,l,j),l = 1,Massoc(j)),j=1,NGA)

	    enddo
	  enddo
2220	  format(3X,'|',<NST>(X,F7.5),' |')

	endif
!-----------------------------------------------------------------------


100	FORMAT(20I3)
101	FORMAT(5F10.2)
102	FORMAT(I3,5F10.2)
111	FORMAT(8F7.2)
116	FORMAT(3X,'|',<NG>F8.3,' |')
115	FORMAT(/,1X,'Binary interaction parameter at T* matrix, [k*]',/)
117	FORMAT(/,1X,'Binary interaction parameter temperature dependence matrix , [k`]',/)
118	FORMAT(/,1X,'Binary dumping factors matrix, [alpha]',/)
119	FORMAT(2I3,4F10.2)
120	FORMAT(1X,90("-"))
130	FORMAT(/,1X,'Los par�metros de interacci�n binaria entre los ','grupos',I2,' y ',I2,' no est�n disponibles.',/,1X,'Se usar�n', &
     &        'valores est�ndar a menos que estos sean ingresador por el usuario')
	RETURN
ENDsubroutine
!-----------------------------------------------------------------------
!
!	La subrutinas PARMOL lee los par�metros moleculares de cada compo-
!	nente:
!	Tc, en K
!	Pc, en atm
!	omega (factor ac�ntrico
!	Tsat, en K
!	Psat, en atm
!	dc, en cm/mol^1/3
!	PM (peso molec), en g/mol o lb/mol-lb
!	c (Peneloux), en cm3/mol
!	A, B, C, D, E, F (coef de Passut&Danner para Hgi, en BTU/mol-lb/�R
!
!
SUBROUTINE PARMOL (NC, NG, NST, inputFile, outputFile)
!
	IMPLICIT real*8 (A-H, O-Z)
	parameter (NCM = 15, NGM = 15, NGAM = 12, NSM = 24)
	character*10 AName
!
	DIMENSION TSAT(NCM), PS(NCM), Z(2), OMEGA(NCM), NYOLD(NGM), EXPOL(NGM)
	DIMENSION ID(NGM),GSTOL(NGM),G1OL(NGM),G2OL(NGM),TSTROL(NGM),      &
     &          TSPLOL(NGM),EPXOL(NGM),QOL(NGM),XKIJOL(NGM,NGM),         &
     &          XK1OL(NGM,NGM),ALFOL(NGM,NGM)
!
	common/versSUB/versGCA
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
!
	common/NAME/AName(NCM)
	common/GCPROM/pmm(ncm),Pen(NCM),HHA(ncm),HHB(ncm),hhc(ncm),HHD(ncm),HHE(NCM),HHF(NCM),     &
&	              HHG(NCM)
	common/ENTALP/iEntalp

	dimension eps_R_bak(NST,NST)
	real*8 kappa_bak(NST,NST),kappa
	integer sigma, sigma_bak(NST,NC), outputFile
!
!  Iniciaci�n de constantes de Antoine
	HA(:NC) = 0D0
	HB(:NC) = 0D0
!
!	write (outputFile, '(//,"Di�metros cr�ticos",/)')
	read (5,*) nUnit    !nUnit 0 = unidades en MPa, nUnit=1 en atm
	DO I=1,NC

! 	  write (*, *) i
	  read (inputFile, *) ANAME(I)
	  read (inputFile, *) TC(I), PC(I), OMEGA(I), TSAT(I), PS(I), DC(I)
	  read (inputFile, *) PMM(I), Pen(i), HHA(I), HHB(I), HHC(I), HHD(I), HHE(i), HHF(i), HHG(i)
	  if (nUnit == 0) then

	    PC(I) = PC(I)/0.101325D0
	    Ps(i) = Ps(i)/0.101325D0

	  endif
	  if (dc(I) == 0) then

	    DC(I)=7.337449944D0*TC(I)/PC(I)
	    DC(I)=DC(I)**.3333333333D0

	  end if
!	  write (outputFile, *) AName(i),' dc = ',dc(i),'cm/mol^(1/3)'

	enddo
1	continue
!
!  Guarda la informaci�n del componente "1", porque para ajustar el dc(i)
!  sobrescribe esta informaci�n, debido al que'l traspase de la misma se
!  hace mediantes "COMMON"s y no mediante entrada de subrutinas.
	DO 8 I=1,NG

	  NYOLD(I)=NY(1,I)

8	enddo
!
!
!  Almacenamiento de variables grupales asociativas
!  Entiendo que es un derroche talvez, almacenar TODO el paquete de
!  variables, pero tambi�n es muy sencillo... y seguro.
!			Francisco, 18/07/2011, v-1.8.5
	sigma_bak(:NST,:NC) = sigma(:NST,:NC)
	kappa_bak(:NST,:NST) = kappa(:NST,:NST)
	eps_R_bak(:NST,:NST) = eps_R(:NST,:NST)
!
!  Almacenamiento de variables moleculares (t�rmino repulsivo)
	TCOLD=TC(1)
	PCOLD=PC(1)
	DCOLD=DC(1)
	DO 3 I=1,NC

! 	write (*, *)  i
!  Ajusta el dc(i) para reproducir el dato de presi�n de vapor
	  NG1=0
	  if (I == 2) DCOLD=DC(1)
	  DO 400 J=1,NG

!  Almacenamiento de variables grupales atractivas
	    if (NY(I,J) > 0) then	!EQ.0) GOTO 400

	      NG1=NG1+1
	      NY(1,NG1)=NY(I,J)
	      ID(NG1)=J
	      GSTOL(NG1)=GSTR(NG1)
	      G1OL(NG1)=G1(NG1)
	      G2OL(NG1)=G2(NG1)
	      TSTROL(NG1)=TSTR(NG1)
	      TSPLOL(NG1)=TSPL(NG1)
	      EXPOL(NG1)=EPX(NG1)
	      QOL(NG1)=Q(NG1)
	      GSTR(NG1)=GSTR(J)
	      G1(NG1)=G1(J)
	      G2(NG1)=G2(J)
	      TSTR(NG1)=TSTR(J)
	      TSPL(NG1)=TSPL(J)
	      EPX(NG1)=EPX(J)
	      Q(NG1)=Q(J)

	    endif

400	  enddo
        DO 401 J=1,NG1

	    IJ=ID(J)
	    DO 4010 K=1,NG1

	      IK=ID(K)
	      XKIJOL(J,K)=XKIJ(J,K)
	      XKIJ(J,K)=XKIJ(IJ,IK)
	      XK1OL(J,K)=AKIJ(J,K)
	      AKIJ(J,K)=AKIJ(IJ,IK)
	      ALFOL(J,K)=ALFA(J,K)
	      ALFA(J,K)=ALFA(IJ,IK)

4010	    enddo

401	  enddo
!
!
!  Traspase a la posici�n "1" de par�metros asociativos
	  m = 0
	  do k = 1,NST

	    if (sigma_bak(k,i) /= 0) then

	      n = m
	      m = m + 1
	      sigma(m,1) = sigma_bak(k,i)

	      do l = k,NST

	        if (sigma_bak(l,i) /= 0) then

	          n = n + 1
	          kappa(m,n) = kappa_bak(k,l)
	          kappa(n,m) = kappa(m,n)
	          eps_R(m,n) = eps_R_bak(k,l)
	          eps_R(n,m) = eps_R(m,n)

	        endif

	      enddo

	    endif

	  enddo
	  NST1 = m !?????   Hasta ac�... 01:00 hrs. v1.8.2
!
!
!  Carga las propiedades del compuesto "i" en la primera posici�n.
!  Par�metros moleculares y repulsivo:
	  TC(1)=TC(I)
	  PC(1)=PC(I)
	  DC(1)=DC(I)
	  if ((TSAT(I) > 1.D-1) .OR. (omega(i) > 1.d-8)) then ! GOTO 7

	    T=TSAT(I)
	    P1=PS(I)
	    PFIT=P1
	    if (Tsat(i) < 1D-7) then
!	    if (TSAT(I) > 1.D-1) GOTO 7


	      T=TC(I)*.7D0
	      P1=PC(I)*DEXP(-(1.D0+OMEGA(I))*2.3026D0)
	      PFIT=P1

	    endif
7	    CALL PARAGC(T,1,NG1,NST1,0)
	    IGP=1
	    IGZ=0
!---------Lazo de ajuste del dc----------------------------------------------------
	    fCrt = 1.0
	    do while(fCrt >= 1.E-8)

5	      D1=D(1)
	      CALL PSAT(T,P1,Z,1,NG1,NST1,1,IGP,IGZ,IER)
	      IGZ=5
	      F1=P1-PFIT
	      FCRT=DABS(F1)

!	    if (FCRT < 1.D-8) GOTO 4

	      P2=P1
	      D(1)=D1+1.D-4
	      CALL PSAT(T,P2,Z,1,NG1,NST1,1,IGP,IGZ,IER)
	      F2=P2-PFIT
	      DFDD=(F2-F1)/1.D-4
	      DPDD=(P2-P1)/1.D-4
! 	      write (*, *) d(1),dFdd
	      DEL=-F1/DFDD
	      DLC=DEL/D(1)
	      if (DLC > .05D0) DEL=.05D0*D(1)
	      if (DLC < -.05D0) DEL=-.05D0*D(1)
	      D(1)=D1+DEL
	      P2=P1+DEL*DPDD
	      Z(2)=Z(2)/P1*P2
	      P1=P2

!	    GOTO 5
	    enddo
!-----------------------------------------------------------------------
!
!	    Nuevo valor del di�metro cr��tico del compuesto "i":
4	    DC(I)=D(1)/(1.D0-.12D0*DEXP(-2.D0*TC(I)/3.D0/T))/1.065655D0
	    TIN=TSAT(I)
	    PIN=PS(I)

	  endif
	  if (TSAT(I) < 1.D-1) TIN=.7D0*TC(I)
	  if (TSAT(I) < 1.D-1) PIN=PC(I)*DEXP(-(1.D0+OMEGA(I))*2.3026D0)
	  HB(I)=DLOG(PC(I)/PIN)/(1.D0/TC(I)-1.D0/TIN)
	  HA(I)=DLOG(PC(I))-HB(I)/TC(I)
!
!	  Reasignaci�n de par�metros grupales atractivos:
	  DO 402 J=1,NG1

	    GSTR(J)=GSTOL(J)
	    G1(J)=G1OL(J)
	    G2(J)=G2OL(J)
	    TSTR(J)=TSTROL(J)
	    TSPL(J)=TSPLOL(J)
	    EPX(J)=EXPOL(J)
	    Q(J)=QOL(J)
	    DO 4020 K=1,NG1

	      XKIJ(J,K)=XKIJOL(J,K)
	      AKIJ(J,K)=XK1OL(J,K)
	      ALFA(J,K)=ALFOL(J,K)

4020	    enddo

402	  enddo

3	enddo	!CONTINUE
!
!  Reasignaci�n de configuracion asociativa y par�metros
	sigma(:NST,:NC) = sigma_bak(:NST,:NC)
	kappa(:NST,:NST) = kappa_bak(:NST,:NST)
	eps_R(:NST,:NST) = eps_R_bak(:NST,:NST)
!
!  Reasignaci�n de par�metros moleculares del comp. 1:
      TC(1)=TCOLD
      PC(1)=PCOLD
      DC(1)=DCOLD

!  Reasignaci�n de configuraci�n grupal dispersiva del comp. 1
!      DO 9 I=1,NG
      do i = 1,NG
!    9 NY(1,I)=NYOLD(I)
	  ny(1,i) = nyOld(i)

	enddo

      if (nUnit == 1) then

	  write (outputFile, 10)

	elseif (nUnit == 0) then

	  write (outputFile, 12)

	endif
	write (outputFile, 13)
	do i = 1,NC

	  write (outputFile, 11) AName(i),TC(I),PC(I),OMEGA(I),TSAT(I),PS(I),DC(I),PMM(i),Pen(i)

	enddo
!
!  Constantes de Passut-Danner: impresi�nn s�lo si al menos se le in-
!  gresa valores a un compuesto:
	iEntalp = NC
	do i = 1,NC
	  if ((HHA(i) == 0) .AND. (HHB(i) == 0) .AND. (HHC(i) == 0) .AND. (HHD(i) == 0) .AND. (HHE(i) == 0) .AND. (HHF(i) == 0)) then

	    iEntalp = iEntalp - 1

	  endif
	enddo
	if (iEntalp /= 0) then

	  write (outputFile, '(/," Passut-Danner constant for ideal gas enthalpy calculation:",/,"   Compound   A(BTU/lb)  B(BTU/(lb �R))  ", &
     &             "C(BTU/(lb �R2)�1E3  D(BTU/(lb �R3))�1E6  E(BTU/(lb �R4))�1E10  F(BTU/(lb �R5))�1E14",/,128("-"))')
	  do i = 1,NC

	    write (outputFile, 14)AName(i),HHA(i),HHB(i),1.D3*HHC(i),1.D6*HHD(i),1D10*HHE(i),1D14*HHF(i)

	  enddo

	endif

   91 FORMAT(8F10.4)
   10 FORMAT(//,X,'Pure component properties',//,3X,'Compound',6X,'Tc(K)',1X,'Pc(atm)',3X,'omega',2X,'Tsat(K)',2X,'Psat(atm)',2X, &
     &       'dc(cm/mol^1/3)',3X,"M(g/mol)",3X,"c(cm3/mol)")
   11 FORMAT(4X,A10,1X,2F7.1,3X,F6.4,F9.2,2X,F8.4,5X,F8.4,7X,F6.2,8X,F6.2)
   12 FORMAT(//,2X,'Pure component properties',/,3X,'Nro.',2X,'Tc(K)',2X,'Pc(MPa)',3X,'omega',2X,'Tsat(K)',2X,'Psat(MPa)',2X,  &
     &       'dc(cm/mol^1/3)',X,"M(g/mol)",X,"c(cm3/mol)")
13	format("  ",100("-"))
14	format(3X,A10,X,F9.6,3X,F9.6,10X,F9.6,12X,F9.6,13X,F9.6,13X,F9.6)
	RETURN
ENDSUBROUTINE PARMOL
!-----------------------------------------------------------------------
!
!
!
SUBROUTINE PSAT (T, PS, Z, NC, NG, NST, IC, IGP, IGZ, IER)
!
!  Esta subrutina provee las presiones de vapor de compuestos puros.
!  Los volumenes de las fases son obtenidos llamando a ZMAX.
!
!  Mensajes de error (Valores de iEr)
!      1    : ning�n error
!     -1    : soluci�n trivial (vL = vV)
!
	IMPLICIT real*8(A-H,O-Z)
	parameter(NCM = 15, NGM = 15, NSM = 24)
	integer sigma
	real*8 LAMBDA
	DIMENSION X(NCM),Z(2)
	COMMON/SAT/PHI3
	COMMON/COORD/ZZ
	COMMON/EXTREM/IEXT(2)
!
!  Variables espec��ficas del t�rmino dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)

!  Versi�n asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!

!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),rNyT(NGM),QT,    &
     &           TETA(NGM),HELP12(2,NGM),HELP15(2,NGM),HELP7(2,NGM),     &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),XLAM1,XLAM2,XLAM3

!	Versi�n asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
	logical::fallo
!
!
	IGOLD=IGZ
	if (IGZ > 1) IGZ=-1
	EPSP = 1.D-10
	PI = 3.1415926536D0
	IER = 1
	x(:NC) = 0D0
	X(IC) = 1.D0
	if (IGP == 0) PS = DEXP(HA(IC) + HB(IC)/T)
!
!-----------------------------------------------------------------------
!
!  C�lculo de iterativo de la presion de vapor.
!  Utiliza el criterio de igualdad de energia libre (Gvap = Gl��q)
	Pcrt = 1.
!	do while(PCRT > EPSP)

3	  CALL ZMAX(Z,1,T,PS,X,IGZ,NC,NG,NST,IF,fallo)
	  if (IEXT(1) /= 0) then	!EQ.0) GOTO 20

	    PS=PS*.75D0
	    GOTO 3

	  endif
20	  if (IEXT(2) /= 0) then	!EQ.0) GOTO 21

	    PS=PS*1.25D0
	    GOTO 3

	  endif
21	  ZCR=DABS(Z(1)-Z(2))
!
!  Fase supercritica:
	  if (ZCR < 1.D-6) IER=-1
	  if (IER < 0) GOTO 99
!	  if (iEr >= 0) then
!
!  Contribuci�n repulsiva de la energia de Helm. (fase vapor y l��quida)
	    XSIV = PS/T/R/Z(1)*PI/6.D0*PHI3
	    XSIL = PS/T/R/Z(2)*PI/6.D0*PHI3
	    CSV=3.D0*XSIV/(1.D0-XSIV)+XSIV/(1.D0-XSIV)**2
	    CSL=3.D0*XSIL/(1.D0-XSIL)+XSIL/(1.D0-XSIL)**2
!
!  Contribuci�n atractiva de la energia de Helm. (fase vapor y l��quida)
	    ATTV=0.D0
	    ATTL=0.D0
	    DO 2 K=1,NG

	      XNY=DFLOAT(NY(IC,K))
	      if (XNY >= 1D-10) then	!LT.1.D-10) GOTO 2

	        ATTV=ATTV+XNY*Q(K)*HELP2(1,K)
	        ATTL=ATTL+XNY*Q(K)*HELP2(2,K)

	      endif

2	    enddo
	    ATTV=ATTV*ZZ/2.D0
	    ATTL=ATTL*ZZ/2.D0
!
!  Contribuci�n asociativa de la energia de Helm. (fase vapor y l��quida)
	    AasocL = 0D0
	    AasocV = 0D0
	    if (NST > 0) then

!	      Vapor:
	      AasocV =  dot_product(sm(:NST),dlog(Xs(1,:NST))) + (Sum_S - sum(sm_Xs(1,:NST)))*.5D0
!	      L��quido:
	      AasocL =  dot_product(sm(:NST),dlog(Xs(2,:NST))) + (Sum_S - sum(sm_Xs(2,:NST)))*.5D0

	    endif
!
	    F = CSV - ATTV - CSL + ATTL - DLOG(Z(1)/Z(2)) - AasocL + AasocV
	    PNY=-F/(Z(1)-Z(2))*PS
	    PCRT=DABS(PNY-PS)
	    PLIM=(PS-PNY)/PS
	    if (PLIM > .1D0) PNY=.9D0*PS
	    if (PLIM < -.1D0) PNY=1.1D0*PS
	    Z(2)=Z(2)*PNY/PS
	    PS=PNY

	if (PCRT > EPSP) GOTO 3
!	  else
!
!  Si la fase es supercritica sale del lazo.
!	    exit

!	  endif
99	  CONTINUE

!	enddo
	IGZ=IGOLD
	RETURN

ENDsubroutine PSAT
!-----------------------------------------------------------------------
!
!
SUBROUTINE ZMAX (Z, iTyp, T, P, XN, iGues, NC, NG, NST, iC,fallo)
!
! 	Esta subrutina encuentra el(los) factor(es) de compresibilidad del
!	gas y/o la fase l��quida a determinada composici�n, temperatura y
!	presi�n.
!	El(los) factor(es) de compresibilidad es(son) encontrado por solu-
!	ci�n iterativa de
!
!           P(EXP)=-RT(D/DV(Ares/RT)-rho)
!
!     La soluci�n es devuelta en Z.
!
!     Si iGuez > � = 1, el usuario debe proporcionar un buen valor ini-
!	cial del factor de compresibilidad deseado (S). De lo contrario se
!	genera autom�ticamente otra estimaci�n inicial.
!
!	Las opciones siguientes especifican la soluci�n deseada
!
!	Opciones
!     iTyp   iGues
!      -1  cualquiera   Fase gaseosa
!       1    0 OR 1     Fases l��quida y gaseosa
!       1      >1       Fase l��quida
!       0       0       Fases l��quida y gaseosa verificando energ��a libre
!                       de Gibbs m��nima
!
!     Mensajes de error (valores de iC)
!      1    : fase correcta
!     -1    : fase incorrecta
!     Si iType = 0 el valor devuelto es el que de menor energ��a libre de
!	Gibbs
!
!-----------------------------------------------------------------------
	IMPLICIT real*8 (A-H, O-Z)

	parameter(NCM = 15, NGM = 15, NSM = 24)
	parameter(maxit=200)
!
	DIMENSION Z(2), XN(NC), X(NCM)
	DIMENSION AA(NGM,NGM), IER1(2)
	logical calc

	real*8 LAMBDA6, LAMBDA
!
!  Variables auxiliares de asociaci�n
	dimension PSI_V(NSM)
	integer sigma
!
!  Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM), HELP4(2,NGM), HELP5(2,NGM), rNyT(NGM), QT,    &
     &           TETA(NGM), HELP12(2,NGM), HELP15(2,NGM), HELP7(2,NGM),     &
     &           E(2,NGM,NGM), PREP(2), DPDV(2), XLAM1, XLAM2, XLAM3
!
!  Variables espec��ficas del t�rmino dispersivo
	COMMON/GROUP2/Q(NGM), A(NGM,NGM), DADT(NGM,NGM), ALFA(NGM,NGM), R,     &
     &              NY(NCM,NGM)
!
!  Propiedades moleculares
	COMMON/MOL/DC(NCM), D(NCM), DT(NCM), HA(NCM), HB(NCM)
	COMMON/CRIT/TC(NCM), PC(NCM)
!
	COMMON/COORD/ZZ
	COMMON/SAT/PHI3
	COMMON/EXTREM/IEXT(2)
!
!  Versi�n asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM), sm(NSM), dXs_dV(2,NSM)

!  Versi�n asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM), dDeldT(NSM,NSM)

	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!
	common/itervolumen/nroIter(2)
	common/GCPROM/PMM(NCM), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM)
	logical::fallo 
	integer::contfallo
	contfallo=0
	fallo = .false.
!
	IOPT=ITYP
	IMIN=0
	if (ITYP == 0) IMIN=1
	if (ITYP == 0) IGUES=0
	if (ITYP == 0) IOPT=1
	IEXT(1)=0
	IEXT(2)=0
	IER1(1)=0
	IER1(2)=0
	ICYC=0
	EPS=1.D-16
	PI=3.1415926536D0
	RT=R*T
!-----------------------------------------------------------------------
!
!	Cantidades independientes del volumen
!
!	Las siguientes variables sirven a la subrutina VIRIAL:
	xNtot = sum(xN(:NC))
	x(:NC) = xN(:NC)/xNtot	!fraccion molar componente i
!
!	Contribuci�n atractiva
!	Fraccion molar grupo k
	rNyT(:NG) = matmul(x(:NC),dfloat(ny(:NC,:NG)))
!	�rea total grupos
	qT = dot_product(rNyT(:NG),q(:NG))
!	Fraccion area grupo k
	teta(:NG) = q(:NG)/qT*rNyT(:NG)
!
!	Contribuci�n asociativa
	if (NST > 0) then
!
!	  matmul: subrutina intr��nseca de FORTRAN para multiplicar matrices:
	  sm(:NST) = matmul(dfloat(sigma(:NST,:NC)),x(:NC))
	  if (NST == 2) then
	    if (sm(2) > sm(1)) then
	  
	      major = 2	!this is only for the Analytical solution of 2 cross associating sites
	        
	    else
	
	      major = 1
	      
	    endif
	  endif
!
	endif
	Sum_S = sum(sm(:NST))
!
!     Contribucion repulsiva
	XLAM1=0.D0
	XLAM2=0.D0
	XLAM3=0.D0
	DO 3 I=1,NC

	  DIA=D(I)
	  DLA=X(I)*DIA
	  XLAM1=XLAM1+DLA        !lambda 1
	  DLA=DLA*DIA
	  XLAM2=XLAM2+DLA        !lambda 2
	  XLAM3=XLAM3+DLA*DIA    !lambda 3

3	enddo
	PHI3=XLAM3
!
!  Coeficiente virial: iniciaci�n del volumen del gas, y algunas va-
!	riables que necesita ZMAX. 
	calc = .TRUE.
	call GCAvirial(NC, NG, NST, T, xn(:NC), x(:NC), B, Bvl, Bat, Bas, calc)
!	write (*, *) B,Bvl,Bat,Bas
!  Inicializacion de las fracciones no asociadas y deriv respecto de V
	Xs(:2,:NST) = 1D0
	dXs_dV(:2,:NST) = 0.D0
	nroIter = 0
!-------------------------------------------------------------------------------------
!
!  Inicialici�n de los volumenes de cada fase
	if (iGues >= 0 .AND. iGues <= 1) then

	  if (iOpt /= -1 .OR. iGues /= 1) then

	    Z(1)=1.D0	           		   !Gas ideal para la fase vapor
	    if (B > 0 .OR. (P < -RT/4D0/B)) then	!m�ximo en presi�n para B < 0

!  Virial para la fase vapor: Z = 1 + B/v (m�s exacto que "1 + B*P/RT"... supuestamente)
	      Z(1) = 1D0 + 2D0*B*P/(R*T + dsqrt(RT**2 + 4D0*RT*B*P))
!
	    endif
!	    write (*, *) B,Z(1)

	  endif
6 	  if (iOpt /= 1 .OR. iGues /= 1) then

	    VTOT=0.D0
	    DO 9 I=1,NC

	      if (X(I) >= 1E-2) then

	        TR=T/TC(I)
	        TCRIT=.99D0*TC(I)
	        if (T > TCRIT) TR=1.D0
	        FAC=1.D0
	        if (T < TCRIT) EXPO=(1.D0-TR)**.2857D0
	        if (T < TCRIT) FAC=.29D0**EXPO
!  Racket para el vol. molar de cada comp.
	        V=.29D0*TC(I)*R/PC(I)*FAC
	        VTOT=VTOT+V*X(I)

	      endif

9	    enddo
!  Vtot = Vtot - dot_product(x(:NC),Pen(:NC))
	    Z(2)=P*VTOT/RT

	  endif

	endif
!-------------------------------------------------------------------------------------
!
!  C�lculo de los volumenes
4	ICYC=ICYC+1	 	! icyc=1 vapor,icyc=2 liquido

	IER=0
	if (IOPT == 1 .AND. IGUES > 1) ICYC=2
	NTRIAL=0
	NHELP=1
	NIT=0
	ROLD=P/Z(ICYC)/RT
	RO=ROLD
	zeta3 = phi3*pi*ro/6
	if (zeta3 > 0.9999) then
	  
!  El volumen inicial es menor que el "covolumen" de la ecuaci�n de Carnahan y Starling
	  ro = 0.9999*6/phi3/pi

	endif
	if (NST > 2 .OR. NST == 2 .AND. (Delta(1,1) /= 0 .OR. Delta(2,2) /= 0)) then

!  Inicio de las fracciones no asociadas por 3 iteraciones de sustici�n directa
	  do i = 1,3
	    do j = 1,NST

	      Xs(iCyc,j) = 1D0/(1D0 + sum(ro*sm(:NST)*Xs(iCyc,:NST)*Delta(j,:NST)) )

	    enddo
	  enddo

	endif
!
!
8	CONTINUE

	NIT=NIT+1
24	CONTINUE
    contfallo = contfallo+1
    if(contfallo > 10000) then
  !      fallo = .true.
       ! return
    endif
!  Contador de iteraciones para l��quido y vapor:
	nroIter(iCyc) = nroIter(iCyc) + 1
	DNOM=RO*QT/RT
	DFVDV=0.D0
	DFV=0.D0
!
	dFVas(iCyc) = 0.D0
	dFVVas = 0.D0
!
!  Contribuci�n atractiva
	DO 35 K=1,NG
		            !Seg�n paper J�rgensen, estas son algunas de las variables:
	  HEL2K=0.D0	!H2(k) = Sum(theta(j)*tau(j,k)*Qt*g(k,j)*rho/TR, j=1,NG)
	  HEL4K=0.D0	!H4(k) = Sum(theta(j)*tau(j,k), j=1,NG)
	  HEL5K=0.D0	!H5(k) = Sum(theta(j)*tau(j,k)*(g(j,k)*Qt*rho/RT)*alfa(j,k)*Dg(j,k)*Qt*rho/RT
	  HEL15K=0.D0
	  HEL12K=0.D0
	  HEL7K=0.D0	!H7(k) = SUM(
	  AKK=A(K,K)
	  DO 5 J=1,NG

	    ARGU=ALFA(J,K)*(A(J,K)-AKK)*DNOM
!	    write (*, *) A(J,K),AKK,DNOM
	    AAH=A(J,K)*DNOM
	    AA(J,K)=AAH
	    E(ICYC,J,K)=DEXP(ARGU)
	    IF (TETA(J) >= 1.D-16) then

	      TETAE=TETA(J)*E(ICYC,J,K)
	      HEL4K=HEL4K+TETAE
	      TETAA=TETAE*AAH
	      HEL2K=HEL2K+TETAA
	      TETAA=TETAA*ARGU
	      HEL5K=HEL5K+TETAA
	      TETAA=TETAA*ARGU
	      HEL15K=HEL15K+TETAA
	      TETAE=TETAE*ARGU
	      HEL12K=HEL12K+TETAE
	      TETAE=TETAE*ARGU
	      HEL7K=HEL7K+TETAE
	      
	    endif

5	  enddo
	  HEL2K=HEL2K/HEL4K
	  HEL5K=HEL5K/HEL4K
	  HEL7K=HEL7K/HEL4K
	  HEL12K=HEL12K/HEL4K
	  HEL15K=HEL15K/HEL4K
	  IF (RNYT(K) >= 1.D-16) then

	    DFV=DFV+RNYT(K)*Q(K)*(HEL5K+HEL2K-HEL2K*HEL12K)
	    DFVDV = DFVDV - RO**2*RNYT(K)*Q(K)*((2.D0*HEL2K*HEL12K**2 - HEL2K*(HEL7K + 2.D0*HEL12K) - 2.D0*HEL5K*HEL12K) + HEL15K + 2.D0*HEL5K)*ZZ/2.D0
	    
	  endif
400	  HELP2(ICYC,K)=HEL2K              !h2k/h4k
	  HELP4(ICYC,K)=HEL4K              !h4k
	  HELP5(ICYC,K)=HEL5K              !h5k/h4k
	  HELP7(ICYC,K)=HEL7K
	  HELP12(ICYC,K)=HEL12K            !h6k/h4k
	  HELP15(ICYC,K)=HEL15K
	  
35	enddo
!
!     DFV y DFVDV son la primera y segunda derivada de Ar/RT con respec al volumen
!
	DFV=RO*DFV*ZZ/2.D0		! (dFdV)att
	DFVDV=DFVDV-RO*DFV*2.D0		! (d2FdV2)att
!
!	Contribuci�n repulsiva
	PI6=PI/6.D0*XLAM3
	Y=1.D0/(1.D0-PI6*RO)
	DYDV=-Y**2*PI6*RO**2
	DYDVV=2.D0/Y*DYDV**2-2.D0*RO*DYDV
!
!     DFVR es la contribuci�n repulsiva de DFV
	DFVR=(3.D0*XLAM1*XLAM2/XLAM3+XLAM2**3/XLAM3**2*(2.D0*Y-1.D0-1.D0/Y)+1.D0/Y)*DYDV
	PREP(ICYC)=-RT*DFVR		!P_rep
!
	DFV=DFV+DFVR			! (dFdV)att+rep
!
!	Ac� se agrega la contribuci�n respulsiva.
	DFVDV=DFVDV+DFVR*DYDVV/DYDV+(XLAM2**3/XLAM3**2*(2.D0+1.D0/Y**2)-1.D0/Y**2)*DYDV**2		 ! (d2FdV2)att+rep
!
!	Contribuci�n asociativa
	if (NST > 0) then
!
	  if (NST == 1) then
	    
!	    Analytical solution for 1 self-associating site (1A)
	    if (Delta(1,1) > 0) then
	    
	      s_Delta_V = sm(1)*Delta(1,1)*ro
	      root(iCyc) = dsqrt(4D0*s_Delta_V + 1D0)
	
	      Xs(iCyc,1) = 2D0/(1D0 + root(iCyc))
! 	      dXs_aux(iCyc,1) = -Xs(iCyc,1)/(1D0 + root)/root  !auxiliary variable for derivative calculations
	      dXs_dV(iCyc,1) = -2D0*Xs(iCyc,1)**2 /root(iCyc)*s_Delta_V*ro

	    endif
	    sm_Xs(iCyc,1) = sm(1)*Xs(iCyc,1)
	    
	  elseif (NST == 2 .AND. Delta(1,1) <= 0 .AND. Delta(2,2) <= 0) then
	    
!         Two cross associating sites, with same or different mole amounts.
	    if (Delta(1,2) > 0) then
	
	      s1_Delta_V = sm(3-major)*Delta(1,2)*ro
	      s2_Delta_V = sm(major)*Delta(1,2)*ro
	      b_aux(iCyc) = 1D0 + s1_Delta_V - s2_Delta_V
	      root(iCyc) = dsqrt( b_aux(iCyc)*b_aux(iCyc) + 4D0*s2_Delta_V )
	
	      Xs(iCyc,major) = 2D0/(b_aux(iCyc) + root(iCyc))

! 	      dXs_aux(iCyc,major) = -Xs(iCyc,major)/(b_aux + root)   !auxiliary variable for derivative calculations	      
	      dXs_dV(iCyc,major) = -Xs(iCyc,major)**2 * (1D0 - b_aux(iCyc) + (b_aux(iCyc)*(1D0 - b_aux(iCyc)) - 2D0*s2_Delta_V)/root(iCyc))*ro/2D0
	      if (sm(major) == sm(3-major)) then
	  
	        Xs(iCyc,3-major) = Xs(iCyc,major)
	        dXs_dV(iCyc,3-major) = dXs_dV(iCyc,major)
	  
	      else
	  	      
	        Xs(iCyc,3-major) = 1D0/(1D0 + s2_Delta_V*Xs(iCyc,major))
! 	        Xs(iCyc,3-major) = (1d0 - Xs(iCyc,major))/Xs(iCyc,major)/s1_Delta_V	
	        dXs_dV(iCyc,3-major) = Xs(iCyc,3-major)**2 * s2_Delta_V*(dXs_dV(iCyc,major) - Xs(iCyc,3-major)*ro)
	      
	      endif
	      sm_Xs(iCyc,major) = sm(major)*Xs(iCyc,major)		
	      sm_Xs(iCyc,3-major) = sm(3-major)*Xs(iCyc,3-major)		
	      
	    endif
	  
	  else
	    
! 	    General case for NST sites with any associating scheme:
	    
!	    Iniciaci�n de las fracciones no asociadas por continuacion. La extrapolaci�n
!	    est� acotada ente 0 y 1:
	    do i = 1,NST

	      auxXs = Xs(iCyc,i) + dXs_dV(iCyc,i)*(1.D0/RO - 1.D0/ROld)
	      if (auxXs > 1D0) then
         
	        Xs(iCyc,i) = 1D0
         
	      elseif (auxXS <= 0.D0) then
         
	        Xs(iCyc,i)= Xs(iCyc,i)/5d0
         
	      else
         
	        Xs(iCyc,i) = auxXs
         
	      endif
         
	    enddo
!  Llamado a calculo en para obtener fracciones no asociadas.
	    call OptiNewton (maxIt, NST, sm(:NST), Delta(:NST,:NST), RO, Xs(iCyc,:NST), H(iCyc,:NST,:NST), indx(iCyc,:NST), in)
! 	  
!  C�lculo de la derivada de la fracci�n no asociada respecto VOLUMEN
!
!	    Number of moles of non-bonded sites
	    sm_Xs(iCyc,:NST) = sm(:NST)*Xs(iCyc,:NST)
	    
!	    Vector -dg/dv provisionally stored in dXs_dV:
	    dXs_dV(iCyc,:NST) = -sm(:NST)*matmul( Delta(:NST,:NST), sm_Xs(iCyc,:NST) )*ro*ro

! 	    do i = 1, NST	    
! 	      sm_Xs(iCyc,i) = 0.D0
! 	      dXs_dV(iCyc,i) = 0.D0
! 	      if (sm(i) >= 1D-16) then
! 	  	
! 	        sm_Xs(iCyc,i) = sm(i)*Xs(iCyc,i)
! 	        do j = 1, NST	  
! 	          if (sm(j) >= 1D-16 .AND. Delta(i,j) > 1D-16) then
! 	        
! 	            dXs_dV(iCyc,i) = dXs_dV(iCyc,i) - sm(i)*sm(j)*Xs(iCyc,j)*Delta(i,j)
! 	          
! 	          endif
! 	        enddo
! 	        dXs_dV(iCyc,i) = dXs_dV(iCyc,i)*ro*ro
! 	        
! 	      endif
! 	    enddo
	    
! 	    call LUDcmp (H(iCyc,:NST,:NST), NST, NST, indx(iCyc,:NST), d1)
	    call LUBksb (H(iCyc,:NST,:NST), NST, NST, indx(iCyc,:NST), dXs_dV(iCyc,:NST))
!
	  endif

!	  Association contribution to d(Ar/RT)/dV y d2(Ar/RT)/dV2:
	  dFVas(iCyc) = 0.5D0*ro*(Sum_S - sum(sm_Xs(iCyc,:NST)))
	  dFVVas = -(dFVas(iCyc) + 0.5D0*dot_product(sm(:NST), dXs_dV(iCyc,:NST)))*ro
!
	endif
!
!	Agregado de la contribuci�n asociativa
	dFV = dFV + dFVas(iCyc)
	dFVdV = dFVdV + dFVVas
!-------------------------------------------------------------------------------------
!
	DPDV(ICYC) = RT*(1.D0 + DFVDV/RO/RO)	! -V**2 * dP/dV
	G = RO*R*T - R*T*DFV - P
	DGRO = R*T*(1.D0 + DFVDV/RO**2)
!
!  C�lculo de (ro)j+1 por el metodo de Newton
!
	IF (DGRO > 0.D0 .OR. NIT <= 1) GOTO 26

!	  At present iteration, dP/drho < 0, which means mechanical unstability. Rho will be restored and
!	  calculation will proceed with a half Newton step.
	  rOld = ro
	  RO   = RO-DRO
	  DRO  = .5D0*DRO
	  GOTO 25

26	IF (NHELP > 1) GOTO 20

!	  Si la fase es l��quida salta
	  IF (ICYC == 2) GOTO 21

	  IF (DGRO > 0.D0) GOTO 20

	    NTRIAL = NTRIAL + 1
	    if (NTRIAL  >  25) IER=100

	  IF (IER /= 0) GOTO 10

	    SRC=RO-G/DGRO

	  if (SRC < 0.D0) GOTO 10

	    RO=RO*.90D0
	    GOTO 24

21	if (IGUES > 2 .OR. IGUES < 0) GOTO 20
   
        if (NTRIAL == 0) GOTO 22
        
        Gcrit=G*GOld
        
        if (G > GOLD .AND. DGRO < 0.D0) GOTO 10
        
        if (GCRIT < 0.D0) GOTO 23

22	GOLD=G
	ROLD=RO
	if (IGUES /= 0 .OR. G < 0.D0) GOTO 27

	  SRO = Z(1) - P/RT/(RO - G/DGRO)

      if (SRO < 0.D0) GOTO 10
27	GC = DABS(G)
	if (GC < EPS) then

      goto 200 !exit Newton loop if g = 0, since the solution has been found...

    endif
	rOld = ro
      RO = RO - (.04D0*G/GC*RO) !rho - 0.04g/|g|�rho = rho�(1 - 0.04 g/|g|) = si g > 0, rho disminuye un 4%; si g < 0, rho aumenta un 4%.
      NTRIAL=NTRIAL+1
      
      GOTO 24 !arriba hasta el comienzo, sin sumar iteraci�n.

   23 NHELP=NHELP+1
	rOld2 = rOld
	rOld  = ro
!  Extrapolation of rho:
!	RO    = -G*(RO - ROLD)/(G - GOLD) + RO
!
!  Modifiqu� esto para  que ro - rOld sea <> 0
	RO    = -G*(RO - ROLD2)/(G - GOLD) + RO
!  Creo que deber�a ser de esta manera:	
! 	RO    = -Gold*(RO - ROLD2)/(G - GOLD) + RO	
      GOTO 24

      
20	NHELP = NHELP + 1
!  Paso del m�todo de Newton tradicional:
	DRO=-G/DGRO
!  Paso m�ximo permitido = 5% del valor de rho para gas:      
	DLTRO=.05D0*RO
!  20% del valor de rho para l�quido:      
	if (ICYC < 2) DLTRO=4.D0*DLTRO
	if (DRO > DLTRO) DRO=DLTRO
	if (DRO < -DLTRO) DRO=-DLTRO

25	RCR=DABS(DRO)

!  Criterio de salida:
	if (RCR < EPS .OR. NIT > 20) GOTO 200
!
!	Almacenamiento del valor anterior de la densidad...   y posterior actualizacion
	rOld = ro
	RO=RO+DRO

	GOTO 8
!
!	Culmin� el calculo del volumen de la fase iCyc
200	if (NIT < 20) GOTO 10
	IER=100*ICYC
10	if (IGUES > 1) GOTO 110

	  PCRIT=1.D3*DABS(G)
	  if (PCRIT > P) IEXT(ICYC)=1
	  if (DGRO < 0.D0 .OR. IEXT(ICYC) == 1) IER = 100*(ICYC - 1) + 99
	  
110	IER1(ICYC)=IER
	Z(ICYC)=P/RO/R/T
	if (IGUES > 1) GOTO 113
	if (IOPT == 1 .AND. ICYC == 1) GOTO 4
	if (IOPT == -1 .AND. ICYC == 2) GOTO 115
	if (IOPT == -1 .AND. IER /= 0) GOTO 4
	if (IOPT == 1) GOTO 112
	IF (IER == 0) GOTO 113
115	IER=99
	Z(1)=Z(2)
	GOTO 113
112	if (IER /= 0) GOTO 114
	if (IER1(1) /= 0) GOTO 113
	ZCRT = DABS(Z(1) - Z(2))
!	if (ZCRT > 1.D-10) GOTO 210
	if (Zcrt > 1.D-8) goto 210
	IER1(2)=199
114	Z(2)=Z(1)

210	IF (IMIN == 0 .OR. IER1(2) /= 0) GOTO 113
!-----------------------------------------------------------------------
!	C�lculo de la energia residual de Gibbs para ambas fases.
!	  r     r
!	 G     A
!	--- = --- + Z - 1
!	R T   R T
!
      GIB1=0.D0
      GIB2=0.D0
!
!	Contribucion atractiva: A_att
	DO 202 K=1,NG

	  GIB1=GIB1+RNYT(K)*Q(K)*HELP2(1,K)
	  GIB2=GIB2+RNYT(K)*Q(K)*HELP2(2,K)

202	enddo

	GIB1=-GIB1*ZZ/2.D0       ! (A)att vapor
	GIB2=-GIB2*ZZ/2.D0       ! (A)att liq
!
!	Contribuci�n asociativa: A_assoc
	if (NST > 0) then

!	  Vapor:
	  Gib1 = Gib1 + dot_product(sm(:NST), dlog(Xs(1,:NST))) + (Sum_S - sum(sm_Xs(1,:NST)))*.5D0	!(A)att + (A)ass vapor
!	  L��quido:
	  Gib2 = Gib2 + dot_product(sm(:NST), dlog(Xs(2,:NST))) + (Sum_S - sum(sm_Xs(2,:NST)))*.5D0	!(A)att + (A)ass l��quido

	endif
!
!	Contribuci�n repulsiva
	Y1 = 1.D0/(1.D0 - PI6*(P/Z(1)/RT))
	R1 = XLAM1*XLAM2/XLAM3*3.D0
	R2 = XLAM2**3/XLAM3**2
!
	GIB1 = GIB1 + R1*(Y1 - 1.D0) + R2*(Y1**2 - Y1 - DLOG(Y1)) + DLOG(Y1) + Z(1) - DLOG(Z(1))	         ! Gres/RT+1=Ares/RT+Z-lnZ vapor
	GIB2 = GIB2 + R1*(Y  - 1.D0) + R2*(Y**2  - Y  - DLOG(Y))  + DLOG(Y)  + Z(2) - DLOG(Z(2))	         ! Gres/RT+1=Ares/RT+Z-lnZ l��quido
!
	IC=1
	IF (GIB1 < GIB2) IC=-1
	
	IPHAS = 1 + (IC + 1)/2
	GOTO 201
113	CONTINUE
	IC=1
	if (IER1(1) /= 0 .AND. IOPT == -1) IC=-1
	if (IER1(2) /= 0 .AND. IOPT == 1) IC=-1
201	CONTINUE
!----------------------------
	RETURN
ENDsubroutine ZMAX
!---------------------------------------------------------------------------------------------
!
!	Esta subrutina obtiene los coeficientes de actividad de una mezcla
!	(preferentemente l��quida) de NC componentes.C
!	Se basa en la f�rmula
!
!	                 ------
!	ln gamma(i) = ln phi(i) - SUM( x(i)*ln phi(i), i=1,NC )
!
!	donde:
!	       ------
!	       phi(i) = coeficiente de fugacidad del comp. i en la mezcla
!
!	       x(i)   = fracci�n molar de i
!
!	       phi(i) = coeficiente de fugacidad de i puro a la misma T y P
!
subroutine GAMMAGC (NC, NG, NST, iTyp, iGz, iTemp, n, T, P, lnGamma, dlnGammadT, Z)
!
	implicit real*8 (A-H, O-Z)
!
	parameter (NCM=15,NGM=15,NSM=24)
	real*8 lnGamma(NCM), n(NCM), Nt
	dimension phi(NCM),dLPhi(NCM,NCM),dLPhiT(NCM),dLPhiP(NCM),dlnGammadT(NCM),x(NCM),aux(NCM)
	logical::fallo
!
	iDer = 0
!	write(3,*)ityp,igz
!
!	Normalizaci�n
	Nt = sum(n(:NC))
!	x(:NC) = n(:NC)/Nt
!
!	Coeficientes de fugacidad de la mezcla:
	call GCEOS (NC, NG, NST, iDer, iTemp, T, P, n, phi, dLPhi, dLPhiT, DLPHIP, Z, IGz, iTyp, IC,fallo)
!
!	Iniciaci�n para c�lculo:
	lnGamma(:NC) = phi(:NC)
	dlnGammadT(:NC) = dLPhiT(:NC)
!
!	Lazo de c�lculo
	do i = 1, NC
!
!	  C�lculo de los coeficientes de fugacidad de los puros:
	  x(:NC) = 1.D-20
	  x(i) = 1D0
	  call GCEOS (NC, NG, NST, iDer, iTemp, T, P, x, phi, dLPhi, dLPhiT, DLPHIP, Zpure, IGz, iTyp, IC,fallo)
	  lnGamma(i) = lnGamma(i) - phi(i)
	  dlnGammadT(i) = dlnGammadT(i) - dLPhiT(i)
!
	enddo
!
!
	return
endsubroutine GAMMAGC
!
!---------------------------------------------------------------------------------------------
!
! 
! 
!	La subrutina VIRIAL halla el 2� coeficiente virial de una mezcla
!
!
subroutine GCAvirial (NC, NG, NST, T, n, x, B, Bvl, Bat, Bas, calc)

	implicit real*8 (a-h, o-z)
	parameter(NCM = 15, NGM = 15, NSM = 24, PI = 3.1415926536D0)
	real*8 n(NC),nt,lambda1,lambda2,lambda3
	dimension x(NC),aux(NSM)
	logical calc
	integer sigma
!
!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),xg(NGM),QT,      &
     &           theta(NGM),HELP6(2,NGM),HELP11(2,NGM),HELP12(2,NGM),    &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),lambda1,lambda2,lambda3
!
	COMMON/COORD/zz
!
!	Variables espec��ficas del t�rmino dispersivo
	COMMON/GROUP2/Q(NGM),g(NGM,NGM),DgDT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)
!
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
!
!	Versi�n asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Versi�n asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!	�Se calcularon variables que aprovecha ZMAX?
	if (.NOT. calc) then
!
	  nt = sum(n(:NC))
	  x(:NC) = n(:NC)/nt
	  xg(:NG) = matmul(x(:NC),dfloat(ny(:NC,:NG)))
	  QT = dot_product(q(:NG),xg(:NG))
	  theta(:NG) = q(:NG)/QT*xg(:NG)

	  lambda1 = dot_product(x(:NC),d(:NC))
	  lambda2 = dot_product(x(:NC),d(:NC)**2)
	  lambda3 = dot_product(x(:NC),d(:NC)**3)

	  if (NST > 0) sm(:NST) = matmul(dfloat(sigma(:NST,:NC)),x(:NC))

	endif
!
!	Contribuci�n atractiva:
	aux(:NG) = matmul(g(:NG,:NG),theta(:NG))
	Bat = -zz/2D0*QT**2*dot_product(theta(:NG),aux(:NG))/R/T
!
!	Contribuci�n repulsiva
	Bvl = PI*(lambda1*lambda2/2D0 + lambda3/6D0)
!
!	Contribuci�n asociativa
	Bas = 0D0
	if (NST > 0) then

	  aux(:NST) = matmul(sm(:NST),Delta(:NST,:NST))
	  Bas = -dot_product(aux(:NST),sm(:NST))/2D0/sum(sm(:NST)*sm(:NST))

	endif
!
	B = Bvl + Bat + Bas
!	write(*,'(<NG>("| ",<NG>(F12.0,2X)," |",/))')((g(i,j),j=1,NG),i=1,NG)
!	write (*, *) QT
!	write (*, *) theta(:NG)
!	write (*, *) B,Bvl,Bat,Bas
!
endsubroutine GCAvirial
!
!---------------------------------------------------------------------------------------------
!
!
!	OptiNewton halla el valor m�ximo de la funci�n Q :
!
!	                          Qmax = A/RT|assoc,eq.
! 
!---------------------Francisco, 14/07/2011------durante v-1.8.0---------
!	C�lculo de la fraccion no asociada, como se muestra en el libro
!	"Thermodynamic models: Fundamentals & computational aspects" de Mollerup & Michelsen y
!	en M.Michelsen Ind.Eng.Chem.Res. 45(2006)8449-8453.
!
!	Siendo Lm�x = Aas/RT|eq:
!		cuando grad(Lmax,X) = 0, se obtiene la formula general de c�l-
!		culo de las fracciones no asociadas.
!
!	Siendo que se trata de una optimizacion de L iterando sobre X, se re-
!	suelve por Newton-Raphson.
!           H (X(u + 1) - X(u)) + g = 0
!
!	Pasos:
!     calcular/estimar X
!     Calcula g
!     calcula H
!     resolver DX: X = X + DX
!   1 si L no aumenta,
!           X = X + alfa�DX    #alfa = 0.5
!           ir a 1
!     fin
!
!	El c�lculo de las derivadas con respecto a V de la fracci�n no aso-
!	ciada Xs, est�n basadas en la modificaci�n de A.E.Andreatta (ver Tesis
!	PhD y T.M.Soria y col. Fluid Phase Equilib. 302(2011)1-9) al traba-
!	jo de S.P.Tan y col. Ind.Eng.Chem.Res. 43(2004)203-208.
subroutine OptiNewton (maxit, NTS, sm, Delta, rho, X, H, indx, in)
!
	implicit real*8 (a-h, o-z)
!
!  N�mero m�ximo de iteraci�nes
!  parameter(maxit = 1000)
!
!  Par�metros de entrada de asociaci�n
	dimension delta(NTS,NTS),sm(NTS)
!
!  Par�metros de salida e intermedios
	dimension X(NTS), Xnew(NTS), g(NTS), H(NTS,NTS), minXindx(1), indx(NTS)
!
!  Par�metros auxiliares
	dimension Xold(NTS), DX(NTS), A(NTS,NTS), gnew(NTS), Hnew(NTS,NTS)
!
!  C�lculo inicial.
	call Qfunction(NTS,X,sm,rho,Delta,Q,g,H)
	im = 0
!
!  Comienzo del lazo de convergencia
	do in = 1,maxit
!
!  Valores de variables auxiliares
	  im = 1 + im
! 	  A = H
	  DX = -g
!
!  Resolviendo sistema de ecuaciones lineales
	  call LUDcmp (H, NTS, NTS, indx, d)
	  call LUBksb (H, NTS, NTS, indx, DX)
!	  Reiniciando alfa
	  alpha = 1.D0
!
!	  Calculando nuevo X
1	  do i = 1,NTS
	    Xnew(i) = X(i) + alpha*DX(i)
!
!	    Si una fracci�n es negativa de disminuye pero no se hace 0
	    if (Xnew(i) <= 0D0) then
	      Xnew(i) = 2.D-1*X(i)
	    endif
	  enddo
!
!	  Salida del lazo de convergencia si el error es menos que una tolerancia.
	  if (maxval(abs(alpha*DX)) <= 1.d-16)exit
!
!	  Nuevo valor de Q:
2	  call Qfunction (NTS, Xnew, sm, rho, Delta, Qnew, gnew, Hnew)
!
!	  Aument� Q?
	  if (Qnew > Q*(1D0 + 1.D-14)) then

3	    Q = Qnew
	    X = Xnew
	    g = gnew
	    H = Hnew

	  else

!	    Si no lo hizo, disminuye a la 1/3 del paso.
	    alpha = alpha/3.
	    im = im + 1
	    goto 1

	  endif
	enddo
!
	return
	stop
100	format(5(F10.8,2X))
	endsubroutine
!---------------------------------------------------------------------------------------------
!
!
!	La subrutina Qfunction calcula el valor de la funci�n Q, su gra-
!	diente y Hessiano para el c�lculo de segundo orden de la fracci�n
!	no asociada, seg�n lo sugerido por M. L. Michelsen, IECR 2002, 45,
!	8449-8453.
!
subroutine Qfunction(NTS,X,sm,rho,Delta,Q,g,H)
!
!	Q = Sum(sm(k)�(ln X(k) - X(k) + 1) - 1/2�Sum(Sum(sm(k)�sm(l)�Delta(k,l)�X(l)/V),l=1,NST),k=1,NST)
!
!	  = Q1 + Q2
!
!	g(k) = sm(k)/X(k) - sm(k) - Sum(sm(k)*sm(l)*Delta(k,l)*X(l)/V),l=1,NST)
!
!	  = sm(k)/X(k) - sm(k) - SUMA
!
!	"H"(k,l) = -(sm(k) + SUMA)/X(k)*d(k,l) - sm(k)*sm(l)*Delta(k,l)/V
!
!	Donde: sm(k)      = es la cantidad de moles de sitios "k"
!	       X(k)       = es la fracci�n de sitio "k" no asociada
!	       Delta(k,l  = fuerza de asociaci�n entre sitios "k" y "l"
!	       NST        = n�mero de sitios totales
!	       d(k,l)	= delta de Kronocker, = 1 si "l" = "k", sino = 0.
!	       g(k)       = gradiente de Q en la direcci�n "k"
!	      "H"(k,l)   = hessiano k,l
!	nota1: tomo a rho = 1/V donde V = vol total. Asume que son iguales porque los moles est�n normalizados (creo).
!
!	Francisco, febrero de 2010.
!
!	---------------------------------------
!	Siendo que cada vez que se corta un c�lculo, el programa corta casi siempre en esta subrutina.
!	Eso indica que esta es un cuello de botella, por eso, comprim�� el c�lculo de las variables
!	Q, g y H en la menor cantidad de lazos posible. v-1.9.18
!
	implicit real*8 (a-h,o-z)
!	implicit real*16 (a-h,o-z)
	dimension X(NTS),sm(NTS),Delta(NTS,NTS),g(NTS),H(NTS,NTS)
!
!	Iniciaci�n
	Q1 = 0.0D0
	Q2 = 0.0D0
	g = 0.0D0
	H = 0.0D0
!
!	C�lculo de Q1, gradiente y hessiano:
	do k = 1,NTS

	  do l = 1,NTS

	    if (Delta(k,l) > 0) then
	      H(k,l) = -sm(k)*sm(l)*Delta(k,l)*rho
!	      Nota: g(k) a�n NO ES el gradiente_k. Fran 18/02/2010.
	      g(k) = g(k) - H(k,l)*X(l)
	      Q2 = Q2 + 5D-1*H(k,l)*X(k)*X(l)
	    endif

	  enddo
	  H(k,k) = H(k,k) - (sm(k) + g(k))/X(k)
!	  Para evitar elementos nulos en la diagonal.
	  if (H(k,k) == 0.0D0)H(k,k) = 1.0D-20
	  Q1 = Q1 + sm(k)*(dLOG(X(k)) - X(k) + 1.D0)
!
!	  Ahora s��, g es el gradiente_k:
	  g(k) = sm(k)/X(k) - sm(k) - g(k)

	enddo
!
!	C�lculo de Q:
	Q = Q1 + Q2
!
	return
endsubroutine
!--------------------------------------------------------------------------------
!
!
!	Halla la soluci�n del sistema
!	                                 U � x = b
!	                                 =   -   -
!
!	por sustitici�n hacia atr�s, complementando la subrutina LUDcmp que
!	descompone a A en U*L.
!	indx es el vector de pivoteo para b.
!
SUBROUTINE LUBksb (A, n, np, indx, b)
	INTEGER n, np, indx(n)
	real*8 a(np,np),b(n)
	INTEGER i,ii,j,ll
	real*8 sum
	ii=0
	do i = 1, n

	  ll=indx(i)
	  sum=b(ll)
	  b(ll)=b(i)
	  if (ii /= 0) then

	    do j=ii,i-1

	      sum=sum-a(i,j)*b(j)

	    enddo

	  else if (sum /= 0.D0) then

	    ii=i

	  endif
	  b(i)=sum

	enddo
	do i = n, 1 , -1

	  sum=b(i)
	  do j=i+1,n
	    sum=sum-a(i,j)*b(j)
	  enddo
	  b(i)=sum/a(i,i)

	enddo
	return
ENDsubroutine
!--------------------------------------------------------------------------------
!
!	LUDcmp halla la matriz U tal que
!
!	                         A = U � L
!	                         =   =   =
!	Donde U es una matriz "diagonal superior" y L "diagonal inferior
!	con sus elementos diagonales unitarios".
!
SUBROUTINE LUDcmp(a, n, np, indx, d)
	INTEGER n, np, indx(n), NMAX
	real*8 d,a(np,np),TINY
	PARAMETER (NMAX=500,TINY=1.0e-20)
	INTEGER i,imax,j,k
	real*8 aamax,dum,sum,vv(NMAX)
      
	d = 1.D0
	do i=1,n

	  aamax=0.D0
	  do j=1,n

	    if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))

	  enddo
	  if (aamax == 0.D0) then

	    write(*,'(/, "Matriz singular en LUDcmp", /)')
! 	    stop

	  endif
        vv(i) = 1.D0/aamax
	enddo
      do j = 1, n

	  do i = 1, j - 1

	    sum=a(i,j)
	    do k=1,i-1

	      sum=sum-a(i,k)*a(k,j)

	    enddo
	    a(i,j)=sum

	  enddo
	  aamax = 0.D0
	  do i = j, n

	    sum = a(i,j)
	    do k=1,j-1

	      sum = sum - a(i,k)*a(k,j)

	    enddo
	    a(i,j)=sum
	    dum=vv(i)*abs(sum)
	    if (dum >= aamax) then

	      imax=i
	      aamax=dum

	    endif

	  enddo
	  if (j /= imax) then

	    do k=1,n

	      dum=a(imax,k)
	      a(imax,k)=a(j,k)
	      a(j,k)=dum
	      
	    enddo
	    d=-d
	    vv(imax)=vv(j) 

	  endif
	  indx(j)=imax
	  if (a(j,j) == 0.D0)a(j,j)=TINY
	  if (j /= n) then

	    dum=1.D0/a(j,j)
	    do i=j+1,n

	      a(i,j)=a(i,j)*dum

	    enddo
        
	  endif

	enddo
!
      return
ENDsubroutine LUDcmp
!
!--------------------------------------------------------------------------------
!===================================================================  
    subroutine GC_Adapt(compuesto,comp,  NG, nu, q, Tstr, gstr, gdot, gddot,kstrT,kdotT,alphaT,gruposgc)
!------------------------------------------------------------------- 
!-------------------------------------------------------------------  
    use Input 
    use CONSTANTES  
    use SubGrupos
    implicit none
    integer,dimension(DiffStructGroups,2),intent(in)::compuesto
    integer,intent(in)::comp
    integer,dimension(15,15),intent(inout)::nu
    real*8,dimension(15),intent(inout)::q, Tstr, gstr, gdot, gddot
    real*8,dimension(15,15),intent(inout)::kstrT,kdotT,alphaT
    integer,intent(inout)::NG,gruposgc(15)
    integer::i,j,k,l
    logical::VL
!COMMONS
    COMMON/INTER/A,kstr,kdot,alpha
    real*8::A(NMG,NMG),kstr(NMG,NMG),kdot(NMG,NMG),alpha(NMG,NMG)    
    

    i = 0
    k=1
    do while (compuesto(k,1) /= 0) 
        call BUSCARAS (compuesto(k,1),gruposgc,15,VL)
         if (.not.VL)then
            i = i+1
            gruposgc(ng+i) = compuesto(k,1)
            q(ng+i) = Obtain_Q(compuesto(i,1))
            Tstr(ng+i) = Obtain_Tstr(compuesto(i,1))
            gstr(ng+i) = Obtain_gstr(compuesto(i,1))
            gdot(ng+i) = Obtain_gdot(compuesto(i,1))
            gddot(ng+i) = Obtain_gddot(compuesto(i,1))

        endif
        do l=1,15
            if(compuesto(k,1)==gruposgc(l))exit
        enddo
        nu(comp,l) = compuesto(k,2)

        
        k = k+1
    enddo
    ng = ng + i

endsubroutine GC_Adapt  