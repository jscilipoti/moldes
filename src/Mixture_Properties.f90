SUBROUTINE MI_SOLV_BREAK (IPAREQ,TEMP1,SLV1,SOLACE,SLOST)
    use constantes
      PARAMETER (NSCM=10,NA=150)
      implicit real*8 (A-H,O-Z)
      LOGICAL SOLACE,NOACE
      COMMON/US/MS(NCOM,DiffStructGroups,2),nms(NCOM)
      DIMENSION X(NCOM),ACT(NCOM),DACT(NCOM,NCOM),TACT(NCOM)
      COMMON/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
 
      SOLACE = .FALSE.
      NDIF = 0
      NACT = 0
      IOUT = 5
      IF (TAZEO.NE.1) THEN
         NC=2
         CALL UNIPAR (NC,NG,TEMP1,MODEL,IOUT,NOACE,MS)
         X(2) = 0.0
         X(1) = 1.0
         CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
         ALFA2 = ACT(2)*PSAT2/(ACT(1)*PSAT1)
         X(2) = 1.0
         X(1) = 0.0
         CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
         ALFA1 = ACT(2)*PSAT2/(ACT(1)*PSAT1)
         IF (ALFA2.GT.1.AND.ALFA1.LT.1) THEN
            I0 = 1
            I1 = 2
            XP = X1AZEO
         ELSE
            I0 = 2
            I1 = 1
            XP = 1 - X1AZEO
         END IF
         NC = 3
         CALL UNIPAR (NC,NG,TEMP1,MODEL,IOUT,NOACE,MS)   
         X(I0) = 0.0
         X(I1) = 1.
         X(3) = 0.
         CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
         FUNA = PSAT2*ACT(2)/(PSAT1*ACT(1)) - 1.
         X(I1) = 0.
         X(3) = 1.
         CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
         FUNB = PSAT2*ACT(2)/(PSAT1*ACT(1)) - 1.
         BE = 1.
         IF (FUNA*FUNB.GT.0.) THEN
            X(3) = 0.1
            X(I1) = 1. - X(3)
 310        CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
            FUNB = PSAT2*ACT(2)/(PSAT1*ACT(1)) - 1.
            IF (DABS(1.-X(3)).LT.1.E-03) THEN
               GO TO 10000
            ELSE
               IF (FUNB.LT.0.) THEN
                  X(3) = X(3) + 0.1
                  X(I1) = 1. - X(3)
                  GO TO 310
               ELSE
                  BE = X(3)
               END IF
            END IF
         END IF
         AL = 0.
         X(3) = BE
5001     X3A = X(3)
         X(3) = (AL + BE)/2.
         X(I1) = 1.0 - X(3)
         CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
         FUN3 = PSAT2*ACT(2)/(PSAT1*ACT(1)) - 1.
         IF ((DABS(FUN3).GT.ERROR).OR.(DABS(X(3)-X3A).GT.ERROR)) THEN
            IF (FUN3*FUNA.GT.0.) THEN
               AL = X(3)
               FUNA = FUN3
               GO TO 5001
            ELSE
               BE = X(3)
               GO TO 5001
            END IF
         ELSE
            SLOST = X(3)
            Z = X(3) + X(I1) + X(I1)*XP/(1 - XP)
	    X(3) = X(3)/Z
	    X(I1)= X(I1)/Z
	    X(I0)= 1.0 - X(I1) - X(3)
            CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
	    IF (ACT(3)*X(3).GT.1) THEN
               GO TO 10000
            ELSE
               SOLACE = .TRUE.
	       GO TO 10000
            END IF
	 END IF
      ELSE
         NC = 3
         CALL UNIPAR (NC,NG,TEMP1,MODEL,IOUT,NOACE,MS)
         IF (SLV1.GT..5) THEN
	    X(2) = .4
	 ELSE
	    X(2) = SLV1
	 END IF
         X(1) = 0.0
	 X(3) = 1.0 - X(2)
         CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
	 IF (ACT(2)*X(2).GT.1.0) GO TO 10000
	 IF (SLV1.GT..5) THEN
	    X(1) = .4
	 ELSE
	    X(1) = SLV1
	 END IF
         X(2) = 0.0
         X(3) = 1 - X(1)
         CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
         IF (ACT(1)*X(1).GT.1.0) GO TO 10000
         SOLACE = .TRUE.
         SLOST = 1.
      END IF
10000 RETURN
ENDsubroutine MI_SOLV_BREAK 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Solvent_Properties(MS,MOP,TEMP1,SLC,SLVP,SLOST,CMDIST,RELVOL,NOACE)
    use Input
    use constantes
    implicit none
    integer,PARAMETER::NCM=15
    integer::i
    integer::ndif,nact,iout,mop,NC,NG
    real*8::SLC,SLVP,SLOST,CMDIST,temp1,x13,x12,relvol,psat2,psat1
    real*8:: X(NCOM),ACT(NCOM),DACT(NCOM,NCOM),TACT(NCOM),ACTAS(NCOM)
    integer::MS(NCOM,DiffStructGroups,2)
    integer::IC,phase_type,guess
    real*8::CZ
    real(8), dimension(15)     :: DLPHIP,phi,dLPhiT !dc, omega, Pc, Psat, Tb, Tc, Tsat, x
    real(8), dimension(15,15) :: dLPhi
    !real(8), dimension(2)       :: CZ !dP_dT, dP_dV, v, 
     
    COMMON/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
    COMMON/MOL/DC,D,DT,HA,HB
    real*8::DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
    real*8::tazeo,x1azeo,error
    LOGICAL NOACE,fallo
!C
      NDIF = 0
      NACT = 0
      IOUT = 5
      noace = .False.
      IF (MOP.EQ.1) THEN
!C
!C         LIQUID-LIQUID EXTRACTION
!C
        SLC = 0.	!selectividad
        SLVP = 0.	!poder solvente 
        SLOST = 0.	!pérdida de solvente
        CMDIST = 0.	!coeficiente de distribución
        NC = 3
        X(1) = 0.0
        X(2) = 0.0
        X(3) = 1.0
        if(model /= 3)then
            CALL UNIPAR (NC,NG,TEMP1,MODEL,IOUT,NOACE,MS)
            IF (NOACE) return
            call UNIFAC (MS,x,NC,Temp1,act)
        else
            call Charge_Commons_GC(ms,NG)
            if(isnan(DC(3)))then
                noace = .true.
                return
            endif
            CALL PARAGC(InputProblem01%T,3,NG,0,0)
            call GCEOS (3, NG, 0, 0, 0, InputProblem01%T, InputProblem01%P, X, phi, dLPhi, dLPhiT, DLPHIP, CZ, guess, 1, IC,fallo)
            act(1:3) = exp(phi(1:3))
        endif
         ! CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
        SLC = ACT(2)/ACT(1)
        SLVP = 1./ACT(1)
        X13 = ACT(1)
        X(1) = 0.
        X(2) = 1.
        X(3) = 0.
        if(model /= 3)then
            call UNIFAC (MS,x,NC,Temp1,act)
        else
            call Charge_Commons_GC(ms,NG)
            CALL PARAGC(InputProblem01%T,3,NG,0,0)
            call GCEOS (3, NG, 0, 0, 0, InputProblem01%T, InputProblem01%P, X, phi, dLPhi, dLPhiT, DLPHIP, CZ, guess, 1, IC,fallo)
            act(1:3) = exp(phi(1:3))
        endif
        !CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
        X12 = ACT(1)
        CMDIST = X12/X13 
        SLOST = 1./ACT(3)
    ELSEIF(MOP==2)THEN 
!C
!C      EXTRACTIVE DISTILLATION
!C
        RELVOL = 0.
        SLVP = 0.
        NC = 3
        CALL UNIPAR (NC,NG,TEMP1,MODEL,IOUT,NOACE,MS)
        IF (NOACE) return
        X(1) = 0.0
        X(2) = 0.0
        X(3) = 1.0
        CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
        RELVOL = ACT(2)*PSAT2/(ACT(1)*PSAT1)
        SLVP = 1./ACT(1)
    elseif(mop==3)then
        SLC = 0.	!selectividad
        SLVP = 0.	!poder solvente 
        SLOST = 0.	!pérdida de solvente
        CMDIST = 0.	!coeficiente de distribución
        NC = 3
        CALL UNIPAR (NC,NG,TEMP1,MODEL,IOUT,NOACE,MS)
        IF (NOACE) return
        X(1) = 0.0
        X(2) = 0.0
        X(3) = 1.0
        call UNIFAC (MS,x,NC,Temp1,act)
        ! CALL UNIFA (NDIF,NACT,NC,NG,TEMP1,X,ACT,DACT,TACT)
        SLC = ACT(2)/ACT(1)
        SLVP = 1./ACT(1)
    END IF
    RETURN
endsubroutine Solvent_Properties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Charge_Commons_GC(ms,NG)
    use PureProp
    use Input
    use subgrupos
    use constantes
    implicit none




!Variable de ENTRADA
    integer,dimension(NCOM,DiffStructGroups,2),intent(in)::MS
!Variable de SALIDA
    integer,intent(out)::NG
!Variables INTERNAS
    integer::i,j,k
    integer::NC
    real*8::TC,PC,Rk
    integer,dimension(15,15)::nu      
    real*8,dimension(15)::q, Tstr, gstr, gdot, gddot  ,Tb ,  TcGC, PcGC, BpGC, dcGC
    real*8::kstrT(15,15),kdotT(15,15),alphaT(15,15) 
    integer:: gruposgc(15)
!COMMONS
    integer, dimension(15,15)     :: ny
    real(8)                       :: dPV, dPDT, R, ZZ
    real(8), dimension(15)        :: d, dcrit, dt, dPDn, HA, HB, Pcrit, Tcrit
    real(8), dimension(15)        :: epx, gs, g1,g2, qArea, Ts, Tspl
    real(8), dimension(15,15)     :: a, dadt, akij, alfa, xkij
  
    COMMON /GROUP1/                GS, G1, G2, TS, TSPL, XKIJ, AKIJ, EPX
    COMMON /GROUP2/                Qarea, A, DADT, ALFA, R, NY
    COMMON /COORD/                 ZZ  
    COMMON /MOL/                   DCrit, D, DT, HA, HB
    COMMON /CRIT/                  TCrit, Pcrit
    common /ScndDer/               dPV, dPDT, dPDn  
    COMMON/INTER/Ai,kstr,kdot,alpha
    real*8::Ai(NMG,NMG),kstr(NMG,NMG),kdot(NMG,NMG),alpha(NMG,NMG)    


    NC = 3
    NG = 0
    nu = 0
    Tstr = 0
    gstr = 0
    gdot = 0
    gddot = 0
    kstrT = 0.
    kdotT = 0.
    alphaT = 0.      
    q = 0.
    R = 82.05D0 !atm·cm3/mol K
    ZZ = 10D0   !Número de coordinación  
    gruposgc=0
    do i=1,nc
        call GC_Adapt(ms(i,:,:),i, NG, nu, q, Tstr, gstr, gdot, gddot,kstrT,kdotT,alphaT,gruposgc)
        call Calc_Prop_Pure(MS(i,:,:),InputProblem01%ifam,InputProblem01%T,TCT = Tc, PCT = Pc)    
        Tcrit(i) = Tc
        Pcrit(i) = Pc
        ! dcrit(i) = (0.08943 * R * Tc / Pc)**(1./3.)   
        j=1 
        Rk=0.0
        do while (ms(i,j,1)/=0)
            Rk=Rk+ Obtain_R(ms(i,j,1))*ms(i,j,2)
            j=j+1
        enddo     
        dcrit(i) = 10 ** (0.4152 + 0.4128*log10(Rk)) !correlación de Espinosa et al. 2002
    enddo 
    i=1
    do while(gruposgc(i) /= 0)
        j = 1
        do while(gruposgc(j) /= 0)
            kstrT (i,j) = kstr(gruposgc(i),gruposgc(j))
            kstrT (j,i) = kstrT (i,j)

            kdotT (i,j) = kdot(gruposgc(i),gruposgc(j))
            kdotT (j,i) = kdotT (i,j)

            alphaT (i,j) = alpha(gruposgc(i),gruposgc(j))
            alphaT (j,i) = alpha(gruposgc(j),gruposgc(i))
            j = j+1
        enddo
        i=i+1
    enddo
! Asignación de los vectores COMMON:  
    do j = 1, NG
    
        gs(j) = gstr(j)
        g1(j) = gdot(j)
        g2(j) = gddot(j)
        Ts(j) = Tstr(j)
        qArea(j) = q(j)
        Tspl(j) = 1D3
        epx(j) = 0D0
        do k = j + 1, NG
      
            xkij(k,j) = kstrT(k,j)
            xkij(j,k) = kstrT(k,j)
            akij(k,j) = kdotT(k,j)
            akij(j,k) = kdotT(k,j)
            alfa(k,j) = alphaT(k,j)
            alfa(j,k) = alphaT(j,k)
        enddo
        
        ny(:NC,j) = nu(:NC,j)
    
    enddo
    

    
    do i = 1, NC

    
    enddo

endsubroutine Charge_Commons_GC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine OW_Partition_Coefficient(POW,COMPOUND)
    use constantes
    PARAMETER(NSCM=10,NG=70)
    IMPLICIT real*8 (A-H,O-Z)
    DIMENSION X(NCOM),ACT(NCOM)
    real*8 lnGamma(NCOM)
    integer maingr(NMG),grupo1(NMG),COMPOUND(DiffStructGroups,2)
    character car*3
    character*8 fs(NMG)
    character*37 nomco
    COMMON/PAREQ/IPAREQ
    COMMON/NOM/FS
    ! COMMON/US/MS(NCOM,NSCM,2),nms(NCOM)
    integer::MS(NCOM,DiffStructGroups,2),nms(NCOM)
    COMMON/PUNSUB/NPUNT(NMG),NGRUP(NMG),NUM
    COMMON/PUNGRU/NPINT(NG),NINTT(NG),NUMINT
 
	nms(1)=3
	nms(2)=1
	nms(3)=MSOL
	i=1
	do while (compound(i,1)/=0)
	    MS(3,I,1)=COMPOUND(I,1)
	    MS(3,I,2)=COMPOUND(I,2)
	    i=i+1
	ENDDO
	ms(1,1,1) = 1
	ms(1,1,2) = 1
	ms(1,2,1) = 2
	ms(1,2,2) = 7
	ms(1,3,1) = 15
	ms(1,3,2) = 1
	if(ipareq.eq.1) ms(1,3,1) = 14
	ms(2,1,1) = 17
	ms(2,1,2) = 1
	T=290.15
	X(1)=.725
	X(2)=.275
	X(3)=0.
	!call Gamma(IPAREQ,3,T,X,ACT)
	call unifac(ms,x,3,T,Act)
	ACT1=ACT(3)
	X(1)=0.
	X(2)=1.
	X(3)=0.
	!call Gamma(IPAREQ,3,T,X,ACT)
	call unifac(ms,x,3,T,Act)
	ACT2=ACT(3)
	Pow=.151*ACT2/ACT1
!
endsubroutine

SUBROUTINE SOLBIN(NG,T)
	parameter(NMG=150)
      implicit real*8(A-H,O-Z)
      DIMENSION DMAT(2,3),DACT1(10,10),DACT2(10,0010),PACT(2,2)       
      DIMENSION X1(10),X2(10),ACT1(10),ACT2(10)
	common/solub/X1,X2
	WRITE(6,*)'Initial guess for the concentrations (in mole percent)'
	write(6,*)'Y12(component 2 in phase 1):'
	read(5,*)Y12
	write(6,*)'Y21(component 1 in phase 2):'
	read(5,*)Y21
	iout=1
      NITER=0                     
      X1(1)=1.D0-Y12/100.D0                                             
      X2(1)=Y21/100.D0                                                  
   10 NITER=NITER+1               
      IF(NITER.GT.10) GOTO 50       
      IF(X1(1).LT.0.D0) X1(1)=0.D0 
      IF(X2(1).LT.0.D0) X2(1)=0.D0
      X1(2)=1.D0-X1(1)   
      X2(2)=1.D0-X2(1)                                                  
      CALL UNIFA3(NG,T,3,X1,ACT1,DACT1,PACT) !actividad en fase 1
      !CALL UNIFA(0,0,3,NG,T,X1,ACT1,DACT1,PACT)                        
      CALL UNIFA3(NG,T,3,X2,ACT2,DACT2,PACT) !actividad en fase 2   
      DO 20 I=1,2                                                       
      DMAT(I,1)=DACT1(I,1)-DACT1(I,2)                                   
      DMAT(I,2)=DACT2(I,2)-DACT2(I,1)                                   
   20 DMAT(I,3)=ACT1(I)-ACT2(I)                                         
      CALL GAUSL(2,3,2,1,DMAT)                                          
      RES=DMAT(1,3)**2+DMAT(2,3)**2                                     
      X1(1)=X1(1)-DMAT(1,3)                                             
      X2(1)=X2(1)-DMAT(2,3)                                             
      IF(RES.GT.1.D-20) GOTO 10                                         
   50 CONTINUE                                                          
      WRITE(6,603)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,603)                                     
      IF(IOUT.NE.6) WRITE(IOUT,604) X1(1),X2(1),X1(2),X2(2)             
      WRITE(6,604) X1(1),X2(1),X1(2),X2(2)                              
      CALL GCON(2,X1,ACT1,DACT1,ICVEX)                                  
      IF(IOUT.NE.6.AND.ICVEX.EQ.-1) WRITE(IOUT,601)                     
      IF(ICVEX.EQ.-1) WRITE(6,601)                                      
      CALL GCON(2,X2,ACT2,DACT2,ICVEX)                                  
      IF(IOUT.NE.6.AND.ICVEX.EQ.-1) WRITE(IOUT,602)                     
      IF(ICVEX.EQ.-1) WRITE(6,602)                                      
	CALL PAUSA
  601 FORMAT(' FALSE SOLUTION IN PHASE 1')                              
  602 FORMAT(' FALSE SOLUTION IN PHASE 2')                              
  603 FORMAT(///,5X,'** BINARY SOLUBILITIES IN MOLE FRACTIONS **',//,11X,'COMPONENT 1',15X,'COMPONENT 2',/)                               
  604 FORMAT(2(2X,2P2D12.2)//)                                          
 100  RETURN                                                            
      ENDsubroutine SOLBIN                                                      

      SUBROUTINE GCON(NK,X,ACT,DACT,ICVEX)                              
      implicit real*8(A-H,O-Z)                                 
      DIMENSION X(3),DG(2),DDG(2,2),ACT(3),DACT(10,10)                  
      ICVEX=1                                                           
      DO 1 I=1,NK                                                       
    1 IF(ACT(I).LT.1.D-10) ACT(I)=1.D-10                                
      DO 5 I=1,NK                                                       
      DO 5 J=1,NK                                                       
    5 DACT(I,J)=DACT(I,J)/ACT(I)                                        
      IF(NK.EQ.3) GOTO 9                                                
      DDG(2,2)=DACT(2,2)-DACT(1,2)-DACT(2,1)+DACT(1,1)                  
      GOTO 30                                                           
9     DO 20 I=2,NK                                                      
      II=I-1                                                            
      DO 20 J=2,NK                                                      
      JJ=J-1                                                            
   20 DDG(II,JJ)=DACT(I,J)-DACT(1,J)-DACT(I,1)+DACT(1,1)                
      IF(X(1).LE.1.D-12.OR.X(2).LE.1.D-12) GOTO 30                      
      DET=DDG(1,1)*DDG(2,2)-DDG(2,1)*DDG(2,1)                           
      IF(DET.LE.0.D0.OR.DDG(1,1).LE.0.D0.OR.DDG(2,2).LE.0.D0) ICVEX=-1  
      GOTO 100                                                          
   30 CONTINUE                                                          
      IF(DDG(2,2).LE.0.D0) ICVEX=-1                                     
  100 CONTINUE                                                          
      RETURN                                                            
      ENDSUBROUTINE GCON         
      
      SUBROUTINE RELVOLPRO (RELVOL,BPOINT)
      use Input
    use constantes
      use input_data, only:pvomegaybp
      IMPLICIT real*8 (A-H,O-Z)
      type(Compound),pointer::ptr_pcp
      DIMENSION X(NCOM),ACT(NCOM),DACT(NCOM,NCOM),TACT(NCOM)
      COMMON/nr/nrcpr
      COMMON/US/MS(NCOM,DiffStructGroups,2),nms(3)
      logical::NOACE
      !COMMON/US/MS(NCOM,NSCM,2)
      

      NDIF = 0
      NACT = 0
      IOUT = 5
      NC = 3
      CALL UNIPAR (NC,NG,BPOINT,MODEL,IOUT,NOACE,MS)
      X(1) = 0.
      X(2) = 1.
      X(3) = 0.
      CALL UNIFA (NDIF,NACT,NC,NG,BPOINT,X,ACT,DACT,TACT)
      call pvomegaybp (nrcpr,ptr_pcp)
      PSAT = EXP(ptr_pcp%a(1) -ptr_pcp%a(2)/(ptr_pcp%BoilingPoint + ptr_pcp%a(3)))
      RELVOL=ACT(3)*760/PSAT
ENDSUBROUTINE RELVOLPRO


subroutine Gamma(IPAREQ,NC,T,X,ACT)
    use constantes
	parameter (NMODEL=3,NGA=70)
      implicit real*8 (A-H,O-Z)
      character*17 tabla(nmodel)
      DIMENSION X(NCOM),ACT(NCOM),TACT(NCOM),DACT(NCOM,NCOM)
      COMMON/PROPIEDADES/PROPER
      COMMON/PUNSUB/NPUNT(NMG),NGRUP(NMG),NUM
      COMMON/PUNGRU/NPINT(NGA),NINTT(NGA),NUMINT
      COMMON/US/MS(NCOM,DiffStructGroups,2),nms(3)
      LOGICAL NOACE
      LOGICAL PROPER
	IF (PROPER) THEN
	  IDEV=0
	ELSE
	  IDEV=6
	ENDIF
      tabla(1)='    liquid-liquid'
	tabla(2)='     liquid-vapor'
	tabla(3)=' ifinite dilution'
	call PARIN (NC,NG,IDEV,NOACE,MS)
      if (NOACE) then
       write (idev,1490) tabla(ipareq)
         ! CALL PAUSA 
	 GO TO 100
	end if
	call PARAM (NC,NG,T)
      NDIF = 0
      NACT = 0
      call unifa(NDIF,NACT,NC,NG,T,X,ACT,DACT,TACT)
	!call unifa2(NDIF,NACT,NC,NG,T,X,ACT,DACT,TACT)
!	formatos
1490  format (' ',/,' * No interaction parameters available for ',&
              'this solution',/,'   in the ',a17,&
              ' UNIFAC table')
 100	return
	end