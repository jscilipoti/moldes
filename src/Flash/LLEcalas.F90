subroutine llecalas(T,PP,compuestos,N,xtemporal,act)
!C  *******************************************************************  
!C  *                                                                 *  
!C  *           PROGRAM   L L E C A L A S  (asociacion incorporada    *  
!C  *                                      para el calculo de flash y *  
!C  *                                      curva binodal (icalc 0 y 1)*
!C  *                                                                 *
!C  *                          ASOCIACI�N CRUZADA					   *		  
!C  *                       VERSI�N GENERALIZADA                      *        
!C  *              FEBRERO 2006 MODIFICADA POR                        *
!C  *                   ALFONSINA  ESTER ANDREATTA                    *        
!C  *        BASADA EN LAS SIMPLLIFICACIONES DE LOS PAPERS:           *
!c  *       Revisada en Octubre del 2007 en el chequeo de estabilidad *
!C  *																   *   	
!c  * Michelsen, et al. (Fluid Phase Equilibria, 180(2001)165-174 )   *		
!C  * Tan, et al.  (Ind. Eng. Chem. Res, 2004,43,203-208).			   *	   	
!C  *																   *	   	
!C  *        Esto permiti�  que todos los casos particulares          *         
!c  *       de asociaci�n se puedan simplificar a un �nico c�lculo. 
!c
!c   V�lido para un m�ximo n�mero grupo asociativo de 12
!c   Con la implementaci�n en el c�lculo de la fracci�n no asociada en el componente puro 
!c   por  el metodo iterativo aqu� implementado se permite que una mol�cula
!c   tenga m�s de un grupo asociativo 14/07/06
!C  El c�lculo se limita a que el n�mero m�ximo de sitios sea 2(por razones matem�ticas)
!c                                                       
!C  *******************************************************************  
!C  *                                           DATE: 24/3 - 1982 /TJ *  
      use Input
      use constantes
      IMPLICIT REAL*8(A-H,O-Z)    
      PARAMETER(NGPM=30)



    integer,intent(in)::N 
    integer,intent(in)::Compuestos(NCOM,DiffStructGroups,2)
    real*8,intent(in)::T,PP,xtemporal(10)   
    real*8,intent(out)::act(NCOM)

!Internas
   ! integer::GAM(N,N)        
                                
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)      
                           
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
!      COMMON/CUFAC/N,NG,P(10,10),T                                      
      COMMON/CY/Y13,Y21,STEP                                            
      COMMON/CA/XC(5),GE(5,2),GC(5,2)                                   
      COMMON/CIPR/IPR                                                   
  !    COMMON/CQT/QT(10,10),Q(10),R(10)                                  
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
!      COMMON/CMODEL/MODEL                                               
      COMMON/COUT/IOUT                                                  
      common/nga/nga,mass(12)
      common/ig/ig
      common/unif/QT(NGPM,NCOM),TAU(NGPM,NGPM),ST(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
      real*8,dimension(10)::DLX,Y  
      real*8,dimension(30)::YVAL,GRAD
      DIMENSION NTEXT(36),X(2),ANT(10,3),XMAT(30,30),WORK(30,5)         
      integer::ICALC,IPR,IOUT,NOVAP,ig            
      character(len=36)::name, name1 
      integer:: parameters, output
        common/NumGrup/NG
      dimension xmj(10),actgam(10),agam(10,4),de(10,10),pe(2,2)
      logical NOACE
!c-----
        output=5673
	OPEN (UNIT=output,FILE='llecalasnew.OUT',FORM='FORMATTED')
	write(output,"(20I3)") (compuestos(5,i,1),compuestos(5,i,2),i=1,10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    z(:)=xtemporal(:) !composici�n temporal para probar llecalas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    act(:) = 0.0
    NC=N
   
    icalc=0
!   icalc:  0-' **** FLASH CALCULATION ****'                            
!           1-' **** BINODAL CURVE CALCULATION ****'
!           2-' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** '
    
!Adapta la variable model de MolDeS a la de Ilecalas
    if (model==1) model=0
!   En Ilecalas:
!       0-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'     
!       1-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'
!   En MolDes:
!        1-' UNIFAC '
!        2-' A-UNIFAC:'
!        3-'GC-EOS'

    ipr=0
!   ipr:    1-' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'
    iout=0
!   iout:   1-'open 'lleasoccuzada.OUT''
    novap=0
!   novap:  0-'VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS'
    ig=1 
!   ig:     0-'write compositions'
!           1-'write compositions and activities'
  
    IF (IOUT.EQ.1) OPEN (UNIT=1,FILE='lleasoccuzada.OUT',FORM='FORMATTED')                                                     
   
    IF(NOVAP/=0) then  !   novap:  0-'VAPOR PHASE NOT INCLUDED IN FLASH-CALCULATIONS'
        write(*,*) "ver comentarios en esta parte del c�digo"
        pause
        !adaptar los siguientes bucles
        !DO J=1,N                                                        
        !    READ(2,*) (ANT(K,J),K=1,3)                                        
        !enddo
        !DO 7 J=1,N                                                        
        !    ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                             
        !    ANT(2,J)=2.302585*ANT(2,J)                                        
        !enddo
    endif
    
    T1=0.                                                             
    NN=0                                                              

    IF(T.EQ.0.) GOTO 10000 !return                                           
      
    IF(PP/=0..and.NOVAP/=0) then                             
        DO I=1,N                                                        
            PRAT(I)=DLOG(PP)-ANT(1,I)+ANT(2,I)/(T-273.15+ANT(3,I))       
        enddo                                        
    endif
                                           
    ZSUM=0.                                                           
    ZMAX=0.                                                           
    DO I=1,N                                                       
        ZSUM=ZSUM+Z(I)                                                    
        IF(Z(I).LT.ZMAX) cycle                                        
        ZMAX=Z(I)                                                         
        MAXZ=I                                                            
    enddo                                                          
    
    IF(T.EQ.T1) GOTO 30   
    
    CALL PARIN (NC,NG,IOUT,NOACE,compuestos)                                            
    CALL param(N,NG,T) !Calcula los coeficientes de actividad residuales grupales                                               
      
    
!================================================================
    !Ac� deber�a ir el if que hace el c�lculo para icalc == 2 
!================================================================
    DO I=1,N                                                       
        DO J=1,N                                                       
            GAM(I,J)=0.D0                                                     
            IF(J.EQ.I) cycle                                             
            CALL GAMINF(NG,I,J,G)                                                
            GAM(I,J)=G                                                        
        enddo    
    enddo
    
30  T1=T                                                              
    NN=NN+1                                                           
    WRITE(6,"(///,' * * * FLASH NUMBER',I3,' * * *',//)   ") NN                                                   
   
    DO I=1,N                                                       
        Z(I)=Z(I)/ZSUM                                                    
    enddo  
    
    WRITE(6,"(' TEMPERATURE =',F10.4,' K, PRESSURE =',F7.3,' ATM, FEED =',&
                F10.2,' MOLES',/,' FEED COMPOSITION (MOLE PERCENT):',/,&
                1X,15(2PF7.3))  ") T,PP,ZSUM,(Z(I),I=1,N)                               
    
    CALL unifaclle(N,1,Z,AL,DA,PACT)                                       
    
    SFAS(1)=1.                                                        
    GNUL=0.                                                           
    
    DO I=1,N                                                       
        XVL(I,1)=1.                                                       
        Z(I)=Z(I)+1.D-20                                                  
        DLX(I)=DLOG(Z(I))                                                 
        A(I)=AL(I)+DLX(I)                                                 
        GNUL=GNUL+Z(I)*AL(I)
    enddo
    
    NF=1  
    !esta l�nea es llamada desde el final de las instrucciones
50  CALL STIG(compuestos,n,Y,S)                                                    
    
    IF(S <= -1.D-7) then                                           
        write(6,"(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')  ")                                                                                       
        do I=1,N                                                       
            YVAL(I)=1.D-5*Y(I)/Z(I)                                           
        enddo
      GOTO 100   
    endif      
   
    DO I=1,N                                                       
        YVAL(I)=DLOG(Y(I))
    enddo
    

    XLAM=1.                    
    ! ipr:    1-' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'
    IF(NF.EQ.1.AND.IPR.GT.0) WRITE(6,"(//,' DIFFERENTIAL STABILITY TEST FOR FEED MIXTURE:')  ")                              
    IF(NF.GT.1.AND.IPR.GT.0) WRITE(6,"(//,' DIFFERENTIAL STABILITY TEST FOR',I2,'-PHASE SYSTEM')  ") NF                             
   
    CALL TMSSJ(30,N,IPR,15,XLAM,1.D-12,FUN,YVAL,GRAD,XMAT,WORK,1)     
   
    IF(FUN >= -1.D-7) then                                         
        
        WRITE(6,"(/,' * SYSTEM IS STABLE *',/)   ")         
        !write(output,*) 1
      !  write(output,"(5(2x,f12.8))") (Z(j),J=1,N)
      !  write(output,"(5(2x,f12.8))") (AL(j),j=1,N)        
       ! write(output,*) "SYSTEM IS STABLE"                                                   

	    write(output,46) T,  (xM(l,1),l=1,N)    !Alfonsina
	    write(output,46) T,  (xM(l,2),l=1,N)!Alfonsina
	    write(output,*)                    !Alfonsina
        46	FORMAT (2X,F12.2, 8X,F12.6, 8X,F12.6 , 8X,F12.6, 8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6) !alfonsina                   
        GOTO 10000  !return      
    elseif(isnan(fun))then
        act(:) = 0.0
        return
    endif      
    
    WRITE(6,"(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')  ")                                                      
                                 
    DO I=1,N                                                       
        YVAL(I)=1.D-5*DEXP(YVAL(I))/Z(I)                                  
    enddo
    
100 NF=NF+1 
    
    do i=1,N
        do while (yval(i) > 1.0) 
            yval(:) = YVAL(:)/10.0
        enddo
    enddo
                                                      
    SFAS(NF)=1.                                                       
    XLAM=.2                                                           
    IF(NF.EQ.2) XLAM=.5                                               
    M=(NF-1)*N                                                        
    IF(IPR.GT.0) WRITE(6,607) NF                                      
       
    CALL TMSSJ(30,M,IPR,60,XLAM,1.D-16,FUN,YVAL,GRAD,XMAT,WORK,2)     
    
    NT=NF*N                                                           
    NB=NT-N                                                           
    
    DO I=1,NB                                                     
        YVAL(NT+1-I)=YVAL(NB+1-I)    
    enddo
      
    WRITE(6,"(//,' RESULT OF',I2,'-PHASE CALCULATION:')  ") NF                                                   
    NVAP=0                                                            
    DO J=1,NF                                                     
        IF(IDUM(J).EQ.1) NVAP=J                                           
    enddo
    
    IF(NVAP.EQ.0) WRITE(6,"(' NO VAPOR PHASE') ")                                        
    IF(NVAP.NE.0) WRITE(6,"(' PHASE',I2,' IS A VAPOR PHASE')  ") NVAP                                   
    
    WRITE(6,"(/,'  PHASE FRACTIONS (PERCENT):',4(5X,I3,2PF7.3,5X))    ") (J,SFAS(J),J=1,NF)                                   
    WRITE(6,"(/,'  COMPOSITION  ',10X,4(8X,I3,9X)) ") (J,J=1,NF)                                           
                                   
    SUM=0.                                                            
    DO I=1,N                                                      
        DLX(I)=XVL(I,NF)*Z(I)/SFAS(NF)                                    
        SUM=SUM+DLX(I)
    enddo
    
    SUM=DLOG(SUM)                                                     
    CALL unifaclle(N,1,DLX,A,DA,PACT)        !                               
      
    DO I=1,N                                                      
        DLX(I)=DLOG(DLX(I))                                               
        A(I)=A(I)+DLX(I)-SUM  
    enddo

      
    do j=1,nf
        do i=1,n
            xmj(i)=xm(i,j)
        enddo
        call unifaclle(n,1,xmj,actgam,de,pe) 
        do i=1,n
            agam(i,j)=actgam(i)
        enddo
    enddo
      
    write (output,*)"Nro fases", NF
      
    DO  I=1,N                                                      
        WRITE(6,613) I,(XM(I,J),J=1,NF)     !composition        
        write(6,1613) i,(agam(i,j),j=1,nf) !Ln(gamma)
        act(i) = xm(i,1)*exp(agam(i,1))
    enddo
        
    GOTO 50                                                           


  501 FORMAT(36A2)                                                      
  502 FORMAT(8F10.2)                                                    
  503 FORMAT(20I3)                                                      
  607 FORMAT(/,' PHASE SPLIT CALCULATION,',I2,' PHASES:')               
  608 FORMAT(1H1)                                                   
  610 FORMAT(///)                                                                  
  613 FORMAT('   X(',I2,')            ',5(8X,F12.8))                    
 1613 format('  ln(G',i2,')            ',5(8x,f12.8))       
  616 FORMAT(//,' * WRONG INPUT SPECIFICATION *',//)                    
  617 FORMAT(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)                         
  618 FORMAT(//,' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'//)                                       
  619 FORMAT(10F12.5)                                                   
  620 FORMAT(' **** FLASH CALCULATION ****')                            
  621 FORMAT(' **** BINODAL CURVE CALCULATION ****',//)                 
  622 FORMAT(' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** ',//)                                                              
  623 FORMAT(1X,'COMPONENTS : ',40A2,//)                                
  624 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'//)     
  625 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'//)    
  626 FORMAT(I5,2F15.4)                                                 
  627 FORMAT(//,' SPECIFIED UNIQUAC R AND Q',/)                         
  628 FORMAT(/,' IOUT = ',I2,/' IF IOUT = 0: OUTPUT ONLY ON UNIT 6',/,  ' IF IOUT = 1: OUTPUT ON BOTH UNIT 6 AND 1')                      
  629 FORMAT(/,' VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS',//)                     
633   FORMAT(///,'   TEMPERATURE =',F10.2,' DEG K')                     

10000 continue
      write(output,*)"*********************************"
      write(output,*)
      
      IF (IOUT.EQ.1) CLOSE (UNIT=1)
      close (unit=3)
    
      if (model==0)model=1
    
      return                                              
    endsubroutine llecalas  

!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
    
    
SUBROUTINE GAMINF(NG,N1,N2,GAM)                                      
       
    use constantes   
    IMPLICIT REAL*8(A-H,O-Z)   
     
     PARAMETER(NGPM=30)                   
   ! COMMON/CUFAC/NK,NG,P(10,10),T                                     
   ! COMMON/CPAR/TAU(10,10),S(10,10),F(10)                             
  !  COMMON/CQT/QT(10,10),Q(10),R(10)     
    common/unif/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),PTemp(NGPM,NGPM)                             
!c------
    common/asoc/nktt,igamt(20,12),nytt(20,12)  
    common/nga/nga,mass(12)
    common/grupas1/rkass(6,12,6,12),enass(6,12,6,12), deloh(6,12,6,12) !Alfonsina
!c------
    dimension xohi0(10),xnohi0(10),tgt(10),xgamk(20)
	common/ioh2/rngoh(12,12)
	common/ioh2sis/rngoht(3,2)

!c------
    dk=1.381e-23
    deloh=0.0
    xnoh=0.0
    xnoh1=0.0
	xoh=0.0
    xgam=0.0
    xgamt=0.0
    do 7777 i=1,10
        xnohi0(i)=0.0
        tgt(i)=0.0
7777 continue
    do 8888 k=1,20
        xgamk(k)=0.0
8888 continue

!c------
    Q1=Q(N2)/Q(N1)                                                    
    R1=R(N2)/R(N1)                                                    
    QR=R1/Q1                                                          
    GAM=F(N2)+Q(N2)*(1.-DLOG(Q1))-R1+DLOG(R1)-5.D0*Q(N2)*(1.-QR+DLOG(R1)-DLOG(Q1))                                                      
    DO 10 I=1,NG                                                      
10      GAM=GAM-S(I,N2)/S(I,N1)*QT(I,N1)-QT(I,N2)*DLOG(S(I,N1))           
!c------

!c------
 
endsubroutine gaminf    


 SUBROUTINE SPLIT(ND,N,IDI,BETA,DEL,G,W)                           
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION G(ND,ND),W(ND,5)                                        
      IDI=0                                                             
      DO 20 I=1,N                                                       
      I1=I-1                                                            
      I2=I+1                                                            
      ID=0                                                              
      TV=G(I,I)                                                         
      SV=TV                                                             
      IF (I1.EQ.0) GO TO 25                                             
      DO 30 J=1,I1                                                      
   30 SV=SV-G(I,J)**2                                                   
   25 IF (SV.LT. DEL*DEL) ID=1                                          
      SVR=DSQRT(DABS(SV))                                               
      IF (SVR.LT. DEL) SVR=DEL                                          
      XM=0.                                                             
      IF (I.EQ.N) GO TO 35                                              
      DO 40 J=I2,N                                                      
      S=G(J,I)                                                          
      IF (I1.EQ.0) GO TO 45                                             
      DO 50 K=1,I1                                                      
   50 S=S-G(I,K)*G(J,K)                                                 
   45 S=S/SVR                                                           
      IF (DABS(S).GT. XM) XM=DABS(S)                                    
   40 G(J,I)=S                                                          
   35 IF (XM.LT. BETA) GO TO 55                                         
      ID=1                                                              
      XM=XM/BETA                                                        
      DO 60 J=I,N                                                       
   60 G(J,I)=G(J,I)/XM                                                  
      SVR=SVR*XM                                                        
   55 IF (ID.EQ.1) W(I,1)=SVR**2-SV                                     
      G(I,I)=SVR                                                        
      IDI=IDI+ID                                                        
   20 CONTINUE                                                          
      RETURN                                                            
    endsubroutine split                                                               
 
    
SUBROUTINE LINE(ND,N,XLAM,GD,G,W)         

      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION GD(ND),G(ND,ND),W(ND,5)                                 
   65 W(1,2)=-GD(1)/(G(1,1)+XLAM)                                       
      IF (N.EQ.1) GO TO 75                                              
      DO 70 I=2,N                                                       
      I1=I-1                                                            
      S=-GD(I)                                                          
      DO 80 J=1,I1                                                      
   80 S=S-G(I,J)*W(J,2)                                                 
   70 W(I,2)=S/(G(I,I)+XLAM)                                            
   75 W(N,3)=W(N,2)/(G(N,N)+XLAM)                                       
      IF (N.EQ.1) GO TO 85                                              
      DO 90 II=2,N                                                      
      I=N+1-II                                                          
      I1=I+1                                                            
      S=W(I,2)                                                          
      DO 100  J=I1,N                                                    
  100 S=S-G(J,I)*W(J,3)                                                 
   90 W(I,3)=S/(G(I,I)+XLAM)                                            
   85 RETURN                                                            
 endsubroutine line   
 



SUBROUTINE CHECK(N,YVAL)                                          
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
    COMMON/CIPR/IPR                                                   
    COMMON/COUT/IOUT                                                  
    DIMENSION YVAL(30)                                                
    
    JMAX=NF                                                           
    SMAX=SFAS(NF)                                                     
    NT=NF-1                                                           
    DO J=1,NT                                                       
        IF(SFAS(J).LT.SMAX) cycle                                        
        SMAX = SFAS(J)                                                      
        JMAX = J                                                            
    enddo                                                          
    if(JMAX /= NF) then                                          
        DO I=1,N                                                       
            DS=1.                                                             
            DO J=1,NT                                                      
                DS=DS-YVAL(I+(J-1)*N)                                             
            enddo
            YVAL(I+(JMAX-1)*N)=DS                                             
        enddo
    endif
    
    IF(NF.EQ.2) GOTO 100                                              
        DO 21 I=1,N                                                       
            XVL(I,NF)=1.                                                      
            DO 21 J=1,NT                                                      
                XVL(I,J)=YVAL(I+(J-1)*N)                                          
   21   XVL(I,NF)=XVL(I,NF)-XVL(I,J)                                      
        DO 23 J=1,NF                                                      
            SFAS(J)=0.                                                        
            DO 22 I=1,N                                                       
                AL(I)=XVL(I,J)*Z(I)                                               
   22           SFAS(J)=SFAS(J)+AL(I)                                             
            DO 23 I=1,N                                                       
   23           XM(I,J)=AL(I)/SFAS(J)                                             
            DO 30 I=1,NT                                                      
                JN=I+1                                                            
                DO 30 J=JN,NF                                                     
                    IF(DABS(XM(MAXZ,I)-XM(MAXZ,J)).LT.5.D-3) GOTO 40                  
   30 CONTINUE                                                          
      GOTO 100                                                          
   40 IV=I                                                              
      JV=J                                                              
      DMAX=0.                                                           
      DO 50 I=1,N                                                       
      R1=DABS(XM(I,IV)-XM(I,JV))                                        
      IF(R1.GT.DMAX) DMAX=R1                                            
   50 CONTINUE                                                          
      IF(DMAX.GT.2.5D-2) GOTO 100                                       
      WRITE(6,601)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,601)                                     
      NF=NF-1                                                           
      NT=NT-1                                                           
      DO 60 I=1,N                                                       
      XVL(I,IV)=XVL(I,IV)+XVL(I,JV)                                     
   60 XVL(I,JV)=XVL(I,NF+1)                                             
  100 CONTINUE                                                          
      IF(NF.LT.3) GOTO 250                                              
      MINF=0                                                            
      DO 200 I=1,NF                                                     
      IF(SFAS(I).LT.1.D-12) MINF=I                                      
  200 CONTINUE                                                          
      IF(MINF.EQ.0) GOTO 250                                            
      WRITE(6,602) MINF,SFAS(MINF)                                      
      IF(IOUT.NE.6) WRITE(IOUT,602) MINF,SFAS(MINF)                     
      DO 220 I=1,NF                                                     
      IF(I.EQ.MINF) GOTO 220                                            
      DO 210 J=1,N                                                      
  210 XVL(J,I)=XVL(J,I)+SFAS(I)*XVL(J,MINF)                             
  220 CONTINUE                                                          
      IF(MINF.EQ.NF) GOTO 250                                           
      NNF=NF-1                                                          
      DO 230 I=1,NNF                                                    
      IF(I.LT.MINF) GOTO 230                                            
      DO 240 J=1,N                                                      
  240 XVL(J,I)=XVL(J,I+1)                                               
  230 CONTINUE                                                          
      NF=NF-1                                                           
  250 CONTINUE                                                          
  601 FORMAT(//,' * * * TWO PHASES ARE IDENTICAL * * *'//)              
  602 FORMAT(/,' * THE AMOUNT OF PHASE',I2,' IS',D12.4,'. THE PHASE IS ELIMINATED *',/)                                                   
      RETURN                                                            
endsubroutine check


                                          
SUBROUTINE STIG(compuestos,nc,Y,S)
    use constantes
     IMPLICIT REAL*8(A-H,O-Z)       
      
      integer,intent(in)::nc
      integer,intent(in)::Compuestos(NCOM,10,2)                                   
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
      COMMON/CUFAC/N,NG,P(10,10),T                                      
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AA(10),DA(10,10),XM(10,4)                                       
      real*8,dimension(10)::Y,V,YGEM                                    
      common/nga/nga,mass(12)
     
      JPUR=0                                                            
      AMAX=0.   
      n=nc
      
      DO I=1,N                                                       
        IF(A(I).LT.AMAX) cycle                                          
        JPUR=I                                                            
        AMAX=A(I)                                       
      enddo
      
      RMAX=1.D5                                                         
      NGN=N                                                             
      
      IF(NF.GT.1) NGN=N+NF                                              
      
      NEG=0                                                             
      main:DO KK=1,NGN                                                   
        JM=KK                                                             
        IF(JPUR.NE.0) JM=JPUR                                             
      
        IF(JM.LE.N) then
            SUM=0.                                                            
            DO I=1,N                                                       
                GG=A(I)-GAM(JM,I)                                                 
                IF(GG.LT.-50.D0) GG=-50.D0                                        
                Y(I)=DEXP(GG)                                                     
                SUM=SUM+Y(I)
            enddo            
        else    
            DO I=1,N                                                       
                Y(I)=Z(I)*(2+XVL(I,JM-N)/SFAS(JM-N))/3                            
            enddo
        endif
   
        NA=3                                                              
      
        DO K=1,NA                                                      
            DO I=1,N                                                       
                Y(I)=Y(I)/SUM                                                     
            enddo
            CALL unifaclle(NC,1,Y,AA,DA,PACT)! !(Compuestos,x,NC,T,Gamma)                                       
            IF(K.EQ.NA) exit                                               
            DO I=1,N                                                       
                Y(I)=DEXP(A(I)-AA(I))                                             
            enddo   

            SUM=0.                                                            
      
            DO I=1,N
                SUM=SUM+Y(I)
            enddo
        enddo
                                                         
        YV1 = 0.0
      
        DO J=1,NF                                                      
            V(J) = 0.0
        enddo
        
        DO I=1,N                                                       
            GD = DLOG(Y(I))+AA(I)-A(I)                                          
            YV1 = YV1+Y(I)*GD                                                   
            DO J=1,NF                                                      
                K=J                                                               
                VV=XVL(I,K)*Z(I)/SFAS(K)                                          
                D=GD*(Y(I)-VV)                                                    
                V(J)=V(J)+D
            enddo
        enddo
        
        YV2=V(1)                                                          
      
        DO J=1,NF                                                      
            IF(V(J).LT.YV2) YV2=V(J)                                          
        enddo

        RT1=YV1                                                           
      
        IF(YV2.GT.0.) RT1=RT1-YV2/2                                       
      
        IF(.not.(NEG.EQ.0.AND.YV1.GT.0.)) then
            RT1=YV1                                                           
            IF(NEG.EQ.0) RMAX=0.                                              
            NEG=1
        endif
        
        IF(RT1.GT.RMAX) cycle                                          
      
        S=YV1                                                             
        RMAX=RT1                                                          
        CC=DEXP(-YV1)                                                     
        DO I=1,N                                                       
            YGEM(I)=Y(I)*CC                                                   
        enddo
        IF(JPUR.NE.0) exit                                            
      enddo main                                                      
    
      DO I=1,N                                                      
        Y(I)=YGEM(I)
      enddo
    
      RETURN                                                            

    endsubroutine STIG         
      

                                                           
 SUBROUTINE GMIX(compuestos,NARG,NDIF,FUN,GRAD,XMAT,YVAL)                     
      IMPLICIT REAL*8(A-H,O-Z)       
      integer::narg
      integer,intent(in)::Compuestos(NARG,10,2)                                  
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),XX(10),DA(10,10),XM(10,4)                                       
      DIMENSION YVAL(30),GRAD(30),X(10),TM(10),FG(10),XMAT(30,30)       
      NG=NF-1                                                           
      N=NARG/NG                                                         
      JD=1                                                              
      IF(NDIF.EQ.2) JD=2                                                
      IF(NDIF.EQ.2) CALL CHECK(N,YVAL)                                  
      IF(NF.NE.NG) GOTO 20                                              
      NARG=NARG-N                                                       
      NG=NG-1                                                           
      DO 10 I=1,N                                                       
      DO 10 J=1,NG                                                      
   10 YVAL(I+(J-1)*N)=XVL(I,J)                                          
   20 FUN=-GNUL                                                         
      DO 50 I=1,N                                                       
      XVL(I,NF)=1.                                                      
      DO 30 J=1,NG                                                      
      XVL(I,J)=YVAL(I+(J-1)*N)                                          
   30 XVL(I,NF)=XVL(I,NF)-XVL(I,J)                                      
      ! write(*,*)xvl
      DO 40 J=1,NF                                                      
      IF(XVL(I,J).GT.0.) GOTO 40                                        
      FUN=0.                                                            
      GOTO 1000                                                         
   40 CONTINUE                                                          
   50 CONTINUE                                                          
      DO 200 J=1,NF                                                     
      SFAS(J)=0.                                                        
      DO 60 I=1,N                                                       
      X(I)=XVL(I,J)*Z(I)                                                
   60 SFAS(J)=SFAS(J)+X(I)                                              
      DO 65 I=1,N                                                       
      XX(I)=X(I)/SFAS(J)                                                
   65 XM(I,J)=XX(I)                                                     
      CALL unifaclle(NARG,JD,XX,FG,DA,PACT) !                                     
      IDUM(J)=NDUM                                                      
      DO 70 I=1,N                                                       
      TM(I)=DLOG(XVL(I,J)/SFAS(J))+FG(I)                                
   70 FUN=FUN+X(I)*TM(I)                                                
      IF(NDIF.EQ.0) GOTO 200                                            
      DO 80 I=1,N                                                       
      S=Z(I)*TM(I)                                                      
      IF(J.EQ.NF) GOTO 75                                               
      GRAD(I+(J-1)*N)=S                                                 
      GOTO 80                                                           
   75 DO 76 K=1,NG                                                      
      NK=I+(K-1)*N                                                      
   76 GRAD(NK)=GRAD(NK)-S                                               
   80 CONTINUE                                                          
      IF(NDIF.EQ.1) GOTO 200                                            
      DO 100 I=1,N                                                      
      ST=Z(I)/SFAS(J)                                                   
      DO 100 L=1,N                                                      
      S=ST*(DA(I,L)-1.)*Z(L)                                            
      IF(L.EQ.I)S=S+Z(I)/XVL(I,J)                                       
      IF(J.EQ.NF) GOTO 90                                               
      XMAT(I+(J-1)*N,L+(J-1)*N)=S                                       
      GOTO 95                                                           
   90 DO 92 K=1,NG                                                      
      DO 92 M=1,K                                                       
      NK=I+(K-1)*N                                                      
      NM=L+(M-1)*N                                                      
      IF(K.NE.M) XMAT(NK,NM)=S                                          
      IF(K.EQ.M) XMAT(NK,NM)=XMAT(NK,NM)+S                              
   92 CONTINUE                                                          
   95 CONTINUE                                                          
  100 CONTINUE                                                          
  200 CONTINUE                                                          
 1000 RETURN                                                            
endsubroutine GMIX                                                               
      
      
SUBROUTINE STABIL(N,NDIF,FUN,GRAD,XMAT,Y)                         
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
    COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
    DIMENSION GRAD(30),XMAT(30,30),Y(30),YEX(10),XEX(10)              
    common/nga/nga,mass(12)

    SUM=0.                                                            
    DO I=1,N                                                       
        YEX(I)=0.                                                         
        IF(Y(I).GT.-40.) YEX(I)=DEXP(Y(I))                                
        SUM=SUM+YEX(I)                                                    
    enddo
    DO I=1,N                                                       
        XEX(I)=YEX(I)/SUM                                                 
    enddo  
    JD=1                                                              
    IF(NDIF.EQ.2) JD=2                                                
    
    CALL unifaclle(N,JD,XEX,AL,DA,PACT)!                                    
      
    FUN=1.                                                            
    DO I=1,N                                                       
        S=Y(I)+AL(I)-A(I)                                                 
        IF(NDIF.EQ.0) cycle                                            
        GRAD(I)=YEX(I)*S                                                  
        FUN=FUN+YEX(I)*(S-1)                                              
    enddo  
      
    IF(NDIF.LT.2) return                                             
      
    DO I=1,N                                                       
        S=XEX(I)                                                          
        DO J=1,I                                                       
            XMAT(I,J)=S*YEX(J)*DA(I,J)                                        
            XMAT(J,I)=XMAT(I,J)                                               
        enddo
        XMAT(I,I)=XMAT(I,I)+YEX(I)+GRAD(I)                                
    enddo
    
    RETURN                                                            
endsubroutine STABIL                                            
      
      
                                                              
      
                                                      
      
SUBROUTINE TMSSJ(ND,N,IPR,NMAX,XLAM,GLIM,F   ,X,GD,G,W,ifunc)     
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/COUT/IOUT                                                  
    integer::ifunc
    DIMENSION X(ND),GD(ND),G(ND,ND),W(ND,5)       
    integer::compuestos(ND,10,2)
      
    NTM=0                                                             
    DEL=1.D-7                                                         
150 NTM=NTM+1                                                         
    GNM=0.                                                            
      
    if(ifunc==1)then
        CALL STABIL(N,2,F,GD,G,X)                                       
    elseif(ifunc==2)then
        CALL GMIX(compuestos,N,2,F,GD,G,X)  
    endif
    DO I=1,N                                                        
        GNM=GNM+GD(I)**2                                                  
    enddo
    IF (GLIM <= GNM) then                                     
        BETA=0.                                                           
        DO I=1,N                                                       
            W(I,4)=X(I)                                                       
            W(I,5)=GD(I)                                                      
            DO J=1,I                                                       
                S=DABS(G(I,J))                                                    
                IF (I.EQ.J) S=S*N                                                 
                IF (S.GT. BETA) BETA=S                                            
                W(I,1)=0.0
            enddo
        enddo
        BETA = DSQRT(BETA/N)                                                
        CALL SPLIT(ND,N,IDI,BETA,DEL,G,W)                                 
        XLM=0.                                                            
        NTS=0                                                             
  350   NTS=NTS+1                                                         
        CALL LINE(ND,N,XLM,GD,G,W)                                        
        SMAX=0.                                                           
        GP=0.                                                             
        DP=0.                                                             
        DO I=1,N                                                      
            S=W(I,3)                                                          
            IF (DABS(S).GT. SMAX) SMAX=DABS(S)                                
            GP=GP+S*GD(I)                                                     
            DP=DP+S*S*W(I,1)                                                  
        enddo
        FAC=1.                                                            
        IF (SMAX.GT. XLAM) FAC=XLAM/SMAX                                  
        DER=FAC*GP                                                        
        ALFA=1.                                                           
  210   FF=ALFA*FAC                                                       
        DO I=1,N                                                      
            X(I)=W(I,4)+FF*W(I,3)                                             
        enddo
        if(ifunc==1)then
            CALL STABIL(N,1,FNEW,GD,G,X)                           
        elseif(ifunc==2)then
            CALL GMIX(compuestos,N,1,FNEW,GD,G,X)   
        endif
 
        if(isnan(fnew).or.isnan(x(1))) return
                                                
        IF(FNEW == 0.) then                                           
            ALFA=ALFA/2.                                                      
            GOTO 210                                                          
        endif
        DELS=FNEW-F                                                       
        IF (.not.(FNEW.NE. 0.  .AND. DELS.LT. 1.D-10)) then                
            IF (NTS.GT. 1) GO TO 125                                          
            DE2=(DELS-ALFA*DER)/ALFA**2/2                                     
            GT=-DER/DE2                                                       
            IF (GT.GT. ALFA) ALFA=ALFA/3.                                     
            IF (GT.LE. ALFA) ALFA=GT                                          
            GOTO 210     
        endif
        
        PRED=FF*GP-FF**2/2*(GP+DP)                                        
        CALPRE=DELS/PRED                                                  
        F=FNEW                                                            
        DO I=1,N                                                      
            W(I,4)=X(I)                                                       
        enddo
        IF (NTS.GE.3) GO TO 125                                           
        IF (ALFA.EQ. 1.  .AND. CALPRE.GT. .8) GO TO 350                   
  125   DO I=1,N                                                      
            X(I)=W(I,4)                                                       
        enddo
        
        IF (IPR.NE.0) WRITE(6,130) NTM,F,CALPRE,GNM,ALFA                  
        IF(IOUT.NE.6.AND.IPR.NE.0) WRITE(IOUT,130) NTM,F,CALPRE,GNM,ALFA  
130     FORMAT(1X,I2,':  F = ',D16.8,'  CALPRE = ',F8.3,'  GRAD = ',D16.8,'  ALFA = ',D10.2)                                                
       
        IF (IPR.GT.1) WRITE(6,131) (X(I),I=1,N)                           
        IF(IOUT.NE.6.AND.IPR.GT.1) WRITE(IOUT,131) (X(I),I=1,N)           
131     FORMAT(' X-VECTOR',10F11.5,(/,9X,10F11.5))                        
        
        IF (IPR.GT. 1) WRITE(6,132)                                       
        IF(IOUT.NE.6.AND.IPR.GT.1) WRITE(IOUT,132)                        
  132   FORMAT(' ')                                                       
        IF (NTM.LT. NMAX) GO TO 150                                       
    endif
    
    IF(IPR.GT.0) WRITE(6,133) NTM,GNM,F                               
      IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,133) NTM,GNM,F              
  133 FORMAT(/,'  NUMBER OF ITERATIONS = ',I2,', NORM OF GRADIENT = ',  D12.5,', F = ',D14.6)                                             
      IF(IPR.GT.0) WRITE(6,131) (X(I),I=1,N)                            
      IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,131) (X(I),I=1,N)           
      RETURN                                                            
endsubroutine TMSSJ 
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                 
                                   
      
 

subroutine unifaclle(NK,NDIF,X,ACT,DACT,PACT)   
    use constantes                       
      IMPLICIT REAL*8(A-H,O-Z)           
      PARAMETER(NGPM=30)                               
!c------
      common/asoc/nktt,igamt(20,12),nytt(20,12)   
      common/nga/nga,mass(12)
      common/grupas1/rkass(6,12,6,12),enass(6,12,6,12),deloh(6,12,6,12)!Alfonsin
!c------
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)            
      common/NumGrup/NG               
     ! COMMON/CUFAC/NK,NG,Porig(10,10),T                                     
   !   COMMON/CPAR/TAU(10,10),S(10,10),F(10)
                     
   !   COMMON/CQT/QT(10,10),Q(10),R(10)                                  
      common/unif/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
      
      DIMENSION X(10),GAM(10),ACT(10),DACT(10,10),THETA(10),PHI(10),RI(10),QI(10),ETA(10),QIL(10),RIL(10),QID(10),ETAL(10),TETAR(10) 
      DIMENSION U(10,10),V(10,10),PACT(2,2),DTAU(2,2,2)                 
!c------
      dimension goh(10),xgamk(20),dxohdx(10),dxxdx(10,10),dasdx1(10,10),dasdx2(10,10),dasdx(10,10)
	common/ioh2/rngoh(12,12)

      dimension dif(12,12), dif1(10,12,12) !Alfonsina
      common/zzzas/xoh(6,12),xohi0(12,6,12),xoh_old(6,12),xohi(6,12),xohi_old(6,12), xohi0_old(12,6,12)  !Alfonsina
	dimension m_lambda(nga*2,nga*2),m_lambda1(nga*2,nga*2) !Alfonsina
	dimension psin(12) !Alfonsina
	dimension indx(20)
      double precision  m_lambda,m_lambda1,xoh,xohi0,xoh_old,xohi0_old  !Alfon
	double precision del, del1, dif, dif1, d1,psin, xgam, xnoh1 !Alfonsina
	integer order !Alfonsina
      double precision sum1, sum2, sum3, sum4, SUMA1J, sumaj !Alfonsina
      dimension xnohi0(12,12),tgt(12),dnohdx(12,12),actas(12) !Alfonsina
	dimension xnoh1(12), xnoh(12),das1(3),das3(3),dxkdni(12,6,12), dxkdnic(12,6,12), dgasdx(12)  !Alfonsina
      dimension dgasdxij (12,12), drhodx(12), drhodni(12,6,12)
	


      dk=1.381e-23
      deloh=0.0
      xnoh=0.0
      xnoh1=0.0
	xoh=0.0
      xgam=0.0
      do 7777 i=1,10
	xohi0=0
      xnohi0=0.0
      tgt(i)=0.0
 7777 continue

      THETS=0.                                                          
      PHS=0.                                                            
      DO 10 I=1,NK                                                      
      THETA(I)=X(I)*Q(I)                                                
      PHI(I)=R(I)*X(I)                                                  
      THETS=THETS+THETA(I)                                              
   10 PHS=PHS+PHI(I)                                                    
      DO 20 I=1,NK                                                      
      RI(I)=R(I)/PHS                                                    
      RIL(I)=DLOG(RI(I))                                                
      QI(I)=Q(I)/THETS                                                  
   20 QIL(I)=DLOG(QI(I))                                                

      do 33 i=1,nk
      goh(i)=0.
      tgt(i)=0.0
      xnohi0=0.0
      xgam=0.0
     
!CCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do j=1,nk
		tgt(j)=0.0
		xgam=0.0
	end do
   33 continue
      do i=1,nk
	if(nga.gt.0) then
	
		do k=1,nktt
			tgt(i)=tgt(i)+nytt(k,i)
		end do

		do j=1,nga  
			xnohi0(i,j)=rngoh(i,j)/R(i)  
		end do
		xgam=xgam+R(i)*x(i)
		end if  
	end do
  
      xnoh1=0d0
	do ja=1,nga
	    do i=1,nk
		xnoh1(ja)=xnoh1(ja)+rngoh(i,ja)*x(i)
	  end do
      end do
	
      do ja=1,nga
		xnoh(ja)=xnoh1(ja)/xgam
	end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c------
      DO 40 I=1,NG                                                      
      ETA(I)=0.                                                         
      DO 45 J=1,NK                                                      
   45 ETA(I)=ETA(I)+S(I,J)*X(J)                                         
   40 ETAL(I)=DLOG(ETA(I))                                              
      DO 55 I=1,NG                                                      
      TETAR(I)=0.                                                       
      DO 55 J=1,NK                                                      
   55 TETAR(I)=TETAR(I)+QT(I,J)*X(J)                                    
      DO 60 I=1,NK                                                      
      QID(I)=1.-RI(I)/QI(I)                                             
      XX=F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                             
      XX=XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                              
      ACT(I)=XX                                                         
      DO 661 J=1,NG                                                     
      U(J,I)=S(J,I)/ETA(J)                                              
      V(J,I)=U(J,I)*TETAR(J)                                            
  661 ACT(I)=ACT(I)-V(J,I)-QT(J,I)*ETAL(J)                              
!c------
  
  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!**************************calculo de las fuerzas de asociacion*******************
      if(nga.ne.0) then
	    DO J=1,NGA 
              IF(MASS(J).EQ.0) GO TO 201
                DO m=1,NGA
                  IF(MASS(m).EQ.0) GO TO 101
                        DO L=1,MASS(J)
                                DO K=1,MASS(m)
                                  IF(ENASS(K,m,L,J).EQ.0) THEN
                                         CONTINUE
                                    ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      deloh(k,m,l,j)=(DEXP(ENASS(K,m,L,J)/T) - 1 )*RKASS(K,m,L,J)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                                  END IF
                                END DO
                        END DO

 101            CONTINUE
                END DO
 201          CONTINUE
        END DO
  

	end if  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	!***********************************calculo de xoh**************************
!c C�lculo de la  fracci�n no asociada Paper:Ind. Eng. Chem. Res, 2004,43,203-208   
!c !Inicializaci�
      if(nga.ne.0) then
      xoh=1.0d0 
      del=1.d0
!c Iteraciones con tolerancia de 10-9
      do while (del>1.0d-10)
      xoh_old=xoh
	do m=1, nga
		do j=1,2
      sum1=0.D0
	do k=1, nga
	sum2=0.D0
	do l=1,2

	sum2=sum2+xoh_old(l,k)*deloh(l,k,j,m)
	end do
	sum1=sum1+sum2*xnoh1(k)
      end do
      xoh(j,m)=1.D0/(1.D0+sum1/xgam)          
	dif(j,m)=dabs((xoh(j,m)-xoh_old(j,m))/xoh(j,m))
	end do
	end do
	del=maxval(dif)
      end do
      end if

			write(4,*)"T=", t           
	write(4,*)"xoh (1,1)=", xoh (1,1)  
	write(4,*)"xoh (1,2)=", xoh(2,1) 
 
	 

!cc Fin del C�lculo de la  fracci�n no asociada 
!c	!*****************************calculo de xohi0**************************************
!C	xohi0=1d0
!c	do i=1, nc
!C		do j=1, nga
!C	       do l=1,mass(j)
!C	do k=1,nga
!C	 do m=1,mass(k)
!C			If (rngoh(i,j).eq.0d0) then
!C					xohi0(i,l,j)=1.d0
!C			elseif (deloh(l,j,m,k).gt.0d0.and.rngoh(i,j).ne.0d0.and.
!C     @		mass(j).eq.2)then
!C		xohi0(i,l,j)=(-1d0+dsqrt(1d0+4d0*xnohi0(i,j)*deloh(l,j,m,k)))/
!C     @			(2d0*xnohi0(i,j)*deloh(l,j,m,k))
!C	
!C			end if
!C		 end do
!C	end do
!C	end do
!C		end do
!c	end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c****************************calculo de xohi0(esta implementaci�n permite que una mol�cula
!c tenga m�s de un grupo asociativo 14/07/06**************************************
!c !Inicializaci�
     
      if(nga.ne.0) then
      xohi0=1.D0
      del1=1.D0
!c Iteraciones con tolerancia de 10-12
      do while (del1>1.0d-10)
      xohi0_old=xohi0
	do m=1, nga
	if	(rngoh(i,m).gt.0d0) then
		do j=1,2
      sum3=0.D0
	do k=1, nga
	sum4=0.D0
	do l=1,2 
	sum4=sum4+ xohi0_old(i,l,k)*deloh(l,k,j,m)*xnohi0(i,k)
	end do
	sum3=sum3+sum4
      end do
      xohi0(i,j,m)=1.D0/(1.D0+sum3)    
 	dif1(i,j,m)=dabs((xohi0(i,j,m)-xohi0_old(i,j,m))/xohi0(i,j,m))
	end do
	else
	end if
	end do
	del1=maxval(dif1)
      end do
      end if

!c*****************************fin del calculo de xohi0**************************************
!C�lculo del gama de asociaci�n ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C actas(M) = LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I  

       if(nga.ne.0) then	

      SUMAJ = 0.D0
      DO J=1,NGA 
      IF(MASS(J).NE.0) THEN      
      DO K=1,MASS(J)
	If(XOH(K,J).gt.1d-13)then
      SUMAJ = SUMAJ + RNGOH(i,j)*(dlog(XOH(K,J)/XOHI0(I,K,J))+0.5D0*(XOHi0(i,K,J)-1))+0.5D0*R(i)*xnoh(j)*(1-xoh(k,j))
	end if
      END DO
      ELSE
      CONTINUE
      END IF
      END DO
      actas(I) = SUMAJ
	end if
      act(i)=act(i)+actas(i)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   60 continue
      NDUM=0                                                            
      IF(NOVAP.EQ.0) GOTO 69                                            
      SS=0                                                              
      DO 61 I=1,NK                                                      
   61 SS=SS+X(I)*(PRAT(I)-ACT(I))                                       
      IF(SS.GT.0.) GOTO 69                                              
      NDUM=1                                                            
      DO 62 I=1,NK                                                      
      ACT(I)=PRAT(I)                                                    
      DO 62 J=1,NK                                                      
   62 DACT(I,J)=0.                                                      
      GOTO 100                                                          
   69 CONTINUE                                                          
      IF(NDIF.EQ.4) GOTO 90                                             
      IF(NDIF.LT.2) GOTO 100                                            
      DO 70 I=1,NK                                                      
      DO 70 J=I,NK                                                      
      XX=Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))         
      DO 75 K=1,NG                                                      
   75 XX=XX+U(K,I)*(V(K,J)-QT(K,J))-U(K,J)*QT(K,I)                      

!********************************calculo de dxkdni Alfonsina**************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!cCalcula los elementos de la matriz deltapq para el c�lculo de la derivada de la fracci�n 
!c no asociada respecto a la fracci�n molar del componente
      psin=0.0d0
	if(nga.ne.0) then
	m_lambda1=0.0d0
	m_lambda=0.0d0
      z=0; y=0
	do n=1,2
	do m=1,nga
     	z=z+1
	do l=1, 2
	do k=1, nga
		y=y+1
      m_lambda(z,y)=xnoh(k)*deloh(l,k,n,m)*xoh(n,m)**2 
	   if (z.eq.y)  then
      m_lambda(z,y)=m_lambda(z,y)+ 1.0d0
      end if
	end do
	end do
	y=0
	end do
	end do
      order=nga*2
	end if 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!c Calculo de  los elementos de la matriz [Yp] para el c�lculo de la derivada de la fracci�n 
!c no asociada respecto a la fracci�n molar del componente
      if(nga.ne.0) then
      do k=1,nga

	do ll=1,2
	do m=1,nga
	drhodni(j,ll,m)=(((rngoh(j,m)*xgam-xnoh1(m)*R(j)))/xgam**2)	  
	end do
	end do

      z=0     
	do ll=1,2
	do m=1,nga
	sum3=0.0d0
      do l=1,nga
	sum4=0.0d0


	do kk=1,2
	sum4=sum4+ (xoh(kk,l)*deloh(kk,l,ll,m))*drhodni(j,kk,l)
	end do
	sum3=sum3+sum4
	end do
	z=z+1
	psin(z)=-(xoh(ll,m)**2)*sum3
	end do
	end do


      N=order
	NP=order
      m_lambda1=m_lambda
       call ludcmp(m_lambda1,N,NP,indx,d1)
       call lubksb(m_lambda1,N,NP,indx,psin)
!c colectando las derivadas en su correspondiente sub�ndice
      z=0
	do m=1,2
	do l=1, nga
	z= z+1
	dxkdni(k,m,l)=psin(z)
	end do
	end do 
	end do


	do l=1,nga
	do m=1,2
	do kk=1,nga
	if (rngoh(i,kk).ne.0) then
      dxkdnic(j,m,l)=dxkdni(kk,m,l)   
	end if
	end do
      end do
	end do

!c fin del c�lculo de la derivada de la fracci�n no asociada respecto a la 
!c fracci�n molar del componente
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C dgasdx(M) = derivada LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I
      if(nga.ne.0) then
   
      DO l=1,NGA 
	   SUMA1J = 0.D0

      IF(MASS(l).NE.0) THEN      
      DO K=1,MASS(l)
	If(XOH(K,l).gt.1d-13)then

      SUMA1J = SUMA1J + RNGOH(i,l)*1.D0/XOH(K,l)*dxkdnic(j,k,l)+0.5D0*&
     r(i)*(drhodni(j,k,l)-xnoh(l)*dxkdnic(j,k,l)-drhodni(j,k,l)*&
     XOH(K,l))


	end if
      END DO
      ELSE
      CONTINUE
      END IF
	xx=xx+ SUMA1J
 	end do   
	end if
		end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DACT(I,J)=XX                                                      
   70 DACT(J,I)=XX    
       IF(NDIF.LT.3) GOTO 100                                           
      DO 80 I=1,NK                                                      
      GAM(I)=DEXP(ACT(I))                                               
   80 ACT(I)=GAM(I)*X(I)                                                
      DO 85 I=1,NK                                                      
      DO 85 J=1,NK                                                      
      DACT(I,J)=ACT(I)*(DACT(I,J)-1.D0)                                 
      IF(J.EQ.I)DACT(I,J)=DACT(I,J)+GAM(I)                              
   85 CONTINUE                                                          
      GOTO 100                                                          
   90 CONTINUE                                                          
      DO 91 I=1,2                                                       
      DO 91 K=1,2                                                       
      DTAU(I,K,K)=0.                                                    
      DO 91 L=1,2                                                       
      IF(L.EQ.K) GOTO 91                                                
      H1=TETAR(L)-QT(L,I)*ETA(L)/S(L,I)                                 
      H2=QT(K,I)-S(L,I)*TETAR(K)/ETA(L)                                 
      DTAU(I,K,L)=-H1*H2/ETA(L)                                         
   91 CONTINUE                                                          
      DO 92 I=1,NK                                                      
      PACT(I,1)=-DTAU(I,1,2)*TAU(1,2)/T*300.D0                          
   92 PACT(I,2)=-DTAU(I,2,1)*TAU(2,1)/T*300.D0                          
  100 RETURN                                                            
endsubroutine unifaclle  

     
 


