    subroutine UNIFAC (Compuestos,x,NC,T,Gamma)
!-------------------------------------------------------------------------
!   Esta subrutina calcula el coeficiente de actividad de cada compuesto 
!   en una mezcla utilizando el modelo UNIFAC (asoc == FALSE) o A-UNIFAC
!   (asoc == TRUE).
!   - Variables de entrada
!       Compuestos: vector de dimensión NC que contiene los compuestos 
!           de la mezcla
!       x: vector de dimensión NC que contiene la composición de la mezcla
!       NC: cantidad de compuestos en la mezcla
!       T: temperatura de la mezcla
!   - Variables de salida
!       Gamma: vector de dimensión NC que contiene los coef. de actividad
!-------------------------------------------------------------------------
    use CONSTANTES
    use SubGrupos
    use Input
    implicit none
!ENTRADA/SALIDA
!...Integers
    integer,intent(in)::NC
    integer,intent(in)::Compuestos(NCOM,DiffStructGroups,2)
!...real*8
    real*8, intent(in)::x(NC),T
    real*8, intent(out)::Gamma(NC)
!VARIABLES INTERNAS
    type(Groups),pointer::group
!...Integers
   ! integer::GrupMezc(NMG)
    integer::NPUNT(NMG),NGRUP(NMG),NPINT(NINT),NINTT(NINT)
    integer::i,k
!...real*8
    real*8::r(NC),q(NC),PhiTot,thetatot
    real*8:: xi(NC),J(NC), L(NC),lnGammaComb(NC),LnGammaRes(NC),LnGammaAs(NC),lnGamma(NC)
    real*8:: lnGammaGrupMezcla,lnGammaGrupPuro
    real*8:: RR(NMG),QQ(NMG),A(NMG,NMG)
!...Commons
   ! common /PUNSUB/ NPUNT(NMG),NGRUP(NMG),NCANT
    !COMMON /INTER/A(NMG,NMG)
    common /GEOM/ Rx,Qx,main
    real*8::Rx(NMG),Qx(NMG)
    integer::main(NMG)
    common/AS/asoc
    logical::asoc    
    
    
    
!...Sentencias
    ngrup(:) = 0
    npunt(:) = 0
    npint(:) = 0
    nintt(:) = 0
!...Grupos presentes en la mezcla
    Do i=1,NC
        k=1
        Do while (compuestos(i,k,1)/=0)
            call CR_PUNTF (compuestos(i,k,1),NPUNT,NGRUP,NPINT,NINTT)
            k=k+1
        enddo
    enddo
    call LEER_PRnew (ipareq,ngrup,RR,QQ)
    call Leer_Int_Par (ipareq,nintt,A)
!***Término combinatorial***
    r(:)=0.0
    q(:)=0.0
    do i=1, NC !volumen y superficie de cada componente
        k=1
        do while (Compuestos(i,k,1)/=0)
            r(i) = r(i) + Compuestos(i,k,2) * Obtain_R(Compuestos(i,k,1))   !R(npunt(Compuestos(i,k,1))) !volumen
            q(i) = q(i) + Compuestos(i,k,2) * Obtain_Q(Compuestos(i,k,1))  !Q(npunt(Compuestos(i,k,1))) !superficie
            k=k+1
        enddo
    enddo
    PhiTot=0.0 
    ThetaTot=0.0 !Fi(:)=0.0, ThetaCom(:)=0.0
    do i=1,NC !Sumatoria de volumenes y superficies
        PhiTot = PhiTot + x(i) * r(i)
        ThetaTot = ThetaTot + x(i) * q(i)
    enddo
!    do i=1, NC !Fracciones de volumen (Fi) y de área (ThetaCom)
!        Fi(i)=x(i)*r(i)/PhiTot  
!        ThetaCom(i)=x(i)*q(i)/ThetaTot    
!    enddo
    J(:)=0.0
    L(:)=0.0
    do i=1,NC !Fracciones de volumen (J) y superficie (L) de cada componente
              !No se multiplica por la fracción molar (x) para simplificar la ec. 2) Fredenslund et al. (1975)
        J(i)= r(i)/PhiTot
        L(i)= q(i)/ThetaTot
    enddo
    do i=1,NC !Contribución combinatorial ec.(2) Fredenslund et al. (1975)
        lnGammaComb(i)= 1 - J(i) + log(J(i)) - 5*q(i)*(1-(J(i)/L(i))+(log(J(i)/L(i))))
    enddo
!***Término Residual***
!...Coeficiente de actividad residual grupal
    lnGammaRes(:)=0.0
    Do i=1,nc
        xi(:)=0.0
        xi(i)=1.0
        k=1
        Do while (compuestos(i,k,1)/=0)
            lnGammaGrupPuro=0.0
            lnGammaGrupMezcla=0.0
            call GammaResGrup(compuestos(i,k,1),xi,compuestos,nc,T,NGRUP,NPINT,NINTT,A,RR,QQ,lnGammaGrupPuro) 
            call GammaResGrup(compuestos(i,k,1),x,compuestos,nc,T,NGRUP,NPINT,NINTT,A,RR,QQ,lnGammaGrupMezcla)    
            lnGammaRes(i)= lnGammaRes(i)+(compuestos(i,k,2)*(lnGammaGrupMezcla-lnGammaGrupPuro))
            k=k+1
        enddo
    enddo
!***Término Asociativo***
    lnGammaAs(:)=0.0
    if(asoc) call GammaAsociativo(compuestos,nc,T,x,lnGammaAs,0)

!***Gamma total***
    lnGamma(:) = lnGammaRes(:) + lnGammaComb(:) + lnGammaAs(:)
    Gamma(:)=exp(lnGamma(:))
 1  continue 

endsubroutine UNIFAC
    


Subroutine GammaResGrup (numgrup,x,Compuestos,NC,T,NGRUP,NPINT,NINTT,A,RR,QQ,lnGammaGrup)
!-------------------------------------------------
!   Descripción
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    use CONSTANTES
    use SubGrupos
    implicit real*8(A-H,O-Z)
!...Pointers
    type(Groups),pointer::group
!...Integers
    integer,intent(in)::numgrup
    integer,intent(in)::Compuestos(NCOM,DiffStructGroups,2)
    integer,intent(in)::NGRUP(NMG),NPINT(NINT),NINTT(NINT)
    real*8, intent(in):: RR(NMG),QQ(NMG)
!...real*8
    real*8, intent(in)::x(NC),T,A(NMG,NMG)
    real*8, intent(out):: lnGammaGrup
    real*8:: Xk(NMG),ThetaGrup(NMG),Psi(NMG,NMG)
!...Characters

!...Logical

!...Commons
    !common /PUNSUB/ NPUNT(NMG),NGRUP(NMG),NCANT
    !common /INTER/ A(NMG,NMG)
    common /GEOM/ Rx(NMG),Qx(NMG),MAIN(NMG)
!...Sentencias
    ThetaGrup(:)=0.0
    Psi(:,:)=0.0
    Xk(:)=0.0
!...Fracción molar del grupo k en la mezcla
!...Denominador
    Xtot=0.0
    do i=1,NC 
        k=1
        do while (Compuestos(i,k,1)/=0)
            Xtot=Xtot+x(i)*Compuestos(i,k,2)
            k=k+1
        enddo    
    enddo
!...Numerador y resultado
    Xk(:)=0.0
    k=1
    Do while (ngrup(k)/=0)
        Do i=1,NC
            ki=1
            Do while (Compuestos(i,ki,1)/=0)
                if (compuestos(i,ki,1)==ngrup(k)) then
                    Xk(k)=Xk(k)+compuestos(i,ki,2)*x(i)
                    exit
                endif
                ki=ki+1
            enddo
        enddo
    Xk(k)=Xk(k)/Xtot !Fracción molar grupo k
    k=k+1
    enddo
!...Fracción de superficie del grupo k en la mezcla
    ThetaTot=0.0 !
    k=1
    Do while (ngrup(k)/=0) !Denominador
        ThetaTot=ThetaTot + Obtain_Q(ngrup(k)) * Xk(k)
        k=k+1
    enddo   
    ThetaGrup(:)=0.0
    k=1
    Do while (ngrup(k)/=0)
        ThetaGrup(k)=(Obtain_Q(ngrup(k)) * Xk(k)) / ThetaTot
        k=k+1
    enddo
    m=1
    Do while (ngrup(m)/=0) !Cálculo de Psi
        n=1
        Do while (ngrup(n)/=0)
            nm=npint(Obtain_MainGroup_Number(ngrup(m))) !lugar que ocupa el grupo ppal m
            nn=npint(Obtain_MainGroup_Number(ngrup(n))) !lugar que ocupa el grupo ppal n
            Psi(nm,nn)=exp(-(A(nm,nn)/T))
            n=n+1
        enddo
        m=m+1
    enddo
!   Mezcla
    lnGammaGrup=0.0

    Term1 = 0.0     ! Term1 = Término ln(Sumatoria de grupos(ThetaGrup*Psi))
    Vtemp = 0.0     ! Variable temporal
    m=1
    Do while (ngrup(m)/=0)
        nm = npint(Obtain_MainGroup_Number(ngrup(m))) !lugar que ocupa el grupo ppal m
        nngrup = npint(Obtain_MainGroup_Number(numgrup)) !lugar que ocupa el grupo ppal numgrup
        Vtemp = Vtemp + (ThetaGrup(m)*Psi(nm,nngrup))
        m=m+1
    enddo
    Term1=log(Vtemp)
    Term2 = 0.0     !Term2 = Término Sumatoria (ThetaGrup*Psi/Sumatoria (ThetaGrup*Psi))
    m=1
    Do while (ngrup(m)/=0)
        Vtemp=0.0
        nm=npint(Obtain_MainGroup_Number(ngrup(m))) !lugar que ocupa el grupo ppal m
        n=1
        Do while (ngrup(n)/=0)
            nn=npint(Obtain_MainGroup_Number(ngrup(n))) !lugar que ocupa el grupo ppal n
            Vtemp=Vtemp+ThetaGrup(n)*Psi(nn,nm)
            n=n+1
        enddo
        nngrup=npint(Obtain_MainGroup_Number(numgrup)) !lugar que ocupa el grupo ppal numgrup
        Term2=Term2+((ThetaGrup(m)*Psi(nngrup,nm))/Vtemp)
        m=m+1
    enddo
    lnGammaGrup = Obtain_Q(numgrup) * (1-Term1-Term2)  
endsubroutine

!=================================================================================    
SUBROUTINE UNIFA(NDIF,NACT,NC,NG,T,X,ACT,DACT,TACT)
!---------------------------------------------------------------------------------
!     X(NCOM)       =   fracción molar de cada componente.
!     QT(NGPM,NCOM) =   suma de la superficie de todos los subgrupos de cada
!                       grupo ppal en el componente.
!     RT(NGPM,NCOM) =   suma del volumen de todos los subgrupos de cada
!                       grupo ppal en el componente.
!     Q(NCOM)       =   área del componente
!     R(NCOM)       =   volumen del componente
!     THETA(NCOM)   =   fracción de área del componente en la mezcla
!     PHI(NCOM)     =   fracción de volumen del componente en la mezcla
!     THETS         =   sumatoria área de los componentes en la mezcla
!     PHS           =   sumatoria volúmenes de los componentes en la mezcla
!     RI(NCOM)      =   fracción de volumen/fracción molecular
!---------------------------------------------------------------------------------

      PARAMETER(NMG=150,NGPM=30,NGA=70,NCOM=3,NSCM=10)
      implicit real*8(A-H,O-Z)
      COMMON/UNIF/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
      COMMON/AS/ASOC
      DIMENSION X(NCOM),GAM(NCOM),ACT(NCOM),DACT(NCOM,NCOM),THETA(NCOM),&
               PHI(NCOM),RI(NCOM),QI(NCOM),QIL(NCOM),RIL(NCOM),&
               QID(NCOM),ETAL(NGPM),TACT(NCOM),U(NGPM,NCOM),V(NGPM,NCOM),ACTAS(NCOM)
      DIMENSION  DETA(NGPM),DS(NGPM,NCOM),ETA(NGPM),TETAR(NGPM),H3(NGPM,NCOM)  
      LOGICAL ASOC
      THETS=0.D0                                                              
      PHS=0.D0                                                                
      DO 10 I=1,NC                                                            
		THETA(I)=X(I)*Q(I)   !numerador ec.(5) -área-
		PHI(I)=R(I)*X(I)     !numerador ec.(4) -volumen-
		THETS=THETS+THETA(I) !Denominador ec.(5)
   10		PHS=PHS+PHI(I)       !Denominador ec.(4)
      DO 20 I=1,NC                                                            
		THETA(I)=THETA(I)/THETS !ec.(5)                    
		PHI(I)=PHI(I)/PHS       !ec.(4)
		RI(I)=R(I)/PHS                                                          
		RIL(I)=DLOG(RI(I))    !1er término ec.(1)
		QI(I)=Q(I)/THETS                                                        
		QID(I)=1.-RI(I)/QI(I)                                                   
   20		QIL(I)=DLOG(QI(I))                                                      
      DO 30 I=1,NC                                                            
		XX=F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                                   
		XX=XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                                    
   30		GAM(I)=XX                                                               
      DO 40 I=1,NG                                                            
		TETAR(I)=0.                                                             
		ETA(I)=0.                                                               
		DO 45 J=1,NC                                                            
			ETA(I)=ETA(I)+S(I,J)*X(J)                                               
   45			TETAR(I)=TETAR(I)+QT(I,J)*X(J)                                          
   40		ETAL(I)=DLOG(ETA(I))                                                    
      DO 50 I=1,NC                                                            
		DO 60 J=1,NG                                                            
			U(J,I)=S(J,I)/ETA(J)                                                    
			V(J,I)=U(J,I)*TETAR(J)                                                  
   60			GAM(I)=GAM(I)-V(J,I)-QT(J,I)*ETAL(J)                                    
		ACT(I)=DEXP(GAM(I))                                                     
   50	  IF(NACT.EQ.1) ACT(I)=ACT(I)*X(I)                                        
      IF(NDIF.EQ.0) GOTO 1000                                                 
      IF (NDIF.EQ.2) GO TO 81                                                 
      DO 70 I=1,NC                                                            
		DO 70 J=I,NC                                                            
			XX=Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))               
			DO 75 K=1,NG                                                            
   75				XX=XX+U(K,I)*(V(K,J)-QT(K,J))-U(K,J)*QT(K,I)                           
			DACT(I,J)=XX                                                            
			DACT(J,I)=XX                                                            
			IF(NACT.EQ.1) GOTO 70                                                   
			DACT(I,J)=DACT(I,J)*ACT(I)                                              
			IF (J.EQ.I) GO TO 70                                                    
			DACT(J,I)=DACT(J,I)*ACT(J)                                              
   70 CONTINUE                                                                
      IF(NACT.EQ.0) GOTO 81                                                   
      DO 80 I=1,NC                                                            
		DO 80 J=1,NC                                                            
			DACT(I,J)=ACT(I)*(DACT(I,J)-1.D0)                                       
			IF(J.EQ.I) DACT(I,J)=DACT(I,J)+DEXP(GAM(I))                             
   80 CONTINUE                                                                
   81 CONTINUE                                                                
      IF(NDIF.EQ.1) GOTO 1000                                                 
      DO 150 K=1,NG                                                           
		DETA(K)=0.D0                                                            
		DO 150 I=1,NC                                                           
			DS(K,I)=0.D0                                                            
			DO 151 M=1,NG                                                           
				IF(QT(M,I).EQ.0.D0) GOTO 151                                            
				DS(K,I)=DS(K,I)-QT(M,I)*DLOG(TAU(M,K))*TAU(M,K)/T                       
  151			CONTINUE                                                                
  150 DETA(K)=DETA(K)+DS(K,I)*X(I)                                            
      DO 152 I=1,NC                                                           
		TACT(I)=0.D0                                                            
		DO 153 K=1,NG                                                           
			H3(K,I)=(-S(K,I)*DETA(K)/ETA(K)+DS(K,I))/ETA(K)                         
			HH=H3(K,I)*(TETAR(K)-QT(K,I)*ETA(K)/S(K,I))                             
  153			TACT(I)=TACT(I)-HH                                                      
  152 TACT(I)=TACT(I)*ACT(I)                                                  
 1000 CONTINUE    
      IF (ASOC) THEN
       ! CALL GammaAsociativo (NC,T,X,ACTAS,NACT)
            DO I=1,NC
                ACT(I)=ACT(I)*ACTAS(I)
            ENDDO
      ENDIF                                                            
      RETURN                                                                  
endsubroutine

SUBROUTINE UNIFA2(NDIF,NACT,NC,NG,T,X,ACT,DACT,TACT)                     
      PARAMETER(NMG=150,NGPM=30,NGA=70,NCOM=3,NSCM=10)
      implicit real*8(A-H,O-Z)
      COMMON/UNIF/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
      DIMENSION X(NCOM),GAM(NCOM),ACT(NCOM),DACT(NCOM,NCOM),THETA(NCOM),&
               PHI(NCOM),RI(NCOM),QI(NCOM),QIL(NCOM),RIL(NCOM),&
               QID(NCOM),ETAL(NGPM),TACT(NCOM),U(NGPM,NCOM),V(NGPM,NCOM)                                                   
      DIMENSION  DETA(NGPM),DS(NGPM,NCOM),ETA(NGPM),TETAR(NGPM),H3(NGPM,NCOM)               
      THETS=0.D0                                                              
      PHS=0.D0                                                                
      DO 10 I=1,NC                                                            
      THETA(I)=X(I)*Q(I)                                                      
      PHI(I)=R(I)*X(I)                                                        
      THETS=THETS+THETA(I)                                                    
   10 PHS=PHS+PHI(I)                                                          
      DO 20 I=1,NC                                                            
      THETA(I)=THETA(I)/THETS                                                 
      PHI(I)=PHI(I)/PHS                                                       
      RI(I)=R(I)/PHS                                                          
      RIL(I)=DLOG(RI(I))                                                      
      QI(I)=Q(I)/THETS                                                        
      QID(I)=1.-RI(I)/QI(I)                                                   
   20 QIL(I)=DLOG(QI(I))                                                      
      DO 30 I=1,NC                                                            
      XX=F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                                   
      XX=XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                                    
   30 GAM(I)=XX                                                               
      DO 40 I=1,NG                                                            
      TETAR(I)=0.                                                             
      ETA(I)=0.                                                               
      DO 45 J=1,NC                                                            
      ETA(I)=ETA(I)+S(I,J)*X(J)                                               
   45 TETAR(I)=TETAR(I)+QT(I,J)*X(J)                                          
   40 ETAL(I)=DLOG(ETA(I))                                                    
      DO 50 I=1,NC                                                            
      DO 60 J=1,NG                                                            
      U(J,I)=S(J,I)/ETA(J)                                                    
      V(J,I)=U(J,I)*TETAR(J)                                                  
   60 GAM(I)=GAM(I)-V(J,I)-QT(J,I)*ETAL(J)                                    
      ACT(I)=DEXP(GAM(I))                                                     
   50 IF(NACT.EQ.1) ACT(I)=ACT(I)*X(I)                                        
      IF(NDIF.EQ.0) GOTO 1000                                                 
      IF (NDIF.EQ.2) GO TO 81                                                 
      DO 70 I=1,NC                                                            
      DO 70 J=I,NC                                                            
      XX=Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))               
      DO 75 K=1,NG                                                            
   75 XX=XX+U(K,I)*(V(K,J)-QT(K,J))-U(K,J)*QT(K,I)                            
      DACT(I,J)=XX                                                            
      DACT(J,I)=XX                                                            
      IF(NACT.EQ.1) GOTO 70                                                   
      DACT(I,J)=DACT(I,J)*ACT(I)                                              
      IF (J.EQ.I) GO TO 70                                                    
      DACT(J,I)=DACT(J,I)*ACT(J)                                              
   70 CONTINUE                                                                
      IF(NACT.EQ.0) GOTO 81                                                   
      DO 80 I=1,NC                                                            
      DO 80 J=1,NC                                                            
      DACT(I,J)=ACT(I)*(DACT(I,J)-1.D0)                                       
      IF(J.EQ.I) DACT(I,J)=DACT(I,J)+DEXP(GAM(I))                             
   80 CONTINUE                                                                
   81 CONTINUE                                                                
      IF(NDIF.EQ.1) GOTO 1000                                                 
      DO 150 K=1,NG                                                           
      DETA(K)=0.D0                                                            
      DO 150 I=1,NC                                                           
      DS(K,I)=0.D0                                                            
      DO 151 M=1,NG                                                           
      IF(QT(M,I).EQ.0.D0) GOTO 151                                            
      DS(K,I)=DS(K,I)-QT(M,I)*DLOG(TAU(M,K))*TAU(M,K)/T                       
  151 CONTINUE                                                                
  150 DETA(K)=DETA(K)+DS(K,I)*X(I)                                            
      DO 152 I=1,NC                                                           
      TACT(I)=0.D0                                                            
      DO 153 K=1,NG                                                           
      H3(K,I)=(-S(K,I)*DETA(K)/ETA(K)+DS(K,I))/ETA(K)                         
      HH=H3(K,I)*(TETAR(K)-QT(K,I)*ETA(K)/S(K,I))                             
  153 TACT(I)=TACT(I)-HH                                                      
  152 TACT(I)=TACT(I)*ACT(I)                                                  
 1000 CONTINUE                                                                
      RETURN                                                                  
endsubroutine

subroutine UNIFA3(NG,T,NDIF,X,ACT,DACT,PACT)                      
	PARAMETER(NGPM=30,NCOM=3)
      implicit real*8(A-H,O-Z)                                  
      COMMON/UNIF/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
      DIMENSION X(10),GAM(10),ACT(10),DACT(10,10),THETA(10),PHI(10),RI(10),QI(10),ETA(10),QIL(10),RIL(10),QID(10),ETAL(10),TETAR(10) 
      DIMENSION U(10,10),V(10,10),PACT(2,2),DTAU(2,2,2)                 
      NK=2                                                              
      THETS=0.                                                          
      PHS=0.                                                            
      DO 10 I=1,NK                                                      
        THETA(I)=X(I)*Q(I) !superficie                                  
        PHI(I)=R(I)*X(I)   !volumen                                     
        THETS=THETS+THETA(I)!Denominador para cálculo fracción de area  
   10   PHS=PHS+PHI(I)     !Denominador para cálculo fracción de volumen
      DO 20 I=1,NK                                                      
        RI(I)=R(I)/PHS  !Fraccion de volumen                            
        RIL(I)=DLOG(RI(I))                                              
        QI(I)=Q(I)/THETS !Fracción de superficie                        
   20   QIL(I)=DLOG(QI(I))                                              
      DO 40 I=1,NG                                                      
        ETA(I)=0.                                                       
        DO 45 J=1,NK                                                    
   45       ETA(I)=ETA(I)+S(I,J)*X(J)                                   
   40   ETAL(I)=DLOG(ETA(I))                                            
      DO 55 I=1,NG                                                      
        TETAR(I)=0.                                                     
        DO 55 J=1,NK                                                    
   55       TETAR(I)=TETAR(I)+QT(I,J)*X(J)                              
      DO 60 I=1,NK                                                      
        QID(I)=1.-RI(I)/QI(I)                                           
        XX=F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                           
        XX=XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                            
        ACT(I)=XX                                                       
        DO 60 J=1,NG                                                    
            U(J,I)=S(J,I)/ETA(J)                                        
            V(J,I)=U(J,I)*TETAR(J)                                      
   60       ACT(I)=ACT(I)-V(J,I)-QT(J,I)*ETAL(J)                        
      IF(NDIF.EQ.4) GOTO 90                                             
      IF(NDIF.LT.2) GOTO 100                                            
      DO 70 I=1,NK                                                      
        DO 70 J=I,NK                                                    
            XX=Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))   
            DO 75 K=1,NG                                                
   75           XX=XX+U(K,I)*(V(K,J)-QT(K,J))-U(K,J)*QT(K,I)            
            DACT(I,J)=XX                                                
   70       DACT(J,I)=XX                                                
      IF(NDIF.LT.3) GOTO 100                                            
      DO 80 I=1,NK                                                      
        GAM(I)=DEXP(ACT(I))                                             
   80   ACT(I)=GAM(I)*X(I)                                              
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
   92   PACT(I,2)=-DTAU(I,2,1)*TAU(2,1)/T*300.D0                        
  100 RETURN                                                            
endsubroutine   


SUBROUTINE GammaAsociativo (MS,NC,T,X,ACTAS,NACT)
!C*********************************************************************
!c     Calcula la contribución asociativa al coeficiente de actividad.
!c     [1] Basado en el trabajo de Ferreirab,Macedob,Bottini "Extension 
!c     of the A-UNIFAC model to mixtures of cross- and self-associating 
!c     compounds". Fluid Phase Equilibria 227 (2005) 165–176.
!C*********************************************************************      
      use Input
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER(NMS=2,NSCM=10)
      DIMENSION MS(NC,NSCM,2)
      DIMENSION xnoh1(12),NPUNT(NMG),NGRUP(NMG),&
     RR(NMG),NPUNTA(NMG),NGRUPA(NMG),ENASS(6,12,6,12),&
     RKASS(6,12,6,12),deloh(6,12,6,12),R(NC),xnohi0(12,12),&
     X(NC),xnoh(12),dif(12,12),xoh_old(6,12),xoh(6,12),&
     xohi0(12,6,12),xohi0_old(12,6,12),xxohi0(12,6,12), dif1(10,12,12),&
     act(NC),actas(NC),NGRUPASD(NMG),NPUNTASD(NMG)
      CHARACTER FS*8,TG*3
      INTEGER A,B,NMAINGR(NMG),rngoh(NC,NSCM),NPUNTMG(NMG),GA,&
             CANIL(10),CS(NMG),TS(NMG,NMS)
      LOGICAL VL
!C-----INICIALIZACIÓN VARIABLES
      xnoh1(:)=0.
      NPUNT(:)=0
      NGRUP(:)=0
      RR(:)=0.
      NPUNTA(:)=0
      NGRUPA(:)=0
      ENASS(:,:,:,:)=0.
      RKASS(:,:,:,:)=0.
      DELOH(:,:,:,:)=0.
      R(:)=0.
      xnohi0(:,:)=0.
      xnoh(:)=0.
      dif(:,:)=0.
      xoh_old(:,:)=0.
      xoh(:,:)=0.
      xohi0(:,:,:)=0.
      xohi0_old(:,:,:)=0.
      xxohi0(:,:,:)=0.
      dif1(:,:,:)=0.
      act(:)=0.
      actas(:)=0.
      CS(:)=0
      NGRUPASD(:)=0
      NPUNTASD(:)=0
      TS(:,:)=0
      NMAINGR(:)=0
      rngoh(:,:)=0
      NPUNTMG(:)=0
      CANIL(:)=0
      CS(:)=0
      TS(:,:)=0
!C-----VARIABLES EN ENTRADA
!C.....Averigua cantidad de anillos aromáticos por componente
      DO J=1,NC
        GA=0
        DO I=1,10
            IF(MS(J,I,1).EQ.0)CYCLE
            IREC2=MS(J,I,1)+(ipareq-1)*150
            READ(14,501,REC=IREC2)TG
 501        FORMAT(28X,A3)
            IF(TG.EQ."AR1")GA=GA+MS(J,I,2)
        ENDDO
        CANIL(J)=GA/6
      ENDDO
!C.....Crea el vector NPUNT-NGRUP y NPUNTMG-NMAINGR (p/ grupos ppales)
      NUM = 0
      NMGR = 0
      NGA = 0
	DO J=1,NC
	  DO I=1,10 
	      IF(MS(J,I,1).EQ.0)CYCLE
	      IREC2=MS(J,I,1)+(ipareq-1)*150      
      	    READ(14,500,REC=IREC2)MGR,RRT,ICS,ITS1,ITS2 !500  FORMAT(2x,I2,37X,D15.8,120X,3I2)
            IF (ICS.NE.0) THEN !GRUPOS ASOCIATIVOS
                NG = MS(J,I,1)
                CALL BUSCARAS (NG,NGRUPA,NMG,VL)
                IF(VL)THEN
                    IF(MS(J,I,1).EQ.10) THEN
                        RNGOH(J,NPUNTA(NG))=CANIL(J)
                    ELSE
                        RNGOH(J,NPUNTA(NG))=MS(J,I,2)
                    ENDIF                  
                ELSE
                    NGA=NGA+1
                    NPUNTA(NG) = NGA
                    NGRUPA(NGA) = NG
                    CS(NGA) = ICS
                    TS(NGA,1) = ITS1
                    TS(NGA,2) = ITS2
                    IF(MS(J,I,1).EQ.10) THEN
                        RNGOH(J,NGA)=CANIL(J)
                    ELSE
                        RNGOH(J,NGA)=MS(J,I,2)
                    ENDIF
                ENDIF
            ENDIF      	    
      	    IF (NPUNT(MS(J,I,1)).EQ.0) THEN !Ordena todos los grupos
                NUM = NUM + 1
                NPUNT(MS(J,I,1)) = NUM
                NGRUP(NUM) = MS(J,I,1)
                RR(NUM)=RRT
                CALL BUSCARAS (MGR,NMAINGR,NMG,VL)
                IF(.NOT.VL) THEN !Crea vectores NPUNTMG y NMAINGR
                    NMGR = NMGR+1
                    NPUNTMG(MGR) = NMGR
                    NMAINGR(NMGR) = MGR
                ENDIF
            END IF
        ENDDO
      ENDDO
!C....................................................................................
!C.....Genera las matrices ENASS Y RKASS con los parámetros de energía y volumen 
!c.....de asociación, respectivamente, según los grupos asociativos de los componentes
!c.....del sistema que se esté corriendo.
!C....................................................................................     
      DO J=1,nga
        DO K=1,nga
            DO B=1,CS(J)
                DO A=1,CS(K)
                    IREC1 = MAINSG(NGRUPA(K),ipareq)+(ipareq-1)*70
                    IF(TS(J,B).EQ.1.OR.TS(K,A).EQ.1)THEN !Si alguno de ambos sitios es del tipo 1
                       CALL LEEPAR (J,IREC1,ipareq,NGRUPA,ENASST,RKASST)
                       ENASS(A,K,B,J)=ENASST
                       RKASS(A,K,B,J)=RKASST                        
                    ELSEIF (TS(J,B).NE.TS(K,A)) THEN
                       CALL LEEPAR (J,IREC1,ipareq,NGRUPA,ENASST,RKASST)
                       ENASS(A,K,B,J)=ENASST
                       RKASS(A,K,B,J)=RKASST
                    ELSE !ELSEIF (TS(J,B).EQ.TS(K,A)) THEN
                       ENASS(A,K,B,J)=0.0
                       RKASS(A,K,B,J)=0.0     
                    ENDIF   
                ENDDO
            ENDDO
        ENDDO
      ENDDO
!C
!C.....Cálculo del volumen (R) de cada componente      
      R(:)=0.0
      DO I=1,NC
        DO J=1,10
            IF(MS(I,J,1).EQ.0)CYCLE
            R(I) = R(I)+RR(NPUNT(MS(I,J,1)))*MS(I,J,2)
        ENDDO
      ENDDO
!C...................................................................
!C.....Cálculo de densidades de grupos asociativos en sus respectivos
!c.....componentes puros (EC.5 [1])
!C...................................................................
      DO I=1,NC
        DO J=1,NGA
            xnohi0(i,j)=rngoh(i,j)/R(i)  
        ENDDO
      ENDDO
!C...................................................................
!C.....Cálculo de densidades de grupos asociativos en la mezcla
!c.....(EC.4 [1])
!C...................................................................
      !DENOMINADOR
      XGAM=0.0
      DO I=1,NC
        XGAM=XGAM+R(I)*X(I)     
      ENDDO
      !NUMERADOR
      xnoh1(:)=0d0
      DO J=1,NGA
        DO I=1,NC
            DO K=1,10
                IF(NGRUPA(J).EQ.MS(I,K,1)) THEN
                    xnoh1(J)=XNOH1(J)+MS(I,K,2)*X(I)
                ENDIF
            ENDDO
        ENDDO
      ENDDO
      DO J=1,NGA
        XNOH(J)=XNOH1(J)/XGAM
      ENDDO
!C...................................................................
!C.....Cálculo de fuerzas de asociación (EC.6 [1])
!C...................................................................
 999  if(nga.ne.0) then
	  DO J=1,NGA 
            DO K=1,NGA
                DO B=1,CS(NPUNTA(NGRUPA(J)))
                    DO A=1,CS(NPUNTA(NGRUPA(K)))
      deloh(A,K,B,J)=(DEXP(ENASS(A,K,B,J)/T) - 1 )*RKASS(A,K,B,J)
                    END DO
                END DO
 101            CONTINUE
            END DO
 201        CONTINUE
        END DO
      end if
!C...................................................................
!C.....Cálculo de fracción de sitios no asociados en la mezcla 
!c.....(EC.2 [1])
!C...................................................................
!c !Inicializació
      if(nga.ne.0) then
        xoh=1.0d0
        del=1.0
!c Iteraciones con tolerancia de 10-9
        do while (del>1.0d-9)
            xoh_old=xoh
	      do J=1, nga
		        do B=1,CS(NPUNTA(NGRUPA(J)))
                    sum1=0.
	              do K=1, nga
	                  sum2=0.
	                      do A=1,CS(NPUNTA(NGRUPA(K)))
	                          sum2=sum2+xoh_old(A,K)*deloh(A,K,B,J)*xnoh(K)
	                      end do
	                  sum1=sum1+sum2
                    end do
                    xoh(B,J)=1.0d0/(1.0d0+sum1)          
	              dif(B,J)=dabs((xoh(B,J)-xoh_old(B,J))/xoh(B,J))
	          end do
	      end do
	      del=maxval(dif)
        end do
      end if
!cc Fin del Cálculo de la  fracción no asociada  
      DO 300 I=1,NC
!C.....................................................................
!C.....Cálculo de fracción de sitios no asociados en el componente puro
!c.....(EC.3 [1])
!C.....................................................................
!c !Inicializació
        if(nga.ne.0) then
            xohi0=1.0d0
            del1=1.0
!c Iteraciones con tolerancia de 10-9
            do while (del1>1.0d-9)
                xohi0_old=xohi0
	          do m=1, nga
	              if(rngoh(i,m).gt.0d0) then
		                do j=1,2
                            sum3=0.
	                      do k=1, nga
	                          sum4=0.
	                          do l=1,2 
	                 sum4=sum4+ xohi0_old(i,l,k)*deloh(l,k,j,m)*xnohi0(i,k)
	                          end do
	                          sum3=sum3+sum4
                            end do
                            xohi0(i,j,m)=1.0d0/(1.0d0+sum3)    
	                      xxohi0(i,j,m)=1.0d0/(1.0d0+sum3)          
	                      dif1(i,j,m)=dabs((xohi0(i,j,m)-xohi0_old(i,j,m))/xohi0(i,j,m))
	                  end do
		            else
	              end if
	          end do
	          del1=maxval(dif1)
            end do
        end if
!C.......................................................................
!C.....Cálculo de la contribución asociativa del coeficiente de actividad
!c.....(EC.7 [1])
!C.......................................................................
        if(nga.EQ.0) THEN
	      ACTAS (I)=0
	  END IF
        if(nga.ne.0) then	
            SUMAJ = 0d0
            DO J=1,NGA 
                IF(CS(NPUNTA(NGRUPA(J))).NE.0) THEN      
                    DO K=1,CS(NPUNTA(NGRUPA(J)))
	                  If(XOH(K,J).gt.1d-13)then
	                      SUMAJ1=RNGOH(i,j)*(dlog(XOH(K,J)/XOHI0(I,K,J))+0.5*(XOHi0(i,K,J)-1))
	                      SUMAJ2=0.5*R(i)*xnoh(j)*(1-xoh(k,j))
	                      SUMAJ=SUMAJ+SUMAJ1+SUMAJ2
	                  end if
                    END DO
                ELSE
                    CONTINUE
                END IF
            END DO
            actas(I) = SUMAJ
	      gami=actas(i)
		end if
!c-----
        IF(NACT.EQ.1) then
            actas(i)=dexp(actas(i))*x(i)
        end if
300   CONTINUE
!C-----FORMATOS
 500  FORMAT(2x,I2,37X,D15.8,120X,3I2)
 502  FORMAT(a8,70d12.5)

endsubroutine