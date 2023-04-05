subroutine open_proppdb(oprop)
!-------------------------------------------------------------------
!      Esta subrunina ABRE los bancos de datos siguientes:
!     - prop.pdb      (UNIT=7)
!-------------------------------------------------------------------      
    implicit none
    logical, intent(out)::oprop

!SENTENCES      
    oprop = .False.  
    open(unit=7,file='Database\PROP.pdb',access='direct',status='old',&
         form='formatted',recl=624,err=20)
	goto 10

20	write(6,*) '*ERROR* EN LA APERTURA DE "PROP.PDB"'
	oprop = .True.

10  continue      
endsubroutine open_proppdb
      
      
subroutine ab_ban1(mod)
!-------------------------------------------------------------------
!      Esta subrunina ABRE los bancos de datos siguientes:
!     - INTRCN.MDS      (UNIT=13)
!     - INTRCNAS.MDS    (UNIT=13)
!     - GRUPOSRAM.MDS   (UNIT=14)
!     - PARVOLAS.MDS    (UNIT=15)  
!     - PARENEAS.MDS    (UNIT=16)     
!-------------------------------------------------------------------
!
    use input_data, only:Model_S
      COMMON/AS/ASOC
      LOGICAL ASOC
      integer, intent(out),optional::mod
      
  
    CALL Model_S(mod)
    if(mod /= 3)then
        if(mod==1)then
            open (unit=13,file='Database\intrcn.mds',status='old',&
                  access='direct',form='formatted',recl=850)   
        else     
            open (unit=13,file='Database\intrcnas.mds',status='old',&
                  access='direct',form='formatted',recl=850)        
        endif
        open (unit=14,file='Database\gruposram.mds',status='old',&
     	      access='direct',form='formatted',recl=300)
        open (unit=15,file='Database\parvolas.mds',status='old',&
              access='direct',form='formatted',recl=850)
        open (unit=16,file='Database\pareneas.mds',status='old',&
              access='direct',form='formatted',recl=850)    
    else

        open (unit=14,file='Database\gruposramgc.mds',status='old',&
     	      access='direct',form='formatted',recl=263)
        open (unit=13,file='Database\intrcngcalpha.mds',status='old',&
              access='direct',form='formatted',recl=730)
        open (unit=16,file='Database\intrcngckapa.mds',status='old',&
              access='direct',form='formatted',recl=730)    
        
    endif
      

      return
      endsubroutine ab_ban1
      
      
      
      subroutine ab_ban
!c-------------------------------------------------------------------
!c      Esta subrunina ABRE los bancos de datos siguientes:
!c     - INTRCN.MDS      (UNIT=13)
!c     - GRUPOSRAM.MDS   (UNIT=14)
!c     - PARVOLAS.MDS    (UNIT=15)
!c     - PARENEAS.MDS    (UNIT=16)     
!c-------------------------------------------------------------------
!c
      COMMON/AS/ASOC
      LOGICAL ASOC
      IF(ASOC)THEN
        open (unit=13,file='Database\intrcnas.mds',status='old',&
             access='direct',form='formatted',recl=850)
      ELSE
        open (unit=13,file='Database\intrcn.mds',status='old',&
             access='direct',form='formatted',recl=850)
      ENDIF
      open (unit=14,file='Database\gruposram.mds',status='old',&
     	    access='direct',form='formatted',recl=300)
      open (unit=15,file='Database\parvolas.mds',status='old',&
             access='direct',form='formatted',recl=850)
      open (unit=16,file='Database\pareneas.mds',status='old',&
             access='direct',form='formatted',recl=850)
      return
      end

      subroutine ci_ban1
!c-------------------------------------------------------------------
!c      Esta subrunina CIERRA los bancos de datos siguientes:
!c     - INTRCN.MDS      (UNIT=13)
!c     - GRUPOSRAM.MDS   (UNIT=14)
!c     - PARVOLAS.MDS    (UNIT=15)
!c     - PARENEAS.MDS    (UNIT=16)   
!c-------------------------------------------------------------------

      close (unit=13)
      close (unit=14)
      close (unit=15)
      close (unit=16)
      return
      end



SUBROUTINE Store_Pr (IPAREQ)
!C-----------------------------------------------------------------------
!C      ESTA SUBRUTINA CARGA TODAS LOS PARAMETROS DEL BANCO GRUPOSRAM QUE
!C      SE USAN EN MOLDES. LAS PROPIEDADES SE ALMACENAN EN COMMONS.
!C-----------------------------------------------------------------------
!C
    use SubGrupos
       PARAMETER (NMG=150)
       IMPLICIT real*8 (A-H,O-Z)
       CHARACTER*8 FS(NMG),NAME
       INTEGER H1
       COMMON /PUNSUB/ NPUNT(NMG),NGRUP(NMG),NCANT
       COMMON /CONT/ PMG(NMG),DELTC(NMG),DELV(NMG),DELPC(NMG),&
		    DELVI(NMG),DELN(NMG),DELTCA(NMG)
       COMMON /STAT/ MM(NMG),MJ(NMG),MK(NMG),MI(NMG),MH(NMG)
       COMMON /NOM/ FS
       COMMON /GEOM/ R(NMG),Q(NMG),MAIN(NMG)
        type(groups),pointer::recorre
        
           
       I=1
       recorre => LSubGroups
	 DO WHILE (NGRUP(I).NE.0)
          CALL Leer_Pr (NGRUP(I),IPAREQ,NAME,M1,J1,K1,I1,H1,PM,RPAR,&
                       QPAR,DTC,DPC,DV,dvi,dn,DTCA)
		  FS(I) = recorre%Name
          MM(I) = recorre%attM
          MJ(I) = recorre%attJ
          MK(I) = recorre%attK
          MI(I) = recorre%attI
          MH(I) = recorre%attH
          PMG(I) = recorre%MW
          DELTC(I) = recorre%dTC
          DELV(I) = recorre%dV
          DELPC(I) = recorre%dPC
          delvi(i)=recorre%dVi
     	  deln(i)=recorre%dN
          R(I) = recorre%R
          Q(I) = recorre%Q
          DELTCA(I)=recorre%dTCa
          I=I+1
          recorre => recorre%next
       ENDDO
       RETURN 
       END

SUBROUTINE Store_In ()
!C-----------------------------------------------------------------------
!C      ESTA SUBRUTINA CARGA TODOS LOS PARAMETROS DEL BANCO INTRCN QUE SE
!C      USAN EN MOLDES. LAS PROPIEDADES SE ALMACENAN EN COMMONS.
!C-----------------------------------------------------------------------
!C
    use Input
   !PARAMETER (NG=70,NMG=150)
    implicit none

!INTEGER
    integer::i,j,id_g1,id_g2

!COMMONS    
    COMMON/PUNGRU/NPINT,NINTT,NUMINT
    integer::NPINT(NINT),NINTT(NINT),NUMINT
    COMMON/INTER/A,kstr,kdot,alpha
    real*8::A(NMG,NMG),kstr(NMG,NMG),kdot(NMG,NMG),alpha(NMG,NMG)

!FUNCIONES
    real*8,external::Leer_Alpha, Leer_Kapa,Leer_In

    if(Model /= 3)then
    I=1
    do while (NINTT(I).ne.0)
       J=1
       do while (NINTT(J).ne.0)

          A(I,J) = Leer_In(NINTT(I),NINTT(J),IPAREQ)

          J=J+1
       enddo
       I=I+1
    enddo

    else

    I=1
    do while (NINTT(I).ne.0)
       J=I
       do while (NINTT(J).ne.0)

             if(Nintt(i) > Nintt(j))then
                 id_g1 = Nintt(j)
                 id_g2 = Nintt(i)
             else
                 id_g1 = Nintt(i)
                 id_g2 = Nintt(j)
             endif
             
             kstr(id_g1,id_g2) = Leer_Kapa(id_g1,id_g2) !; kstr(id_g2,id_g1) = param(1)
             if(isnan(kstr(id_g1,id_g2))) kstr(id_g1,id_g2) = 1
             kstr(id_g2,id_g1) = kstr(id_g1,id_g2)
             
             kdot(id_g1,id_g2) = Leer_Kapa(id_g2,id_g1) !; kdot(id_g2,id_g1) = param(2)
             if(isnan(kdot(id_g1,id_g2))) kdot(id_g1,id_g2) = 0
             kdot(id_g2,id_g1) = kdot(id_g1,id_g2) 
             
             alpha(id_g1,id_g2) = Leer_Alpha(id_g1,id_g2)
             if(isnan(alpha(id_g1,id_g2))) alpha(id_g1,id_g2) = 0
             alpha(id_g2,id_g1) = Leer_Alpha(id_g2,id_g1)
             if(isnan(alpha(id_g2,id_g1))) alpha(id_g2,id_g1) = 0

          J=J+1
       enddo
       I=I+1
    enddo
    
    endif    
    
    RETURN
    ENDSUBROUTINE Store_In

SUBROUTINE ESCPAR (NC,MODEL,IDEV)
    use constantes
      IMPLICIT real*8 (A-H,O-Z)
      LOGICAL NOACE
      COMMON/US/MS(NCOM,DiffStructGroups,2),nms(3)

      CALL PARIN (NC,NG,IDEV,NOACE,MS)
      RETURN
END

subroutine parin(NC,NG,IDEV,NOACE,MS)
!
!		Escribe par�metros de Q, R y a para subgrupos
!
    use CONSTANTES
    use Input
    use SubGrupos
    PARAMETER(NGPM=30,NGA=70,NSCM=10,NMAXCOMP=30)
    implicit real*8(A-H,O-Z)
    
    logical,intent(out)::noace
    
    type(Compound),pointer::recorreSolutes
  !  type(FinalStructure),pointer::recorreSolventes
    integer::MS(NCOM,DiffStructGroups,2),nms(NCOM)
    COMMON/PUNSUB/NPUNT(NMG),NGRUP(NMG),NCANT
    COMMON/PUNGRU/NPINT(NGA),NINTT(NGA),NUMINT
    COMMON/UNIF/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
    COMMON/INTER/A,kstr,kdot,alpha
    real*8::A(NMG,NMG),kstr(NMG,NMG),kdot(NMG,NMG),alpha(NMG,NMG)
    COMMON/GEOM/RR(NMG),QQ(NMG),MAIN(NMG)
    DIMENSION RT(NGPM,NCOM),NGM(NGPM)
    DIMENSION NY(NMAXCOMP,NGPM),JH(NMG),IH(NGPM)
    LOGICAL exists

      NOACE = .FALSE.
      NK=NC !cuando se llama de Solvent_Properties, es igual a 3

      QT(:,:)=0.D0 !QT (30,3)
      RT(:,:)=0.D0 !RT (30,3)
      P(:,:) = 0.D0
      JH(:)=0
      
!-----Ordena de < a > , en el vector IH, todos los grupos que
!     participan en la corrida (CAR, CPR y mol�cula a dise�ar o solvente(s) a evaluar)  

 !     nk=0 !Cuenta la cantidad de componentes
 !!     IC=1
 !     recorreSolutes => InputProblem01%MixtureInput%Solutes
 !     do while (associated(recorreSolutes)) !Agrega a ih los grupos de solutos
 !       i=1
 !       do while(recorresolutes%Formula(i,1)/=0)
 !           call insert_group(recorresolutes%Formula(i,1),ih,NGPM)
 !           i=i+1
 !       enddo
 ! !      ic=ic+i-1
 !       nk=nk+1
 !       recorreSolutes => recorresolutes%next
 !     enddo
      
      !i=1
      !do while(InputProblem01%MixtureInput%PCR%Formula(i,1)/=0) !Agrega a ih los grupos del PCR
      !  call insert_group(InputProblem01%MixtureInput%PCR%Formula(i,1),ih,NGPM)
      !  i=i+1
      !enddo
      !nk=nk+1 
      
      !i=1
      !do while(MS(3,i,1)/=0) !Agrega a ih los grupos del solvente
      !  call insert_group(MS(3,i,1),ih,NGPM)
      !  i=i+1
      !enddo      
      !nk=nk+1

      
      
      !AREGLAR ESTO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ngrup(:) = 0
    npunt(:) = 0
    npint(:) = 0
    nintt(:) = 0
    ih(:) = 0
!...Grupos presentes en la mezcla
    Do i=1,NK
        k=1
        Do while (ms(i,k,1)/=0)
            call CR_PUNTF (ms(i,k,1),NPUNT,NGRUP,NPINT,NINTT)
            k=k+1
        enddo
    enddo
    call LEER_PRnew (ipareq,ngrup,RR,QQ)
    call Leer_Int_Par (ipareq,nintt,A)
      do i=1,NK
          j=1
          do while(ms(i,j,1)/=0)
              call insert_group(ms(i,j,1),ih,NGPM)
              j=j+1
          enddo
          
      enddo
      
      i=1
      do while(ih(i)/=0) !cuenta antidad de grupos en ih
        i=i+1
      enddo
      ic=i-1
      
      call ordenar_arreglo (ih,NGPM,ic) !ordena el vector ih
     
      
      DO 73 I=1,IC ! Guarda en la posici�n de JH que corresponde al n� del grupo, la posici�n en IH
73    JH(IH(I))=I
      
!-----En guarda las cantidades de cada subgrupo seg�n el mismo orden de columna que
!     IH, pero en tres filas distintas (dependiendo de si es CAR, CPR o Sol.)
      
      NY(:,:)=0
      
      do i=1, NK
          j=1
          do while(ms(i,j,1)/=0)
            N1=ms(i,j,1)
            N2=ms(i,j,2) 
            N3=JH(N1)
            NY(i,N3)=N2                    
            j=j+1  
          enddo
      enddo
      
      !recorreSolutes => InputProblem01%MixtureInput%Solutes
      !i=0
      !do while (associated(recorreSolutes))!Carga los solutos
      !  j=1
      !  i=i+1        
      !  do while(recorresolutes%Formula(j,1)/=0)
      !      N1=recorresolutes%Formula(j,1)
      !      N2=recorresolutes%Formula(j,2) 
      !      N3=JH(N1)
      !      NY(i,N3)=N2                    
      !      j=j+1
      !  enddo
      !  recorreSolutes => recorresolutes%next
      !enddo   
      !j=1
      !i=i+1
      !do while(InputProblem01%MixtureInput%PCR%Formula(j,1)/=0) !Carga el PCR
      !    N1 = InputProblem01%MixtureInput%PCR%Formula(j,1)
      !    N2 = InputProblem01%MixtureInput%PCR%Formula(j,2) 
      !    N3=JH(N1)
      !    NY(i,N3)=N2                    
      !    j=j+1
      !enddo   
      !j=1
      !i=i+1
      !do while(MS(3,j,1)/=0) !Carga el solvente
      !    N1=MS(3,j,1)
      !    N2=MS(3,j,2) 
      !    N3=JH(N1)
      !    NY(i,N3)=N2                    
      !    j=j+1
      !enddo      

      I=0
      NGMGL=0
      NGM(:)=0
      DO 80 K=1,IC
  !       call buscaras(main(npunt(ih(k))),NGM,NGPM,exists)
         NSG=IH(K)
         !NGMNY=MAIN(NPUNT(NSG))
         NGMNY=mainsg (NSG,ipareq)
         IF(NGMNY.NE.NGMGL) I=I+1
   !     if(exists) cycle
    !     i=i+1
         NGM(I)=NGMNY
         NGMGL=NGMNY
         DO 80 J=1,NK
            RT(I,J)=RT(I,J)+NY(J,K)*Obtain_R(NSG)
            QT(I,J)=QT(I,J)+NY(J,K)*Obtain_Q(NSG)        
80    continue
      NG=I
      IF ((IDEV.EQ.6).OR.(IDEV.EQ.2)) THEN
         WRITE(IDEV,608) (IH(K),K=1,IC)				!(' SUB GROUPS :',20I3)
         WRITE(IDEV,609) (mainsg (IH(K),ipareq),K=1,IC)	!(' MAIN GROUPS:',20I3)
         WRITE(IDEV,610)								!(' COMPONENT')
         DO 191 I=1,NK								
191         WRITE(IDEV,"(6X,I2,5X,20I3)") I,(NY(I,K),K=1,IC)		!(6X,I2,5X,20I3)
      END IF
   85 CONTINUE
!-----Checkea si exiten par�metros de interacci�n
    if(model /= 3)then
      DO 20 I=1,NG !NG=I
         DO 20 J=1,NG
            NI=NGM(I)
            NJ=NGM(J)
            P(I,J)=A(NPINT(NI),NPINT(NJ))
            IF (DABS(P(I,J) - 9000.).LE.0.1) THEN
               NOACE = .TRUE.
               GO TO 1000
            END IF
   20 CONTINUE
      NN=IH(K)
      IF ((IDEV.EQ.6).OR.(IDEV.EQ.2)) THEN
         WRITE(IDEV,612)
         DO 960 K=1,IC
            NN=IH(K)
  960       WRITE(IDEV,613) NN,Obtain_R(NN),Obtain_Q(NN)
         WRITE(IDEV,604)
         DO 127 I=1,NG
  127       WRITE(IDEV,603) (P(I,J),J=1,NG)
      END IF
    endif
      DO 30 I=1,NK
         Q(I)=0.D0
         R(I)=0.D0
         DO 30 K=1,NG
            Q(I)=Q(I)+QT(K,I)
   30       R(I)=R(I)+RT(K,I)    
!c     IF(IOUT.NE.6) GOTO 260
!c     DO 40 I=1,NK
!c  40    WRITE(6,606) I,R(I),Q(I)
!c     DO 140 I=1,NK
!c 140    WRITE(2,606) I,R(I),Q(I)
  260 CONTINUE
  501 FORMAT(20I3)
  502 FORMAT(8F10.2)
  503 FORMAT(I3,2F10.2)
  504 FORMAT(2I3,F10.2)
  603 FORMAT(1X,10F12.3)
  604 FORMAT(' ',/,' INTERACTION PARAMETERS',/)
  605 FORMAT(' UNIFAC MOLECULAR R AND Q',/)
  606 FORMAT(I5,2F15.4)
  607 FORMAT(' ** WARNING: NUMBER OF SUB GROUPS MUST NOT EXCEED 20 **')
  608 FORMAT(' SUB GROUPS :',20I3)
  609 FORMAT(' MAIN GROUPS:',20I3)
  610 FORMAT(' COMPONENT')
!  611 FORMAT(6X,I2,5X,20I3)
  612 FORMAT(' ',/,' GROUP R- AND Q-VALUES',/)
  613 FORMAT(1X,I3,2F10.4)
  627 FORMAT(' SPECIFIED UNIQUAC R AND Q',/)
  699 FORMAT(/)
 1000 RETURN
endsubroutine parin


SUBROUTINE UNIPAR (NC,NG,T,MODEL,IOUT,NOACE,MS) !NG= n� main goup dist
    use constantes
      IMPLICIT real*8 (A-H,O-Z)
      integer::MS(NCOM,DiffStructGroups,2)
      LOGICAL NOACE

      CALL PARIN (NC,NG,IOUT,NOACE,MS)
      CALL PARAM (NC,NG,T)
      RETURN
endsubroutine

SUBROUTINE PARIN2 (NC,NG,MODEL,IPAREQ,IDEV,NOACE)
!
!		Escribe par�metros de Q,R y a para subgrupos
!
      PARAMETER(NGPM=30,NGA=70,NCOM=3,NSCM=10,NMG=150)
      implicit real*8(A-H,O-Z)
      COMMON/US/MS(NCOM,NSCM,2),nms(NCOM)
      COMMON/UNIF/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
      COMMON/PUNSUB/NPUNT(NMG),NGRUP(NMG),NCANT
      COMMON/PUNGRU/NPINT(NGA),NINTT(NGA),NUMINT
    COMMON/INTER/A,kstr,kdot,alpha
    real*8::A(NMG,NMG),kstr(NMG,NMG),kdot(NMG,NMG),alpha(NMG,NMG)
      COMMON/GEOM/RR(NMG),QQ(NMG),MAIN(NMG)
      DIMENSION RT(NGPM,NCOM),NGM(NGPM)
      DIMENSION NY(NCOM,NGPM),JH(NMG),IH(NGPM)
      LOGICAL NOACE
!

      IPAREQ = IPAREQ
      NOACE = .FALSE.
      NK=NC
      DO 15 I=1,NGPM
         DO 15 J=1,NK
            QT(I,J)=0.D0
   15       RT(I,J)=0.D0
      DO 16 I=1,NGPM
         DO 16 J=1,NGPM
            P(I,J) = 0.D0
   16 CONTINUE
      DO 49 I=1,NMG
   49    JH(I)=0
      IC=1
      DO 71 I=1,NK
         DO 70 J=1,NMS(I)
            IF(MS(I,J,1).EQ.0) GOTO 71
            IH(IC)=MS(I,J,1)
            IF(IC.EQ.1) GOTO 69
            IF(IH(IC).EQ.IH(IC-1)) GOTO 70
            IF(IH(IC).GT.IH(IC-1)) GOTO 69
            IF(IC.GT.2) GOTO 55
            IHH=IH(1)
            IH(1)=IH(2)
            IH(2)=IHH
            GO TO 69
   55       I1=IC-1
            DO 65 I2=1,I1
               IF(IH(IC).GT.IH(I2)) GOTO 65
               IF(IH(IC).EQ.IH(I2)) GOTO 70
               I4=IC-I2
               DO 61 I3=1,I4
   61             IH(IC+1-I3)=IH(IC-I3)
               IH(I2)=MS(I,J,1)
   65       CONTINUE
   69       IC=IC+1
   70    CONTINUE
   71 CONTINUE
      IC=IC-1
      DO 73 I=1,IC
   73    JH(IH(I))=I
      DO 72 I=1,NCOM
         DO 72 J=1,NGPM
   72       NY(I,J)=0
      DO 75 I=1,NK
         DO 74 J=1,NMS(I)
            IF(MS(I,J,1).EQ.0) GOTO 75
            N1=MS(I,J,1)
            N2=MS(I,J,2)
            IF(N1.EQ.0) GOTO 75
            N3=JH(N1)
   74       NY(I,N3)=N2
   75 CONTINUE
      I=0
      NGMGL=0
      DO 80 K=1,IC
         NSG=IH(K)
         NGMNY=MAIN(NPUNT(NSG))
         IF(NGMNY.NE.NGMGL) I=I+1
         NGM(I)=NGMNY
         NGMGL=NGMNY
         DO 80 J=1,NK
            RT(I,J)=RT(I,J)+NY(J,K)*RR(NPUNT(NSG))
   80       QT(I,J)=QT(I,J)+NY(J,K)*QQ(NPUNT(NSG))
      NG=I
      IF ((IDEV.EQ.6).OR.(IDEV.EQ.2)) THEN
         WRITE(IDEV,608) (IH(K),K=1,IC)
         WRITE(IDEV,609) (MAIN(NPUNT(IH(K))),K=1,IC)
         WRITE(IDEV,610)
         DO 191 I=1,NK
191         WRITE(IDEV,611) I,(NY(I,K),K=1,IC)
      END IF
   85 CONTINUE
      DO 20 I=1,NG
         DO 20 J=1,NG
            NI=NGM(I)
            NJ=NGM(J)
            P(I,J)=A(NPINT(NI),NPINT(NJ))
            IF (DABS(P(I,J) - 9000.).LE.0.1) THEN
               NOACE = .TRUE.
               GO TO 1000
            END IF
   20 CONTINUE
      NN=IH(K)
      IF ((IDEV.EQ.6).OR.(IDEV.EQ.2)) THEN
         WRITE(IDEV,612)
         DO 960 K=1,IC
            NN=IH(K)
  960       WRITE(IDEV,613) NN,RR(NPUNT(NN)),QQ(NPUNT(NN))
         WRITE(IDEV,604)
         DO 127 I=1,NG
  127       WRITE(IDEV,603) (P(I,J),J=1,NG)
      END IF
      DO 30 I=1,NK
         Q(I)=0.D0
         R(I)=0.D0
         DO 30 K=1,NG
            Q(I)=Q(I)+QT(K,I)
   30       R(I)=R(I)+RT(K,I)
!     IF(IOUT.NE.6) GOTO 260
!     DO 40 I=1,NK
!  40    WRITE(6,606) I,R(I),Q(I)
!     DO 140 I=1,NK
! 140    WRITE(2,606) I,R(I),Q(I)
  260 CONTINUE
  501 FORMAT(20I3)
  502 FORMAT(8F10.2)
  503 FORMAT(I3,2F10.2)
  504 FORMAT(2I3,F10.2)
  603 FORMAT(1X,10F12.3)
  604 FORMAT(' ',/,' INTERACTION PARAMETERS',/)
  605 FORMAT(' UNIFAC MOLECULAR R AND Q',/)
  606 FORMAT(I5,2F15.4)
  607 FORMAT(' ** WARNING: NUMBER OF SUB GROUPS MUST NOT EXCEED 20 **')
  608 FORMAT(' SUB GROUPS :',20I3)
  609 FORMAT(' MAIN GROUPS:',20I3)
  610 FORMAT(' COMPONENT')
  611 FORMAT(6X,I2,5X,20I3)
  612 FORMAT(' ',/,' GROUP R- AND Q-VALUES',/)
  613 FORMAT(1X,I3,2F10.4)
  627 FORMAT(' SPECIFIED UNIQUAC R AND Q',/)
  699 FORMAT(/)
 1000 RETURN
endsubroutine


SUBROUTINE PARAM(NC,NG,T)                                               
!
!	Calcula los coeficientes de actividad residuales grupales
!
    use input
      PARAMETER(NGPM=30,NGA=70,NSCM=10)
      implicit real*8(A-H,O-Z)
      COMMON/UNIF/QT(NGPM,NCOM),TAU(NGPM,NGPM),S(NGPM,NCOM),F(NCOM),Q(NCOM),R(NCOM),P(NGPM,NGPM)
      DO 30 I=1,NG                                                            
         DO 30 J=1,NG                                                          
   30       TAU(I,J)=DEXP(-P(I,J)/T) !P: par�metro de interacci�n grupal. Ec. 9 Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975)
                                                             
      DO 50 I=1,NC                                                            
         DO 50 K=1,NG                                                          
		  S(K,I)=0.D0                                                        
            DO 50 M=1,NG
   50          S(K,I)=S(K,I)+QT(M,I)*TAU(M,K) !QT: fracci�n de �rea del grupo k. S: Ec. 7 Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975)                        
    
      DO 60 I=1,NC                                                            
         F(I)=1.D0                                                             
         DO 60 J=1,NG 
   60       F(I)=F(I)+QT(J,I)*DLOG(S(J,I)) !Segundo t�rmino de Ec. 7 Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975)                             
      RETURN                                                                  
endsubroutine                                                                   

  



!
SUBROUTINE TABLA_PARAMETROS (ipareq)
!-----------------------------------------------------------
!
!-----------------------------------------------------------
    
    integer,intent(out)::ipareq
!
    if(ipareq/=3)then
        do while (ipareq<1 .or. ipareq>3)
            write (6,520)
            read (5,*) ipareq
        enddo
	else
	    ipareq = 4
    endif

	call SubGroups_Characterisation()
	
	
520  format (/,' Choose Unifac parameters: ',/,28x,&
                ' liquid-liquid:      1',&
             /,28x,' liquid-vapor:       2',&
         	  /,28x,' infinite dilution:  3',&
      	    //,52x,'> ',$)
!521  format (/,' Option not available in the actual version ')
endsubroutine tabla_parametros


       SUBROUTINE LEER_PRnew (IPAREQ,ngrup,RR,QQ)
!-----------------------------------------------------------------------
!      ESTA SUBRUTINA CARGA TODAS LOS PARAMETROS DEL BANCO GRUPOSRAM QUE
!      SE USAN EN MOLDES. LAS PROPIEDADES SE ALMACENAN EN COMMONS.
!-----------------------------------------------------------------------
!
       PARAMETER (NMG=150)
       IMPLICIT real*8 (A-H,O-Z)
       integer,intent(in)::ipareq,ngrup(NMG)
       real*8, intent(out)::RR(NMG), QQ(NMG)
       CHARACTER*8 FS(NMG),NAME
       INTEGER H1
       
     
       I=1
	    DO WHILE (NGRUP(I).NE.0)
          CALL Leer_Pr (NGRUP(I),IPAREQ,NAME,M1,J1,K1,I1,H1,PM,RPAR,&
                       QPAR,DTC,DPC,DV,dvi,dn,DTCA)
          RR(I) = RPAR
          QQ(I) = QPAR
          I=I+1
       ENDDO
       RETURN 
       END
 
        SUBROUTINE Leer_Int_Par (IPAREQ,nintt,A)
!-----------------------------------------------------------------------
!      ESTA SUBRUTINA CARGA TODOS LOS PARAMETROS DEL BANCO INTRCN QUE SE
!      USAN EN MOLDES. LAS PROPIEDADES SE ALMACENAN EN COMMONS.
!-----------------------------------------------------------------------
!
       PARAMETER (NG=70,NMG=150)
       IMPLICIT real*8 (A-H,O-Z)
       EXTERNAL Leer_In
       real*8 Leer_In,a(nmg,nmg)
       integer nintt(ng)

!
       I=1
       do while (NINTT(I).ne.0)
          J=1
          do while (NINTT(J).ne.0)
             A(I,J) = Leer_In(NINTT(I),NINTT(J),IPAREQ)
             J=J+1
          enddo
          I=I+1
       enddo
       RETURN
       END
       
       

                                                    


      subroutine inter1 (imprim,kgrup,mgr,mint,iint)
!c-----------------------------------------------------------------------------
!c     Esta subrutina  controla  que  cada  grupo  almacenado  en  mgr, tenga 
!c     parametros de interaccion con todos los subgrupos almacenados en kgrup.
!c     Los subgrupos de mgr que satisfacen esa condicion se incorporan a mint.
!c
!c     Variables de entrada:
!c     --------------------
!c                   ipareq: modelo de parametros de equilibrio:
!c                           1:liquido-liquido
!c                           2:liquido-vapor
!c                           3:liquido-vapor a dilucion infinita.
!c                   imprim: 0: sin salida.
!c                           1: salida por pantalla.
!c                           2: salida por pantalla y por archivo (unidad 2).
!c                    kgrup: vector de  numeros  de  identificacion  de  grupos
!c                           principales base.
!c                      mgr: vector de  numeros  de  identificacion  de  grupos
!c                           principales a ser analizado.  
!c
!c     Variables de salida:
!c     -------------------
!c                     mint: vector de  numeros  de  identificacion  de  grupos
!c                           principales que satisfacen la condicion  explicada
!c                           en el encabezamiento.  Su dimension es nsel.
!c                     iint: numero de elementos en el vector mint.
!c
!c     Las dimensiones de  los  vectores  quedan  fijadas  por  los  parametros
!c     siguientes:
!c
      use CONSTANTES
      use Input
!c
!c-----------------------------------------------------------------------------
      implicit none
      external nom_gru1,Leer_In
      integer:: kgrup(nsel),mgr(nsel),mint(nsel)
      character*8 fs1,fs2,nom_gru1
      logical inter
      
      integer::iinit,i,j,imprim,iint,k1,k2
      real*8::dabs,apar,alphapar,kstrpar,kdotpar 
!c     common ipant     
    !c                  
    real*8::Leer_In,Leer_Alpha,Leer_Kapa

    iint = 0
    inter = .true.
    i=1
	do while(mgr(i)/=0)!20 i=1,igr
        j=1 
        do while(kgrup(j)/=0)!30 j=1,ifin
            k1 = kgrup(j)
            k2 = mgr(i)
            fs1 = nom_gru1 (k1,ipareq)
            fs2 = nom_gru1 (k2,ipareq)
            if(model /= 3)then
                apar = Leer_In (k1,k2,ipareq)  
            else
                alphapar = Leer_Alpha(k1,k2)
                kstrpar = Leer_Kapa(k1,k2)
                kdotpar = Leer_Kapa(k2,k1)
            endif
            if(model /= 3)then    
                if (dabs(apar-9000.).le.0.1) then
                    inter = .false.
                    if (imprim.gt.0) then
                        write (6,"(' ',/,' ** Interaction parameters not available for: ',2(a8,2x))")  fs1,fs2
                        if (imprim.eq.2) write (2,"(' ',/,' ** Interaction parameters not available for: ',2(a8,2x))")  fs1,fs2
                    end if
                end if
            else
                !Completar para la GC
            endif
            j=j+1
        enddo !30      continue
        if (inter) then
            iint = iint + 1
            mint(iint) = mgr(i)
        end if
        inter = .true.
        i=i+1
    enddo !20   continue

    return
    endsubroutine inter1


      subroutine inter2 (ipareq,imprim,kgrup,ifin,mgr,igr,inter)
!c-----------------------------------------------------------------------------
!c     Esta subrutina controla si existe por lo menos un grupo del vector kgrup
!c     que tenga parametro de interaccion con algun grupo del vector mgr.
!c
!c     Variables de entrada:
!c     --------------------
!c                   ipareq: modelo de parametros de equilibrio:
!c                           1:liquido-liquido
!c                           2:liquido-vapor
!c                           3:liquido-vapor a dilucion infinita.
!c                   imprim: 0: sin salida.
!c                           1: salida por pantalla.
!c                           2: salida por pantalla y por archivo (unidad 2).
!c                    kgrup: vector de  numeros  de  identificacion  de  grupos
!c                           principales base.  Su dimension es nsel.
!c                     ifin: numero de elementos en el vector kgrup.
!c                      mgr: vector de  numeros  de  identificacion  de  grupos
!c                           principales a ser analizado.  Su dimension es nsel.
!c                      igr: numero de elementos en el vector mgr.
!c
!c     Variables de salida:
!c     -------------------
!c                    inter: true: se satisface la condicion  explicada  en  el
!c                                 encabezamiento.
!c                          false: no se satisface la condicion explicada en el
!c                                 encabezamiento.
!c
      use CONSTANTES
!c
!c-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8,external:: Leer_In
      character*8,external::nom_gru1
      dimension kgrup(nsel),mgr(nsel)
      character*8 fs1,fs2
      logical inter
!c     common ipant
!c
      inter = .false.
      i=1
      do while (mgr(i)/=0) !30 i=1,igr
         if (inter) exit
         j=1
         do while(kgrup(j)/=0)!20 j=1,ifin
            k1 = kgrup(j)
            k2 = mgr(i)
            fs1 = nom_gru1 (k1,ipareq)
            fs2 = nom_gru1 (k2,ipareq)
            apar = Leer_In (k1,k2,ipareq)
            if (dabs(apar-9000.).gt.0.1) then
               inter = .true.
            else
               if (imprim.gt.0) then
                  write (6,10)  fs1,fs2
                  if (imprim.eq.2) then
                     write (2,10)  fs1,fs2
                  end if
               end if
            end if
            j=j+1
         enddo !20      continue
         i=i+1
      enddo !30   continue
!40   continue
!c
!c---- formatos
 10   format (' ',/,' ** Interaction parameters not available for: ',2(a8,2x))
      return
      end


      subroutine cambiar_status (grupos,ngrup,ngm1v,ngm,mst,JST1,KST1,IST1,hst)
!c-----------------------------------------------------------------------------
!c     Esta subrutina permite cambiar las propiedades de combinacion para los
!c     subgrupos que participan en el diseno molecular de solventes.
!c
!c     Variables de entrada:
!c     --------------------
!c                       fs: vector  de  nombres  de  identificacion  de  los
!c                           subgrupos UNIFAC. Su dimension es NMG.
!c                   grupos: vector  de  numeros  de  identificacion  de  los
!c                           subgrupos  cuyas   propidedades  de  combinacion 
!c                           pueden ser modificadas. Su dimension es nsel.
!c                    ngrup: cantidad de elementos en el vector grupos.
!c
!c     Variables de salida:
!c     -------------------
!c                     ngm1V: vector  de  numeros  de  identificacion  de   los
!c                           subgrupos cuyas propiedades  de  combinacion  son
!c                           modificadas. Su dimension es nsel.
!c                      ngm: numero de elementos en el vector ngm1V.
!c                      mst: vector con las propiedades de combinacion  m para
!c                           los subgrupos de ngm1V. Su dimension es nsel.
!c                      JST1: idem j.
!c                      KST1: idem k.
!c                      IST1: idem i.
!c                      hst: idem h.
!c
!c     Observacion:  Esta   subrutina  esta   preparada   para   modificar  20
!c                   componentes.
!c
!c     Las dimensiones de los  vectores  quedan  fijadas  por  los  siguientes
!c     parametros:
!c  
      use SubGrupos
      use CONSTANTES
!c
!c-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer grupos(NMG),ngm1V(nsel),mst(nsel),JST1(nsel),KST1(nsel),IST1(nsel),hst(nsel)
      character fs(NMG)*8,letra*1,return*1
!c      character coma*1
      logical enc
!c
      idev = 6
150   do 160 i=1,nsel
         mst(i) = 0
         JST1(i) = 0
         KST1(i) = 0
         IST1(i) = 0
         hst(i) = 0
 160  continue
!c     call limp (idev)
  60  write (6,10)
      do 170 i=1,nsel
         ngm1V(i) = 0
 170  continue
3010  call Groups_Present ()
      write (6,20)
      read (5,30,err=60) (ngm1V(i),i=1,20)
 30   format (20i3)
      ngm = 0
      do 230 i=1,nsel
         if (ngm1V(i).le.0) go to 240
         enc = .false.
         do 210 j=1,ngrup 
            if (ngm1V(i).eq.grupos(j)) then
               enc = .true.
               go to 220
            end if
 210     continue
 220     if (.not.enc) then
            write (6,40) ngm1V(i)
            read (5,50) return
           go to 60
         else
            ngm = ngm + 1
         end if
 230  continue
!c     call limp (idev)
 240  do 180 i=1,ngm
         write (6,70)
         write (6,80)
 110     write (6,190) fs(ngm1V(i))
         read (5,90) letra
3030     write (6,250) letra
         read (5,260,err=3030) num
         if (num.le.0) go to 180 
         if ((letra.eq.'M').or.(letra.eq.'m')) then
             mst(i) = num
         else if ((letra.eq.'J').or.(letra.eq.'j')) then
             JST1(i) = num
         else if ((letra.eq.'K').or.(letra.eq.'k')) then
             KST1(i) = num
         else if ((letra.eq.'I').or.(letra.eq.'i')) then
             IST1(i) = num
         else if ((letra.eq.'H').or.(letra.eq.'h')) then
             hst(i) = num
         else
             write (6,100) letra
             read (5,50) return
             go to 110
         end if
!c     call limp (idev)
 180  continue
      write (6,120)
      do 200 i=1,ngm
         write (6,130)fs(ngm1V(i)),mst(i),JST1(i),KST1(i),IST1(i),hst(i)
 200  continue
      write (6,140)
      read (5,50) return
      if (return.eq.'1') then
         go to 150
      end if
140   format (' ',/,' If it is OK. <ret>.  If not 1.   ',$)
120   format (' ',//,' Combination Properties Specified',//,' Component    M    J    K    I    H')
130   format (/,2x,a8,5i5)
100   format (' ',/,' * Incorrect Attachment Name: ',a1,'   <ret> ',$)
 90   format (a1)
260   format (i3)
190   format (1x,/,12x,'Subgroup:            ',a8,/,12x,'-status letter:           ',$)
250   format (12x,'-number of ',a1,' attachments: ',$)
 70   format (' ',/,' Combination Properties Specification',&
                 //,'     M,J,K,I,H: Number of attachments of type ',&
             'M,J,K,I,H,',/,'                present in a group',&
                /,'     M: Unrestricted Carbon Attachment',&
             ' (i.e. "CH3"),',&
                /,'     J: Unrestricted Radial Carbon Attachment',&
             ' (i.e. "CH2"),',&
                /,'     K: Restricted Carbon Attachment,',&
                /,'     I: Aromatic Carbon Attachment in a not',&
             ' substituded',/,&
                   '        Aromatic Group (i.e. "ACH"),',/,&
                   '     H: Aromatic Carbon with the Hidrogen',&
             ' substituded for',/,&
                   '        other group (i.e. "ACH2","ACNO2").')
 80   format (' ',/,' Give the Combination Properties for each',&
             ' subgroup.',/,' Example:   Subgroup:               (OH)',&
                          /,'            -status letter:           M',&
                          /,'            -number of M attachments: 1')
 10   format (' ',/,' Choose the groups that you want to specifie the',' Combination Properties.')
 20   format (' ',/,' Write the chosen group numbers separated with ','commas (maximum 20 groups).')
 40   format (' ',/,' * Incorrect chosen group: ',i3,'   <ret> ',$)
 50   format (a)
      return
      end


      subroutine no_isomeros (igrupi,ifig,mgr,igr,maingr,jgrup,numj)
!c-----------------------------------------------------------------------------
!c     Esta subrutina selecciona los subgrupos  de  valencia  simple  que  no 
!c     formaran estructuras isomeras.
!c
!c     Variables de entrada:
!c     --------------------
!c                   igrupi: vector de numeros de identificacion de subgrupos
!c                           de valencia simple  posibles  para  el  problema
!c                           corriente.  Su dimension es nsel.
!c                     ifig: numero de elementos en el vector igrupi.
!c                      mgr: vector de numeros de identificacion de subgrupos
!c                           de valencia dual seleccionados. Su dimension  es
!c                           nsel.
!c                      igr: numero de elmentos en el vector mgr.
!c                   maingr: vector de numeros de identificacion de los grupos
!c                           principales   para   cada  subgrupo  UNIFAC.   Su 
!c                           dimension es NMG.
!c
!c     Variables de salida:
!c     -------------------
!c                    jgrup: vector de numeros de identificacion de subgrupos
!c                           que   no   formaran  estructuras  isomeras.   Su 
!c                           dimension es nsel.
!c                     numj: numero de elementos en el vector jgrup.
!c
!c     Las dimensiones de los parametros  estan  fijadas  por  los  parametros
!c     siguientes
!c
      use CONSTANTES
!c
!c-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension igrupi(nsel),mgr(nsel),jgrup(nsel),maingr(NMG)
      logical esisom

      numj = 0
      do 10 i=1,ifig
         esisom = .false.
         do 260 j=1,igr
            if (igrupi(i).eq.maingr(mgr(j))) then
               if (maingr(mgr(j)).ne.1) then
                  esisom = .true.
                  go to 270
               end if
            end if
 260     continue
 270     if (.not.(esisom)) then
            numj = numj + 1
            jgrup(numj) = igrupi(i)
         end if
 10   continue
      return
      end
      
      
       subroutine car_car (ipareq,isvt,isdt,idvt,i3vt,i4vt,icyct,iart,igrpt,&
                                    nsvt,nsdt,ndvt,n3vt,n4vt,ncyct,nart,ngrpt)
!c-------------------------------------------------------------------------------
!c
!c      Esta subrutina carga los numeros de los subgrupos segun su valencia
!c      para  el  dise�o  molecular  de solventes y los nombres de todos los
!c      subgrupos del banco de datos correspondiente.
!c
!c      Datos de entrada:
!c      ----------------
!c                       ipareq: banco de datos de equilibrio:
!c                             1: liquido-liquido
!c                             2: liquido-vapor
!c                             3: liquido-vapor a dilucion infinita.
!c
!c      Variables de salida:
!c      -------------------
!c                            fs: vector  que  contiene  los  nombres  de 
!c                                todos los subgrupos del banco de  datos
!c                                de equilibrio ipareq.
!c                           isv: subgrupos terminales de valencia simple.
!c                           isd: subgrupos terminales de valencia simple
!c                                sin subgrupos de valencia  dual  en  el
!c                                mismo grupo principal.
!c                           idv: subgrupos   de   valencia   dual   para
!c                                compuestos alifaticos.
!c                           i3v: subgrupos   tri   valentes   para
!c                                compuestos alifaticos.
!c                           i4v: subgrupos  tetra  valentes  para
!c                                compuestos alifaticos.
!c                          icyc: subgrupos   de   valencia   dual   para
!c                                compuestos ciclicos.
!c                           iar: subgrupos para compuestos aromaticos.
!c                          igrp: subgrupos que son compuestos moleculares.
!c                           nsv: numero de subgrupos en el vector isv.
!c                           nsd: numero de subgrupos en el vector isd.
!c                           ndv: numero de subgrupos en el vector idv.
!c                           n3v: numero de subgrupos en el vector i3v.
!c                           n4v: numero de subgrupos en el vector i4v.
!c                          ncyc: numero de subgrupos en el vector icyc.
!c                           nar: numero de subgrupos en el vector iar.
!c                          ngrp: numero de subgrupos en el vector igrp.
!c
!c      Funciones y Subrutinas Intervinientes:
!c      -------------------------------------
!c      Del paquete MANUNI:
!c                              CARAC
!c                              MAX_SUB
!c
!c-----------------------------------------------------------------------
!c
!c
      use CONSTANTES
      use SubGrupos
      implicit real*8 (a-h,o-z)
!     Variables de entrada
      INTEGER,INTENT(IN)::IPAREQ
!     Variables de salida
      type(groups),pointer::recorre
      !CHARACTER*8,INTENT(OUT),OPTIONAL::FS(NMG)
      INTEGER,INTENT(inout)::isvt(NMG),isdt(NMG),&
     idvt(NMG),icyct(NMG),iart(NMG),igrpt(NMG),i3vt(NMG),i4vt(NMG),nsvt,&
     nsdt,ndvt,n3vt,n4vt,ncyct,nart,ngrpt
!     Otras variables
      CHARACTER*8 FST(NMG)
      character*3 car

      numsgr = Size_LSubGroups()
      nsvT = 0
      nsdT = 0
      ndvT = 0
      n3vT = 0
      n4vT = 0
      ncycT = 0
      narT = 0
      ngrpT = 0
    recorre=>LSubGroups
    i=1
	do while(associated(recorre))
		!call carac (i,ipareq,FST(i),recorre%Charact) !CARACTERIZACIONES SUBGRUPALES
		!if(.not.associated(recorre%next)) call pausa
		if ((recorre%Charact.eq.'SV1').or.(recorre%Charact.eq.'SV2')) then
			nsvT = nsvT + 1
			isvT(nsvT) = i
			if (recorre%Charact.eq.'SV2') then
				nsdT = nsdT + 1
				isdT(nsdT) = i
			end if
		else if ((recorre%Charact.eq.'DV1').or.(recorre%Charact.eq.'DV2')) then
			ndvT = ndvT + 1
			idvT(ndvT) = i
			if (recorre%Charact.eq.'DV2') then
				ncycT = ncycT + 1
				icycT(ncycT) = i
			end if
		else if (recorre%Charact.eq.'3V1') then
			n3vT = n3vT + 1
			i3vT(n3vT) = i
		else if (recorre%Charact.eq.'4V1') then
			n4vT = n4vT + 1
			i4vT(n4vT) = i
		else if (recorre%Charact.eq.'AR1') then
			narT = narT + 1
			iarT(narT) = i
		else if (recorre%Charact.eq.'GR1') then
			ngrpT = ngrpT + 1
			igrpT(ngrpT) = i
		end if
		recorre=>recorre%next
		i=i+1
      enddo		
      return
      end

      subroutine sel_gru_fam (ipareq,family,ifam,iar,imaar,&
                             imprim,kgrup,ifin2,mfvi,ndv,entro,mgr,igr,&
                             ngru1,mgru1,prob)
!c-----------------------------------------------------------------------------
!c     Esta subrutina lee los subgrupos de la familia de componentes a generar
!c     por pantalla, seleccionando aquellos para los cuales existe informacion
!c     de los parametros de interaccion.  Si la familia  es  de  subgrupos  de 
!c     valencia simple  se realiza  una  preseleccion  para  eliminar aquellos 
!c     grupos que forman estructuras isomeras con los subgrupos ingresados  en
!c     mfvi.
!c
!c     Variables de entrada:
!c     --------------------
!c                       fs: vector de nombres de identificacion de  subgrupos
!c                           UNIFAC.  Su dimension es NMG.
!c                   maingr: vector de numeros de grupos principales para cada
!c                           subgrupo UNIFAC.  Su dimension es NMG.
!c                   ipareq: modelo de parametros de equilibrio:
!c                           1: liquido-liquido
!c                           2: liquido-vapor
!c                           3: liquido-vapor a dilucion infinita.
!c                   family: cadena de 21 caracteres  con  el  nombre  de  la 
!c                           familia seleccionada.  
!c                     ifam: tipo de familia de solventes seleccionada.
!c                      iar: vector de numeros de identificacion de subgrupos
!c                           de la familia seleccionada.
!c                    imaar: numero de elementos en el vector iar.
!c                   imprim: 0: sin salida.
!c                           1: salida por pantalla.
!c                           2: salida por pantalla y archivo (unidad 2).
!c                    kgrup: vector de numeros  de  identificacion  de  grupos 
!c                           pertenecientes al CAR y CPR para seleccionar  los  
!c                           subgrupos de iar que tengan parametros de 
!c                           interaccion.
!c                    ifin2: cantidad de elementos en el vector kgrup.
!c                      mfvi: vector de numeros de identificacion de  subgrupos
!c                           de valencia dual para  determinar  los  subgrupos
!c                           de valencia  simple  que  no  forman  estructuras 
!c                           isomeras.  Este vector es necesario  solo  en  el
!c                           caso de que family sea igual a    'single valence
!c                           groups'.
!c                      ndv: cantidad de elementos en el vector mfvi.
!c                    entro: false: es la primera vez que  entra  en  el  caso
!c                                  corriente.
!c                           true: entra nuevamente con un caso  anterior  con 
!c                                 el mismo kgrup y familia.
!c                      
!c     Variables de salida:
!c     -------------------
!c                      mgr: vector de numeros de identificacion de subgrupos
!c                           de la familia iar seleccionados.
!c                      igr: cantidad de elementos en el vector mgr.
!c                    ngru1: vector de numeros de identificacion de subgrupos
!c                           seleccionables de la familia iar,  obtenido  por
!c                           esta subrutina la primera vez que entra.
!c                    mgru1: cantidad de elementos en el vector ngru1.
!c                     prob: false: existen subgrupos seleccionables.
!c                           true: no existen subgrupos seleccionables.
!c   
!c     Las dimensiones de los vectores quedan fijadas por los siguientes 
!c     parametros:
!c
      use CONSTANTES
      use SubGrupos
!c
!c-----------------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      integer grupos(nsel),lgrup(nsel),igrupi(nsel),maingr(NMG),&
             iar(nsel),kgrup(nsel),mgr(nsel),mfvi(nsel),jgrup(nsel),&
             ngru1(nsel)
      integer::igr
      character*8 fs(NMG)
      character*21 family
      logical prob,entro
!c
      idev = 6
      idevr= 5
      prob = .false.
      if (.not.(entro)) then !viene de afuera como falso
!c         entro = .true.
		call asignar_grupos_principales (iar,imaar,grupos)
		call armar_grupos_distintos (grupos,imaar,lgrup,0,ifinl)
		call inter1 (imprim,kgrup,lgrup,igrupi,ifig)
		if (imprim.gt.0) then
			CALL PAUSA 
		end if
		if (ifig.eq.0) then
!c			call limp (idev)         
			write (idev,1030)
			write (idev,1040)
			read (idevr,510) return
			prob = .true.
		else
!c			call limp (idev)         
			if (family.eq.'single valence groups') then
!c				Se anul� la restricci�n de grupos terminales 	
!c			    if ((ifam.ge.3).and.(ifam.le.5)) then
!c			    call no_isomeros (igrupi,ifig,mfvi,ndv,maingr,jgrup,
!c     *		                      numj)
!c			    else
				do 310 j=1,ifig
					jgrup(j) = igrupi(j)
 310				continue
				numj = ifig
!c			    end if
				write (idev,610)
				write (idev,1270) family
				call Sel_Groups (iar,imaar,jgrup,numj,mgr,igr)
				do 320 i=1,numj
					ngru1(i) = jgrup(i)
 320				continue
				mgru1 = numj
			else
				write (idev,600) family
				write (idev,1270) family
				call Sel_Groups (iar,imaar,igrupi,ifig,mgr,igr)
				do 330 i=1,ifig
					ngru1(i) = igrupi(i)
 330				continue
				mgru1 = ifig
			end if
		end if
      else
!c		call limp (idev)         
		if (family.eq.'single valence groups') then
			write (idev,610)
		else
			write (idev,600) family
		end if
		write (idev,1270) family
		call Sel_Groups (iar,imaar,ngru1,mgru1,mgr,igr)
      end if
!c
!c---- formatos
 510  format (a)
 980  format (i2)
 600  format (/,' Choose the groups that participate in the',&
             ' design  of',/,' structures (maximum 10',a21,').',&
             //,10x,'If you choose the group:',&
             ' 1',/,10x,'If not                 : <ret>.')
 610  format (/,' Choose the single valence groups (maximum 10 ','single valence groups).')
 680  format (20i3)
1030  format (' ',/,' ','* No interaction parameters for this family ','of groups',/,'   and the two main components.')
1040  format ('   Change the family of groups:  <ret>. ',$)
1270  format (' ',/,' ',29x,'* ',a21,' *',/)
      return
endsubroutine sel_gru_fam

SUBROUTINE LEEPAR (J,IREC1,IPAREQ,NGRUPA,ENASST,RKASST)
!c-----------------------------------------------------------------
!c     Lee los par�metros de las bases de datos PARVOLAS.MDS 
!c     (volumen de asociaci�n, UNIT=15) y PARENEAS.MDS (energ�a de 
!c     asociaci�n, UNIT=16)  
!c-----------------------------------------------------------------
IMPLICIT real*8 (A-H,O-Z)
PARAMETER (NMG=150)
INTEGER J,IREC1,IPAREQ
DIMENSION NGRUPA(NMG)
      READ(16,502,REC=IREC1)FS,(ENASST,I=1,MAINSG(NGRUPA(J),IPAREQ))
      READ(15,502,REC=IREC1)FS,(RKASST,I=1,MAINSG(NGRUPA(J),IPAREQ))     
 502  FORMAT(a8,70d12.5)      
END


real*8 function Leer_In(i1,i2,ipareq)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el parametro de interaccion entre los 
!c      grupos i1 y i2 de la tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      character*8 fs

      if (i1.eq.i2) then
          aint1 = 0.0
      else
	    irec4=i1+(ipareq-1)*70
	    read (13,20,rec=irec4) fs,(aint1,i=1,i2)
      end if
      
      Leer_In = aint1

  20  format (a8,70d12.5)
      return
end function Leer_In

real*8 function Leer_Alpha(i1,i2)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el parametro de interaccion entre los 
!c      grupos i1 y i2 de la tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      character*8 fs
      real*8::aint1

      if (i1.eq.i2) then
          aint1 = 0.0
      else
       ! do i4=1,60
       ! do i3=1,60
        read (13,20,rec=i1) fs,(aint1,i=1,i2)
       ! enddo
       ! enddo
      end if
      
      Leer_Alpha = aint1

  20  format (a8,60d12.5)
!20  format (a8,708x,d12.5)
      return
end function Leer_Alpha

real*8 function Leer_Kapa(i1,i2)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el parametro de interaccion entre los 
!c      grupos i1 y i2 de la tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      character*8 fs

      if (i1.eq.i2) then
          aint1 = 0.0
      else
	    irec4=i1
	    read (16,20,rec=irec4) fs,(aint1,i=1,i2)
      end if
      
      Leer_Kapa = aint1

  20  format (a8,60d12.5)
      return
end function Leer_Kapa


      integer function mainsg (i1,ipareq)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el numero de grupo principal correspon-
!c      diente al subgrupo i1 de la tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------
!c
	irec1=i1+(ipareq-1)*150
      read (14,10,rec=irec1) main
      mainsg = main
  10  format (i4)
      return
endfunction mainsg

character*8 function nom_gru1 (i1,ipareq)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el nombre o formula quimica del grupo
!c      Unifac i1 y tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      character fs*8

	irec4=i1+(ipareq-1)*70
      read (13,20,rec=irec4) fs
      nom_gru1 = fs

  20	format (a8)
      return
end

integer function max_sub (ipareq)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el numero maximo de subgrupos en la ta-
!c      bla Unifac de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
	character fs*8

	num=(ipareq-1)*150
	fs= ''
      do while (fs .ne. 'fin')
	  num=num+1
	  read (14,10,rec=num) fs
	enddo

  	max_sub = num-(ipareq-1)*150-1

  10  format (4x,a8)
      return
end

integer function max_sub_int (ipareq)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el numero maximo de subgrupos en la ta-
!c      bla Unifac de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
	character mg*8

	num=(ipareq-1)*70
	mg= ''
      do while (mg .ne. 'fin')
	  num=num+1
	  read (13,10,rec=num) mg
	enddo

  	max_sub_int = num-(ipareq-1)*70-1

  10  format (a8)
      return
end

subroutine Leer_Pr (i1,ipareq,fs,mm,mj,mk,mi,mh,pmg,rpar,qpar,&
     			   deltc,delpc,delv,delvi,deln,DTCA)
!c-------------------------------------------------------------------
!c      Esta subrutina devuelve las contribuciones grupales para la 
!c      prediccion de propiedades fisicas del grupo Unifac i1 y tabla 
!c      de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c
!c      Las contribuciones grupales son las siguientes:
!c                   pmg: peso molecular
!c                 deltc: temperatura
!c                 delpc: presion
!c                  delv: densidad
!c-------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      character fs*8, citem*3

	irec2=i1+(ipareq-1)*150

	read(14,30,rec=irec2)fs,item,item,citem,mm,&
                           mj,mk,mi,mh,rpar,qpar,&
     			    pmg,deltc,delpc,delv,&
     			    delvi,deln,DTCA

   30	format (12x,a8,2i4,a3,5i2,9d15.8)
      return
end

subroutine SubGroups_Characterisation()  
    use Input     
    use SubGrupos
    implicit none

!Variables de ENTRADA/SALIDA    
!Variables internas
    type(groups),pointer::nuevo
    integer::i,j,irec2
    integer::comod
!Sentencias
    i=1
    do while (.True.)
        call Create_Group(nuevo)
        nuevo%Number=i
        irec2=i+(ipareq-1)*150

        if(model /= 3) then

   30	    format (i4,2a8,8x,a3,5i2,9d15.8,4i2)
	        read(14,30,rec=irec2)nuevo%NumberMainrGroup,&
	                         nuevo%NameMainGroup,&
	                         nuevo%Name,&
	                         nuevo%Charact,&
	                         nuevo%attM,&
	                         nuevo%attJ,&
	                         nuevo%attK,&
	                         nuevo%attI,&
	                         nuevo%attH,&
	                         nuevo%R,&
	                         nuevo%Q,&
	                         nuevo%MW,&
	                         nuevo%dTC,&
	                         nuevo%dPC,&
	                         nuevo%dV,&
	                         nuevo%dVi,&
	                         nuevo%dN,&
	                         nuevo%dTCa,&
	                         (nuevo%AssoSites(j),j=1,4)
	     else

   31	    format (i4,2a8,8x,a3,5i2,9d15.8,5i2,5d15.8)
	        read(14,31,rec=i)nuevo%NumberMainrGroup,&
	                         nuevo%NameMainGroup,&
	                         nuevo%Name,&
	                         nuevo%Charact,&
	                         nuevo%attM,&
	                         nuevo%attJ,&
	                         nuevo%attK,&
	                         nuevo%attI,&
	                         nuevo%attH,&
	                         nuevo%R,&
	                         nuevo%Q,&
	                         nuevo%MW,&
	                         nuevo%dTC,&
	                         nuevo%dPC,&
	                         nuevo%dV,&
	                         nuevo%dVi,&
	                         nuevo%dN,&
	                         nuevo%dTCa,&
	                         (nuevo%AssoSites(j),j=1,4),&
	                         comod,&
	                         nuevo%Tstr,&
	                         nuevo%Q,&
	                         nuevo%gstr,&
	                         nuevo%gdot,&
	                         nuevo%gddot    
	     endif
	     if(nuevo%NameMainGroup/="fin     ") then
	        call Incorporate_Group (LSubGroups,nuevo)
	     else
	        exit
	     endif
	     i=i+1
    enddo
endsubroutine SubGroups_Characterisation

subroutine carac (i1,ipareq,fs,caract)
!c-------------------------------------------------------------------
!c      Esta  subrutina  devuelve  la  caracterizacion  del subgrupo 
!c      Unifac i1 y tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c
!c      Las caracterizaciones subgrupales son las siguientes:
!c                   sv1: subgrupos de  valencia  simple  excluyendo
!c                        los que en su grupo  principal  no  tienen
!c                        subgrupos de valencia dual. 
!c                   sv2: subgrupos de valencia  simple  que  en  su
!c                        grupo  principal  no  tienen  subgrupos de
!c                        valencia dual.  
!c                   dv1: subgrupos de valencia dual excluyendo  los
!c                        que pueden formar compuestos ciclicos.
!c                   dv2: subgrupos  de  valencia  dual  que  pueden
!c                        formar compuestos ciclicos.
!c                   ar1: subgrupos aromaticos.
!c                   gr1: subgrupos   que    conforman    compuestos
!c                        moleculares. 
!c				     3v1: Trivalentes
!c				     4v1: Tetravalentes
!c-------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      character fs*8, caract*3

	irec2=i1+(ipareq-1)*150

	read(14,30,rec=irec2)fs,item,item,caract,mm,&
                           mj,mk,mi,mh,rpar,qpar,&
     			    pmg,deltc,delpc,delv,&
     			    delvi,deln
   30	format (12x,a8,2i4,a3,5i2,8d15.8)
      return
end

subroutine car_combprop (i1,ipareq,mm,mj,mk,mi,mh)
!c-------------------------------------------------------------------
!c      Esta subrutina devuelve las propiedades de combinaci�n  
!c      del subgrupo Unifac i1 y tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c
!c-------------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      character fs*8, caract*3

	irec2=i1+(ipareq-1)*150

	read(14,30,rec=irec2)fs,item,item,caract,mm,&
                           mj,mk,mi,mh,rpar,qpar,&
     			    pmg,deltc,delpc,delv,delvi,deln
 30	format (12x,a8,2i4,a3,5i2,8d15.8)
      return
end


