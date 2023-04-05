subroutine score (nsol)
!--------------------------------------------------------------
!   Esta subroutine clasifica las estructuras de acuerdo a 
!   alguna de sus propiedades. Una vez realizada la clasifica-
!   ción, se ordenan las estructuras de manera ascendente en
!   función de la cantidad de grupos funcionales (para cade-
!   nas alifáticas) o carbonos aromáticos sustituídos (para 
!   estructuras aromáticas).
!--------------------------------------------------------------
use StructuresDesign
use Input
use input_data,only:Search_Isomers
implicit none
!INTERFACES
interface
    subroutine Sort_by_Property(punt)
    use StructuresDesign
    implicit none
    type(FinalStructure),pointer,intent(inout)::punt
    endsubroutine Sort_by_Property
endinterface
interface
    subroutine ordenar_FunctionalGroups(punt)
    use StructuresDesign
    implicit none
    type(FinalStructure),pointer,intent(inout)::punt
    endsubroutine ordenar_FunctionalGroups
    endinterface
!..Variables de ENTRADA/SALIDA
integer,intent(in)::nsol
!Variables INTERNAS/COMMONS
type(FinalStructure),pointer::recorre
integer::i

!Sentencias----------------------------------------------------   
    call Sort_by_Property(FMSs)
    i=1
    recorre=>FMSs
    do while(.True.)
        recorre%position=i
        if(.not.associated(recorre%next))exit
        recorre=>recorre%next
        i=i+1
    enddo
    if (ipareq.eq.2) call Search_Isomers(nsol) !Adaptar para las otras tablas de parámetros
    call ordenar_FunctionalGroups (FMSs)
endsubroutine score

!==============================================================
subroutine Sort_by_Property(punt)
!--------------------------------------------------------------
!   Esta subrutina clasifica estructuras en función de la 
!   propiedad seleccionada
!--------------------------------------------------------------
use StructuresDesign
use Input
implicit none
!INTERFACES
interface
    subroutine acomodar (punt,anterior,recorre)
    use StructuresDesign
    implicit none
    type(FinalStructure),pointer,intent(inout)::punt,anterior,recorre
    endsubroutine acomodar
endinterface
!Variables de ENTRADA/SALIDA
type(FinalStructure),pointer,intent(inout)::punt
!Variables INTERNAS
type(FinalStructure),pointer::recorre
type(FinalStructure),pointer::tope
type(FinalStructure),pointer::anterior
logical::orden
!Sentencias----------------------------------------------------
    recorre=>punt
    orden = .False.
    if(InputProblem01%mop/=4)then !PROVISORIO para reacciones !!!!
    do while(associated(recorre%next).and..not.associated(recorre%next,tope))
        do while(associated(recorre%next).and..not.associated(recorre%next,tope))
          !  if(.True..and.(recorre%DistCoefficient < recorre%next%DistCoefficient))orden = .True.
            if(.True..and.(recorre%SolventPower < recorre%next%SolventPower))orden = .True.
            if(orden)then
                call acomodar(punt,anterior,recorre)
                orden = .False.
                cycle
            endif
            anterior=>recorre
            recorre=>recorre%next
        enddo
        anterior=>recorre
        tope=>recorre
        recorre=>punt
    enddo
    else
    do while(associated(recorre%next).and..not.associated(recorre%next,tope))
        do while(associated(recorre%next).and..not.associated(recorre%next,tope))
          !  if(.True..and.(recorre%DistCoefficient < recorre%next%DistCoefficient))orden = .True.
            if(.True..and.(recorre%conversion < recorre%next%conversion))orden = .True.
            if(orden)then
                call acomodar(punt,anterior,recorre)
                orden = .False.
                cycle
            endif
            anterior=>recorre
            recorre=>recorre%next
        enddo
        anterior=>recorre
        tope=>recorre
        recorre=>punt
    enddo        
    endif
endsubroutine Sort_by_Property

!======================================================================
subroutine ordenar_FunctionalGroups(punt)
!----------------------------------------------------------------------
!   Esta subrutina clasifica estructuras en función de la cantidad de 
!   grupos funcionales (alifáticas)o carbonos sustituídos (aromáticas)
!   en su estructura.
!----------------------------------------------------------------------
use StructuresDesign
implicit none
!INTERFACES
interface
    subroutine acomodar (punt,anterior,recorre)
    use StructuresDesign
    implicit none
    type(FinalStructure),pointer,intent(inout)::punt,anterior,recorre
    endsubroutine acomodar
endinterface
!Variables de ENTRADA/SALIDA
type(FinalStructure),pointer,intent(inout)::punt
!Variables INTERNAS
type(FinalStructure),pointer::recorre
type(FinalStructure),pointer::tope
type(FinalStructure),pointer::anterior
logical::orden
!Sentencias------------------------------------------------------------
    recorre=>punt
    orden = .False.
    do while(associated(recorre%next).and..not.associated(recorre%next,tope))
        do while(associated(recorre%next).and..not.associated(recorre%next,tope))
            if(.True..and.(recorre%FunctionalGroups > recorre%next%FunctionalGroups))orden = .True.
            if(orden)then
                call acomodar(punt,anterior,recorre)
                orden = .False.
                cycle
            endif
            anterior=>recorre
            recorre=>recorre%next
        enddo
        anterior=>recorre
        tope=>recorre
        recorre=>punt
    enddo
endsubroutine

!======================================================================
subroutine acomodar (punt,anterior,recorre)
!----------------------------------------------------------------------
!   Esta subrutina intercambia las posiciones de dos nodos adyascentes
!   (recorre y anterior) de la lista a la que apunta "punt"
!----------------------------------------------------------------------
use StructuresDesign
implicit none
!Variables de ENTRADA/SALIDA
type(FinalStructure),pointer,intent(inout)::punt,anterior,recorre
!Variables INTERNAS
type(FinalStructure),pointer::aux
!Sentencias----------------------------------------------------
  if(associated(recorre,punt))then
      aux=>punt
      punt=>recorre%next
      recorre%next=>recorre%next%next
      punt%next=>aux
      anterior=>punt
  elseif(.not.associated(recorre%next%next))then
      aux=>recorre%next
      nullify(recorre%next)
      anterior%next=>aux
      aux%next=>recorre
      recorre=>aux
  else
      aux=>recorre%next
      recorre%next=>recorre%next%next
      anterior%next=>aux
      aux%next=>recorre
      anterior=>aux
  endif
endsubroutine

!==============================================================
subroutine MoldesInvertido (MS)
!--------------------------------------------------------------
!     Descripción:
!     Esta subrutine invierte las posiciones entre el soluto y
!     el solvente en la matriz MS para poder Evaluate las 
!     propiedades de mezcla según el enfoque "Moldes Invertido"
!--------------------------------------------------------------          
implicit none
integer, intent(inout)::MS(3,10,2)
integer::MSIN(10,2),L
    MSIN(:,:)=0
    do L=1,10
        MSIN(L,1) = MS(3,L,1)
		MSIN(L,2) = MS(3,L,2)
        MS(3,L,1) = MS(1,L,1)
        MS(3,L,2) = MS(1,L,2) 
        MS(1,L,1) = MSIN(L,1)
        MS(1,L,2) = MSIN(L,2)
    enddo
endsubroutine
 
!==============================================================
recursive subroutine COMPARAR (IWANT,ISUS,A,IFIND)
!--------------------------------------------------------------

!-------------------------------------------------------------- 
      INTEGER IWANT(10,2),ISUS(10,2)
      INTEGER A
      LOGICAL IFIND
      IF (IWANT(A,1).NE.ISUS(A,1).OR.IWANT(A,2).NE.ISUS(A,2)) THEN
        IFIND=.FALSE.
        return
      ELSE IF (A.EQ.10) THEN
        IFIND=.TRUE.
        return
      ELSE
        CALL COMPARAR (IWANT,ISUS,A+1,IFIND)
      ENDIF      
      return
endsubroutine
!==============================================================
SUBROUTINE ORDENAR (ORD)
!--------------------------------------------------------------

!-------------------------------------------------------------- 
      IMPLICIT NONE
      INTEGER ORD (10,2)
      INTEGER I,J,AUX,AUXB
      DO I=10,1,-1
        DO J=10,12-I,-1
            IF (ORD(J,1).LT.ORD(J-1,1).AND.ORD(J,1).GT.0) THEN
                AUX=ORD(J,1)
                AUXB=ORD(J,2)
                ORD(J,1)=ORD(J-1,1)
                ORD(J,2)=ORD(J-1,2)
                ORD(J-1,1)=AUX
                ORD(J-1,2)=AUXB
            ENDIF
        ENDDO
      ENDDO
      RETURN
endsubroutine

        SUBROUTINE CR_PUNT (MS,MSOL,MRAF,NGDV,MDV,NGSV,MSV,NGSV1,MSV1,IPAREQ,IFAM,IALAR,NALAR)
!C--------------------------------------------------------------------------
!C       ESTA SUBRUTINA CREA EL VECTOR NGRUP, QUE CONTIENE TODOS LOS GRUPOS 
!C       PRESENTES EN EL CASO QUE SE ESTA CORRIENDO.
!C	  NPUNT es el vector que contiene las posiciones que los distintos
!c	  grupos ocupan en el vector NGRUP. ej. NPUNT(9)=1 ; NPUNT(12)=2 ...
!C--------------------------------------------------------------------------
!C
        PARAMETER (NSCM=10,NCOM=3,NA=150,NG=70,NMOD=3,NMG=150,NSVA=20,NAA=4)
        implicit real*8 (A-H,O-Z)
        EXTERNAL MAINSG
        DIMENSION MS(NCOM,NSCM,2),NGDV(0:NA),NGSV(NSVA),NGSV1(NSVA),IALAR(NAA)
        COMMON/PUNSUB/NPUNT(NA),NGRUP(NA),NUM
        COMMON/GRUPAR/IACH(NMOD),IACCH2(NMOD),IACCL(NMOD),IACNH2(NMOD),IACNO2(NMOD)
        COMMON/PUNGRU/NPINT(NG),NINTT(NG),NUMINT
        COMMON /GEOM/ R(NMG),Q(NMG),MAIN(NMG)

	  DO 10 I=1,NA !NA = 150
           NPUNT(I) = 0
  10    CONTINUE
        DO 70 I=1,NG
           NPINT(I) = 0 !NG = 70
  70    CONTINUE
        NUM = 0
        DO 20 I=1,MSOL ! MSOL = cantidad de grupos distintos del CAR
           NUM = NUM + 1
           NPUNT(MS(1,I,1)) = NUM
           NGRUP(NUM) = MS(1,I,1)
  20    CONTINUE
	  DO 30 I=1,MRAF ! MRAF = cantidad de grupos distintos del CPR
           IF (NPUNT(MS(2,I,1)).EQ.0) THEN
              NUM = NUM + 1
              NPUNT(MS(2,I,1)) = NUM
              NGRUP(NUM) = MS(2,I,1)
           END IF
  30    CONTINUE
        DO 40 I=1,MDV ! MDV= cant de subgrupos intermedios inresados por el user
           IF (NPUNT(NGDV(I)).EQ.0) THEN !NGDV = cont los núm ident de c/ grupo
              NUM = NUM + 1
              NPUNT(NGDV(I)) = NUM
              NGRUP(NUM) = NGDV(I)
           END IF
  40    CONTINUE
        DO 50 I=1,MSV ! MSV= cant de subgrup terminales sel por el usuario
           IF (NPUNT(NGSV(I)).EQ.0) THEN
              NUM = NUM + 1
              NPUNT(NGSV(I)) = NUM
              NGRUP(NUM) = NGSV(I)
           END IF
  50    CONTINUE
        DO 60 I=1,MSV1 !MSV1= MSV que no comparten un (...) valencia dual
           IF (NPUNT(NGSV1(I)).EQ.0) THEN
              NUM = NUM + 1
              NPUNT(NGSV1(I)) = NUM
              NGRUP(NUM) = NGSV1(I)
           END IF
  60    CONTINUE
        IF (IFAM.EQ.4) THEN
           DO 90 I=1,NALAR
              IF (IALAR(I).EQ.1) THEN
                 IF (NPUNT(IACH(IPAREQ)).EQ.0) THEN
                    NUM = NUM + 1
                    NPUNT(IACH(IPAREQ)) = NUM
                    NGRUP(NUM) = IACH(IPAREQ)
                 END IF
                 IF (NPUNT(IACCH2(IPAREQ)).EQ.0) THEN
                    NUM = NUM + 1
                    NPUNT(IACCH2(IPAREQ)) = NUM
                    NGRUP(NUM) = IACCH2(IPAREQ)
                 END IF
              ELSE IF ((IALAR(I).EQ.2).AND.(NPUNT(IACCL(IPAREQ)).EQ.0)) THEN
                 NUM = NUM + 1
                 NPUNT(IACCL(IPAREQ)) = NUM
                 NGRUP(NUM) = IACCL(IPAREQ)
              ELSE IF ((IALAR(I).EQ.3).AND.(NPUNT(IACNH2(IPAREQ)).EQ.0)) THEN
                 NUM = NUM + 1
                 NPUNT(IACNH2(IPAREQ)) = NUM
                 NGRUP(NUM) = IACNH2(IPAREQ)
              ELSE IF ((IALAR(I).EQ.4).AND.(NPUNT(IACNO2(IPAREQ)).EQ.0))THEN
                 NUM = NUM + 1
                 NPUNT(IACNO2(IPAREQ)) = NUM
                 NGRUP(NUM) = IACNO2(IPAREQ)
              END IF
  90       CONTINUE
        END IF
        NUMINT = 0
        DO 80 I=1,NUM
           MAIN(I) = MAINSG(NGRUP(I),IPAREQ)
	   IF (NPINT(MAIN(I)).EQ.0) THEN
              NUMINT = NUMINT + 1
              NPINT(MAIN(I)) = NUMINT
              NINTT(NUMINT) = MAIN(I)
           END IF
  80    CONTINUE
        RETURN
        ENDsubroutine

        SUBROUTINE CR_PUNT2 (NC,IPAREQ)
!--------------------------------------------------------------------------
!       ESTA SUBRUTINA CREA EL VECTOR NGRUP, QUE CONTIENE TODOS LOS GRUPOS 
!       PRESENTES EN EL CASO QUE SE ESTA CORRIENDO.
!	  NPUNT es el vector que contiene las posiciones que los distintos
!	  grupos ocupan en el vector NGRUP. ej. NPUNT(9)=1 ; NPUNT(12)=2 ...
!--------------------------------------------------------------------------
!
        PARAMETER (NSCM=10,NCOM=3,NA=150,NG=70,NMOD=3,NMG=150,NSVA=20,NAA=4)
        implicit real*8 (A-H,O-Z)
        EXTERNAL MAINSG
      COMMON/US/MS(NCOM,NSCM,2),nms(NCOM)
        COMMON/PUNSUB/NPUNT(NA),NGRUP(NA),NUM
        COMMON/PUNGRU/NPINT(NG),NINTT(NG),NUMINT
        COMMON /GEOM/ R(NMG),Q(NMG),MAIN(NMG)
!
	  DO 10 I=1,NA
           NPUNT(I) = 0
  10    CONTINUE
        DO 70 I=1,NG
           NPINT(I) = 0
  70    CONTINUE
        NUM = 0
        DO 20 I=1,Nms(1)
           NUM = NUM + 1
           NPUNT(ms(1,I,1)) = NUM
           NGRUP(NUM) = ms(1,I,1)
  20    CONTINUE
	DO J=2,NC
	  DO 30 I=1,Nms(J)
           IF (NPUNT(ms(J,I,1)).EQ.0) THEN
              NUM = NUM + 1
              NPUNT(ms(J,I,1)) = NUM
              NGRUP(NUM) = ms(J,I,1)
           END IF
  30    CONTINUE
	END DO
        NUMINT = 0
        DO 80 I=1,NUM
           MAIN(I) = MAINSG(NGRUP(I),IPAREQ)
	   IF (NPINT(MAIN(I)).EQ.0) THEN
              NUMINT = NUMINT + 1
              NPINT(MAIN(I)) = NUMINT
              NINTT(NUMINT) = MAIN(I)
           END IF
  80    CONTINUE
        RETURN
endsubroutine

!==========================================================================
SUBROUTINE CR_PUNTF (IGRUP,NPUNT,NGRUP,NPINT,NINTT)
!--------------------------------------------------------------------------
!   ESTA SUBRUTINA CREA EL VECTOR NGRUP, QUE CONTIENE TODOS LOS GRUPOS 
!   PRESENTES EN EL CASO QUE SE ESTA CORRIENDO.
!	NPUNT es el vector que contiene las posiciones que los distintos
!	grupos ocupan en el vector NGRUP. ej. NPUNT(9)=1 ; NPUNT(12)=2 ...
!--------------------------------------------------------------------------
!   
    use Input
    !use CONSTANTES, only:NINT,NMG
    implicit none

!Variables de ENTRADA
    integer,intent(in)::igrup
!Variables de ENTRADA/SALIDA      
    integer,dimension(NMG)::NPUNT,NGRUP
    integer,dimension(NINT):: NPINT,NINTT
!Variables INTERNAS      
    integer::j,k,maingrup
!COMMONS      
 !   common/PAREQ/IPAREQ
 !   integer::ipareq
    common/GEOM/R,Q,MAIN
    integer,dimension(NMG)::MAIN
    real*8,dimension(NMG)::R,Q
!FUNCTIONS    
    integer,external::mainsg

!SENTENCIAS
    if(npunt(igrup)==0)then             !PREGUNTA SI ESTÁ CARGADO EL SUBGRUPO
      do j=1,NMG                        !BUSCA PRIMER POSICIÓN VACÍA
          if(ngrup(J)==0)exit
      enddo
      npunt(igrup) = J
      ngrup(J) = igrup
      maingrup = mainsg(ngrup(J),ipareq)!BUSCA A QUÉ GRUPO PPAL PERTENECE
      if(npint(maingrup)==0)then        !PREGUNTA SI ESTÁ CARGADO EL GRUPO PPAL
          do K=1,NINT                   !BUSCA PRIMER POSICIÓN VACÍA
              if(nintt(K)==0)exit
          enddo
          npint(maingrup) = k
          nintt(K) = maingrup
          main(k) = maingrup
      endif      
    endif
    
endsubroutine CR_PUNTF


      subroutine ordenar_arreglo (arreglo,nelem,elem)
!c-----------------------------------------------------------------------------
!c     Esta subrutina ordena los primeros "elem" elementos del vector "arreglo" 
!      por medio  del  metodo  de  la burbuja. 
!c 
!c     Variables de entrada:
!c     -------------------
!c                 arreglo: vector de numeros enteros desordenado.  
!c                          Su dimension es nelem.
!c                   nelem: dimension del vector arreglo. 
!                     elem: número de elementos a ordenar
!c
!c     Variables de salida:
!c     -------------------
!c                 arreglo: vector de numeros enteros ordenado
!c
!c-----------------------------------------------------------------------------
!c
      implicit none
      integer::nelem,i,j,temp,elem
      integer,dimension(nelem)::arreglo

      do 20 i=1,elem
         do 30 j=elem,i+1,-1
            if (arreglo(j-1).gt.arreglo(j)) then
               temp = arreglo(j-1)
               arreglo(j-1) = arreglo(j)
               arreglo(j) = temp
            end if
 30      continue
 20   continue
      return
      end


      subroutine asignar_grupos_principales (subgrup,nsub,grupos)
!c-----------------------------------------------------------------------------
!c     Esta subrutina genera el  vector  de  grupos  principales  que 
!c     corresponde al vector de subgrupos subgrup
!c
!c     Variables de entrada:
!c     --------------------
!c                   maingr: puntero que apunta a la lista con todos los 
!                            subgrupos de la tabla ipareq
!c                  subgrup: vector de numeros de identificacion de  subgrupos.
!c                           Su dimension es nsel.
!c                     nsub: cantidad de elementos en subgrup
!c
!c     Variables de salida:
!c     -------------------
!c                   grupos: vectpr de numeros de identificacion de los grupos
!c                           principales del vector de subgrupos subgrup.
!c
!c     Los parametros que fijan la dimension de los vectores son:
!c
      use CONSTANTES
      use SubGrupos
!c
!c-----------------------------------------------------------------------------
!c
      implicit none
      integer::grupos(nsel),subgrup(nsel),nsub !,maingr(NMG)
      integer::i
      
      i=1
      grupos(:)=0
      do while (subgrup(i)/=0)
            grupos(i)= Obtain_MainGroup_Number(subgrup(i))
            i=i+1
      enddo
      return
      end

subroutine armar_grupos_distintos (grupos,ngrup,kgrup,ifin0,ifin)
!c-----------------------------------------------------------------------------
!c     Esta subrutina ingresa al vector kgrup los grupos principales del vector
!c     grupos que no estaban previamente en kgrup.
!c
!c     Variables de entrada:
!c     --------------------
!c                   grupos: vector de  numeros  de  identificacion  de  grupos 
!c                           principales  que  pueden  ingresar  a  kgrup.   Su 
!c                           dimension es nsel.
!c                    ngrup: cantidad de elementos en el vector grupos.
!c                    kgrup: vector de  numeros  de  identificacion  de  grupos
!c                           principales de entrada.  Su dimension es nsel.
!c                    ifin0: cantidad  de  elementos  de  entrada  en el vector
!c                           kgrup.
!c
!c     Variables de salida:
!c     -------------------
!c                    kgrup: vector de  numeros  de  identificacion  de  grupos
!c                           principales de salida.
!c                     ifin: cantidad de elementos de salida en el vector kgrup.
!c
!c     Los parametros que fijan las dimensiones de los vectores son:
!c
      use CONSTANTES
!c
!c-----------------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      logical hacer
      integer grupos(nsel),kgrup(nsel)

      ifin = ifin0
      ICERO=IFIN0+1
!      DO 15 I=ICERO,NSEL
!		KGRUP(I)=0
!15    CONTINUE



      do 10 i=1,ngrup
		hacer = .true.
		if (ifin.gt.0) then
			do 360 j=1,ifin
				if (.not.(hacer)) go to 370
				if (kgrup(j).eq.grupos(i)) then
					hacer = .false.
				end if
 360			continue
 370		end if
		if (hacer) then
			kgrup(ifin+1) = grupos(i)
			ifin = ifin + 1
		end if
 10   continue

      return
      end
      
!C ============================================================================
      subroutine Sel_groups (igrp,max,kgrup,iult,mgr,igr)
!c-----------------------------------------------------------------------------
!c     Esta subrutina permite seleccionar los subgrupos por tipo de moleculas
!c     a generar, almacenados en el vector igrp.
!c
!c     Variables de entrada:
!c     --------------------
!c                  igrp: vector de numeros de identificacion de los subgrupos 
!c                        seleccionables.  Su dimension es nsel.
!c                   max: cantidad de elementos en el vector igrp.
!c                 kgrup: 
!c                  iult: cantidad de elementos en el vector kgrup.
!             LSubGroups: Lista con todos los subgrupos de la tabla ipareq
!c
!c     Variables de salida:
!c     -------------------
!c     
!c                   mgr: vector de  numeros  de  identificacion   de subgrupos
!c                        seleccionados.  Su dimension es nsel.
!c                   igr: cantidad de subgrupos seleccionados.
!c
!c     Los parametros que fijan las dimensiones tienen los siguientes valores:
!c
      use CONSTANTES
      use SubGrupos
!c
!c-----------------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      dimension igrp(nsel),mgr(nsel),kgrup(nsel)
      INTEGER grusel(nsel)
      logical enc
      COMMON/GSELECCIONADOS/GRUSEL

      enc = .false.
      igr = 0
      K=1
      DO WHILE (GRUSEL(K).NE.0)
        K=K+1
      ENDDO 
      do 10 i = 1,max
         do 40 j = 1,iult
            if (Obtain_MainGroup_Number(igrp(i)).eq.kgrup(j)) then
               enc = .true.
            end if
 40      continue
         if (enc) then 
3010        write (6,20) Obtain_SubGroup_Name(igrp(i))
            read (5,30,err=3010) is
            if (is.eq.1) then
               igr = igr + 1
               mgr(igr) = igrp(i)
               GRUSEL(K+IGR-1)=IGRP(I)
            end if
         end if
         enc = .false.
 10   continue
!c
!c---- formatos
 20   format (43x,a8,' > ',$)
 30   format (i2)
      return
      end

      SUBROUTINE CR_PUNT3(NCANTCOM,COMPUESTOS,NCANT,NPUNT,NGRUP,CANTGRUP)
!--------------------------------------------------------------------------
!       ESTA SUBRUTINA CREA EL VECTOR NGRUP, QUE CONTIENE TODOS LOS GRUPOS 
!       PRESENTES EN EL CASO QUE SE ESTA CORRIENDO.
!	  NPUNT es el vector que contiene las posiciones que los distintos
!	  grupos ocupan en el vector NGRUP. ej. NPUNT(9)=1 ; NPUNT(12)=2 ...
!--------------------------------------------------------------------------
!
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (NG=10,NMG=150)
      INTEGER COMPUESTOS(NCANTCOM,NG,2)
      DIMENSION NPUNT(NMG),NGRUP(NMG)
      INTEGER CANTGRUP(NCANTCOM)
      LOGICAL LOG
!
	DO 10 I=1,NG
        NPUNT(I) = 0
  10  CONTINUE

      NCANT = 0
      DO 20 I=1,NCANTCOM
        NCANTGRUP=0
        DO 30 J=1,10
            IF (COMPUESTOS(I,J,1)/=0) THEN
                NCANTGRUP=NCANTGRUP+1
                CALL BUSCARAS (COMPUESTOS(I,J,1),NGRUP,NMG,LOG)
                IF (.NOT.LOG) THEN
                    NCANT = NCANT + 1
                    NPUNT(COMPUESTOS(I,J,1)) = NCANT
                    NGRUP(NCANT) = COMPUESTOS(I,J,1)
                ENDIF
            ENDIF
  30    CONTINUE
        CANTGRUP(I)=NCANTGRUP
  20  CONTINUE

      RETURN
      END

      SUBROUTINE BUSCARCADENA (NGRUP,NCANT,GRUPO,LOG)
!--------------------------------------------------------------------------
!     Esta subrutina busca en el vector NGRUP los grupos que ya fueron    
!     cargados. Si encuentra un grupo cargado log=true.
!--------------------------------------------------------------------------
      
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (NA=30)
      CHARACTER GRUPO*35, NGRUP(NA)*35
      LOGICAL LOG
      
      LOG=.FALSE.
      DO I=1,NCANT
        IF (NGRUP(I)==GRUPO)THEN
            LOG=.TRUE.
            EXIT
        ENDIF
      ENDDO
      RETURN
      END

SUBROUTINE BUSCARAS (N,VECTOR,tam,VL)
!c---------------------------------------------------------------------
!c     Comprueba la existencia del grupo N en VECTOR
!c--------------------------------------------------------------------
    implicit none
    integer::tam,n,i
    INTEGER VECTOR(tam)
    LOGICAL VL

    VL=.FALSE.
    I=0
    DO WHILE (.NOT.VL)
        I=I+1
	  IF (N.EQ.VECTOR(I)) VL=.TRUE.
	  IF (I.EQ.tam) EXIT
	ENDDO
endsubroutine

