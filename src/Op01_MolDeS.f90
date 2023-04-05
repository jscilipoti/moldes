subroutine Moldes(MolecularDesign)
!-------------------------------------------------------------------------------
!     Subrutinas y Funciones Utilizadas
!     ---------------------------------
!     Programa modificado para que calcule viscosidad
!     Del paquete MANUNI:
!                         NOM_GRU1
!                         Size_LSubGroups()
!                         MAINSG
!                         AINT
!-------------------------------------------------------------------------------
!USE
    use blockdatas
    use CONSTANTES
    use design
    use Evaluation, only:Constraints_Selection,Evaluate_mixture,Evaluate_Pure
    use GRPS
    use Input
    use input_data , only:seleccion_grupos,enter_mixture,rec_inf,leer_comp,checkdatos,elegir_CXR,pvomegaybp
   ! use PROP_PDB, only:checkdatos,elegir_CXR,pvomegaybp
    use PropertiesData
    use PureProp 
    use StructuresDesign
    use SubGrupos
    

    implicit none
!Variables de entrada
    logical,intent(in)::MolecularDesign
!Variables INTERNAS    
    !type(Compound),pointer::recorreSolutes    
    integer::i,j,k,l,n,i0,InputMethod,i2,qLimits,ierror,icomp,isize,i1,tsolv,nousari,SolventsNumber
    integer::idev,idevr,lenfile,imprim,imacy,imaar,imagr,ngrmax,mainsg,ifin0,ifin1,ifin2,imp,imasv
    integer::imasd,imadv,imaev,ima3v,ima4v
    integer::n3v,n4v,jist
    integer,dimension(NSEL):: mgr,maux,lgrup,mgrup,igrupi,iaux,mdum,ngr,jgrup,nsubgr,ngru1,ngru2,grupos
    integer,dimension(NMG)::iar,igrp,idv,icyc,isv,isd,i3v,i4v,grupo1
    integer,dimension(DiffStructGroups)::ngr1,ngr2,ngr3,ngr4,ngr5
    integer,dimension(NCOM,DiffStructGroups,2)::MS
    integer,dimension(NMG,5)::NGSDV, NGSSV
    integer,dimension(NMAX)::nmst,nicomp
    integer,dimension(NMAX,2)::ISOL
    integer,dimension(NMAX,DiffStructGroups,2)::mstemp
    real*8::temp,temp1,dtc2s,solv,slsup1,slsupl,dist,pmma,sclli,select
    real*8::pcch,pcch2,tcch2
    real*8,dimension(20)::pmg
    real*8,dimension(NMG)::deltc,delpc
    real*8,dimension(NMAX,8)::VLIQ
                             real*8::xtemporal(10),act(NCOM),conversion,kexp
                             real*8,external::newton
    character*1::opt,iopt,optinv,sigue
    character*8::nom_gru1
    character*8,dimension(NMG)::fs
    logical::prob,inter,esc,entro,entro1,entro2,entro3,mezcla,solace,salir,puro,nousarl
    character*35::nomb    
    character*37::nomco
    type(Compound),pointer::recorreSolutes  
    type(Compound),pointer::ptr_pcp
    type(FinalStructure),pointer::ptr_FMS,recorre_FMS
    type(CalculatedProperties),pointer::ListProperties
    type(CalculatedPerformance),pointer::ListPerformance 
!COMMONS
    common/EVAL/inv
    logical::INV      
    common/int/ifin2,kgrup,maingr
    integer,dimension(NSEL)::kgrup
    integer,dimension(NMG)::maingr
    common /molcom/arch !arch=false cuando molde2 o rec_inf no encontraron el archivo .mdi
    logical::arch
    COMMON/PUNSUB/npunt,ngrup,ncant
    integer::ncant 
    integer,dimension(NMG)::NPUNT,NGRUP
    COMMON/PUNGRU/npint,nintt,numint
    integer::numint
    integer,dimension(nint)::npint,nintt
    COMMON/GRUPES/ICH,IACOH,IACCH,IAC,IOH,ICOOH,IACCOO
    integer,dimension(NCOM)::ICH,IACOH,IACCH,IAC,IOH,ICOOH,IACCOO
    COMMON/GRUPAL/ICH3,ICH2
    integer,dimension(NCOM)::ICH3,ICH2
    COMMON/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
    real*8::PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
    COMMON/SOLVISOM/NISOMSOLV,IDENTSOLV,NOMBRESOLV,FORMULASOLV,TB,DENST,VRELSOLV
    integer,dimension(NMAX)::nisomsolv
    integer,dimension(NMAX,8)::identsolv
    character*35,dimension(NMAX,8)::nombresolv,formulasolv
    real*8,dimension(NMAX,8)::TB,DENST,VRELSOLV
!EXTERNAL
    external max_sub,mainsg,nom_gru1,Leer_In   

!SENTENCIAS
    idev = 6
	idevr = 5
    MS(:,:,:)=0
    allocate(InputProblem01)    
!---Pregunta si se corre en modo MOLDES Invertido
	InputProblem01%inv = .FALSE.
    if(MolecularDesign)then
        i1 = 0
        do while (i1==0)
            write (6,"(////' ','Select one option :',/&
                   //,11x,' Run a solvent molecular design problem     : 1',&
                   //,11x,' Run a solute molecular design problem      : 2',&
                   //,60x,'> ',$)")
            read (5,'(a)') optinv
            i1 = index ('12',optinv)            
        enddo
        if (i1==2) InputProblem01%inv = .TRUE.
    endif 

!---Opciones de corrida
    InputMethod=0
    do while(InputMethod==0)
        write (idev,"(' ','Options:',//,11x,'Create a new input file : 1',/,11x,'Run a created input file: 2',//,42x,'> ',$)") 
        read (idevr,"(a)") opt
        InputMethod = index ('12',opt)
    enddo
 
    if (InputMethod==1) then !Create a new input file : 1
    
        ! Pide el título del problema    
2110    write (idev,"('Give problem title (maximum 70 characters)')") 
        read (idevr,510,err=2110) titl
        InputProblem01%ProblemTitle = titl
        
        ! Pide el nombre del archivo .mdi        
3100    write (idev,"(/,' Give input file name (maximum 35 characters).')")
        read (idevr,510,err=3100) nomb
        InputProblem01%FileName = nomb
        
        entro2 = .false.
        entro3 = .false.
        arch = .true.
  	    call enter_mixture()! (MOP,temp
        ms(1,:,:)=InputProblem01%MixtureInput%Solutes%Formula(:,:) !Esto es para un solo soluto MODIFICAR!
        ms(2,:,:)=InputProblem01%MixtureInput%PCR%Formula(:,:)  	  
          
        !Carga de subgrupos de solutos en NPUNT y NPINT
        recorreSolutes => InputProblem01%MixtureInput%Solutes
        do while(associated(recorreSolutes))
            k=1
            do while(recorreSolutes%Formula(k,1)/=0)
                call CR_PUNTF(recorreSolutes%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
                k=k+1
            enddo
            recorreSolutes => recorreSolutes%next
        enddo

        !Carga de subgrupos del PCR en NPUNT y NPINT    
        k=1
        do while (InputProblem01%MixtureInput%PCR%Formula(k,1)/=0)
            call CR_PUNTF(InputProblem01%MixtureInput%PCR%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
            k=k+1
        enddo
    
105     entro1 = .false. 

        !Selección de grupos para el diseño molecular
        if(MolecularDesign)then 
            ifam=0 !Seleccion de la familia de componentes a generar
            do while ((ifam.lt.1).or.(ifam.gt.4))
                write (idev,590) !Give the family of solvents to be generated: ....
                read (idevr,*) ifam
            enddo
        
            imprim = 0
            if(inv)then
                write (idev,1601) !now you have to choose
            else      
3050            write (idev,1600) !now you have to choose......screen(1);screen and output
            endif
            read (idevr,980,err=3050) imprim
            write (idev,1610) family(ifam)
            if (imprim.eq.2) then
                lenfile = index(nomb,' ')-1
                open (unit=2,file=nomb(1:lenfile)//'.par',form='formatted')
		        write (2,5000)
		        write (2,1610) family(ifam)
            end if
            entro = .false.
        
1320        CALL SELECCION_GRUPOS(IMPRIM,KGRUP,IFIN2,NGR1,NGR2,NGR3,NGR4,InputProblem01%MainGroups,family)

300         continue

            call Constraints_Selection()
            !Carga de subgrupos intermedios en NPUNT y NPINT
            i=1
            do while(ngdv(i)/=0)
                call CR_PUNTF(ngdv(i),NPUNT,NGRUP,NPINT,NINTT)
                i=i+1
	        enddo

            !Carga de subgrupos terminales en NPUNT y NPINT
            i=1
            do while(ngsv(i)/=0)
                call CR_PUNTF(ngsv(i),NPUNT,NGRUP,NPINT,NINTT)
                i=i+1
	        enddo

            !Carga de subgrupos terminales que no comparten un mismo grupo ppal con subgrupos de valencia dual. 	
            i=1
            do while(ngsv1(i)/=0)
                call CR_PUNTF(ngsv1(i),NPUNT,NGRUP,NPINT,NINTT)
                i=i+1
	        enddo

	        if(InputMethod==2)call Constraints_Selection()
	        ERROR = 1.E-3
            
            !Tipos de enlaces de cada subgrupo elegido para el diseño molecular
	        DO 6 I=1,MDV	!MDV: número de subgrupos intermedios
	        	NGSDV(I,1)=Obtain_attM(NGDV(I))
	        	NGSDV(I,2)=Obtain_attJ(NGDV(I))
	        	NGSDV(I,3)=Obtain_attK(NGDV(I))
	        	NGSDV(I,4)=Obtain_attI(NGDV(I))
	        	NGSDV(I,5)=Obtain_attH(NGDV(I))
   6        CONTINUE
	        DO 7 I=1,MSV	!MSV número de grupos terminales
	        	NGSSV(I,1)=Obtain_attM(NGSV(I))
	        	NGSSV(I,2)=Obtain_attJ(NGSV(I))
	        	NGSSV(I,3)=Obtain_attK(NGSV(I))
	        	NGSSV(I,4)=Obtain_attI(NGSV(I))
	        	NGSSV(I,5)=Obtain_attH(NGSV(I))
   7        CONTINUE
            if (imprim.eq.2) close (unit=2)
        else !if not MolecularDesign; ingreso de solvente a evaluar
            nsol=100 !compuestos por página
            nullify(FMSs)
            nomco = 'Solvent                 '
            write(idev,"(//,' How many solvents do you want to evaluate?:',6X,' > ',$)")
            read(5,*) SolventsNumber
            do i=1,SolventsNumber
                nullify(ptr_FMS)
                allocate(ptr_FMS)
                nomco(9:)= '#'//char(48+i)//':'
                call leer_comp (ipareq,nomco,ptr_FMS%Formula)
                write(6,"(' ','Give the structure type of the solvent you have ',&
     		    'entered:',//,&
                10x,'-Aromatic                           : 1',//&
                10x,'-Single substance groups            : 2',//&
                10x,'-Aliphatic                          : 3',//&
                10x,'-Cyclic                             : 4',///&
     	        51x,' > ',$)")
     	  
     	        read(5,*)ptr_FMS%tsolv
     	        tsolv=ptr_FMS%tsolv
       !!!!CHECKEAR PARÁMETROS DE INTERACCIÓN!!!!!

!Carga de subgrupos del SOLVENTE en NPUNT y NPINT     
                k=1
                do while (ptr_FMS%Formula(k,1)/=0)
                    call CR_PUNTF(ptr_FMS%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
                    k=k+1
                enddo
                call Evaluate_Pure(.True.,ptr_FMS%Formula,InputProblem01%T,tsolv,nousari,nousarl,ListProperties)
                ms(3,:,:)=ptr_FMS%Formula
            
                recorreSolutes => InputProblem01%MixtureInput%Solutes
                nullify(ListPerformance)            
                do while(associated(recorreSolutes))
                    ms(1,:,:) = recorreSolutes%Formula            
                    call Evaluate_Mixture(MolecularDesign,.True.,ms,NCOM,InputProblem01%T,recorreSolutes%BoilingPoint,&
                        tsolv,nousari,nousari,nousarl,ListPerformance,nousari)
                    recorreSolutes => recorreSolutes%next
                enddo
                call Incorporate_Solvent (InputProblem01%mop,ptr_FMS%Formula,tsolv,ListProperties,ListPerformance)
            enddo
        endif

!ESCRITURA INPUT FILE         
        lenfile = index(nomb,' ')-1    
        open(unit=1,file=nomb(1:lenfile)//'.mdi',form='formatted')
        call Write_Input_File(MolecularDesign)
        close (unit=1)
						 
    else !Run a created input file: InputMethod==2
    !Lectura del input file  
        nsol=100
375     write (idev,2060)
		read (idevr,510,err=375) nomb    
        call rec_inf(nomb,MolecularDesign,ptr_FMS,arch)	    
        do while(.not.arch)
            write(6,*) '*****',nomb(1:lenfile)//'.mdi',': file not found *****'
376         write (idev,2060)
		    read (idevr,510,err=376) nomb
		    call rec_inf(nomb,MolecularDesign,ptr_FMS,arch)	
        enddo	
        
        nullify(FMSs)
!        recorre_FMS => ptr_FMS
!        do while(associated(recorre_FMS))
!!Carga de subgrupos del SOLVENTE en NPUNT y NPINT     
!            k=1
!            do while (recorre_FMS%Formula(k,1)/=0)
!                call CR_PUNTF(recorre_FMS%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
!                k=k+1
!            enddo
!            call Evaluate_Pure(.True.,recorre_FMS%Formula,InputProblem01%T,recorre_FMS%tsolv,nousari,nousarl,ListProperties)
!            ms(3,:,:)=recorre_FMS%Formula       
!            
!            ptr_pcp => InputProblem01%MixtureInput%Solutes
!            nullify(ListPerformance)            
!            do while(associated(ptr_pcp))
!                ms(1,:,:) = ptr_pcp%Formula                
!                call Evaluate_Mixture(MolecularDesign,.True.,ms,NCOM,InputProblem01%T,ptr_pcp%BoilingPoint,&
!                    recorre_FMS%tsolv,nousari,nousari,nousarl,ListPerformance,nousari)
!                    ptr_pcp => ptr_pcp%next
!            enddo
!            call Incorporate_Structure (recorre_FMS%Formula,recorre_FMS%tsolv,FMSs,ListProperties,ListPerformance)                             
!            recorre_FMS => recorre_FMS%next
!        enddo
    endif

    CALL Store_Pr (ipareq)
    CALL Store_In ()
    
!Abre archivo de salida
    lenfile = index(nomb,' ')-1
	OPEN (UNIT=2,file=nomb(1:lenfile)//'.MDO',FORM='FORMATTED')

81  CONTINUE 
    SOLV = SOLV/100.	 !SOLV  : Minimum solvent power required (wt)
    SLSUP1 = SLSUP1/100. !SLSUP1: SolventLoss of an intermediate structure      
    SLSUPL = SLSUPL/100. !SLSUPL: SolventLoss of a final structure 

    IF(IFAM.EQ.5)THEN !IFAM = 5 --> Cyclic Solvents
       DELTC(NPUNT(ICH2(IPAREQ)))=.013
       DELPC(NPUNT(ICH2(IPAREQ)))=.184
       DELPC(NPUNT(ICH(IPAREQ)))=.192
    END IF
    IF (IFAM.EQ.2) THEN
       ICOMP = 2
    ELSE IF ((IFAM.EQ.1).OR.(IFAM.EQ.4)) THEN
       ICOMP = 1
    ELSE IF (IFAM.EQ.5) THEN
       ICOMP = 3
    ELSE !IF (IFAM .EQ. 3) THEN
       ICOMP = 0
    END IF

    IF (MOP.EQ.2) THEN !Destilación extractiva
	    IF (TAZEO.EQ.1.0) THEN
		    TEMP1 = BP1		!CAR normal boiling point
	    ELSE
		    TEMP1 = TAZEO	!normal azeotropic temperature
	    END IF
	    PSAT1 = EXP(A1 -A2/(TEMP1 + A3))
	    PSAT2 = EXP(B1 -B2/(TEMP1 + B3))
    endif

!C   PRINTING OF GROUP PARAMETERS FOR THE KEYS			   ! ?
    CALL ESCLIS (2,MolecularDesign,MSOL,NGDV,NGSDV,MDV,NGSV,NGSSV,MSV,A1,A2,A3,B1,B2,B3,TAZEO,X1AZEO,IFAM)
    CALL ESCLIS (6,MolecularDesign,MSOL,NGDV,NGSDV,MDV,NGSV,NGSSV,MSV,A1,A2,A3,B1,B2,B3,TAZEO,X1AZEO,IFAM)


!Diseño Molecular
    
    ms(2,:,:)=InputProblem01%MixtureInput%PCR%Formula(:,:)  	
    salir=.False.
    if(MolecularDesign)CALL STRUCTURE_GENERATOR (JIST,MS,3,NGSSV,NGSDV,SALIR,.FALSE.)
    if(salir)then
        close(unit=2)
        return
    endif

    N=0
!    DO 710 J=1,JIST
!		ISOL(J,1) = J
!		ISOL(J,2) = J
! 710	CONTINUE
    PRINT*,CHAR(7)

    ISIZE=J-1
    call score(nsol)
803 CONTINUE    

    call write_results (mop,6,ifam,nsol)
    call write_results (mop,2,ifam,nsol)

2080 if (arch) then
2068    write(idev,2069)
        read (5,510) opt
		i0 = index ('12',opt)
		if (i0.eq.1) then
			!user=.true.
			!call molde2 (nomb,opt)
		else 
			if (i0.eq.0) then 
				goto 2068
			end if
		end if
		if (mop.eq.2) then
			write (idev,2070) nomb
		else if (mop.eq.1) then
			write(idev,2073) nomb
		end if
		read (5,510) opt
		i1 = index ('123456789',opt)
!c		call limp (idev)
		if (i1.eq.0) then
			go to 2080
		else if (i1.eq.1) then
			!go to 400 CORREGIR
		else if ((i1.eq.2).or.(i1.eq.3).or.(i1.ge.7)) then
			if (.not.(entro2)) then
!c				call limp (idev)
				write (idev,*) '*** Reading data from ',nomb(1:lenfile)//'.mdi'
				call limp (idev)
				call rec_inf (nomb,MolecularDesign,ptr_FMS,arch)
				if (mop.eq.2) then
					dato2(1) = SCLLI 
					dato2(2) = SELECT
					dato2(3) = SOLV  
					dato2(4) = PMMA
				end if
			end if
			if (arch) then
				write (idev,340)
  390				read (idevr,510) opt
				i2 = index ('12',opt)
				if (i2.eq.0) then
					go to 390
				else if (i2.eq.1) then
  380					write (idev,360)
					read (idevr,510,err=380) nomb
				end if
				if ((i1.eq.2).or.(i1.ge.8)) then
					entro2 = .true.
					go to 300
				else if (i1.eq.7) then
					if (.not.(entro2)) then
						call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,imasv,imasd,imadv,ima3v,ima4v,imacy,imaar,imagr)
					end if
					!goto 1860
				else if (i1.eq.3) then
					if (.not.(entro3)) then
!c
!c---Lectura de los datos segun la seleccion de parametros
						write (idev,240) 
						ngrmax = Size_LSubGroups()
						do 315 i=1,ngrmax
							grupo1(i) = i
							maingr(i) = mainsg (i,ipareq)
 315				    continue
						call car_car (ipareq,isv,isd,idv,i3v,i4v,icyc,iar,igrp,imasv,imasd,imadv,ima3v,ima4v,imacy,imaar,imagr)
						ifin0 = 0
						call asignar_grupos_principales (idr,msol,grupos)
						call armar_grupos_distintos (grupos,msol,kgrup,ifin0,ifin1)
						ifin0 = ifin1
						call asignar_grupos_principales (idf,mraf,grupos)
						call armar_grupos_distintos (grupos,mraf,kgrup,ifin0,ifin2)
					end if
					entro2 = .true.
					entro3 = .true.
               		if (i1.eq.3) then
						go to 105
					else
						goto 1320
					end if
				end if
			else
 314				write (idev,2060)
				read (idevr,510,err=314) nomb
				entro2 = .false.
				entro3 = .false.
				arch = .true.
				!go to 400 CORREGIR
			end if
		else if (i1.eq.4) then
			call ci_ban1
			CLOSE(UNIT=IMP)
			go to 2110
		else if (i1.eq.5) then
 310			write (idev,2060)
			read (idevr,510,err=310) nomb
			entro2 = .false.
			entro3 = .false.
			!go to 400 CORREGIR
		end if
	else
 316		write (idev,2060)
		read (idevr,510,err=316) nomb
		entro2 = .false.
		entro3 = .false.
		arch = .true.
		!go to 400 CORREGIR
	end if
10000 call ci_ban1
    !  CLOSE(UNIT=IMP)
      1111    CLOSE (UNIT=2)
!C	CLOSE (UNIT=4)
 2020	continue

!c---- formatos
2069	format ('  Do you want to evaluate other specific ', 'solvents?',//,50x,'yes: 1',/,51x,'no: 2',/,60x,'>',$)
   1	format ('  Choose the restriction level for the Molecular Design',' of Solvents:',//,&
     		'  To use the original Brignole Feasibility Criteria ',&
     		/,'  (better UNIFAC predictions):',50x,'1',/,&
     		'  To perform a "less restricted" synthesis by using the',&
     		' modified BFC (more ',&
     		/,'  solvent structures but possibly reactive and',&
     		' less reliable predictions):',6x,'2')

 360  format (' ',//,' Enter the new input file name ','(maximum 35 characters): ',$)
 340  format (' Do you want to keep the last results file?',//,&
             '                                          Yes: 1',/,&
             '                                          No : 2',/,&
             60x,'> ',$)
5000  format (' This file has UNIFAC pairs of groups which have no interaction parameters.',//)
2070  format (' ','Options:',/&
            /,11x,'Run the ',a6,' separation problem with',&
            /,11x,'                               the same data: 1',&
            /,11x,'            different restriction parameters: 2',&
            /,11x,'    a new solvent family or different groups: 3',&
            /,11x,'            different combination properties: 7',&
            /,11x,'different evaluation of pre-final structures: 8',&
            /,11x,'Create a new input file                     : 4',&
            /,11x,'Run another file                            : 5',&
            /,11x,'Exit                                        : 6',&
    	      //,60x,'> ',$)
2073  format (' ','Options:',/&
             /,11x,'Run the ',a6,' separation problem with',&
             /,11x,'                               the same data: 1',&
             /,11x,'            different restriction parameters: 2',&
             /,11x,'    a new solvent family or different groups: 3',&
             /,11x,'            different combination properties: 7',&
             /,11x,'different evaluation of pre-final structures: 8',&
!c    *        /,11x,'            different restriction level: 8',&
     	    /,11x,'                        other temperature   : 9',&
             /,11x,'Create a new input file                     : 4',&
             /,11x,'Run another file                            : 5',&
             /,11x,'Exit                                        : 6',&
     	      //,60x,'> ',$)

2060  format (' ',//,' Enter the input file name (maximum 35 characters): ',$)
1910  format (' ','*** Checking interaction parameters ***',/)
 410  format (' ',/,' Do you want to change some group combination property ? (Y/N)  > ',$)
1600  format (' ',/,' * Now you have to choose the groups that take ',&
             'part in the Molecular',/,&
                   '   Design of Solvents.  The program performs a ',&
             'preliminary selection',/,&
                   '   in  order to  eliminate  those groups which ',&
             'have  no  interaction',/,&
                   '   parameters  with  respect  to  the others ',&
             'main components.',/,&
                   '   If you want to  see  the  pairings  without ',&
             'interaction parameters',/,&
                   '   enter the following code:',/,&
              22x,'screen:                              1',/,&
              22x,'screen and output file (nointer):    2',/,&
              22x,'no printing:                      <ret>',/,&
              63x,'> ',$)
1601  format (' ',/,' * Now you have to choose the groups that take ',&
             'part in the Molecular',/,&
                   '   Design of Solutes.  The program performs a ',&
             'preliminary selection',/,&
                   '   in  order to  eliminate  those groups which ',&
             'have  no  interaction',/,&
                   '   parameters  with  respect  to  the others ',&
             'main components.',/,&
                   '   If you want to  see  the  pairings  without ',&
             'interaction parameters',/,&
                   '   enter the following code:',/,&
              22x,'screen:                              1',/,&
              22x,'screen and output file (nointer):    2',/,&
              22x,'no printing:                      <ret>',/,&
              63x,'> ',$)
1610  format (' ',/,' *** ',a21,' preliminary selection ***')
1500  format (' ',/,' * No interaction parameters between the aromatic',&
             ' part of the',/,'   structure of the solvent and the ',&
             'dual valence groups already',/,'   entered.')
1510  format (' ',/,' * No interaction parameters between the aromatic',&
            ' part of the',/,'   structure of the solvent and the ',&
            'single valence groups',/,'   already entered.')
1490  format (' ',/,' * No interaction parameters available for ',&
             'this family of solvents',/,'   in the ',a17,&
             ' UNIFAC table')
1040  format ('   Change the family of groups:  <ret>. ',$)
1080  format (' ',/,' ','* No interaction parameters for some of the ',&
            'subgroups',/,'   already entered and the two main ',&
            'components.',/,'   Please change the subgroups: <ret>. '&
            ,$)
1090  format ('   Do you want to change some subgroups?',/,&
             '                         First subgroups entered:  1',/,&
             '                         Single valence subgroups: 2',&
             /,60x,'> ',$)
 970  format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 800  format (' ',a63,f10.4)
 801  format (' 12 Tem : Operation temperature                   ','          (K)',f8.2)
 802  format (' ',/,' Give new operation temperature                 ','                    > ',$)
 810  format (' ',a64,1x,i8)
 990  format (' ',a63,2x,i8)
 820  format ('  ')
 830  format (' ','Code',3x,'Parameter',43x,'Default Value',/)
 840  format (' ',/,' If you want to change a parameter value give code',' number.','  If not 0.',/,71x,'> ',$)
 860  format (' ',a63,7x,'> ',$)
 870  format (' ',a64,6x,'> ',$)
 890  format (' ',a63,7x,'> ',$)
 880  format (' ','Code',3x,'Parameter',47x,'New Value',/)
 !610  format (a70)
! 620  format (20i3)
! 630  format (6i3)
! 625  format (11i3)
 !626  format (3i3)
 500  format (2X,'****************************************************',&
             '************'&
             ,//,27X,'MOLDES 3.0',&
             //,6X,' Based on the algorithm of Gani et al. (1983),',/,&
             6X,' Brignole et al. (1986) and Pretel et al. (1994).',//,&
             6x,' Present version prepared by Patricio Araya, ',&
             'Eduardo Pretel',/,6X,&
             ' and completed by Martin Cismondi from the',&
     		' original version',/,6x,' of E. A. Brignole'&
             ,//,6X,&
             ' July, 2000.',//,6X,&
             ' Planta Piloto de Ingenieria Quimica, 8000 Bahia Blanca,'&
             ,/,6X,&
             ' Argentina.',//,2X,&
             '***************************************************',&
             '************',//,&
             '  MOLDES 3.0 is a Group Contribution Methodology to ',&
             'Computer Aided',/,&
     	      '  Molecular  Design.    It   selects  solvents  for ',&
             ' Liquid-Liquid',/,& 
             '  Extraction  or  Extractive  Distillation  Separati',&
             'ons of  binary',/,& 
             '  mixtures.   The selection  is  based on  UNIFAC   ',&
             'prediction  of',/,&
             '  solvent properties and on physical and molecular  ',&
             'constraints.')
 510  format (a) 
 240  format (/,'*** Loading data ***')
 980  format (i2)
 680  format (20i3)
 590  format (/,' Give the family of solvents to be generated:',///,&
             10x,'-Aromatic structures                : 1',//&
             10x,'-Single substance groups            : 2',//&
             10x,'-No-aromatics with up to 12 groups in',/&
             10x,' the final structure                : 3',//&
             10x,'-Cyclic structures                  : 4',///&
     	      51x,' > ',$)

    return
endsubroutine Moldes

