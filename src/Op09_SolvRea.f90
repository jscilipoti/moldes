subroutine SolvRea ()
!-------------------------------------------------
!   Descripci�n
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    use Input
    use GRPS
    use StructuresDesign
    use CONSTANTES
    use design
    use input_data,only:Leer_Comp,seleccion_grupos
    use Evaluation, only:Constraints_Selection
    use SubGrupos
    
    implicit none
!INTERFACES
interface
    subroutine ab_ban1(mod)
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface  

    integer,parameter::NSVA=20,NSCM=10,NMAX=60000,NG=70
    integer,external::max_sub
    
!...Integers
    integer::maingr(NMG),grupo1(NMG),ids(nsel),nys(nsel),MS(NCOM,DiffStructGroups,2),nr,np,MGR(nsel),comp(DiffStructGroups,2)
    integer::NGSSV(NSVA,5),NGSDV(NSVA,5),MM(NMG),MJ(NMG),MK(NMG),MI(NMG),MH(NMG)
    integer::MST(NMAX,NSCM,2), iter
    integer::reagents(10,10,2),products(10,10,2) 
    integer::ngrmax,i,j,k,l,icomp,imprim,ifin2,jist,inp
    integer::npunt(NMG),ngrup(NMG),npint(NG),nintt(NG),num,numint
    integer::SolventsNumber, idev, idevr, InputMethod, lenfile, CantReact, CantProd
!...real*8
    real*8:: Temp, ChemConst, Conversion, conversionorg, delta
    real*8::xtemporal(10)
!...Characters
    character*1::opt
    character*8 FS(NMG)
    character*3 CAR
    character*35::nomb 
    character*37 nomco
!...Logical
    logical::SALIR,MolecularDesign
!...
    type(compound),pointer:: ptr_Comp
    type(FinalStructure),pointer::ptr_FMS, recorre_FMS
    type(CalculatedProperties),pointer::ListProperties   !Provisorio para llamar a incorporate_solvent
    type(CalculatedPerformance),pointer::ListPerformance !Provisorio para llamar a incorporate_solvent    
!...Commons
    COMMON /STAT/ MM,MJ,MK,MI,MH
    COMMON/US/MS
    common/NOM/FS
    common/PUNSUB/NPUNT,NGRUP,NUM
    common/PUNGRU/NPINT,NINTT,NUMINT
!-----PRUEBAS
    integer::kgrup(nsel),ngr1(10),ngr2(10),ngr3(10),ngr4(10)
    CHARACTER*21 family(8)
    real*8,external:: NEWTON

!...Sentencias
    idev = 6
	idevr = 5    
    ifin2 = 0
!   inicializaci�n de variables GRPS


    allocate(InputProblem01)    
    InputProblem01%mop = 4
    
!---Opciones de corrida (crear archivo o leer archivo)
    InputMethod=0
    do while(InputMethod==0)
        write (idev,"(' ','Options:',//,11x,'Create a new input file : 1',/,11x,'Run a created input file: 2',//,42x,'> ',$)") 
        read (idevr,"(a)") opt
        InputMethod = index ('12',opt)
    enddo    
!Selecci�n de modelo termodin�mico y apertura de bases de datos 
    call ab_ban1(model)
!---Selecci�n de tabla de par�metros UNIFAC
    call tabla_parametros(ipareq) 	


    if (InputMethod==1)then !Create a new input file : 1
        
        ! Pide el t�tulo del problema    
2110    write (idev,"('Give problem title (maximum 70 characters)')") 
        read (idevr,510,err=2110) titl
        InputProblem01%ProblemTitle = titl
        
        ! Pide el nombre del archivo .mdi        
3100    write (idev,"(/,' Give input file name (maximum 35 characters).')")
        read (idevr,510,err=3100) nomb
        InputProblem01%FileName = nomb     
510     format (a)
        
!-----Ingresar REACTIVOS
        i=1
        inp=0
        nullify(InputProblem01%Reagents)
	    do while (.False.)
	        nullify(ptr_Comp)
	        allocate(ptr_Comp)
	        nomco = 'Reagent                   '
            nomco(9:)= '#'//char(48+i)//':'
            call leer_comp (ipareq,nomco,ptr_Comp%Formula)         !se presentan y seleccionan grupos

        !Carga de subgrupos del REACTIVO en NPUNT y NPINT     
            k=1
            do while (ptr_Comp%Formula(k,1)/=0)
                call CR_PUNTF(ptr_Comp%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
                k=k+1
            enddo
        
            call incorporate_compound (ptr_Comp,InputProblem01%reagents)
  	    
  	        write(6,"('/,Another reagent?'/,' yes: <ret>.  not: 1.      ',$)")
            read (5,"(i2)") inp
            if(inp==1)exit
  	        i=i+1
  	    enddo
  	    CantReact = i
!-----Ingresar PRODUCTOS
        i=1
        inp=0
        nullify(InputProblem01%Products)
	    do while (.False.)
	        nullify(ptr_Comp)
	        allocate(ptr_Comp)
	        nomco = 'Product                   '
            nomco(9:)= '#'//char(48+i)//':'
            call leer_comp (ipareq,nomco,ptr_Comp%Formula)         !se presentan y seleccionan grupos

        !Carga de subgrupos del REACTIVO en NPUNT y NPINT     
            k=1
            do while (ptr_Comp%Formula(k,1)/=0)
                call CR_PUNTF(ptr_Comp%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
                k=k+1
            enddo
        
            call incorporate_compound (ptr_Comp,InputProblem01%Products)
  	    
  	        write(6,"('/,Another Product?'/,' yes: <ret>.  not: 1.      ',$)")
            read (5,"(i2)") inp
            if(inp==1)exit
  	        i=i+1
        enddo
        CantProd = i
        
       ! write(6,"(/,' Give reaction temperature (K)',19X,' > ',$)")
        !read (5,*)InputProblem01%T
        InputProblem01%T = 298.15
        call Constraints_Selection()
     !   write(6,"(' Give chemical equilibrium constant: ',$)")
     !   read(5,"(d)") InputProblem01%ChemEqConst
    !ESCRITURA INPUT FILE         
        lenfile = index(nomb,' ')-1    
        open(unit=1,file=nomb(1:lenfile)//'.mdi',form='formatted')
        !
        write(1,*) CantReact
        
        !
        close (unit=1)        
        
        
        
    else !si se lee desde archivos
    
    endif
!================================================================================================
!   PRUEBA SOLVENTES PARA REACCIONES
!================================================================================================
    
    
    
    MolecularDesign = .True.
    if(MolecularDesign)then
1200    write (6,590) !Give the family of solvents to be generated: ....
        read (5,*,err=1200) ifam
        if ((ifam.lt.1).or.(ifam.gt.4))go to 1200    

3050    write (6,1600) !now you have to choose......screen(1);screen and output
        read (5,"(i2)",err=3050) imprim   
        CALL SELECCION_GRUPOS (IMPRIM,KGRUP,IFIN2,NGR1,NGR2,NGR3,NGR4,MGR,family) !family
        call Store_Pr (ipareq)
        call Store_In ()
        DO 6 I=1,MDV	!MDV: n�mero de subgrupos intermedios
			NGSDV(I,1)=MM(NPUNT(NGDV(I)))
			NGSDV(I,2)=MJ(NPUNT(NGDV(I)))
			NGSDV(I,3)=MK(NPUNT(NGDV(I)))
			NGSDV(I,4)=MI(NPUNT(NGDV(I)))
			NGSDV(I,5)=MH(NPUNT(NGDV(I)))
6		CONTINUE
	    DO 7 I=1,MSV	!MSV n�mero de grupos terminales
			NGSSV(I,1)=MM(NPUNT(NGSV(I)))
			NGSSV(I,2)=MJ(NPUNT(NGSV(I)))
			NGSSV(I,3)=MK(NPUNT(NGSV(I)))
			NGSSV(I,4)=MI(NPUNT(NGSV(I)))
			NGSSV(I,5)=MH(NPUNT(NGSV(I)))
7       CONTINUE
        allocate(InputProblem01%MixtureInput)
      !  nullify(InputProblem01%MixtureInput%Solutes)    
        salir=.False.
        CALL STRUCTURE_GENERATOR (JIST,MS,3,NGSSV,NGSDV,SALIR,.FALSE.)    
 
        !composici�n mezcla. Temporal para probar llecalas
        recorre_FMS => FMSs
        !cargo el oleico 1  1  2 14  6  1 43  1
        ms(1,:,:) = reshape((/1,2,6,43,0,0,0,0,0,0,1,14,1,1,0,0,0,0,0,0/),(/10,2/))
        !butanol a la segunda posici�n
        ms(2,:,:) = reshape((/1,2,15,0,0,0,0,0,0,0,1,3,1,0,0,0,0,0,0,0/),(/10,2/))        
        !pongo el agua a la tercera posici�n
        ms(3,:,:) = reshape((/17,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0/),(/10,2/))
        !cargo el butyl oleate 1 2 2 17 6 1 23 1
        ms(4,:,:) = reshape((/1,2,6,23,0,0,0,0,0,0,2,17,1,1,0,0,0,0,0,0/),(/10,2/))    
    
        i=0
        do while(associated(recorre_FMS))
        

            !cargo el solvente
            ms(5,:,:) = recorre_FMS%Formula(:,:)
        
            xtemporal(1) = 0.0007
            xtemporal(2) = 0.0011
            xtemporal(3) = 0.9848
            xtemporal(4) = 0.0
            xtemporal(5) = 0.0135
        
            ChemConst = 78.0
            conversionorg=0.5
            conversion = conversionorg
            !InputProblem01%T=298.0
            InputProblem01%P = 1
            iter = 0
            do while(.True.)
                conversion = NEWTON (conversion,ChemConst,xtemporal,ms,InputProblem01%T,InputProblem01%P,5)
                if(iter == 0)then
                    if(conversion <= 1)then
                        exit
                    else
                        iter = iter+1
                        conversion = conversionorg + 0.1
                        cycle
                    endif
                elseif(conversion<1)then
                    exit
                elseif(conversion>1 .and. iter<8)then
                    continue
                elseif(iter==8)then
                    exit
                else !deber�a entrar cuando conversion == NaN
                    conversion = conversionorg + 0.1*iter
                endif
                iter = iter+1
                conversion = conversionorg + 0.1
            enddo
            
            recorre_FMS%conversion = conversion 
        !  call llecalas(InputProblem01%T,InputProblem01%P,ms,3,xtemporal,act)
        
            recorre_FMS=>recorre_FMS%next
            i=i+1
        enddo
        nsol=100
        call score(nsol)
        call write_results (mop,6,ifam,nsol)
    else !Para solventes ingresados por el usuario
        
            nullify(FMSs)
            nomco = 'Solvent                 '
            write(6,"(//,' How many solvents do you want to evaluate?:',6X,' > ',$)")
            read(5,*) SolventsNumber
            do i=1,SolventsNumber
                nullify(ptr_FMS)
                allocate(ptr_FMS)
                nomco(9:)= '#'//char(48+i)//':'
                call leer_comp (ipareq,nomco,ptr_FMS%Formula)
!                write(6,"(' ','Give the structure type of the solvent you have ',&
!     		    'entered:',//,&
!                10x,'-Aromatic                           : 1',//&
!                10x,'-Single substance groups            : 2',//&
!                10x,'-Aliphatic                          : 3',//&
!                10x,'-Cyclic                             : 4',///&
!     	        51x,' > ',$)")
!     	  
!     	        read(5,*)ptr_FMS%tsolv
!     	        tsolv=ptr_FMS%tsolv
       !!!!CHECKEAR PAR�METROS DE INTERACCI�N!!!!!

!Carga de subgrupos del SOLVENTE en NPUNT y NPINT     
                k=1
                do while (ptr_FMS%Formula(k,1)/=0)
                    call CR_PUNTF(ptr_FMS%Formula(k,1),NPUNT,NGRUP,NPINT,NINTT)
                    k=k+1
                enddo

                call Incorporate_Solvent (InputProblem01%mop,ptr_FMS%Formula,3,ListProperties,ListPerformance)
            enddo
        

    endif

    pause

!...Formats
601 format(' Give the group composition of the compound: ',&
     	   ' id1,ny1,id2,ny2,etc.',/,                      &
           10x,'id1,id2: subgroup identification number',/,&
           10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
970 format (' ',/,' ','If it is OK: <ret>.  If not: 1.      ',$)
 590  format (/,' Give the family of solvents to be generated:',///,&
             10x,'-Aromatic solvents                  : 1',//&
             10x,'-Single substance groups            : 2',//&
             10x,'-No-aromatics with up to 12 groups in',/&
             10x,' the final structure                : 3',//&
             10x,'-Cyclic solvents                    : 4',///&
     	      51x,' > ',$)
 1600  format (' ',/,' * Now you have to choose the groups that take ',&
             'part in the Molecular',/,&
                   '   Design of Solvents.  The program performs a ',&
             'preliminary selection',/,&
                   '   in  order to  eliminate  those groups which ',&
             'have  no  interaction',/,&
                   '   parameters  with  respect  to  the two main ',&
             'components.',/,&
                   '   If you want to  see  the  pairings  without ',&
             'interaction parameters',/,&
                   '   enter the following code:',/,&
              22x,'screen:                              1',/,&
              22x,'screen and output file (nointer):    2',/,&
              22x,'no printing:                      <ret>',/,&
              63x,'> ',$)
    endsubroutine
    
    subroutine SeleccionarCompuesto(compuestos,n,ms)
!-------------------------------------------------
!   Descripci�n
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    parameter(NMG=150)
    implicit real*8(A-H,O-Z)
!...VARIABLES DE ENTRADA
    integer::compuestos(n,10,2),n
!...VARIABLES DE SALIDA
    integer,intent(out)::ms(10,2)
!...VARIABLES INTERNAS
!   Integers
    integer::Sel
!...Commons
    character*8 FS(NMG)
    common/NOM/FS
!...Sentencias
    do i=1, n
        write(6,"(/,I2,A1,$)") i,"-"
        j=1
        do while (compuestos(i,j,1)/=0)
            write(6,"(A9,I2,$)") FS(compuestos(I,J,1)), compuestos(I,J,2)
            j=j+1
        enddo
        !write(6,"(/)")
    enddo
    write (6,"(//,'Select the compound. > ',$)")
    read  (5,"(i2)") Sel
    write (6,"(/,'Compound selected:',$)")
    j=1
    do while (compuestos(Sel,j,1)/=0)
        write(6,"(A9,I2,$)") FS(compuestos(Sel,J,1)), compuestos(Sel,J,2)
        j=j+1
    enddo
    write(*,*)''
    j=1
    do while (compuestos(Sel,j,1)/=0)
        ms(j,1)=compuestos(Sel,j,1)
        ms(j,2)=compuestos(Sel,j,2)
        j=j+1
    enddo
    call pausa
    endsubroutine
!
!
