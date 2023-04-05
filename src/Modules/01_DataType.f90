module constantes
    save
    integer,parameter::NCOM=10
    integer,parameter::NMG = 150 !(dimensión máxima para vectores que guardan inf. sobre subgrupos)
    integer,parameter::NINT = 70 !Dimensión para vectores que guardan inf. sobre parámetros de interacción
    integer,parameter::DiffStructGroups = 10 !Cantidad máxima de grupos distintos en una estructura
    integer,parameter::NSEL = 150
    integer,parameter::COMPUESTOSxPAGINA_CONSOLA = 6, COMPUESTOSxPAGINA_ARCHIVO = 100
endmodule

module Boundaries
    integer::lugar(100,2)
    type bound
        character(len=8)::ID
        character(len=56)::cartel
        character(len=7)::unit
        logical::ExistLowerBound,ExistUpperBound,ExistRunValue
        real*8:: LowerBound,Upperbound,RunValue
    endtype bound
endmodule Boundaries


module SubGrupos
    type groups
        character*8::Name
        integer::Number
        character*8::NameMainGroup
        integer::NumberMainrGroup
        integer::attM,attJ,attK,attI,attH
        character*3::Charact
        real*8::R,Q,MW,dTC,dPC,dV,dVi,dN,dTCa
        real*8::Tstr,gstr,gdot,gddot
        integer::AssoSites(4)
        real*8::IntPar(70)
        type(Groups),pointer::next
    endtype groups
    type(Groups),pointer::LSubGroups
  
  contains
    subroutine Create_Group(nuevo)
        implicit none
        type(Groups),pointer,intent(out)::nuevo
        allocate(nuevo)
        nullify(nuevo%next)
    endsubroutine Create_Group   
    
    
!===================================================== 
    subroutine Incorporate_Group (puntero,nuevo)
!-----------------------------------------------------
!   Esta subroutine incorpora una estructura
!   aceptada junto con todas sus propiedades estimadas
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-----------------------------------------------------
 !   use constantes
    implicit none
!Variables de ENTRADA
    type(Groups),pointer,intent(in)::nuevo
!Variables de ENTRADA/SALIDA
    type(Groups),pointer,intent(inout)::puntero   
!Variables INTERNAS
    type(Groups),pointer::recorre
!Sentencias
    if (.not.associated(puntero)) then 
        allocate(puntero)
        !puntero%Formula=compound
        puntero=>nuevo
    else
        recorre=>puntero
        do while (associated(recorre%next))
            recorre=>recorre%next
        enddo
        recorre%next=>nuevo
    end if
      ! RETURN
  endsubroutine Incorporate_Group
!=======================================================
  subroutine Search_Group (Number,puntero,group)
!-------------------------------------------------------
!   Se ingresa la lista que contiene todos los grupos en 
!   la variable "Puntero". Devuelve el vector "group" 
!   que apunta al nodo donde se aloja el grupo "Number".
!-------------------------------------------------------
 !   use constantes
    implicit none
!Variables de ENTRADA
    integer,intent(in)::number
    type(Groups),pointer,intent(in)::puntero
!Variables de SALIDA
    type(Groups),pointer,intent(out)::group
!Variables INTERNAS
    type(Groups),pointer::recorre
!Sentencias
    recorre=>puntero
    do while (.True.)
        if (recorre%Number == Number)then
            group=>recorre
            exit
        endif
        recorre=>recorre%next
    enddo
  endsubroutine Search_Group
!===========================================================
  character*8 function Obtain_SubGroup_Name(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_SubGroup_Name = group%Name
  endfunction Obtain_SubGroup_Name
  
!===========================================================
  character*8 function Obtain_MainGroup_Name(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_MainGroup_Name = group%NameMainGroup
  endfunction Obtain_MainGroup_Name
!===========================================================
  integer function Obtain_MainGroup_Number(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_MainGroup_Number = group%NumberMainrGroup
  endfunction Obtain_MainGroup_Number

!===========================================================
  integer function Obtain_attM(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_attM = group%attM
  endfunction Obtain_attM
  
!===========================================================
  integer function Obtain_attJ(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_attJ = group%attJ
  endfunction Obtain_attJ

!===========================================================
  integer function Obtain_attK(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_attK = group%attK
  endfunction Obtain_attK 

!===========================================================
  integer function Obtain_attI(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_attI = group%attI
  endfunction Obtain_attI

!===========================================================
  integer function Obtain_attH(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_attH = group%attH
  endfunction Obtain_attH
!===========================================================
  real*8 function Obtain_Q(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_Q = group%Q
  endfunction Obtain_Q
!===========================================================
  real*8 function Obtain_R(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_R = group%R
  endfunction Obtain_R
!===========================================================
  real*8 function Obtain_MW(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_MW = group%MW
  endfunction Obtain_MW
  !===========================================================
  real*8 function Obtain_dV(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_dV = group%dV
  endfunction Obtain_dV
  !===========================================================
  real*8 function Obtain_dPC(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_dPC = group%dPC
  endfunction Obtain_dPC
  
!===========================================================
  real*8 function Obtain_dVi(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_dVi = group%dVi
  endfunction Obtain_dVi
  
!===========================================================
  real*8 function Obtain_dN(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_dN = group%dN
  endfunction Obtain_dN
    
  !===========================================================
  real*8 function Obtain_dTC(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_dTC = group%dTC
  endfunction Obtain_dTC
  !===========================================================
  real*8 function Obtain_dTCa(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_dTCa = group%dTCa
  endfunction Obtain_dTCa

!===========================================================
  real*8 function Obtain_Tstr(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_Tstr = group%Tstr
  endfunction Obtain_Tstr 
    
!===========================================================
  real*8 function Obtain_gstr(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_gstr = group%gstr
  endfunction Obtain_gstr 
  
!===========================================================
  real*8 function Obtain_gdot(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_gdot = group%gdot
  endfunction Obtain_gdot
  
!===========================================================
  real*8 function Obtain_gddot(Number)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    integer,intent(in)::Number
!Variables INTERNAS
    type(Groups),pointer::group
!SENTENCIAS
    call Search_Group (Number,LSubGroups,group)
    Obtain_gddot = group%gddot
  endfunction Obtain_gddot  

!===========================================================
  integer function Size_LSubGroups()
!-----------------------------------------------------------
!   Devuelve la cantidad de grupos (nodos) en LSubGroups
!-----------------------------------------------------------  
    implicit none
!Variables INTERNAS
    integer::i
!SENTENCIAS
    type(Groups),pointer::recorre
    i=0
    recorre => LSubGroups
    do while (associated(recorre))
        i = i + 1
        recorre => recorre%next
    enddo
    nullify(recorre)
    Size_LSubGroups = i

  endfunction Size_LSubGroups


endmodule SubGrupos

!****************************************************
module StructuresDesign
    use constantes
    type FinalStructure
        integer::Formula(DiffStructGroups,2)
        real*8::MolecularWeight
        real*8::BoilingPoint
        real*8::CriticalPressure
        real*8::CriticalTemperature
        real*8::CriticalVolume
        real*8::Density
        real*8::VapourPressure
        real*8::Hvap
        real*8::DifferenceBP 
        real*8::RDER
        real*8::Viscosity
        real*8::Selectivity
        real*8::SolventPower
        real*8::RelVolatility
        real*8::SolventLost 
        real*8::DistCoefficient
        real*8::Kow
        real*8::conversion !conversión de reacción
        integer::tsolv !Tipo de estructura molecular
        integer::GroupsNumber        
        integer::FunctionalGroups
        integer::position
        integer::NumberIsomers
        integer::SolutesNumber
        type(CalculatedPerformance),pointer::Performance
        type(FinalStructure),pointer::next
        type(isomers),pointer::nextisomer
    endtype FinalStructure
    !Creación del puntero FMSs
    type(FinalStructure),pointer::FMSs
    
    type isomers
        integer::index !número identificatorio
        character(len=35)::name !nombre del compuesto
        character(len=35)::FormChem !fórmula química
        integer::Formula(DiffStructGroups,2) !Configuración UNIFAC
        real*8::BoilingPoint
        real*8::LiquidMolarVolume
        type(isomers),pointer::next
    endtype
    type CalculatedProperties
        integer::Formula(DiffStructGroups,2)
        real*8::MolecularWeight
        real*8::BoilingPoint
        real*8::CriticalPressure
        real*8::CriticalTemperature
        real*8::CriticalVolume
        real*8::Density
        real*8::VapourPressure        
        real*8::Hvap
        real*8::RDER
        real*8::Viscosity
        real*8::Kow
        integer::GroupsNumber        
        integer::FunctionalGroups
        integer::position
        integer::NumberIsomers
    endtype CalculatedProperties
    type CalculatedPerformance
        integer::SoluteFormula(DiffStructGroups,2)
        real*8::Selectivity
        real*8::SolventPower
        real*8::RelVolatility
        real*8::SolventLost 
        real*8::DistCoefficient
        real*8::DifferenceBP
        type(CalculatedPerformance),pointer::next        
    endtype CalculatedPerformance
  
  contains
!==========================================================================
  subroutine Incorporate_Structure (compound,tsolv,puntero,ListProperties,ListPerformance)
!-------------------------------------------------
!   Esta subroutine incorpora una estructura
!   aceptada junto con todas sus propiedades estimadas
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    use constantes
    implicit none
!Variables de ENTRADA
    integer,intent(in)::tsolv
    integer,dimension(DiffStructGroups,2),intent(in)::compound
    type(CalculatedPerformance),pointer,intent(in)::ListPerformance
    type(CalculatedProperties),pointer,intent(in)::ListProperties
!Variables de ENTRADA/SALIDA
    type(FinalStructure),pointer,intent(inout)::puntero    
!Variables INTERNAS
    type(FinalStructure),pointer::interno
    type(CalculatedPerformance),pointer::recorrePerformance
    real*8::SumSel,SumSolPow,SumRelVol,SumSOlvLost,SumDistCoeff
    integer::i
!Sentencias
    nullify(interno)
    if (.not.associated(puntero)) then 
        allocate(puntero)
        puntero%Formula=compound
        nullify(puntero%next)
    else
        allocate(interno)
        interno%Formula=compound
        interno%next=>puntero
        puntero=>interno
        nullify(interno)
    end if
    puntero%tsolv = tsolv
    puntero%Performance => ListPerformance
    puntero%MolecularWeight = ListProperties%MolecularWeight
    puntero%BoilingPoint = ListProperties%BoilingPoint
    puntero%Hvap = ListProperties%Hvap
    puntero%RDER = ListProperties%RDER !DensT/DENS2 !!!!!!!!!!!!!!! ver qué pasa con esto
	puntero%Viscosity = ListProperties%Viscosity
	puntero%GroupsNumber = ListProperties%GroupsNumber
	puntero%FunctionalGroups = ListProperties%FunctionalGroups  
	puntero%Kow = ListProperties%Kow
    puntero%CriticalTemperature = ListProperties%CriticalTemperature
    puntero%CriticalVolume = ListProperties%CriticalVOlume
    puntero%CriticalPressure = ListProperties%CriticalPressure
    recorrePerformance => ListPerformance
    i = 0
    SumSel = 0
    SumSolPow = 0
    SumRelVol = 0
    SumSOlvLost = 0
    SumDistCoeff = 0
    do while(associated(recorrePerformance))
        SumSel = SumSel + recorrePerformance%Selectivity
        SumSolPow = SumSolPow + recorrePerformance%SolventPower
        SumRelVol = SumRelVol + recorrePerformance%RelVolatility
        SumSolvLost = SumSolvLost + recorrePerformance%SolventLost
        SumDistCoeff = SumDistCoeff + recorrePerformance%DistCoefficient
        i = i+1
        recorrePerformance => recorrePerformance%next     
    enddo
    puntero%Selectivity = SumSel/i
    puntero%SolventPower = SumSOlPow/i
    puntero%RelVolatility = SumRelVol/i
    puntero%SolventLost = SumSolvLost/i
    puntero%DistCoefficient = SumDistCoeff/i
    puntero%SolutesNumber = i

  endsubroutine Incorporate_Structure   

!==========================================================================
  subroutine Incorporate_Performance (puntero,ListPerformance)
!-------------------------------------------------
!   Esta subroutine incorpora una estructura
!   aceptada junto con todas sus propiedades estimadas
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    use constantes
    implicit none
!Variables de ENTRADA/SALIDA
    type(CalculatedPerformance),pointer,intent(inout)::puntero 
    type(CalculatedPerformance),pointer,intent(inout)::ListPerformance
!Variables INTERNAS
    type(CalculatedPerformance),pointer::interno
!Sentencias
    allocate(interno)
    nullify(interno)
    if (.not.associated(ListPerformance)) then 
        allocate(ListPerformance)
        ListPerformance=>puntero
        nullify(puntero%next)
    else
        allocate(interno)
        interno=>ListPerformance
        do while(associated(interno%next))
            interno=>interno%next
        enddo
        interno%next=>puntero
        nullify(puntero%next)
    end if
      ! RETURN
  endsubroutine Incorporate_Performance
  
!===========================================================
  integer function Size_LFMSs(FMSs)
!-----------------------------------------------------------
!   Devuelve el nombre del subgrupo Number
!-----------------------------------------------------------  
    implicit none
!Variables de ENTRADA
    type(FinalStructure),pointer,intent(in)::FMSs
!Variables INTERNAS
    integer::i
!SENTENCIAS
    type(FinalStructure),pointer::recorre
    i=0
    recorre => FMSs
    do while (associated(recorre))
        i = i + 1
        recorre => recorre%next
    enddo
    nullify(recorre)
    Size_LFMSs = i

  endfunction Size_LFMSs  
endmodule StructuresDesign

!***************************************************************
module Input
    use CONSTANTES
    integer::contador
    integer::model
    integer::ipareq
    type DataProblem01
        real*8::T,P
        real*8::ChemEqConst 
        integer::kilout
        integer::MainGroups(NMG) !reemplaza al vector MGR
        integer::SubGroups(NMG)
        integer::mop !separation operation
        integer::ifam !familia de solventes a generar
        character*70::ProblemTitle
        character*35::FileName
        logical::inv ! verdadero para MOLDES invertido
        type(Mixture),pointer::MixtureInput
        type(compound),pointer::Reagents
        type(compound),pointer::Products          
    endtype DataProblem01
    type(DataProblem01),pointer::InputProblem01
        
    !type DataProblem02
    !    real*8::T,P
    !    integer::MainGroups(NMG) !reemplaza al vector MGR
    !    integer::SubGroups(NMG)        
    !    integer::ifam !familia de solventes a generar
    !    integer::mop
    !    real*8::ChemEqConst        
    !    character*70::ProblemTitle
    !    character*35::FileName
    !    type(compound),pointer::Reagents
    !    type(compound),pointer::Products        
    !endtype DataProblem02
    !type(DataProblem02),pointer::InputProblem01
    
    type Mixture
        integer::SolutesNumber
        integer::MainGroups(NMG)
        integer::SubGroups(NMG)
        real*8::Tazeo, Xazeo
        type(Compound),pointer::Solutes
        type(Compound),pointer::PCR
        type(Mixture),pointer::next
    endtype Mixture
    
    type Compound
        character*35::Name    
        integer::Formula(DiffStructGroups,2)
        integer::GroupsNumber
        real*8::MW
        real*8::TC
        real*8::PC
        real*8::VC
        real*8::BoilingPoint
        real*8::vliq
        real*8::dens
        real*8::a(3)
        type(Compound),pointer::next
    endtype Compound
    
contains    
!---Subroutines    
  subroutine incorporate_compound (ptr_pcp,puntero)
!-------------------------------------------------
!   Esta subroutine incorpora una estructura
!   aceptada junto con todas sus propiedades estimadas
!   - Variables de entrada
!       vars:
!   - Variables de salida
!       vars:
!-------------------------------------------------
    use constantes
    implicit none
!Variables de ENTRADA/SALIDA
    type(Compound),pointer,intent(inout)::ptr_pcp
    type(Compound),pointer,intent(inout)::puntero    
!Variables INTERNAS
    type(Compound),pointer::recorre
    integer::i
!Sentencias
    nullify(recorre)
    if (.not.associated(puntero))then 
        allocate(puntero)
        i=1
        do while(ptr_pcp%Formula(i,1)/=0) 
            i=i+1
        enddo
        ptr_pcp%GroupsNumber = i-1        
        puntero => ptr_pcp
        
        nullify(puntero%next)
    else
        allocate(recorre)
        i=1
        do while(ptr_pcp%Formula(i,1)/=0) 
            i=i+1
        enddo
        ptr_pcp%GroupsNumber = i-1
        recorre => puntero
        do while(associated(recorre%next))
            recorre=>recorre%next
        enddo
        
        recorre%next => ptr_pcp
        nullify(ptr_pcp%next)
        nullify(recorre)
    end if
      ! RETURN
  endsubroutine incorporate_compound 
  




    
endmodule Input
module blockdatas
      parameter (nmodel=3)
      implicit real*8 (a-h,o-z)
      dimension idato(3),dato2(5),dato3(14)
      INTEGER iach(nmodel),iacch2(nmodel),iaccl(nmodel),&
             iacnh2(nmodel),iacno2(nmodel)
      character*63 cartel(8),concar(2),carte2(6),cartelmop3(4)
      character*64 icartl(3),icart2(3)
      character*17 tabla(nmodel)
      character*21 family(8)
      CHARACTER*36 cartel3(7)

      data family /'      aromatic groups','     molecular groups',&
                 3*'  intermediate groups','single valence groups',&
     			' three valence groups','  four valence groups'/
      data tabla /'    liquid-liquid','     liquid-vapor',&
     			' ifinite dilution'/
      data cartel&
             /'1 Seli: Minimum selectivity of an intermediate structure &
       (wt) ','2 Sele: Minimum selectivity of a final structure        &
       (wt) ','3 Solv: Minimum solvent power required                  &
     (wt %) ','4 Slsi: Solvent loss of an intermediate structure       &
     (wt %) ','5 SLost: Solvent loss of a final structure               &
     (wt %) ','6 Dist: Minimum distribution coefficient required       &
       (wt) ','7 Pmma: Maximum molecular weight of the final structure &
            ','8 Den2: Density of raffinate                           (&
     gr/ml) '/
      data carte2 &
             /'1 Voli: Minimum volatility obtained with an intermediate &
      struct','2 Vole: Minimum volatility obtained with a final structu &
     re     ','3 Solv: Minimum solvent power required                   &
     (wt %) ','4 Pmma: Maximum molecular weight of the final structure  &
            ','5 Dtcs: Minimun boiling temperature difference between   &
     solvent','        and less volatile component                     &
        (K) '/
      data cartel3&
             /'1 Pmma: Molecular weight            '&
             ,'2 Bpnt: Boiling point            [K]'&
             ,'3 Tcrt: Critical temperture      [K]'&
             ,'4 Pcrt: Critical pressure      [atm]'&
             ,'5 Vcrt: Critical volume    [cm3/mol]'&
             ,'6 Visc: Viscosity            [mPa*s]'&
             ,'7 Dens: Density                     '/ 
      data cartelmop3&
             /'1 Seli: Minimum selectivity of an intermediate structure &
       (wt) ','2 Sele: Minimum selectivity of a final structure        &
     (wt %) ','3 Dist: Minimum distribution coefficient required       &
       (wt) ','4 Pmma: Maximum molecular weight of the final structure &
            '/
       data icartl                                           &
     /'9 Nsol: Maximum number of solvents to be listed       &           
     ','10  Is: Number of solvents per page                  &         
      ','11 out: Type of output                              &         
       '/
       data icart2                                          &
     /'6 Nsol: Maximum number of solvents to be listed      &          
     ','7   Is: Number of solvents per page                 &          
      ','8  out: Type of output                             &          
       '/
       data concar /'        0 : Screen and output file     &             
                ','        1 : Output file                  &          
                '/
      data iach /9,10,10/
      data iacch2 /12,13,13/
      data iaccl /40,54,50/
      data iacnh2 /43,37,34/
      data iacno2 /47,58,55/
      data idato /100,21,0/
!      data dato /1.,8.,5.,0.5,0.05,0.2,240.,1./
!      data dato /1.,4.,1.,25.,10.0,0.2,240.,1./
      !data dato /.25,1.,.25,150.,150.0,0.02,440.,1./
      data dato2 /1.,1.2,15.,240,30./
      data dato3 /0.,9000.,0.,9000.,0.,9000.,0.,9000.,0.,9000.,0.,9000.,0.,9000./
endmodule
 
BLOCK DATA BANCO
       PARAMETER (NMOD=3)
       implicit real*8(A-H,O-Z)
       COMMON/GRUPES/ICH(NMOD),IACOH(NMOD),IACCH(NMOD),IAC(NMOD),IOH(NMOD),ICOOH(NMOD),IACCOO(NMOD)
       COMMON/GRUPAR/IACH(NMOD),IACCH2(NMOD),IACCL(NMOD),IACNH2(NMOD),IACNO2(NMOD)
       COMMON/GRUPAL/ICH3(NMOD),ICH2(NMOD)
!C
!C---DATAS PARA GRUPOS QUE SE UTILIZAN EN EL CUERPO DEL PROGRAMA EN
!C---FORMA ESPECIAL. LOS 3 NUMEROS DE CADA DATA CORRESPONDEN AL NUMERO
!C---DE SUBGRUPO EN CADA TABLA DE PARAMETROS (LIQ-LIQ,LIQ-VAP,DIL.INF.)
	  DATA ICH2/2,2,2/
        DATA ICH/3,3,3/
        DATA IACOH/18,18,18/
        DATA IACCH2/12,13,13/
        DATA IACCH/13,14,14/
        DATA IAC/10,11,11/
        DATA IACH/9,10,10/
        DATA IACCL/40,54,50/
        DATA IACNH2/43,37,34/
        DATA IACNO2/47,58,55/
        DATA ICH3/1,1,1/
        DATA IOH/14,15,15/
        DATA ICOOH/23,43,0/
        DATA IACCOO/0,0,65/
END