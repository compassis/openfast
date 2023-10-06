MODULE SeaFEM
    
   USE ISO_C_BINDING    
   USE SeaFEM_Types  
   USE NWTC_Library
   
   IMPLICIT NONE
   
   !DEC$ ATTRIBUTES C, EXTERN, DLLIMPORT :: Fast_waves_global
   TYPE(C_FUNPTR) :: proc
   INTEGER(C_INTPTR_T) :: module_handle
   PROCEDURE(EXCHANGE_FAST_DATA), pointer :: EXCHANGE_DATA
   PROCEDURE(RUNNING_FAST_UPDATE), pointer :: UPDATE_SEAFEM
   PROCEDURE(END_FAST_COUPLING), pointer :: END_TIMELOOP

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: SeaFEM_Ver = ProgDesc( 'SeaFEM', 'v1.00.00', '22-September-2023' )
   
   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SeaFEM_Init                           ! Initialization routine
   
   PUBLIC :: SeaFEM_CalcOutput                     ! Routine for computing outputs 
   
   ABSTRACT INTERFACE
        SUBROUTINE EXCHANGE_FAST_DATA(q,qdot,qdotdot,SeaFEM_Return_Forces,flag) BIND(C)
        USE ISO_C_BINDING
        IMPLICIT NONE
        REAL(C_FLOAT), INTENT(OUT), DIMENSION(*) :: q,qdot,qdotdot
        REAL(C_FLOAT), INTENT(IN), DIMENSION(*)  :: SeaFEM_Return_Forces
        INTEGER(C_INT), INTENT(OUT)              :: flag
        END SUBROUTINE EXCHANGE_FAST_DATA
    END INTERFACE
   
   ABSTRACT INTERFACE
        SUBROUTINE RUNNING_FAST_UPDATE() BIND(C)
        USE ISO_C_BINDING
        IMPLICIT NONE
        END SUBROUTINE RUNNING_FAST_UPDATE
   END INTERFACE
   
   ABSTRACT INTERFACE
        SUBROUTINE END_FAST_COUPLING() BIND(C)
        USE ISO_C_BINDING
        IMPLICIT NONE
        END SUBROUTINE END_FAST_COUPLING
    END INTERFACE   
   
    INTERFACE 
        FUNCTION LoadLibrary(lpFileName) BIND(C,NAME='LoadLibraryA')
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INTPTR_T, C_CHAR
            IMPLICIT NONE 
            CHARACTER(KIND=C_CHAR) :: lpFileName(*) 
            !GCC$ ATTRIBUTES STDCALL :: LoadLibrary 
            !DEC$ ATTRIBUTES STDCALL :: LoadLibrary
            INTEGER(C_INTPTR_T) :: LoadLibrary 
        END FUNCTION LoadLibrary 
    
        FUNCTION GetProcAddress(hModule, lpProcName) BIND(C, NAME='GetProcAddress')
            USE, INTRINSIC :: ISO_C_BINDING, ONLY:  & 
                C_FUNPTR, C_INTPTR_T, C_CHAR
            IMPLICIT NONE
            !GCC$ ATTRIBUTES STDCALL :: GetProcAddress
            !DEC$ ATTRIBUTES STDCALL :: GetProcAddress
            TYPE(C_FUNPTR) :: GetProcAddress
            INTEGER(C_INTPTR_T), VALUE :: hModule
            CHARACTER(KIND=C_CHAR) :: lpProcName(*)
        END FUNCTION GetProcAddress  
   END INTERFACE
   
   CONTAINS
   
   SUBROUTINE SeaFEM_Init( InitInp, u, p, OtherState, y )
        ! This routine is called at the start of the simulation to perform initialization steps.
        ! The parameters are set here and not changed during the simulation.
        ! The initial states and initial guess for the input are defined.
        !..................................................................................................................................
   
        TYPE(SeaFEM_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
        TYPE(SeaFEM_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
        TYPE(SeaFEM_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
        TYPE(SeaFEM_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
        TYPE(SeaFEM_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
        
        ! local variables                                 
        INTEGER(IntKi)                                    :: ErrStat2                            ! local error status
        CHARACTER(1024)                                   :: ErrMsg2                             ! local error message
        
        ! Initialize the NWTC Subroutine Library
        CALL NWTC_Init( )
        
        ! Display the module information
        CALL DispNVD( SeaFEM_Ver )     
        
        ! Define initial system states here:
        OtherState%T               = 0.0
        OtherState%Out_Flag        = 0
      
        ! Create the input mesh to store platform motions
        CALL MeshCreate( BlankMesh         = u%SeaFEMMesh      &
                        ,IOS               = COMPONENT_INPUT   &
                        ,Nnodes            = 1                 &
                        ,ErrStat           = ErrStat2          &
                        ,ErrMess           = ErrMsg2           &
                        ,TranslationDisp   = .TRUE.            &
                        ,Orientation       = .TRUE.            &
                        ,TranslationVel    = .TRUE.            &
                        ,RotationVel       = .TRUE.            &
                        ,TranslationAcc    = .TRUE.            &
                        ,RotationAcc       = .TRUE.            )
         
        ! Create the node on the mesh
        CALL MeshPositionNode (u%SeaFEMMesh, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/), ErrStat2, ErrMsg2)
     
        ! Create the mesh element
        CALL MeshConstructElement (u%SeaFEMMesh, ELEMENT_POINT, ErrStat2, ErrMsg2, 1)

        ! Commits the mesh
        CALL MeshCommit (u%SeaFEMMesh, ErrStat2, ErrMsg2) 

        ! Creates an output mesh to store SeaFEM loads
        CALL MeshCopy ( SrcMesh   = u%SeaFEMMesh        &
                       ,DestMesh  = y%SeaFEMMesh        &
                       ,CtrlCode  = MESH_SIBLING        &
                       ,IOS       = COMPONENT_OUTPUT    &
                       ,ErrStat   = ErrStat2            &
                       ,ErrMess   = ErrMsg2             &
                       ,Force     = .TRUE.              &
                       ,Moment    = .TRUE.              )       
      
   END SUBROUTINE SeaFEM_Init

   SUBROUTINE SeaFEM_CalcOutput( t, u, p, OtherState, y, ErrStat, ErrMsg )
        ! Routine for computing outputs, used in both loose and tight coupling.
        !..................................................................................................................................
   
        REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
        TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
        TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
        TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
        TYPE(SeaFEM_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
        INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
        CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
        
        ! Local variables
        REAL(ReKi)                           :: q(6), qdot(6), qdotdot(6)    ! Platform motions
        REAL(ReKi)                           :: rotdisp(3)                   ! Small angle rotational displacements
        REAL(ReKi)                           :: SeaFEM_Return_Forces(6)      ! SeaFEM loads
        INTEGER(IntKi)                       :: I 
        
        ! Load exported procedures from SeaFEM
        module_handle=LoadLibrary(C_NULL_CHAR)
        proc=GetProcAddress(module_handle,"Exchange_Fast_Data"C)
        CALL C_F_PROCPOINTER(proc,EXCHANGE_DATA)
        proc=GetProcAddress(module_handle,"Running_Fast_Update"C)
        CALL C_F_PROCPOINTER(proc,UPDATE_SEAFEM)
         proc=GetProcAddress(module_handle,"End_Fast_Coupling"C)
        CALL C_F_PROCPOINTER(proc,END_TIMELOOP)
        
        ! Determine the rotational angles from the direction-cosine matrix
        rotdisp = GetSmllRotAngs( u%SeaFEMMesh%Orientation(:,:,1), ErrStat, ErrMsg )              
              
        ! Displacements, velocities and accelerations are obteined from the input mesh (12 iterations for time step increment)
        q       = reshape((/real(u%SeaFEMMesh%TranslationDisp(:,1),ReKi),rotdisp(:)/),(/6/))
        qdot    = reshape((/u%SeaFEMMesh%TranslationVel(:,1),u%SeaFEMMesh%RotationVel(:,1)/),(/6/))
        qdotdot = reshape((/u%SeaFEMMesh%TranslationAcc(:,1),u%SeaFEMMesh%RotationAcc(:,1)/),(/6/))   
        
        ! Updates SeaFEMs time step
        IF(OtherState%T==t)THEN
        !     WRITE(*,*) "Simulation time = ",t
        ELSE
            CALL UPDATE_SEAFEM() 
            ! WRITE(*,*) "Simulation time = ",t
            OtherState%T=t
        END IF
        
        ! Data exchange between SeaFEM and OpenFAST (motions sent and loads received) 
        CALL EXCHANGE_DATA(q,qdot,qdotdot,SeaFEM_Return_Forces,OtherState%flag_SeaFEM)
        
        ! Ends SeaFEM computation
        IF (t>=p%TMax) THEN
            IF(OtherState%Out_flag==(1+2*p%Iterations))THEN
                CALL END_TIMELOOP() 
            ELSE
                OtherState%Out_flag=OtherState%Out_Flag+1
            END IF
        END IF       
        
        ! SeaFEM loads are stored in the output mesh
        DO I=1,3
            y%SeaFEMMesh%Force(I,1)=SeaFEM_Return_Forces(I)
            ! WRITE(*,'(A,I1,A,E)') "Returned Forces Value SF[",I,"] = ",SeaFEM_Return_Forces(I)
        END DO
        DO I=1,3
            y%SeaFEMMesh%Moment(I,1)=SeaFEM_Return_Forces(I+3)
            ! WRITE(*,'(A,I1,A,E)') "Returned Forces Value SF[",I+3,"] = ",SeaFEM_Return_Forces(I+3)
        END DO
                
   END SUBROUTINE SeaFEM_CalcOutput   
   
END MODULE SeaFEM