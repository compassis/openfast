MODULE SeaFEM
    
   USE ISO_C_BINDING    
   USE SeaFEM_Types  
   USE NWTC_Library
   
   IMPLICIT NONE
   
   !DEC$ ATTRIBUTES C, EXTERN, DLLIMPORT :: Fast_waves_global
   !TYPE(C_PTR) :: Fast_waves_global 
   TYPE(C_FUNPTR) :: proc
   PROCEDURE(EXCHANGE_FAST_DATA), pointer :: EXCHANGE_DATA
   PROCEDURE(RUNNING_FAST_UPDATE), pointer :: UPDATE_SEAFEM
   PROCEDURE(END_SF_TIMELOOP), pointer :: END_TIMELOOP
   INTEGER(C_INTPTR_T) :: module_handle
   TYPE(C_PTR) :: linux_handle=C_NULL_PTR
   INTEGER, PARAMETER, PUBLIC :: RTLD_LAZY=1, RTLD_NOW=2, RTLD_GLOBAL=256, RTLD_LOCAL=0
   CHARACTER(C_CHAR), DIMENSION(1), SAVE, TARGET, PRIVATE :: dummy_string="?"
   
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: SeaFEM_Ver = ProgDesc( 'SeaFEM', 'v1.00.00', '01-July-2021' )
   
   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SeaFEM_Init                           ! Initialization routine
   
   PUBLIC :: SeaFEM_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: SeaFEM_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   
   PUBLIC :: SeaFEM_CalcOutput                     ! Routine for computing outputs   
   
   PUBLIC :: SeaFEM_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   
   PUBLIC :: SeaFEM_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   
   PUBLIC :: SeaFEM_UpdateDiscState                ! Tight coupling routine for updating discrete states
   
   PUBLIC :: EXCHANGE_FAST_DATA
    
   PUBLIC :: RUNNING_FAST_UPDATE
   
   PUBLIC :: DLOpen, DLSym
   
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
        SUBROUTINE END_SF_TIMELOOP() BIND(C)
        USE ISO_C_BINDING
        IMPLICIT NONE
        END SUBROUTINE END_SF_TIMELOOP
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
         
         FUNCTION GetModuleHandle(lpModuleName) BIND(C, NAME='GetModuleHandle')
            USE, INTRINSIC :: ISO_C_BINDING, ONLY:  & 
                C_FUNPTR, C_CHAR
            IMPLICIT NONE
            !GCC$ ATTRIBUTES STDCALL :: GetModuleHandle
            !DEC$ ATTRIBUTES STDCALL :: GetModuleHandle
            TYPE(C_FUNPTR) :: GetModuleHandle
            CHARACTER(KIND=C_CHAR) :: lpModuleName(*)
         END FUNCTION GetModuleHandle 
    END INTERFACE
    
   INTERFACE ! All we need is interfaces for the prototypes in <dlfcn.h>
      FUNCTION DLOpen(file,mode) RESULT(handle) BIND(C,NAME="dlopen")
         USE ISO_C_BINDING
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: file
         INTEGER(C_INT), VALUE :: mode
         TYPE(C_PTR) :: handle
      END FUNCTION
      
      FUNCTION DLSym(handle,name) RESULT(funptr) BIND(C,NAME="dlsym")
         USE ISO_C_BINDING
         TYPE(C_PTR), VALUE :: handle
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: name
         TYPE(C_FUNPTR) :: funptr ! A function pointer
      END FUNCTION
   
      FUNCTION DLError() RESULT(error) BIND(C,NAME="dlerror")
         USE ISO_C_BINDING
         TYPE(C_PTR) :: error
      END FUNCTION         
   END INTERFACE
   
    CONTAINS
    
    FUNCTION C_F_STRING(CPTR) RESULT(FPTR)
      ! Convert a null-terminated C string into a Fortran character array pointer
      TYPE(C_PTR), INTENT(IN) :: CPTR ! The C address
      CHARACTER(KIND=C_CHAR), DIMENSION(:), POINTER :: FPTR
      
      INTERFACE ! strlen is a standard C function from <string.h>
         FUNCTION strlen(string) RESULT(len) BIND(C,NAME="strlen")
            USE ISO_C_BINDING
            TYPE(C_PTR), VALUE :: string ! A C pointer
         END FUNCTION
      END INTERFACE   
      
      IF(C_ASSOCIATED(CPTR)) THEN
         CALL C_F_POINTER(FPTR=FPTR, CPTR=CPTR, SHAPE=[strlen(CPTR)])
      ELSE
         ! To avoid segfaults, associate FPTR with a dummy target:
         FPTR=>dummy_string
      END IF
            
   END FUNCTION
    
   SUBROUTINE SeaFEM_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
        ! This routine is called at the start of the simulation to perform initialization steps.
        ! The parameters are set here and not changed during the simulation.
        ! The initial states and initial guess for the input are defined.
        !..................................................................................................................................

        TYPE(SeaFEM_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
        TYPE(SeaFEM_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
        TYPE(SeaFEM_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
        TYPE(SeaFEM_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
        TYPE(SeaFEM_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
        TYPE(SeaFEM_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
        TYPE(SeaFEM_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
        TYPE(SeaFEM_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                        !   only the output mesh is initialized)
        TYPE(SeaFEM_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables 
        REAL(DbKi),                       INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                        !   (1) SeaFEM_UpdateStates() is called in loose coupling &
                                                                        !   (2) SeaFEM_UpdateDiscState() is called in tight coupling.
                                                                        !   Input is the suggested time from the glue code;
                                                                        !   Output is the actual coupling interval that will be used
                                                                        !   by the glue code.
        TYPE(SeaFEM_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
        INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
        CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
        
        ! local variables
         
        INTEGER(IntKi)                                    :: NumOuts                           
        INTEGER(IntKi)                                    :: ErrStat2                            ! local error status
        CHARACTER(1024)                                   :: ErrMsg2                             ! local error message
      
        ! Initialize variables

        ErrStat = ErrID_None
        ErrMsg  = ""
        NumOuts = 2
        
#ifdef _WIN32 
      module_handle=LoadLibrary(C_NULL_CHAR)
      proc=GetProcAddress(module_handle,"Exchange_Fast_Data"C)
      CALL C_F_PROCPOINTER(proc,EXCHANGE_DATA)
      proc=GetProcAddress(module_handle,"Running_Fast_Update"C)
      CALL C_F_PROCPOINTER(proc,UPDATE_SEAFEM)
      proc=GetProcAddress(module_handle,"End_Fast_Coupling"C)
      CALL C_F_PROCPOINTER(proc,END_TIMELOOP)
#else
      linux_handle=DLOpen(C_NULL_CHAR,IOR(RTLD_NOW, RTLD_GLOBAL))
      IF(.NOT.C_ASSOCIATED(linux_handle)) THEN
        WRITE(*,*) "Error in dlopen: ", C_F_STRING(DLError())
        STOP
      END IF
      proc=DLSym(linux_handle,"Exchange_Fast_Data"//C_NULL_CHAR)
      IF(.NOT.C_ASSOCIATED(funptr)) THEN
        WRITE(*,*) "Error in dlsym (Exchange_Fast_Data): ", C_F_STRING(DLError())
        STOP
      END IF
      CALL C_F_PROCPOINTER(proc,EXCHANGE_DATA)
      proc=DLSym(linux_handle,"Running_Fast_Update"//C_NULL_CHAR)
      IF(.NOT.C_ASSOCIATED(funptr)) THEN
        WRITE(*,*) "Error in dlsym (Running_Fast_Update): ", C_F_STRING(DLError())
        STOP
      END IF
      CALL C_F_PROCPOINTER(proc,UPDATE_SEAFEM)
      proc=DLSym(linux_handle,"End_Fast_Coupling"//C_NULL_CHAR)
      IF(.NOT.C_ASSOCIATED(funptr)) THEN
        WRITE(*,*) "Error in dlsym (End_Fast_Coupling): ", C_F_STRING(DLError())
        STOP
      END IF
      CALL C_F_PROCPOINTER(proc,END_TIMELOOP)
#endif
        
        ! Initialize the NWTC Subroutine Library

        CALL NWTC_Init( )
        
        ! Display the module information

        CALL DispNVD( SeaFEM_Ver )

        ! Call SeaFEM to obtain initial values (Time, Gravity and TMax)
      
        ! Define parameters here:

        p%DT  = Interval
        
        ! Define initial system states here:

        x%DummyContState           = 0
        xd%DummyDiscState          = 0
        z%DummyConstrState         = 0
        OtherState%T               = 0.0
        OtherState%perDOF          = 0
        OtherState%Out_Flag        = 1

        ! Define initial guess for the system inputs here:

        ! Define system output initializations (set up mesh) here:
      
        ! Create the input and output meshes associated with lumped load at the WAMIT reference point (WRP)
      
         call MeshCreate( BlankMesh        = u%PRPMesh            &
                        ,IOS               = COMPONENT_INPUT   &
                        ,Nnodes            = 1           &
                        ,ErrStat           = ErrStat2          &
                        ,ErrMess           = ErrMsg2           &
                        ,TranslationDisp   = .TRUE.            &
                        ,Orientation       = .TRUE.            &
                        ,TranslationVel    = .TRUE.            &
                        ,RotationVel       = .TRUE.            &
                        ,TranslationAcc    = .TRUE.            &
                        ,RotationAcc       = .TRUE.)
         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WAMIT_Init')
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL Cleanup()
               RETURN
            END IF
         
            ! Create the node on the mesh
  
         CALL MeshPositionNode (u%PRPMesh                                &
                                 , 1                              &
                                 , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)    &  
                                 , ErrStat2                           &
                                 , ErrMsg2                            )
      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WAMIT_Init')

      
            ! Create the mesh element
         CALL MeshConstructElement (  u%PRPMesh              &
                                     , ELEMENT_POINT      &                         
                                     , ErrStat2           &
                                     , ErrMsg2            &
                                     , 1              &
                                                 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WAMIT_Init')

      CALL MeshCommit ( u%PRPMesh              &
                        , ErrStat2            &
                        , ErrMsg2             )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WAMIT_Init')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF      

        call MeshCopy ( SrcMesh   = u%PRPMesh           &
                       ,DestMesh  = y%PRPMesh           &
                       ,CtrlCode  = MESH_SIBLING     &
                       ,IOS       = COMPONENT_OUTPUT &
                       ,ErrStat   = ErrStat2         &
                       ,ErrMess   = ErrMsg2          &
                       ,Force     = .TRUE.           &
                       ,Moment    = .TRUE.           )
        
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WAMIT_Init')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF    
      u%PRPMesh%RemapFlag  = .TRUE.
      y%PRPMesh%RemapFlag  = .TRUE.
        
        CONTAINS

        SUBROUTINE CleanUp()
      
        END SUBROUTINE CleanUp      
      
   END SUBROUTINE SeaFEM_Init
   
   SUBROUTINE SeaFEM_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
        ! This routine is called at the end of the simulation.
        !..................................................................................................................................

              TYPE(SeaFEM_InputType),           INTENT(INOUT)  :: u           ! System inputs
              TYPE(SeaFEM_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
              TYPE(SeaFEM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
              TYPE(SeaFEM_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_OutputType),          INTENT(INOUT)  :: y           ! System outputs
              TYPE(SeaFEM_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables    
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
              
             ! Initialize ErrStat

              ErrStat = ErrID_None
              ErrMsg  = ""

             ! Place any last minute operations or calculations here:

             ! Close files here:

             ! Destroy the input data:

          CALL SeaFEM_DestroyInput( u, ErrStat, ErrMsg )

             ! Destroy the parameter data:

          CALL SeaFEM_DestroyParam( p, ErrStat, ErrMsg )

             ! Destroy the state data:

          CALL SeaFEM_DestroyContState(   x,           ErrStat, ErrMsg )
          CALL SeaFEM_DestroyDiscState(   xd,          ErrStat, ErrMsg )
          CALL SeaFEM_DestroyConstrState( z,           ErrStat, ErrMsg )
          CALL SeaFEM_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
          
             ! Destroy misc variables:
      
          CALL SeaFEM_DestroyMisc( m, ErrStat, ErrMsg )

             ! Destroy the output data:

          CALL SeaFEM_DestroyOutput( y, ErrStat, ErrMsg )              

   END SUBROUTINE SeaFEM_End
   
   SUBROUTINE SeaFEM_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
        ! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
        ! Continuous, constraint, and discrete states are updated to values at t + Interval.
        !..................................................................................................................................
   
              REAL(DbKi),                        INTENT(IN   ) :: t            ! Current simulation time in seconds
              INTEGER(IntKi),                    INTENT(IN   ) :: n               ! Current step of the simulation: t = n*Interval
              TYPE(SeaFEM_InputType),            INTENT(IN   ) :: Inputs(:)       ! Inputs at InputTimes
              REAL(DbKi),                        INTENT(IN   ) :: InputTimes(:)   ! Times in seconds associated with Inputs
              TYPE(SeaFEM_ParameterType),        INTENT(IN   ) :: p               ! Parameters
              TYPE(SeaFEM_ContinuousStateType),  INTENT(INOUT) :: x               ! Input: Continuous states at t;
                                                                                  ! Output: Continuous states at t + Interval
              TYPE(SeaFEM_DiscreteStateType),    INTENT(INOUT) :: xd              ! Input: Discrete states at t;
                                                                                  ! Output: Discrete states at t + Interval
              TYPE(SeaFEM_ConstraintStateType),  INTENT(INOUT) :: z               ! Input: Constraint states at t;
                                                                                  ! Output: Constraint states at t + Interval
              TYPE(SeaFEM_OtherStateType),       INTENT(INOUT) :: OtherState      ! Other/optimization states
              TYPE(SeaFEM_MiscVarType),          INTENT(INOUT) :: m               !< Initial misc/optimization variables  
              INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat         ! Error status of the operation
              CHARACTER(*),                      INTENT(  OUT) :: ErrMsg          ! Error message if ErrStat /= ErrID_None
              
              ! Local variables

              TYPE(SeaFEM_ContinuousStateType)                 :: dxdt            ! Continuous state derivatives at t
              TYPE(SeaFEM_DiscreteStateType)                   :: xd_t            ! Discrete states at t (copy)
              TYPE(SeaFEM_ConstraintStateType)                 :: z_Residual      ! Residual of the constraint state functions (Z)
              TYPE(SeaFEM_InputType)                           :: u               ! Instantaneous inputs
              INTEGER(IntKi)                                   :: ErrStat2        ! Error status of the operation (secondary error)
              CHARACTER(LEN(ErrMsg))                           :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
              
              ! Initialize variables

              ErrStat   = ErrID_None           ! no error has occurred
              ErrMsg    = ""
              
              ! Get first time derivatives of continuous states (dxdt):

              CALL SeaFEM_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
              IF ( ErrStat >= AbortErrLev ) THEN
                 CALL SeaFEM_DestroyContState( dxdt, ErrStat2, ErrMsg2)
                 RETURN
              END IF
              
              ! Update discrete states:
              ! Note that xd [discrete state] is changed in SeaFEM_UpdateDiscState() so xd will now contain values at t+Interval
              ! We'll first make a copy that contains xd at time t, which will be used in computing the constraint states
              
              CALL SeaFEM_CopyDiscState( xd, xd_t, MESH_NEWCOPY, ErrStat, ErrMsg )

              CALL SeaFEM_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
              IF ( ErrStat >= AbortErrLev ) THEN
                 CALL SeaFEM_DestroyConstrState( Z_Residual, ErrStat2, ErrMsg2)
                 CALL SeaFEM_DestroyContState(   dxdt,       ErrStat2, ErrMsg2)
                 CALL SeaFEM_DestroyDiscState(   xd_t,       ErrStat2, ErrMsg2) 
                 RETURN
              END IF
              
              ! Solve for the constraint states (z) here:
              
              CALL SeaFEM_CalcConstrStateResidual( t, u, p, x, xd_t, z, OtherState, m, Z_Residual, ErrStat, ErrMsg )
              IF ( ErrStat >= AbortErrLev ) THEN
                 CALL SeaFEM_DestroyConstrState( Z_Residual, ErrStat2, ErrMsg2)
                 CALL SeaFEM_DestroyContState(   dxdt,       ErrStat2, ErrMsg2)
                 CALL SeaFEM_DestroyDiscState(   xd_t,       ErrStat2, ErrMsg2) 
                 RETURN
              END IF
              
              ! Destroy local variables before returning
         
              CALL SeaFEM_DestroyConstrState( Z_Residual, ErrStat2, ErrMsg2)
              CALL SeaFEM_DestroyContState(   dxdt,       ErrStat2, ErrMsg2)
              CALL SeaFEM_DestroyDiscState(   xd_t,       ErrStat2, ErrMsg2) 
       
   END SUBROUTINE SeaFEM_UpdateStates
   
   SUBROUTINE SeaFEM_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg)
        ! Routine for computing outputs, used in both loose and tight coupling.
        !..................................................................................................................................
        use, intrinsic :: iso_c_binding
   
              REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
              TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
              TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
              TYPE(SeaFEM_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
              TYPE(SeaFEM_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                              !   nectivity information does not have to be recalculated)
              TYPE(SeaFEM_MiscVarType),         INTENT(INOUT)  :: m           ! Other/optimization states
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
              INTEGER(IntKi)                                   :: ErrStat2        ! Error status of the operation (secondary error)
              CHARACTER(LEN(ErrMsg))                           :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
              
              REAL(ReKi)                           :: q(6), qdot(6), qdotsq(6), qdotdot(6)
              REAL(ReKi)                           :: rotdisp(3)                              ! small angle rotational displacements
              REAL(ReKi)                           :: SeaFEM_Return_Forces(6)
              INTEGER(IntKi)                       :: I
              
              ! Initialize ErrStat
              
              ErrStat = ErrID_None
              ErrMsg  = ""
              
                   ! Determine the rotational angles from the direction-cosine matrix
                rotdisp = GetSmllRotAngs ( u%PRPMesh%Orientation(:,:,1), ErrStat2, ErrMsg2 )
                   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
              
                !BORJA: Aqu\ED se obtienen los movimientos, velocidades y aceleraciones desde los datos de mesh.
                q      = reshape((/real(u%PRPMesh%TranslationDisp(:,1),ReKi),rotdisp(:)/),(/6/))
                qdot   = reshape((/u%PRPMesh%TranslationVel(:,1),u%PRPMesh%RotationVel(:,1)/),(/6/))
                qdotsq   = abs(qdot)*qdot
                qdotdot   = reshape((/u%PRPMesh%TranslationAcc(:,1),u%PRPMesh%RotationAcc(:,1)/),(/6/)) 
              
              IF(OtherState%calcJacobian .AND. OtherState%perDOF.NE.0) THEN
                  ! Jacobian is being calculated but the velocities and positions are not being updated...
                  !Perturbation of acceleration is set to constant as 1.0 for the moment!! JCC: CHANGE THIS!!
                  !q(OtherState%perDOF)=q(OtherState%perDOF)+p%DT*p%DT/4*1.0e-0  !perturbacion
                  !qdot(OtherState%perDOF)=qdot(OtherState%perDOF)+p%DT/2*1.0e-0 !perturbacion
              END IF
              IF(OtherState%perDOF.EQ.6)THEN
                  OtherState%perDOF=0
              END IF
              
              IF(OtherState%T==t)THEN
                 ! WRITE(*,*) "Simulation time = ",t
              ELSE
                  CALL UPDATE_SEAFEM() ! BORJA: Update seafem
                 ! WRITE(*,*) "Simulation time = ",t
                  OtherState%T=t
              END IF
              
              !BORJA: Aqu\ED se intercambia la informaci\F3n directamente con el ejecutable SeaFEM. Mandamos Movimientos y recibimos fuerzas.
              CALL EXCHANGE_DATA(q,qdot,qdotdot,SeaFEM_Return_Forces,OtherState%flag_SeaFEM)
              
              IF (t>=p%TMax) THEN
                  IF(OtherState%Out_flag==(2+2*p%Iterations))THEN
                      CALL END_TIMELOOP()
                  ELSE
                      OtherState%Out_flag=OtherState%Out_Flag+1
                  END IF
              END IF
              
              DO I=1,3
                 y%PRPMesh%Force(I,1)=SeaFEM_Return_Forces(I)
             !    WRITE(*,'(A,I1,A,E)') "Returned Forces Value SF[",I,"] = ",SeaFEM_Return_Forces(I)
              END DO
              DO I=1,3
                y%PRPMesh%Moment(I,1)=SeaFEM_Return_Forces(I+3)
            !     WRITE(*,'(A,I1,A,E)') "Returned Forces Value SF[",I+3,"] = ",SeaFEM_Return_Forces(I+3)
              END DO
      
                 ! Compute outputs here:
              y%DummyOutput    = 2.0_ReKi                
   
   END SUBROUTINE SeaFEM_CalcOutput
   
   SUBROUTINE SeaFEM_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, Z_residual, ErrStat, ErrMsg )
        ! Tight coupling routine for solving for the residual of the constraint state functions
        !..................................................................................................................................

              REAL(DbKi),                       INTENT(IN   )  :: Time        ! Current simulation time in seconds
              TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
              TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
              TYPE(SeaFEM_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
              TYPE(SeaFEM_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables 
              TYPE(SeaFEM_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state functions using
                                                                              !     the input values described above
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
              
              ! Initialize ErrStat

              ErrStat = ErrID_None
              ErrMsg  = ""

              ! Solve for the residual of the constraint state functions here:

              Z_residual%DummyConstrState = 0

   END SUBROUTINE SeaFEM_CalcConstrStateResidual
   
   SUBROUTINE SeaFEM_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
        ! Tight coupling routine for computing derivatives of continuous states
        !..................................................................................................................................

              REAL(DbKi),                       INTENT(IN   )  :: Time        ! Current simulation time in seconds
              TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
              TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
              TYPE(SeaFEM_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
              TYPE(SeaFEM_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables 
              TYPE(SeaFEM_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at t
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
              
              ! Initialize ErrStat

              ErrStat = ErrID_None
              ErrMsg  = ""
              
              ! Compute the first time derivatives of the continuous states here:

              dxdt%DummyContState = 0

   END SUBROUTINE SeaFEM_CalcContStateDeriv
   
   SUBROUTINE SeaFEM_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
        ! Tight coupling routine for updating discrete states
        !..................................................................................................................................

              REAL(DbKi),                       INTENT(IN   )  :: Time        ! Current simulation time in seconds
              INTEGER(IntKi),                   INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
              TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
              TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
              TYPE(SeaFEM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                              !   Output: Discrete states at t + Interval
              TYPE(SeaFEM_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables 
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
              
             ! Initialize ErrStat

             ErrStat = ErrID_None
             ErrMsg  = ""

             ! Update discrete states here:

             xd%DummyDiscState = 0.0  

   END SUBROUTINE SeaFEM_UpdateDiscState
   
END MODULE SeaFEM