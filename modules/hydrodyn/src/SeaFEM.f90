MODULE SeaFEM
    
   USE ISO_C_BINDING    
   USE SeaFEM_Types  
   USE NWTC_Library
   
   IMPLICIT NONE
   
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
   
   CONTAINS
    
   SUBROUTINE SeaFEM_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
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
      
        CALL MeshCreate( BlankMesh         = u%PRPMesh         &  
                        ,IOS               = COMPONENT_INPUT   &
                        ,Nnodes            = 1                 &
                        ,ErrStat           = ErrStat2          &
                        ,ErrMess           = ErrMsg2           &
                        ,TranslationDisp   = .TRUE.            &
                        ,Orientation       = .TRUE.            &
                        ,TranslationVel    = .TRUE.            &
                        ,RotationVel       = .TRUE.            &
                        ,TranslationAcc    = .TRUE.            &
                        ,RotationAcc       = .TRUE.)
        
        ! Create the node on the mesh
         
        CALL MeshPositionNode ( u%PRPMesh                          &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
      
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SeaFEM_Init')
         
         ! Create the mesh element
         
         CALL MeshConstructElement ( u%PRPMesh           &
                                  , ELEMENT_POINT        &                         
                                  , ErrStat2             &
                                  , ErrMsg2              &
                                  , 1                    )
         
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SeaFEM_Init')
         
         CALL MeshCommit ( u%PRPMesh        &
                         , ErrStat2         &
                         , ErrMsg2          )
   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SeaFEM_Init')
         
         CALL MeshCopy (   SrcMesh      = u%PRPMesh              &
                          ,DestMesh     = y%PRPMesh              &
                          ,CtrlCode     = MESH_NEWCOPY           &
                          ,IOS          = COMPONENT_OUTPUT       &
                          ,ErrStat      = ErrStat2               &
                          ,ErrMess      = ErrMsg2                &
                          ,Force        = .TRUE.                 &
                          ,Moment       = .TRUE.                 )
     
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SeaFEM_Init:y%Mesh')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF  
         
         u%PRPMesh%RemapFlag  = .TRUE.
         y%PRPMesh%RemapFlag  = .TRUE.
         
        ALLOCATE( y%WriteOutput(NumOuts), STAT = ErrStat )
        IF ( ErrStat/= 0 ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Error allocating output header and units arrays in SeaFEM_Init'
            RETURN
        END IF
        
        y%DummyOutput = 0
        y%WriteOutput = 0
        
        CONTAINS

        SUBROUTINE CleanUp()
      
        END SUBROUTINE CleanUp      
      
   END SUBROUTINE SeaFEM_Init
   
   SUBROUTINE SeaFEM_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
        ! This routine is called at the end of the simulation.
        !..................................................................................................................................

              TYPE(SeaFEM_InputType),           INTENT(INOUT)  :: u           ! System inputs
              TYPE(SeaFEM_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
              TYPE(SeaFEM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
              TYPE(SeaFEM_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_OutputType),          INTENT(INOUT)  :: y           ! System outputs
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   END SUBROUTINE SeaFEM_End
   
   SUBROUTINE SeaFEM_UpdateStates( Time, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
        ! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
        ! Continuous, constraint, and discrete states are updated to values at t + Interval.
        !..................................................................................................................................
   
              REAL(DbKi),                        INTENT(IN   ) :: Time            ! Current simulation time in seconds
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
              INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat         ! Error status of the operation
              CHARACTER(*),                      INTENT(  OUT) :: ErrMsg          ! Error message if ErrStat /= ErrID_None
       
   END SUBROUTINE SeaFEM_UpdateStates
   
   SUBROUTINE SeaFEM_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
        ! Routine for computing outputs, used in both loose and tight coupling.
        !..................................................................................................................................
        use, intrinsic :: iso_c_binding
   
              REAL(DbKi),                       INTENT(IN   )  :: Time        ! Current simulation time in seconds
              TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
              TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
              TYPE(SeaFEM_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
              TYPE(SeaFEM_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                              !   nectivity information does not have to be recalculated)
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   END SUBROUTINE SeaFEM_CalcOutput
   
   SUBROUTINE SeaFEM_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
        ! Tight coupling routine for solving for the residual of the constraint state functions
        !..................................................................................................................................

              REAL(DbKi),                       INTENT(IN   )  :: Time        ! Current simulation time in seconds
              TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
              TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
              TYPE(SeaFEM_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
              TYPE(SeaFEM_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state functions using
                                                                              !     the input values described above
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   END SUBROUTINE SeaFEM_CalcConstrStateResidual
   
   SUBROUTINE SeaFEM_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
        ! Tight coupling routine for computing derivatives of continuous states
        !..................................................................................................................................

              REAL(DbKi),                       INTENT(IN   )  :: Time        ! Current simulation time in seconds
              TYPE(SeaFEM_InputType),           INTENT(IN   )  :: u           ! Inputs at t
              TYPE(SeaFEM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
              TYPE(SeaFEM_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
              TYPE(SeaFEM_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
              TYPE(SeaFEM_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
              TYPE(SeaFEM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
              TYPE(SeaFEM_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at t
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   END SUBROUTINE SeaFEM_CalcContStateDeriv
   
   SUBROUTINE SeaFEM_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
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
              INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
              CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   END SUBROUTINE SeaFEM_UpdateDiscState
   
END MODULE SeaFEM