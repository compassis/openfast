MODULE SeaFEM
    
   USE ISO_C_BINDING    
   USE SeaFEM_Types  
   USE NWTC_Library
   
   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: SeaFEM_Ver = ProgDesc( 'SeaFEM', 'v1.00.00', '22-September-2023' )
   
   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SeaFEM_Init                           ! Initialization routine
   
   PUBLIC :: SeaFEM_CalcOutput                     ! Routine for computing outputs 
   
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
        OtherState%Out_Flag        = 1
      
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
 
        ! Determine the rotational angles from the direction-cosine matrix
        rotdisp = GetSmllRotAngs( u%SeaFEMMesh%Orientation(:,:,1), ErrStat, ErrMsg )              
              
        ! Displacements, velocities and accelerations are obteined from the input mesh (12 iterations for time step increment)
        q       = reshape((/real(u%SeaFEMMesh%TranslationDisp(:,1),ReKi),rotdisp(:)/),(/6/))
        qdot    = reshape((/u%SeaFEMMesh%TranslationVel(:,1),u%SeaFEMMesh%RotationVel(:,1)/),(/6/))
        qdotdot = reshape((/u%SeaFEMMesh%TranslationAcc(:,1),u%SeaFEMMesh%RotationAcc(:,1)/),(/6/))        
                
   END SUBROUTINE SeaFEM_CalcOutput   
   
END MODULE SeaFEM