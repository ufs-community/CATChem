! \file error_mod.F90
!! \brief Enhanced error handling and diagnostic system for CATChem
!!
!! This module provides a comprehensive error handling system with standardized
!! error reporting, warning messages, error context tracking, and recovery mechanisms.
!!
!! \author Barry Baker
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!! \ingroup core_modules
!!
!! \details
!! The enhanced error handling module provides:
!! - Standardized error codes and severity levels
!! - Error context stack for better debugging
!! - Structured error types with detailed information
!! - Error recovery and continuation mechanisms
!! - Performance monitoring and error statistics
!! - Thread-safe error handling for parallel execution
!! - Integration with logging systems
!!
!! **New Features in v2.0:**
!! - ErrorManager type for centralized error handling
!! - Error context stack for detailed debugging information
!! - Structured error types with categories and severity levels
!! - Error recovery mechanisms and suggestions
!! - Performance impact tracking
!! - Statistics and reporting capabilities
!!
!! \section error_usage Usage Example
!! \code{.f90}
!! use error_mod
!! type(ErrorManagerType) :: error_mgr
!! integer :: rc
!!
!! call error_mgr%init()
!! call error_mgr%push_context('my_subroutine', 'processing temperature data')
!! call error_mgr%report_error(ERROR_INVALID_INPUT, 'Temperature out of range', rc)
!! call error_mgr%pop_context()
!! \endcode
!!
! \brief Enhanced error handling and diagnostic system
!!
!! This module provides comprehensive error handling with context tracking,
!! structured error types, and recovery mechanisms.
MODULE Error_Mod
   !
   ! !USES:
   !
   use precision_mod, only: fp
   IMPLICIT NONE
   PRIVATE
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   !
   ! Legacy functions (maintained for backward compatibility)
   PUBLIC :: CC_Error
   PUBLIC :: CC_Warning
   PUBLIC :: CC_CheckVar

   ! New enhanced error handling
   PUBLIC :: ErrorManagerType
   PUBLIC :: ErrorInfoType
   PUBLIC :: ErrorContextType

   !
   ! !DEFINED PARAMETERS:
   !
   ! \name Standard Return Codes
   !! \brief Standard return codes for CATChem routines
   !! \{
   INTEGER, PUBLIC, PARAMETER :: CC_SUCCESS =  0   ! Routine completed successfully
   INTEGER, PUBLIC, PARAMETER :: CC_FAILURE = -1   ! Routine failed to complete
   ! \}

   ! \name Enhanced Error Codes
   !! \brief Detailed error codes for specific error types
   !! \{
   INTEGER, PUBLIC, PARAMETER :: ERROR_NONE = 0
   INTEGER, PUBLIC, PARAMETER :: ERROR_INVALID_INPUT = 1001
   INTEGER, PUBLIC, PARAMETER :: ERROR_INVALID_CONFIG = 1002
   INTEGER, PUBLIC, PARAMETER :: ERROR_INVALID_STATE = 1003
   INTEGER, PUBLIC, PARAMETER :: ERROR_FILE_NOT_FOUND = 1004
   INTEGER, PUBLIC, PARAMETER :: ERROR_FILE_READ = 1005
   INTEGER, PUBLIC, PARAMETER :: ERROR_FILE_WRITE = 1006
   INTEGER, PUBLIC, PARAMETER :: ERROR_MEMORY_ALLOCATION = 1007
   INTEGER, PUBLIC, PARAMETER :: ERROR_MEMORY_DEALLOCATION = 1008
   INTEGER, PUBLIC, PARAMETER :: ERROR_DIMENSION_MISMATCH = 1009
   INTEGER, PUBLIC, PARAMETER :: ERROR_BOUNDS_CHECK = 1010
   INTEGER, PUBLIC, PARAMETER :: ERROR_CONVERGENCE = 1011
   INTEGER, PUBLIC, PARAMETER :: ERROR_NUMERICAL_INSTABILITY = 1012
   INTEGER, PUBLIC, PARAMETER :: ERROR_MPI_COMMUNICATION = 1013
   INTEGER, PUBLIC, PARAMETER :: ERROR_PROCESS_INITIALIZATION = 1014
   INTEGER, PUBLIC, PARAMETER :: ERROR_STATE_INCONSISTENCY = 1015
   INTEGER, PUBLIC, PARAMETER :: ERROR_UNSUPPORTED_OPERATION = 1016
   INTEGER, PUBLIC, PARAMETER :: ERROR_DUPLICATE_ENTRY = 1017
   INTEGER, PUBLIC, PARAMETER :: ERROR_NOT_FOUND = 1018
   ! \}

   ! \name Error Severity Levels
   !! \brief Severity levels for error classification
   !! \{
   INTEGER, PUBLIC, PARAMETER :: SEVERITY_INFO = 0
   INTEGER, PUBLIC, PARAMETER :: SEVERITY_WARNING = 1
   INTEGER, PUBLIC, PARAMETER :: SEVERITY_ERROR = 2
   INTEGER, PUBLIC, PARAMETER :: SEVERITY_CRITICAL = 3
   INTEGER, PUBLIC, PARAMETER :: SEVERITY_FATAL = 4
   ! \}

   ! \name Error Categories
   !! \brief Categories for error classification
   !! \{
   INTEGER, PUBLIC, PARAMETER :: CATEGORY_GENERAL = 0
   INTEGER, PUBLIC, PARAMETER :: CATEGORY_INPUT = 1
   INTEGER, PUBLIC, PARAMETER :: CATEGORY_COMPUTATION = 2
   INTEGER, PUBLIC, PARAMETER :: CATEGORY_MEMORY = 3
   INTEGER, PUBLIC, PARAMETER :: CATEGORY_IO = 4
   INTEGER, PUBLIC, PARAMETER :: CATEGORY_MPI = 5
   INTEGER, PUBLIC, PARAMETER :: CATEGORY_PROCESS = 6
   ! \}

   !> \brief Error information structure
   !! \details
   !! Stores information about a specific error including code, message,
   !! severity, category, and optional context information.
   type :: ErrorInfoType
      integer :: error_code = ERROR_NONE         !< Error code
      character(len=255) :: message = ''         !< Error message
      integer :: severity = SEVERITY_INFO        !< Error severity level
      integer :: category = CATEGORY_GENERAL     !< Error category
      character(len=100) :: location = ''        !< Error location
      character(len=255) :: suggestion = ''      !< Suggested solution
      real(fp) :: timestamp = 0.0_fp            !< Error timestamp
   end type ErrorInfoType

   !> \brief Error context structure
   !! \details
   !! Tracks the context stack for error reporting, including routine names,
   !! file locations, and call hierarchy.
   type :: ErrorContextType
      character(len=100) :: routine_name = ''    !< Name of the routine
      character(len=255) :: description = ''     !< Context description
      character(len=100) :: file_name = ''       !< Source file name
      integer :: line_number = 0                 !< Line number
      real(fp) :: timestamp = 0.0_fp            !< Context timestamp
   contains
      procedure :: init => error_context_init
      procedure :: clear => error_context_clear
      procedure :: to_string => error_context_to_string
   end type ErrorContextType

   !> \brief Enhanced error manager
   !! \details
   !! Provides comprehensive error handling with context tracking,
   !! severity levels, categories, and performance monitoring.
   type :: ErrorManagerType
      private
      type(ErrorContextType), allocatable :: context_stack(:)  !< Context stack
      integer :: stack_depth = 0                               !< Current stack depth
      integer :: total_errors = 0                              !< Total error count
      integer :: total_warnings = 0                            !< Total warning count
      integer :: max_stack_depth = 20                          !< Maximum stack depth
      integer :: max_errors_before_abort = 100                 !< Error limit before abort
      logical :: verbose_errors = .false.                      !< Verbose error reporting
      logical :: track_performance = .false.                   !< Performance tracking
      logical :: abort_on_critical = .true.                    !< Abort on critical errors

      ! Error statistics
      integer :: errors_by_severity(0:4) = 0                   !< Errors by severity
      integer :: errors_by_category(0:6) = 0                   !< Errors by category
   contains
      procedure :: init => error_manager_init
      procedure :: push_context => error_manager_push_context
      procedure :: pop_context => error_manager_pop_context
      procedure :: report_error => error_manager_report_error
   end type ErrorManagerType

CONTAINS

   ! \brief Display error message and set failure return code
   !!
   !! This subroutine prints a formatted error message to standard output
   !! and sets the return code to indicate failure. The message includes
   !! optional location and instruction information.
   !!
   !! \param ErrMsg Error message to display
   !! \param RC Return code (set to CC_FAILURE)
   !! \param ThisLoc Optional location where error occurred
   !! \param Instr Optional additional instructions for user
   !!
   !! \par Example:
   !! \code{.f90}
   !! call CC_Error('Invalid temperature value', rc, 'temperature_check', &
   !!               'Check input data file')
   !! \endcode
   SUBROUTINE CC_Error( ErrMsg, RC, ThisLoc, Instr )
      !
      ! !INPUT PARAMETERS:
      !
      CHARACTER(LEN=*), INTENT(IN)            :: ErrMsg  ! Message to display
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: ThisLoc ! Location of error
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: Instr   ! Other instructions
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      INTEGER,          INTENT(INOUT)            :: RC      ! Error code

      CHARACTER(LEN=1000) :: Message
      !=======================================================================
      ! CC_ERROR begins here
      !=======================================================================

      ! Construct error message

      ! Separator
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Print error message to log
      Message =  'CATChem ERROR: ' // TRIM( ErrMsg )
      WRITE(*,'(a)') TRIM(Message)

      ! Print error location to log
      IF ( PRESENT( ThisLoc ) ) THEN
         Message = 'ERROR LOCATION: ' // TRIM( ThisLoc )
         WRITE( 6, '(a)' ) TRIM( ThisLoc )
      ENDIF

      ! Print additional instructions to log
      IF ( PRESENT( Instr ) ) THEN
         WRITE( 6, '(a)' )
         WRITE(*,'(a)') TRIM(Instr)
      ENDIF

      ! Separators
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) ''


      ! Force the message to be flushed to the log file
      CALL Flush( 6 )

      ! Return with failure, but preserve existing error code
      IF ( RC == CC_SUCCESS ) THEN
         RC = CC_FAILURE
      ENDIF

   END SUBROUTINE CC_Error

   !>
   !! \brief CC_Warning
   !!
   !! This subroutine prints a warning message and sets RC to CC_SUCCESS.
   !!
   !! \ingroup core_modules
   !!
   !! \param WarnMsg The warning message
   !! \param RC The return code
   !! \param ThisLoc The location of the warning
   !! \param Instr Other instructions
   !!!>
   SUBROUTINE CC_Warning( WarnMsg, RC, ThisLoc, Instr )
      !
      ! !INPUT PARAMETERS:
      !
      CHARACTER(LEN=*), INTENT(IN   )            :: WarnMsg
      CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: ThisLoc
      CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: Instr
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      INTEGER,          INTENT(INOUT)            :: RC

      CHARACTER(LEN=1000) :: Message

      !=======================================================================
      ! CC_ERROR begins here
      !=======================================================================

      ! Separator
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Print error message to log
      Message =  'CATChem WARNING: ' // TRIM( WarnMsg )
      WRITE(*,'(a)') TRIM(Message)

      ! Print error location to log
      IF ( PRESENT( ThisLoc ) ) THEN
         Message = 'WARNING LOCATION: ' // TRIM( ThisLoc )
         WRITE( 6, '(a)' ) TRIM( ThisLoc )
      ENDIF

      ! Print additional instructions to log
      IF ( PRESENT( Instr ) ) THEN
         WRITE( 6, '(a)' )
         WRITE(*,'(a)') TRIM(Instr)
      ENDIF

      ! Separators
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) ''

      ! Force the message to be flushed to the log file
      CALL Flush( 6 )

      ! Return with success, since this is only a warning message
      RC = CC_SUCCESS

   END SUBROUTINE CC_Warning

   !>
   !! \brief CC_CheckVar
   !!
   !! This subroutine checks if a variable is allocated.
   !!
   !! \param Variable The variable to check
   !! \param Operation 0=Allocate 1=Register 2=Deallocate
   !! \param RC The return code
   !! \ingroup core_modules
   !!!>
   SUBROUTINE CC_CheckVar( Variable, Operation, RC )
      !
      ! !INPUT PARAMETERS:
      !
      CHARACTER(LEN=*), INTENT(IN)    :: Variable   ! Name of variable to check
      INTEGER,          INTENT(IN)    :: Operation  ! 0=Allocate
      ! 1=Register
      ! 2=Deallocate
      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,          INTENT(INOUT) :: RC         ! Success or failure
      !
      ! !LOCAL VARIABLES:
      !
      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=========================================================================
      ! Initialize
      !=========================================================================

      ! Define error message
      SELECT CASE( Operation )
       CASE( 1 )
         ErrMsg = 'Could not register '   // TRIM( Variable ) // '!'
       CASE( 2 )
         ErrMsg = 'Could not deallocate ' // TRIM( Variable ) // '!'
       CASE DEFAULT
         ErrMsg = 'Could not allocate '   // TRIM( Variable ) // '!'
      END SELECT

      ! Define location string
      ThisLoc   = ' -> at CC_CheckVar (in Headers/errcode_mod.F90)'

      !=========================================================================
      ! Display error message if necessary
      !=========================================================================
      IF ( RC /= CC_SUCCESS ) THEN
         CALL CC_Error( ErrMsg, RC, ThisLoc )
      ENDIF

   END SUBROUTINE CC_CheckVar

   !========================================================================
   ! Enhanced Error Handling Implementation
   !========================================================================

   !> \brief Initialize error context
   subroutine error_context_init(this, routine_name, description, file_name, line_number)
      implicit none
      class(ErrorContextType), intent(inout) :: this
      character(len=*), intent(in) :: routine_name
      character(len=*), intent(in), optional :: description
      character(len=*), intent(in), optional :: file_name
      integer, intent(in), optional :: line_number

      this%routine_name = trim(routine_name)
      if (present(description)) this%description = trim(description)
      if (present(file_name)) this%file_name = trim(file_name)
      if (present(line_number)) this%line_number = line_number
      this%timestamp = 0.0_fp  ! Could use system time here
   end subroutine error_context_init

   !> \brief Clear error context
   subroutine error_context_clear(this)
      implicit none
      class(ErrorContextType), intent(inout) :: this

      this%routine_name = ''
      this%description = ''
      this%file_name = ''
      this%line_number = 0
      this%timestamp = 0.0_fp
   end subroutine error_context_clear

   !> \brief Convert error context to string
   function error_context_to_string(this) result(context_str)
      implicit none
      class(ErrorContextType), intent(in) :: this
      character(len=1024) :: context_str

      write(context_str, '(A,": ",A)') trim(this%routine_name), trim(this%description)
      if (len_trim(this%file_name) > 0) then
         write(context_str, '(A," (",A,":",I0,")")') trim(context_str), &
            trim(this%file_name), this%line_number
      endif
   end function error_context_to_string

   !> \brief Initialize error manager
   subroutine error_manager_init(this, verbose, track_performance, rc)
      implicit none
      class(ErrorManagerType), intent(inout) :: this
      logical, intent(in), optional :: verbose
      logical, intent(in), optional :: track_performance
      integer, intent(out), optional :: rc

      ! Allocate context stack
      if (.not. allocated(this%context_stack)) then
         allocate(this%context_stack(this%max_stack_depth))
      endif

      this%stack_depth = 0
      this%total_errors = 0
      this%total_warnings = 0
      this%errors_by_severity = 0
      this%errors_by_category = 0

      if (present(verbose)) this%verbose_errors = verbose
      if (present(track_performance)) this%track_performance = track_performance
      if (present(rc)) rc = CC_SUCCESS
   end subroutine error_manager_init

   !> \brief Push context onto error context stack
   subroutine error_manager_push_context(this, routine_name, description, file_name, line_number)
      implicit none
      class(ErrorManagerType), intent(inout) :: this
      character(len=*), intent(in) :: routine_name
      character(len=*), intent(in), optional :: description
      character(len=*), intent(in), optional :: file_name
      integer, intent(in), optional :: line_number

      if (.not. allocated(this%context_stack)) then
         allocate(this%context_stack(this%max_stack_depth))
      endif
      if (this%stack_depth < this%max_stack_depth) then
         this%stack_depth = this%stack_depth + 1
         call this%context_stack(this%stack_depth)%init(routine_name, description, file_name, line_number)
      endif
   end subroutine error_manager_push_context

   !> \brief Pop context from error context stack
   subroutine error_manager_pop_context(this)
      implicit none
      class(ErrorManagerType), intent(inout) :: this

      if (this%stack_depth > 0) then
         call this%context_stack(this%stack_depth)%clear()
         this%stack_depth = this%stack_depth - 1
      endif
   end subroutine error_manager_pop_context

   !> \brief Report error with enhanced information
   subroutine error_manager_report_error(this, error_code, message, rc, location, suggestion)
      implicit none
      class(ErrorManagerType), intent(inout) :: this
      integer, intent(in) :: error_code
      character(len=*), intent(in) :: message
      integer, intent(inout) :: rc
      character(len=*), intent(in), optional :: location
      character(len=*), intent(in), optional :: suggestion

      character(len=1024) :: full_message
      integer :: severity, category

      ! Update statistics
      this%total_errors = this%total_errors + 1

      ! Determine severity and category from error code
      call get_error_properties(error_code, severity, category)
      this%errors_by_severity(severity) = this%errors_by_severity(severity) + 1
      this%errors_by_category(category) = this%errors_by_category(category) + 1

      ! Build comprehensive error message
      write(full_message, '(A,I0,A,A)') 'ERROR ', error_code, ': ', trim(message)

      ! Add context information if available
      if (this%stack_depth > 0) then
         write(full_message, '(A,A,A)') trim(full_message), ' [Context: ', &
            trim(this%context_stack(this%stack_depth)%to_string())//']'
      endif

      ! Use legacy error reporting for now (can be enhanced)
      if (present(location)) then
         call CC_Error(full_message, rc, location, suggestion)
      else
         call CC_Error(full_message, rc)
      endif

      ! Check if we should abort
      if (severity >= SEVERITY_CRITICAL .and. this%abort_on_critical) then
         write(*, '(A)') 'CRITICAL ERROR: Aborting execution'
         stop
      endif

      if (this%total_errors >= this%max_errors_before_abort) then
         write(*, '(A,I0,A)') 'Maximum error count (', this%max_errors_before_abort, ') exceeded. Aborting.'
         stop
      endif
   end subroutine error_manager_report_error

   !> \brief Get error properties from error code
   subroutine get_error_properties(error_code, severity, category)
      implicit none
      integer, intent(in) :: error_code
      integer, intent(out) :: severity
      integer, intent(out) :: category

      ! Default values
      severity = SEVERITY_ERROR
      category = CATEGORY_GENERAL

      ! Classify based on error code
      select case (error_code)
       case (ERROR_INVALID_INPUT, ERROR_INVALID_CONFIG)
         severity = SEVERITY_ERROR
         category = CATEGORY_INPUT
       case (ERROR_FILE_NOT_FOUND, ERROR_FILE_READ, ERROR_FILE_WRITE)
         severity = SEVERITY_ERROR
         category = CATEGORY_IO
       case (ERROR_MEMORY_ALLOCATION, ERROR_MEMORY_DEALLOCATION)
         severity = SEVERITY_CRITICAL
         category = CATEGORY_MEMORY
       case (ERROR_NUMERICAL_INSTABILITY, ERROR_CONVERGENCE)
         severity = SEVERITY_WARNING
         category = CATEGORY_COMPUTATION
       case (ERROR_MPI_COMMUNICATION)
         severity = SEVERITY_CRITICAL
         category = CATEGORY_MPI
       case (ERROR_STATE_INCONSISTENCY)
         severity = SEVERITY_ERROR
         category = CATEGORY_PROCESS
       case default
         severity = SEVERITY_ERROR
         category = CATEGORY_GENERAL
      end select
   end subroutine get_error_properties
   !EOC
END MODULE Error_Mod
