!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

  ! This is a module that provides some basic functionality for doing things 
  ! with ESMF and the NASA tile file

module NASATileFile
  use ESMF
  use NetCDF
  
  implicit none
  
#define NTF_SRCLINE __FILE__,__LINE__
#define NTF_CONTEXT file=__FILE__,line=__LINE__
#define NTF_PASSTHRU msg="Internal subroutine call returned Error"

  ! Public interfaces available from this module
  public NTF_Write

  ! Make things private, except that which is explicitly make public
  private
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

    ! Given an ESMF XGrid output a NASA tile file
    subroutine NTF_Write(xgrid, fileName, rc)
#undef  NTF_METHOD
#define NTF_METHOD "NTF_Write()"

      !
      ! !ARGUMENTS:
      type(ESMF_XGrid),          intent(in) :: xgrid
      character (len=*),         intent(in) :: fileName
      integer,                   intent(out), optional :: rc

      ! Local variables
      integer :: localrc

      write(*,*) "Inside NTF_Write() with filename=",fileName
      
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine NTF_Write
               
end module
