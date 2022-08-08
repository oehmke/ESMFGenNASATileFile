!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

program GenNASATileFile
  use ESMF
  !  use NUOPC
  use NASACatchFile
  
  implicit none
  integer :: localrc
  type(ESMF_Field) :: catchField
  
  ! Initialize ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
    defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  

  ! Create Field from file
  catchField=NCF_CreateField("catch_rast_small.nc", rc=localrc)
  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
   
  
  ! Finalize ESMF
  call ESMF_Finalize()
  
end program  
