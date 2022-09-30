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
      integer, parameter :: unitNum=14
      integer :: localrc
      type(ESMF_Field) :: sideBIndField
      type(ESMF_Field) :: areaField
      type(ESMF_Array) :: sideBIndArray
      integer, pointer :: sideBIndPtr(:)
      real(ESMF_KIND_R8), pointer :: areaPtr(:) 
      integer :: i, clbnd(1), cubnd(1)
      integer :: lDE, localDECount
      type(ESMF_VM) :: vm
      integer :: localPet
      real(ESMF_KIND_R8), pointer :: centroids(:) 
      integer :: numElements
      

      ! Get VM
      call ESMF_VMGetCurrent(vm, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return     

      ! Get VM info
      call ESMF_VMGet(vm, localPet=localPet, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return     
      
      write(*,*) "Inside NTF_Write() with filename=",fileName

      ! Create a Field to hold sideBInd on the Xgrid
      sideBIndField=ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_I4, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Get Array
      call ESMF_FieldGet(sideBIndField,  array=sideBIndArray, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Get Info from XGrid
      call ESMF_XGridGet(xgrid, &
           sideBGeomIndArray=sideBIndArray, &
           rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Create a Field to hold Area on the Xgrid
      areaField=ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Get area
      call ESMF_FieldRegridGetArea(areaField, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Get localDECount
      call ESMF_FieldGet(sideBIndField,  localDECount=localDECount, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      
      
      
      ! Open file
      open(unit=unitNum, file=fileName)
      
      ! Loop over Field
      do lDE=0,localDECount-1

         ! Get tile type
         call ESMF_FieldGet(sideBIndField, lDE, &
              farrayPtr=sideBIndPtr, computationalLBound=clbnd, computationalUBound=cubnd, &
              rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      
         
         ! Get area
         call ESMF_FieldGet(areaField, lDE, farrayPtr=areaPtr, rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      
         
         ! Loop writing data to file
         do i=clbnd(1), cubnd(1)
            write(unitNum,*) sideBIndPtr(i),areaPtr(i)
         enddo
         
      enddo
      
      ! Close file
      close(unit=unitNum)
      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine NTF_Write

    subroutine getCentroidIntoFields(xgrid, centroidlonField, centroidLatField, rc)
      type(ESMF_XGrid),          intent(in) :: xgrid
      type(ESMF_Field),          intent(inout) :: centroidLon
      type(ESMF_Field),          intent(inout) :: centroidLat
      integer,                   intent(out), optional :: rc
      
      real(ESMF_KIND_R8), pointer :: lonPtr(:)
      real(ESMF_KIND_R8), pointer :: latPtr(:) 
      integer :: i, clbnd(1), cubnd(1)
      integer :: lDE, localDECount



      ! Get localDECount
      call ESMF_FieldGet(centroidLonField,  localDECount=localDECount, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      
            
      ! Loop over Field
      do lDE=0,localDECount-1

         ! Get Lon pointer
         call ESMF_FieldGet(centroidLonField, lDE, &
              farrayPtr=lonPtr, computationalLBound=clbnd, computationalUBound=cubnd, &
              rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      
         
         ! Get Lat pointere
         call ESMF_FieldGet(centroidLatField, lDE, farrayPtr=latPtr, rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      
         
         ! Loop copying coords to Fields
         do i=clbnd(1), cubnd(1)

         enddo
         
      enddo
      

           
    end subroutine getCentroidIntoFields

end module
