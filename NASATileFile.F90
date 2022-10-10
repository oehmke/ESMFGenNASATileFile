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
      type(ESMF_Field) :: centroidLonField, centroidLatField
      type(ESMF_Array) :: sideBIndArray
      integer, pointer :: sideBIndPtr(:)
      real(ESMF_KIND_R8), pointer :: areaPtr(:) 
      integer :: i, clbnd(1), cubnd(1)
      integer :: lDE, localDECount
      type(ESMF_VM) :: vm
      integer :: localPet
      real(ESMF_KIND_R8), pointer :: centroidLonPtr(:) , centroidLatPtr(:) 
      integer :: numElements
      

      ! Get VM
      call ESMF_VMGetCurrent(vm, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return     

      ! Get VM info
      call ESMF_VMGet(vm, localPet=localPet, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return     
      
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


      ! Create Fields to hold centroid coordinates
      centroidLonField=ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      centroidLatField=ESMF_FieldCreate(xgrid, typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Get coordinates for XGrid centroids
      call getCentroidIntoFields(xgrid, centroidlonField, centroidLatField, rc=localrc)
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

         ! Get centroid
         call ESMF_FieldGet(centroidLonField, lDE, farrayPtr=centroidLonPtr, rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      

         call ESMF_FieldGet(centroidLatField, lDE, farrayPtr=centroidLatPtr, rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      

         
         ! Loop writing data to file
         do i=clbnd(1), cubnd(1)

            ! Write cell types
            if (sideBIndPtr(i) == 1) then 
               ! Write land type
               ! TODO: need to add lake and landice when we have them coming out
               write(unitNum,'(i)',advance="no") 100
            else if (sideBIndPtr(i) == 2) then 
               ! Write ocean type
               write(unitNum,'(i)',advance="no") 0
            endif

            ! Write cell area
            write(unitNum,'(f)',advance="no") areaPtr(i)

            ! Write centroids
            write(unitNum,'(f f)',advance="no") centroidLonPtr(i), centroidLatPtr(i)
 
            ! TODO: fill in second set of info specifc to grid types
            if (sideBIndPtr(i) == 1) then 
               ! Write land info
            else if (sideBIndPtr(i) == 2) then 
               ! Write ocean info
            endif

            ! End line
            write(unitNum,'(A)',advance="yes") " "
         enddo
         
      enddo
      
      ! Close file
      close(unit=unitNum)

      ! Destroy tmp Fields
      call ESMF_FieldDestroy(sideBIndField, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      call ESMF_FieldDestroy(areaField, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      call ESMF_FieldDestroy(centroidLonField, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      call ESMF_FieldDestroy(centroidLatField, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      


      
      ! Return successfully
      if (present(rc)) rc = ESMF_SUCCESS
    end subroutine NTF_Write

    ! Get centroids into Fields
    ! Right now this is just a wrapper to go from XGrids output of a fortran array to 
    ! to an ESMF Field, eventually I'm going to add an option to get ESMF Fields/Arrays
    ! directly out of XGridGet, so this will go away. 
    ! Putting these into Fields (vs just a Fotran array) is 
    ! future proofing to allow to to redist, etc. as necessary to output in parallel.
    ! Also, I think that it gives a more consistent view of the data to have it all 
    ! in Fields on the XGrid DistGrid. 
    subroutine getCentroidIntoFields(xgrid, centroidlonField, centroidLatField, rc)
      type(ESMF_XGrid),          intent(in) :: xgrid
      type(ESMF_Field),          intent(inout) :: centroidLonField
      type(ESMF_Field),          intent(inout) :: centroidLatField
      integer,                   intent(out), optional :: rc
      
      real(ESMF_KIND_R8), allocatable :: centroids(:,:) 
      real(ESMF_KIND_R8), pointer :: lonPtr(:)
      real(ESMF_KIND_R8), pointer :: latPtr(:) 
      integer :: i, clbnd(1), cubnd(1)
      integer :: lDE, localDECount
      integer :: localrc
      integer :: dimCount, elementCount
      integer :: pos


      ! Get localDECount
      call ESMF_FieldGet(centroidLonField,  localDECount=localDECount, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Get element info
      call ESMF_XGridGet(xgrid, elementCount=elementCount, dimCount=dimCount, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      

      ! Allocate space for centroids
      allocate(centroids(elementCount,dimCount)) 
      
      ! Get centroids
      call ESMF_XGridGet(xgrid, centroid=centroids, rc=localrc)
      if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
           NTF_CONTEXT, rcToReturn=rc)) return      
            
      ! Loop over Field putting centroids into Field
      pos=1
      do lDE=0,localDECount-1

         ! Get Lon pointer
         call ESMF_FieldGet(centroidLonField, localDE=lDE, &
              farrayPtr=lonPtr, computationalLBound=clbnd, computationalUBound=cubnd, &
              rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      
         
         ! Get Lat pointere
         call ESMF_FieldGet(centroidLatField, localDE=lDE, farrayPtr=latPtr, rc=localrc)
         if (ESMF_LogFoundError(localrc, NTF_PASSTHRU, &
              NTF_CONTEXT, rcToReturn=rc)) return      
         
         ! Loop copying coords to Fields
         do i=clbnd(1), cubnd(1)
           lonPtr(i)=centroids(pos,1)
           latPtr(i)=centroids(pos,2)
           pos=pos+1
         enddo         
      enddo
     
    end subroutine getCentroidIntoFields

end module
