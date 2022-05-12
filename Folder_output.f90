subroutine Folder_output (k, istep, foldername)


    use FastScapeContext
    implicit none

    integer, intent(in) :: k, istep
    character(len=k), intent(in) :: foldername


    ffoldername = foldername
    kk = k
    iistep = istep


    return
  end subroutine Folder_output
