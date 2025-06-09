        type :: inversion_status 
            character :: yon*5
            real*8, allocatable :: bc(:)
        end type

        type :: search_status 
            character :: yon_ct
            real*8 :: stp 
        end type
