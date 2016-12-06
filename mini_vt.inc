! ==============================================================
!                  Intel Trace Analyzer setup 
! ==============================================================
     common /vtclasses/ main_class, lap_class, non_lap_class, buffering_class
     integer ::         main_class, lap_class, non_lap_class, buffering_class
     common /vtstates/ &
                vt_outside, vt_prebuffer, vt_buffer, vt_diag, vt_mocklap_inner, vt_mocklap_outer
        integer :: &
                vt_outside, vt_prebuffer, vt_buffer, vt_diag, vt_mocklap_inner, vt_mocklap_outer
