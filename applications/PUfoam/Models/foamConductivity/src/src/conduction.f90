!subroutines for calculation of effective conductivity of the foam (conduction-only)
!author: pavel.ferkl@vscht.cz
module conduction
    use constants
    implicit none
    private
    public effcond
contains
!********************************BEGINNING*************************************
!determine effective conductivity of the foam
subroutine effcond
    real(dp) :: xw,xs,f
    character(len=30) :: fmt='(2x,A,1x,es9.3,1x,A)'
    write(*,*) 'Conduction:'
    write(mfi,*) 'Conduction:'
    xw=2*(1+cond1/(2*cond2))/3
    xs=(1+4*cond1/(cond1+cond2))/3
    f=(1-fs)*xw+fs*xs
    kgas=cond1*por/(por+(1-por)*f)
    ksol=cond2*f*(1-por)/(por+(1-por)*f)
    effc=kgas+ksol
    eqc_ross=krad+effc
    gcontr=kgas/eqc_ross
    scontr=ksol/eqc_ross
    rcontr=krad/eqc_ross
    write(*,fmt) 'effective conductivity:',effc*1e3_dp,'mW/m/K'
    write(*,fmt) 'equivalent conductivity:',eqc_ross*1e3_dp,'mW/m/K'
    write(*,fmt) 'contribution of gas:',gcontr*1e2_dp,'%'
    write(*,fmt) 'contribution of solid:',scontr*1e2_dp,'%'
    write(*,fmt) 'contribution of radiation:',rcontr*1e2_dp,'%'
    write(mfi,fmt) 'effective conductivity:',effc*1e3_dp,'mW/m/K'
    write(mfi,fmt) 'equivalent conductivity:',eqc_ross*1e3_dp,'mW/m/K'
    write(mfi,fmt) 'contribution of gas:',gcontr*1e2_dp,'%'
    write(mfi,fmt) 'contribution of solid:',scontr*1e2_dp,'%'
    write(mfi,fmt) 'contribution of radiation:',rcontr*1e2_dp,'%'
end subroutine effcond
!***********************************END****************************************
end module conduction
