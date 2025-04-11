module DarkEnergyThawingPPF
    use DarkEnergyPPF
    use DarkEnergyInterface
    use classes
    implicit none

    private

    type, extends(TDarkEnergyPPF) :: TDarkEnergyThawingPPF
    
    contains
        procedure, nopass :: PythonClass => TDarkEnergyThawingPPF_PythonClass
        procedure, nopass :: SelfPointer => TDarkEnergyThawingPPF_SelfPointer
        procedure :: w_de => TDarkEnergyThawingPPF_w_de
        procedure :: grho_de => TDarkEnergyThawingPPF_grho_de
    end type TDarkEnergyThawingPPF

    public TDarkEnergyThawingPPF
    contains

    subroutine TDarkEnergyThawingPPF_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkEnergyThawingPPF), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType) 
    P => PType

    end subroutine TDarkEnergyThawingPPF_SelfPointer

    function TDarkEnergyThawingPPF_PythonClass()
        character(LEN=:), allocatable :: TDarkEnergyThawingPPF_PythonClass

        TDarkEnergyThawingPPF_PythonClass = 'DarkEnergyThawingPPF'
    end function TDarkEnergyThawingPPF_PythonClass

    function TDarkEnergyThawingPPF_w_de(this, a)
    class(TDarkEnergyThawingPPF) :: this
    real(dl) :: TDarkEnergyThawingPPF_w_de, al
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        TDarkEnergyThawingPPF_w_de= max(-1., this%w_lam+ this%wa*(1._dl-a))
    else
        al=dlog(a)
        if(al <= this%equation_of_state%Xmin_interp) then
            TDarkEnergyThawingPPF_w_de= max(-1._dl, this%equation_of_state%F(1))
        elseif(al >= this%equation_of_state%Xmax_interp) then
            TDarkEnergyThawingPPF_w_de= max(-1._dl, this%equation_of_state%F(this%equation_of_state%n))
        else
            TDarkEnergyThawingPPF_w_de = max(-1._dl, this%equation_of_state%Value(al))
        endif
    endif

    end function TDarkEnergyThawingPPF_w_de  ! equation of state of the PPF DE
    
    function TDarkEnergyThawingPPF_grho_de(this, a) result(grho_de) !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyThawingPPF) :: this
    real(dl) :: grho_de, al, fint, a_piv, a_eff, coeff, pow
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        ! pivot point (DE swaps between lambda and w0wa)
        if (this%wa == 0) then ! account for potential div by 0 problems
            a_piv = 2._dl
        else
            a_piv = (1._dl + this%w_lam + this%wa)/this%wa
        endif
        a_eff = min(a_piv, 1._dl) ! ensures that grho_de=1 @ a=1

        if (this%wa < 0) then ! lambda -> w0wa
            if (a > a_piv) then ! w0wa
                coeff = a
                pow = 1._dl - a
            else ! lambda
                coeff = a_eff
                pow = 1 - a_eff
            endif
        else !  lambda -> w0wa
            if (a > a_piv) then ! lambda
                coeff = 1._dl
                pow = 0.d0
            else ! w0wa
                coeff = a / a_eff
                pow = a_eff - a
            endif
        endif
        grho_de = a**4 * coeff**(-3._dl*(1._dl+this%w_lam+this%wa))
        if (this%wa /= 0) grho_de = grho_de * exp(-3._dl*this%wa*pow)
    else
        if(a == 0.d0)then
            grho_de = 0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
        else
            if (a>=1) then
                fint= 1
            else
                al = dlog(a)
                if(al <= this%logdensity%X(1)) then
                    ! assume here w=w_de(a_min)
                    fint = exp(this%logdensity%F(1) + (1. - 3. * this%equation_of_state%F(1))*(al - this%logdensity%X(1)))
                else
                    fint = exp(this%logdensity%Value(al))
                endif
            end if
            grho_de = fint
        endif
    endif

    end function TDarkEnergyThawingPPF_grho_de
end module DarkEnergyThawingPPF    


