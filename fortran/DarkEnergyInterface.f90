    module DarkEnergyInterface
    use precision
    use interpolation
    use classes
    implicit none

    private

    type, extends(TCambComponent) :: TDarkEnergyModel
        logical :: is_cosmological_constant = .true.
        integer :: num_perturb_equations = 0
        !Whether we are working with the Dark Fluid model (ie dark matter and dark energy are treated ad a single fluid)
        logical :: is_df_model = .false.
        real(dl) :: omch2_eff = 0.d0
        real(dl) :: Omega_DE_eff = 0.d0
        real(dl) :: Omega_c_eff = 0.d0
        real(dl) :: xi = 1._dl
    contains
    procedure :: Init
    procedure :: BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
    procedure :: diff_rhopi_Add_Term
    procedure :: PerturbationInitial
    procedure :: PerturbationEvolve
    procedure :: PrintFeedback
    ! do not have to implement w_de or grho_de if BackgroundDensityAndPressure is inherited directly
    procedure :: w_de
    procedure :: w_de_only
    procedure :: dw_de
    procedure :: cs2_de_a
    procedure :: cs2_de_k
    procedure :: cs2_de_ktau
    procedure :: grho_de
    procedure :: grho_cdm
    procedure :: eval_delta
    procedure :: eval_ddelta
    procedure :: eval_dddelta
    procedure :: Effective_w_wa !Used as approximate values for non-linear corrections
    end type TDarkEnergyModel

    type, extends(TDarkEnergyModel) :: TDarkEnergyEqnOfState
        !Type supporting w, wa or general w(z) table
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (an effective value, used e.g. for halofit)
        real(dl) :: wa = 0._dl !may not be used, just for compatibility with e.g. halofit
        real(dl) :: cs2_lam = 1_dl !rest-frame sound speed, though may not be used
        integer :: n_loga = 0
        integer :: n_sigma = 0
        real(dl), dimension (:), allocatable :: loga 
        real(dl), dimension (:), allocatable :: sigma 
        real(dl), dimension (:), allocatable :: delta 
        logical :: use_tabulated_w = .false.  !Use interpolated table; note this is quite slow.
        logical :: use_tabulated_cs2_a = .false.  !Use interpolated table
        logical :: use_tabulated_cs2_k = .false.  !Use interpolated table
        logical :: use_tabulated_cs2_ktau = .false.  !Use interpolated table
        logical :: no_perturbations = .false. !Don't change this, no perturbations is unphysical
        !Interpolations if use_tabulated_w=.true.
        Type(TCubicSpline) :: equation_of_state, logdensity, equation_of_state_DE_only, sound_speed_a, sound_speed_k, sound_speed_ktau
    contains
    procedure :: ReadParams => TDarkEnergyEqnOfState_ReadParams
    procedure :: Init => TDarkEnergyEqnOfState_Init
    procedure :: SetwTable => TDarkEnergyEqnOfState_SetwTable
    procedure :: SetDeltaTable => TDarkEnergyEqnOfState_SetDeltaTable
    procedure :: SetCs2Table_a => TDarkEnergyEqnOfState_SetCs2Table_a
    procedure :: SetCs2Table_k => TDarkEnergyEqnOfState_SetCs2Table_k
    procedure :: SetCs2Table_ktau => TDarkEnergyEqnOfState_SetCs2Table_ktau
    procedure :: PrintFeedback => TDarkEnergyEqnOfState_PrintFeedback
    procedure :: w_de => TDarkEnergyEqnOfState_w_de
    procedure :: w_de_only => TDarkEnergyEqnOfState_w_de_only
    procedure :: dw_de => TDarkEnergyEqnOfState_dw_de
    procedure :: cs2_de_a => TDarkEnergyEqnOfState_cs2_de_a
    procedure :: cs2_de_k => TDarkEnergyEqnOfState_cs2_de_k
    procedure :: cs2_de_ktau => TDarkEnergyEqnOfState_cs2_de_ktau
    procedure :: grho_de => TDarkEnergyEqnOfState_grho_de
    procedure :: grho_cdm => TDarkEnergyEqnOfState_grho_cdm
    procedure :: Effective_w_wa => TDarkEnergyEqnOfState_Effective_w_wa
    Procedure :: eval_delta => TDarkEnergyEqnOfState_eval_delta
    Procedure :: eval_ddelta => TDarkEnergyEqnOfState_eval_ddelta
    Procedure :: eval_dddelta => TDarkEnergyEqnOfState_eval_dddelta
    end type TDarkEnergyEqnOfState

    public TDarkEnergyModel, TDarkEnergyEqnOfState
    contains

    function w_de(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a

    w_de = -1._dl

    end function w_de  ! equation of state of the PPF DE
    
    function w_de_only(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: w_de_only, al
    real(dl), intent(IN) :: a

    w_de_only = -1._dl

    end function w_de_only
    
    function dw_de(this, a, de_only)
    class(TDarkEnergyModel) :: this
    real(dl) :: dw_de
    real(dl), intent(IN) :: a
    integer, intent(in) :: de_only

    dw_de = 0._dl

    end function dw_de
    
    function eval_delta(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: eval_delta
    real(dl), intent(IN) :: a

    eval_delta = 0._dl

    end function eval_delta
    
    function eval_ddelta(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: eval_ddelta
    real(dl), intent(IN) :: a

    eval_ddelta = 0._dl

    end function eval_ddelta
    
    function eval_dddelta(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: eval_dddelta
    real(dl), intent(IN) :: a

    eval_dddelta = 0._dl

    end function eval_dddelta

    function cs2_de_a(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: cs2_de_a
    real(dl), intent(IN) :: a

    cs2_de_a = 1._dl

    end function cs2_de_a
    
    function cs2_de_k(this, k)
    class(TDarkEnergyModel) :: this
    real(dl) :: cs2_de_k
    real(dl), intent(IN) :: k

    cs2_de_k = 1._dl

    end function cs2_de_k
    
    function cs2_de_ktau(this, ktau)
    class(TDarkEnergyModel) :: this
    real(dl) :: cs2_de_ktau
    real(dl), intent(IN) :: ktau

    cs2_de_ktau = 1._dl

    end function cs2_de_ktau

    function grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyModel) :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a

    grho_de =0._dl

    end function grho_de
    
    function grho_cdm(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: grho_cdm, al, fint
    real(dl), intent(IN) :: a

    grho_cdm =0._dl

    end function grho_cdm

    subroutine PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyModel) :: this
    integer, intent(in) :: FeedbackLevel

    end subroutine PrintFeedback


    subroutine Init(this, State)
    use classes
    class(TDarkEnergyModel), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    end subroutine Init

    subroutine BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1_dl
    else
        ! Ensure a valid result
        if (a > 1e-10) then
            grhov_t = grhov * this%grho_de(a) / (a * a)
        else
            grhov_t = 0._dl
        end if
        if (present(w)) w = this%w_de(a)
    end if

    end subroutine BackgroundDensityAndPressure

    subroutine Effective_w_wa(this, w, wa)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    w = -1
    wa = 0

    end subroutine Effective_w_wa


    subroutine PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe=0
    dgqe=0

    end subroutine PerturbedStressEnergy


    function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    ppiedot = 0._dl

    end function diff_rhopi_Add_Term

    subroutine PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y, cs2_lam)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a,adotoa, k, z, y(:), w, cs2_lam
    integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve

    subroutine PerturbationInitial(this, y, a, tau, k)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k
    !Get intinitial values for perturbations at a (or tau)
    if (.not. this%is_df_model) then
        !For standard adiabatic perturbations can usually just set to zero to good accuracy
        y = 0
    endif
    
    end subroutine PerturbationInitial

    subroutine TDarkEnergyEqnOfState_SetwTable(this, a, w, n)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), w(n)
    real(dl), allocatable :: integral(:)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'w table must end at a=1'

    this%use_tabulated_w = .true.
    call this%equation_of_state%Init(log(a), w)

    allocate(integral(this%equation_of_state%n))
    ! log (rho) =  -3 int dlna (1+w)
    call this%equation_of_state%IntegralArray(integral)
    integral  = -3*( (this%equation_of_state%X-this%equation_of_state%X(1)) + integral) + 4*this%equation_of_state%X
    integral = integral - integral(this%equation_of_state%n) !log(a^4 rho_de)) normalized to 0 at a=1
    call this%logdensity%Init(this%equation_of_state%X, integral)
    !Set w and wa to values today (e.g. as too simple first guess for approx fittings etc).
    this%w_lam = w(size(a))
    this%wa = -this%equation_of_state%Derivative(0._dl)

    end subroutine TDarkEnergyEqnOfState_SetwTable
    
    subroutine TDarkEnergyEqnOfState_SetDeltaTable(this, loga, sigma, delta, n_loga, n_sigma)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: n_loga, n_sigma
    real(dl), intent(in) :: loga(n_loga), sigma(n_sigma), delta(n_loga*n_sigma)
    real(dl) :: step
    integer i
    
    !this%use_tabulated_w = .true.
    
    !allocate(this%loga(n_loga))
    this%loga = loga
    this%sigma = sigma
    this%n_loga = n_loga
    this%n_sigma = n_sigma
    this%delta = delta
    
    end subroutine TDarkEnergyEqnOfState_SetDeltaTable
    
    function TDarkEnergyEqnOfState_eval_delta(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl), intent(in) :: a
    real(dl) :: TDarkEnergyEqnOfState_eval_delta, loga
    integer :: i, j, ind
    
    loga = log10(a)
    TDarkEnergyEqnOfState_eval_delta = 0._dl
    do i=1,this%n_sigma
        do j=1,this%n_loga
            ind = (i-1)*this%n_loga+j
            TDarkEnergyEqnOfState_eval_delta = TDarkEnergyEqnOfState_eval_delta + this%delta(ind)*exp(-0.5_dl*(loga-this%loga(j))**2/this%sigma(i)**2)
        end do
    end do

    end function TDarkEnergyEqnOfState_eval_delta
    
    ! Evaluate d(delta)/dlog(a)
    function TDarkEnergyEqnOfState_eval_ddelta(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl), intent(in) :: a
    real(dl) :: TDarkEnergyEqnOfState_eval_ddelta, loga
    integer :: i, j, ind
    
    loga = log10(a)
    TDarkEnergyEqnOfState_eval_ddelta = 0._dl
    do i=1,this%n_sigma
        do j=1,this%n_loga
            ind = (i-1)*this%n_loga+j
            TDarkEnergyEqnOfState_eval_ddelta = TDarkEnergyEqnOfState_eval_ddelta - this%delta(ind)/this%sigma(i)**2*(loga-this%loga(j))*exp(-0.5_dl*(loga-this%loga(j))**2/this%sigma(i)**2)
        end do
    end do
    
    TDarkEnergyEqnOfState_eval_ddelta = TDarkEnergyEqnOfState_eval_ddelta*log(10._dl)

    end function TDarkEnergyEqnOfState_eval_ddelta
    
    ! Evaluate d^2(delta)/dlog(a)^2
    function TDarkEnergyEqnOfState_eval_dddelta(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl), intent(in) :: a
    real(dl) :: TDarkEnergyEqnOfState_eval_dddelta, loga
    integer :: i, j, ind
    
    loga = log10(a)
    TDarkEnergyEqnOfState_eval_dddelta = 0._dl
    do i=1,this%n_sigma
        do j=1,this%n_loga
            ind = (i-1)*this%n_loga+j
            TDarkEnergyEqnOfState_eval_dddelta = TDarkEnergyEqnOfState_eval_dddelta - this%delta(ind)/this%sigma(i)**2*(1._dl-(loga-this%loga(j))**2/this%sigma(i)**2)*exp(-0.5_dl*(loga-this%loga(j))**2/this%sigma(i)**2)
        end do
    end do
    
    TDarkEnergyEqnOfState_eval_dddelta = TDarkEnergyEqnOfState_eval_dddelta*log(10._dl)**2

    end function TDarkEnergyEqnOfState_eval_dddelta


    subroutine TDarkEnergyEqnOfState_SetCs2Table_a(this, a, cs2_a, n)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), cs2_a(n)
    real(dl), allocatable :: integral(:)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'cs2 table must end at a=1'

    this%use_tabulated_cs2_a = .true.
    call this%sound_speed_a%Init(log(a), cs2_a)

    end subroutine TDarkEnergyEqnOfState_SetCs2Table_a
    
    subroutine TDarkEnergyEqnOfState_SetCs2Table_k(this, k, cs2_k, n)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), cs2_k(n)
    real(dl), allocatable :: integral(:)

    this%use_tabulated_cs2_k = .true.
    call this%sound_speed_k%Init(log(k), cs2_k)

    end subroutine TDarkEnergyEqnOfState_SetCs2Table_k
    
    subroutine TDarkEnergyEqnOfState_SetCs2Table_ktau(this, ktau, cs2_ktau, n)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: ktau(n), cs2_ktau(n)
    real(dl), allocatable :: integral(:)

    this%use_tabulated_cs2_ktau = .true.
    call this%sound_speed_ktau%Init(log(ktau), cs2_ktau)

    end subroutine TDarkEnergyEqnOfState_SetCs2Table_ktau


    function TDarkEnergyEqnOfState_w_de(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_w_de, delta_a, ddelta_a
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w .and. .not. this%is_df_model) then
        TDarkEnergyEqnOfState_w_de= this%w_lam+ this%wa*(1._dl-a)
    else
        delta_a = this%eval_delta(a)
        ddelta_a = this%eval_ddelta(a)
        TDarkEnergyEqnOfState_w_de = -this%Omega_DE_eff/(this%Omega_DE_eff+this%Omega_c_eff/a**3)-1._dl/(3._dl*(1._dl+delta_a))*ddelta_a
    endif

    end function TDarkEnergyEqnOfState_w_de  ! equation of state of the PPF DE

    function TDarkEnergyEqnOfState_w_de_only(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_w_de_only, rho_DE, delta_a, ddelta_a, denom
    real(dl), intent(IN) :: a
    
    delta_a = this%eval_delta(a)
    if (delta_a/this%xi>1e-3) then
        ddelta_a = this%eval_ddelta(a)
        rho_DE = this%Omega_DE_eff+(this%Omega_DE_eff+this%Omega_c_eff/a**3)*delta_a/(1._dl-exp(-2*delta_a/this%xi))
        denom = 1-exp(-2*delta_a/this%xi)
        TDarkEnergyEqnOfState_w_de_only = -1._dl-1._dl/(3._dl*rho_DE)/denom**2*(-3*this%Omega_c_eff/a**3*delta_a*denom+(this%Omega_DE_eff+this%Omega_c_eff/a**3)*ddelta_a*(denom-2._dl/this%xi*delta_a))
    else
        TDarkEnergyEqnOfState_w_de_only = -1._dl+this%xi*this%Omega_c_eff/a**3/(2._dl*this%Omega_DE_eff+(this%Omega_DE_eff+this%Omega_c_eff/a**3)*this%xi)
    end if

    end function TDarkEnergyEqnOfState_w_de_only  ! equation of state of the PPF DE
    
    function TDarkEnergyEqnOfState_dw_de(this, a, de_only)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_dw_de, delta_a, ddelta_a, dddelta_a, w_de, rho_DE, denom, Om_c_a, Om_DF_a
    real(dl), intent(IN) :: a
    integer, intent(in) :: de_only !if 1 then DE only, otherwise DF
    
    delta_a = this%eval_delta(a)
    ddelta_a = this%eval_ddelta(a)
    dddelta_a = this%eval_dddelta(a)
        
    if (de_only == 1) then
        w_de = this%w_de_only(a)
        rho_DE = this%Omega_DE_eff+(this%Omega_DE_eff+this%Omega_c_eff/a**3)*delta_a/(1._dl-exp(-2*delta_a/this%xi))
        if (delta_a/this%xi>1e-3) then
            Om_c_a = this%Omega_c_eff/a**3
            Om_DF_a = this%Omega_DE_eff+Om_c_a
            TDarkEnergyEqnOfState_dw_de = 3*(1._dl+w_de)**2-1._dl/(3._dl*rho_DE*denom)*(12._dl/this%xi*rho_DE*(1._dl+w_de)*ddelta_a*exp(-2*delta_a/this%xi)+9._dl*Om_c_a*delta_a-6._dl*Om_c_a*ddelta_a+6._dl/this%xi*Om_c_a*delta_a*ddelta_a+Om_DF_a*dddelta_a-Om_DF_a*2._dl/this%xi*ddelta_a**2-Om_DF_a*2._dl/this%xi*delta_a*dddelta_a/denom)
        else
            TDarkEnergyEqnOfState_dw_de = -3._dl*this%xi*this%Omega_DE_eff*this%Omega_c_eff/a**3*(2._dl+this%xi)/(2*this%Omega_DE_eff+this%xi*(this%Omega_DE_eff+this%Omega_c_eff/a**3))**2
        end if
    else
        TDarkEnergyEqnOfState_dw_de = -3*this%Omega_DE_eff*this%Omega_c_eff/a**3/(this%Omega_DE_eff+this%Omega_c_eff/a**3)**2-1._dl/3._dl/(1+delta_a)**2*(dddelta_a*(1+delta_a)-ddelta_a**2)
    end if
    end function TDarkEnergyEqnOfState_dw_de

    function TDarkEnergyEqnOfState_cs2_de_a(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_cs2_de_a, al
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_cs2_a) then
        TDarkEnergyEqnOfState_cs2_de_a= this%cs2_lam
    else
        al=dlog(a)
        if(al <= this%sound_speed_a%Xmin_interp) then
            TDarkEnergyEqnOfState_cs2_de_a= this%sound_speed_a%F(1)
        elseif(al >= this%sound_speed_a%Xmax_interp) then
            TDarkEnergyEqnOfState_cs2_de_a= this%sound_speed_a%F(this%sound_speed_a%n)
        else
            TDarkEnergyEqnOfState_cs2_de_a = this%sound_speed_a%Value(al)
        endif
    endif

    end function TDarkEnergyEqnOfState_cs2_de_a
    
    function TDarkEnergyEqnOfState_cs2_de_k(this, k)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_cs2_de_k, kl
    real(dl), intent(IN) :: k

    if(.not. this%use_tabulated_cs2_k) then
        TDarkEnergyEqnOfState_cs2_de_k= 1._dl
    else
        kl=dlog(k)
        if(kl <= this%sound_speed_k%Xmin_interp) then
            TDarkEnergyEqnOfState_cs2_de_k= this%sound_speed_k%F(1)
        elseif(kl >= this%sound_speed_k%Xmax_interp) then
            TDarkEnergyEqnOfState_cs2_de_k= this%sound_speed_k%F(this%sound_speed_k%n)
        else
            TDarkEnergyEqnOfState_cs2_de_k = this%sound_speed_k%Value(kl)
        endif
    endif

    end function TDarkEnergyEqnOfState_cs2_de_k
    
    function TDarkEnergyEqnOfState_cs2_de_ktau(this, ktau)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_cs2_de_ktau, ktaul
    real(dl), intent(IN) :: ktau

    if(.not. this%use_tabulated_cs2_ktau) then
        TDarkEnergyEqnOfState_cs2_de_ktau= 1._dl
    else
        ktaul=dlog(ktau)
        if(ktaul <= this%sound_speed_ktau%Xmin_interp) then
            TDarkEnergyEqnOfState_cs2_de_ktau= this%sound_speed_ktau%F(1)
        elseif(ktaul >= this%sound_speed_ktau%Xmax_interp) then
            TDarkEnergyEqnOfState_cs2_de_ktau= this%sound_speed_ktau%F(this%sound_speed_ktau%n)
        else
            TDarkEnergyEqnOfState_cs2_de_ktau = this%sound_speed_ktau%Value(ktaul)
        endif
    endif

    end function TDarkEnergyEqnOfState_cs2_de_ktau


    subroutine TDarkEnergyEqnOfState_Effective_w_wa(this, w, wa)
    class(TDarkEnergyEqnOfState), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    w = this%w_lam
    wa = this%wa

    end subroutine TDarkEnergyEqnOfState_Effective_w_wa

    function TDarkEnergyEqnOfState_grho_de(this, a) result(grho_de) !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w .and. .not. this%is_df_model) then
        grho_de = a ** (1._dl - 3. * this%w_lam - 3. * this%wa)
        if (this%wa/=0) grho_de=grho_de*exp(-3. * this%wa * (1._dl - a))
    else
        if(a == 0.d0)then
            grho_de = 0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
        else
            if (a>=1) then
                fint= 1
            else
                al = dlog(a)
                if (.not. this%is_df_model) then
                    if(al <= this%logdensity%X(1)) then
                        ! assume here w=w_de(a_min)
                        fint = exp(this%logdensity%F(1) + (1. - 3. * this%equation_of_state%F(1))*(al - this%logdensity%X(1)))
                    else
                        fint = exp(this%logdensity%Value(al))
                    endif
                else
                    fint = (this%Omega_c_eff*a+this%Omega_DE_eff*a**4)*(1+this%eval_delta(a))/((this%Omega_c_eff+this%Omega_DE_eff)*(1+this%eval_delta(1._dl)))
                end if
            end if
            grho_de = fint
        endif
    endif

    end function TDarkEnergyEqnOfState_grho_de
    
    function TDarkEnergyEqnOfState_grho_cdm(this, a) result(grho_cdm)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: grho_cdm, al, fint, delta_a
    real(dl), intent(IN) :: a
    
    delta_a = this%eval_delta(a)
    if (delta_a/this%xi<1e-3) then
        grho_cdm = this%Omega_c_eff/a-this%xi/2._dl*(this%Omega_DE_eff*a**2+this%Omega_c_eff/a)
  
    else
        grho_cdm = this%Omega_c_eff/a+(this%Omega_DE_eff*a**2+this%Omega_c_eff/a)*delta_a/(1-exp(2*delta_a/this%xi))
    end if

    end function TDarkEnergyEqnOfState_grho_cdm

    subroutine TDarkEnergyEqnOfState_PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') &
        &   this%w_lam, this%wa

    end subroutine TDarkEnergyEqnOfState_PrintFeedback

    subroutine TDarkEnergyEqnOfState_ReadParams(this, Ini)
    use IniObjects
    use FileUtils
    class(TDarkEnergyEqnOfState) :: this
    class(TIniFile), intent(in) :: Ini
    real(dl), allocatable :: table(:,:)

    this%use_tabulated_w = Ini%Read_Logical('use_tabulated_w', .false.)
    if(.not. this%use_tabulated_w)then
        this%w_lam = Ini%Read_Double('w', -1.d0)
        this%wa = Ini%Read_Double('wa', 0.d0)
        ! trap dark energy becoming important at high redshift 
        ! (will still work if this test is removed in some cases)
        if (this%w_lam + this%wa > 0) &
             error stop 'w + wa > 0, giving w>0 at high redshift'
    else
        call File%LoadTxt(Ini%Read_String('wafile'), table)
        call this%SetwTable(table(:,1),table(:,2), size(table(:,1)))
    endif

    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)
    this%use_tabulated_cs2_a = Ini%Read_Logical('use_tabulated_cs2_a', .false.)
    if(this%use_tabulated_cs2_a)then
        call File%LoadTxt(Ini%Read_String('cs2file_a'), table)
        call this%SetCs2Table_a(table(:,1),table(:,2), size(table(:,1)))
    endif
    
    this%use_tabulated_cs2_k = Ini%Read_Logical('use_tabulated_cs2_k', .false.)
    if(this%use_tabulated_cs2_k)then
        call File%LoadTxt(Ini%Read_String('cs2file_k'), table)
        call this%SetCs2Table_a(table(:,1),table(:,2), size(table(:,1)))
    endif
    
    this%is_df_model = Ini%Read_Logical('is_df_model', .false.)
    this%omch2_eff = Ini%Read_Double('omch2_eff', 0.d0)
    this%Omega_DE_eff = Ini%Read_Double('Omega_DE_eff', 0.d0)
    this%Omega_c_eff = Ini%Read_Double('Omega_c_eff', 0.d0)

    end subroutine TDarkEnergyEqnOfState_ReadParams


    subroutine TDarkEnergyEqnOfState_Init(this, State)
    use classes
    class(TDarkEnergyEqnOfState), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_cosmological_constant = .not. this%use_tabulated_w .and. &
        &  abs(this%w_lam + 1._dl) < 1.e-6_dl .and. this%wa==0._dl .and. .not. this%is_df_model

    end subroutine TDarkEnergyEqnOfState_Init


    end module DarkEnergyInterface
