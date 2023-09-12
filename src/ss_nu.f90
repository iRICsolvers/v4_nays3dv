module ss_nu_m
	use common_hh
	implicit none
contains
	! ------------------------------------------------
	function ss_nu(usta,hs,xi)
	!-------------------------------------------------
		real(8)::ss_nu,usta,hs,xi

		if(j_snu.eq.0) then
			ss_nu=snu*al_ep
		else if(j_snu.eq.1) then
			ss_nu=kappa*usta*hs*al_ep/6.
		else
			ss_nu=kappa*usta*hs*xi*(1.-xi)*al_ep
		end if
		ss_nu=max(ss_nu,snu)

	end function ss_nu
end module ss_nu_m
