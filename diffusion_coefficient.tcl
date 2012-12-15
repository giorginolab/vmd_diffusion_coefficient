
# Core functions to compute mean square displacements

package provide diffusion_coefficient 1.0

namespace eval ::diffusion_coefficient:: {
    # Variables matching command line options. Note to self: putting
    # them in an array was not a good idea. It makes sense when there
    # is a 1-to-1 correspondence between CLI arguments and GUI
    # elements, but provides no advantage in this case, when the GUI
    # provides additional/different functionality.
    variable arg
    variable arg_defaults {
	selection "water and name OH"
	dt		1
	alongx		1  
	alongy		1 
	alongz		1
	remove_drift	1
	from		-
	to		-
	step		-
	window_from	-
	window_to	-
	window_every	-
    }
    array set arg $arg_defaults

    # List of args in "preferred" order
    variable arg_list {selection dt alongx alongy alongz remove_drift 
	from to step  window_from window_to window_every  }

    # Status text, bound by the GUI, otherwise unused
    variable status_text
}


# User-accessible proc
proc diffusion_coefficient { args } { return [eval ::diffusion_coefficient::diffusion_coefficient $args] }


# Help
proc ::diffusion_coefficient::diffusion_coefficient_usage { } {
    variable arg
    variable arg_list
    puts "VMD Diffusion Coefficient tool. Computes one, two or three-dimensional"
    puts "MSD-based diffusion coefficients of a chosen molecular species. "
    puts " "
    puts "Usage: diffusion_coefficient <args> {msd|diffusion}"
    puts "Args (with defaults):"
    foreach k $arg_list {
	puts "   -$k \"$arg($k)\""
    }
    puts " "
    puts "msd: return mean-squared displacements at each lag time tau"
    puts "diffusion: return MSD(tau)/(2*D*tau)"
    puts " "
    puts "See http://multiscalelab.org/utilities/DiffusionCoefficientTool"
}


# Command line parsing (sets namespace variables). TODO: allow short
# substrings, e.g. -sel
proc ::diffusion_coefficient::parse_args {args} {
    variable dp_args
    foreach {a v} $args {
	if {![regexp {^-} $a]} {
	    error "Argument should start with -: $a"
	} 
	set a [string trimleft $a -]
	if {![info exists dp_args($a)]} {
	    error "Unknown argument: $a"
	} 
	set dp_args($a) $v
    }
}


# Main entry point. 
proc ::diffusion_coefficient::diffusion_coefficient {args} {
    variable arg
    variable arg_defaults
    array set arg $arg_defaults
    set_default_interval
    set_default_lags
    if {[llength $args]==0} {
	diffusion_coefficient_usage
	return
    } 
    eval parse_args $args
    parray arg


    # Compute the bare histogram

}


# Default analysis interval is the whole trajectory, all frames
proc ::diffusion_coefficient::set_default_interval {} {
    variable arg
    if [molinfo num] {
	set nf [molinfo top get numframes]
	set arg(window_from) 0
	set arg(window_to) [expr $nf-1]
	set arg(window_every) 1
    }
}


# Default lags is from 1/10 to 1/2 of the whole trajectory, at steps
# of 1/50
proc ::diffusion_coefficient::set_default_lags {} {
    variable arg
    if [molinfo num] {
	set nf [molinfo top get numframes]
	set arg(from) [expr $nf/10]
	set arg(to)   [expr $nf/2 ]
	set arg(step) [expr $nf/50]
    }
}



# Performs sanity checks on selection
proc diffusion_coefficient::check_selection {} {
    variable arg
    set as [atomselect top $arg(selection)]
    if {[$as num]==0} {
	$as delete
	error "Atom selection is empty" 
    } 
    if { [$as num] != [llength [lsort -uniq [$as get fragment]]] } {
	$as delete
	error "Each selected atom should belong to a separate molecule" 
    }
    $as delete
    if {[molinfo top get numframes]<=2} {
	error "Not enough trajectory frames"
    }
}


# Return a zero-centered version of the input list
proc diffusion_coefficient::veccenter {l} {
    set m [vecmean $l]
    set N [llength $l]
    set mN [lrepeat $N $m]
    set r [vecsub $l $mN]
    return $r
}


# If loaded in gui, update message. Otherwise, print.
proc diffusion_coefficient::status {msg} {
    variable status_text $msg
    update
    puts "$msg"
}



# Uses class-variables xt, yt, zt
proc diffusion_coefficient::msd_between {t0 t1 } {
    variable arg
    set alongx $arg(alongx) 
    set alongy $arg(alongy)
    set alongz $arg(alongz)

    variable xt;   variable yt;   variable zt

    set N [llength [lindex $xt 0]]
    set dx2 0
    set dy2 0
    set dz2 0
    
    if {$alongx==1} {
	set dx [vecsub [lindex $xt $t0] [lindex $xt $t1]]
	set dx2 [vecdot $dx $dx]
    }
    if {$alongy==1} {
	set dy [vecsub [lindex $yt $t0] [lindex $yt $t1]]
	set dy2 [vecdot $dy $dy]
    }
    if {$alongz==1} {
	set dz [vecsub [lindex $zt $t0] [lindex $zt $t1]]
	set dz2 [vecdot $dz $dz]
    }
    set msd [expr ($dx2+$dy2+$dz2)/$N ]
    return $msd
}




# Compute the average MSD. Takes data from the currently-loaded
# molecule, returns MSD (an array of floats) indexed by lag time.
# Drift removal: found in http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1303338/
proc diffusion_coefficient::compute_avg_msd {} {
    variable arg

    set selection $arg(selection)
    set remove_drift $arg(remove_drift)

    set from $arg(from)
    set to $arg(to)
    set step $arg(step)
    set window_from $arg(window_from)
    set window_to $arg(window_to)
    set window_every $arg(window_every)

    set alongx $arg(alongx); set alongy $arg(alongy); set alongz $arg(alongz)
    
    variable xt;   variable yt;   variable zt

    check_selection

    status "Initializing"
    set as [atomselect top $selection]
    set T [molinfo top get numframes]
    set N [$as num];		# TODO check N>0


    if {$to=="last"}		{ set to [expr $N-1] }
    if {$window_to=="last"}	{ set window_to [expr $N-1] }

    # make three monster arrays x/y/z arranged for easy indexing
    # lindex $xt 4   returns the vector of all X's at time 4
    set xt {};    set yt {};     set zt {}
    for {set t 0} {$t<$T} {incr t} {
	$as frame $t
	set xvec [$as get x]
	set yvec [$as get y]
	set zvec [$as get z]
	if {$remove_drift==1} {
	    set xvec [veccenter $xvec]
	    set yvec [veccenter $yvec]
	    set zvec [veccenter $zvec]
	}
	lappend xt $xvec
	lappend yt $yvec
	lappend zt $zvec
    }

    # Form windows of varying sizes
    for {set ws $from} {$ws<=$to} {incr ws $step} {
	set msdavg 0
	set ns 0
	# and slide them
	for {set t0 $window_from} \
	    {$t0<[expr $window_to-$ws]} \
	    {incr t0 $window_every} {
		set t1 [expr $t0+$ws]
		set msd [msd_between  $t0 $t1]
		set msdavg [expr $msdavg+$msd]
		incr ns
	    }
	set msdm($ws) [expr 1.*$msdavg/$ns]

	status [format "Computing: %2.0f%% done" \
		    [expr 100.*($ws-$from)/($to-$from)] ]
    }

    # return
    status "Ready"
    $as delete
    return [array get msdm]
}
