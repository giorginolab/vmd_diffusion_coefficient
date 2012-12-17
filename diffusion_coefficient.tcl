
# Core functions for VMD Diffusion Coefficient Plugin. VMD Diffusion
# Coefficient tool. Computes one, two or three-dimensional MSD-based
# diffusion coefficients of a chosen molecular species. See the
# LICENSE file for copyright information.

# Toni Giorgino, ISIB-Consiglio Nazionale delle Ricerche, 2012
# http://multiscalelab.org/utilities/DiffusionCoefficientTool 

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
	interval_from	-
	interval_to	-
	interval_stride	-
	d               -
	msd             -
    }
    array set arg $arg_defaults

    # List of args in "preferred" order
    variable arg_list {selection dt alongx alongy alongz remove_drift 
	from to step  interval_from interval_to interval_stride  }

    # Status text, bound by the GUI, otherwise unused
    variable status_text
}


# User-accessible proc
proc diffusion_coefficient { args } { return [eval ::diffusion_coefficient::diffusion_coefficient $args] }


# Help
proc ::diffusion_coefficient::diffusion_coefficient_usage { } {
    variable arg
    variable arg_list
    puts "VMD Diffusion Coefficient tool. Computes one, two or three-dimensional
MSD-based diffusion coefficients of a chosen molecular species.

Usage: diffusion_coefficient <options> <command>

Command is one of:
-msd <NN>     Compute mean squared displacement (MSD) at a tau of
              NN frames; equivalent to msd_interval -from NN -to NN.
              Returns a value as Angstrom^2 . This is the recommended 
              way of using the plugin.
-msd range    Compute MSD for taus between -from and -to (mandatory)
              Returns two lists of {tau} {MSD(tau)}
-d range      Compute D(tau)=MSD(tau)/(2*D*tau) between -from and
              -to (mandatory). Returns two lists of {tau} {D(tau)}

See http://multiscalelab.org/utilities/DiffusionCoefficientTool for
definitions. If you don't understand what MSD(tau) is, don't use this
tool.

Toni Giorgino, ISIB-National Research Council of Italy.

Options (with defaults):"
    foreach k $arg_list {
	puts "   -$k \"$arg($k)\""
    }

}


# Command line parsing (sets namespace variables). TODO: allow short
# substrings, e.g. -sel
proc ::diffusion_coefficient::parse_args {args} {
    variable arg

    for {set i 0} {$i<[llength $args]} {incr i} {
	set a [lindex $args $i]
	if [regexp {^-} $a] {
	    set a [string trimleft $a -]
	    if {![info exists arg($a)]} {
		error "Unknown option: $a"
	    } else {
		incr i
		set v [lindex $args $i]
		set arg($a) $v
	    }
	} else {
	    error "Unknown command: $a"
	}
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


    if { ($arg(msd)=="-" && $arg(d)=="-") ||
         ($arg(msd)!="-" && $arg(d)!="-") } {
	error "Exactly one of -msd or -d must be given"
    }

    # Execute
    if { $arg(msd)=="range" } {
	return [compute_avg_msd]; # MSD range
    } elseif [string is integer $arg(msd)] {
	set tau $arg(msd); # MSD integer
	set arg(from) $tau
	set arg(to)   $tau
	set arg(step) 1
	lassign [compute_avg_msd] tlist msdlist
	return [lindex $msdlist 0]
    } elseif { $arg(d)=="range" } { # D range
	    lassign [compute_avg_msd] tlist msdlist
	    set dlist [msd_to_d $tlist $msdlist]
	    return [list $tlist $dlist]
    } else {
	error "Unknown invokation type."
    }

}


# Default analysis interval is the whole trajectory, all frames
proc ::diffusion_coefficient::set_default_interval {} {
    variable arg
    if [molinfo num] {
	set nf [molinfo top get numframes]
	set arg(interval_from) 0
	set arg(interval_to) [expr $nf-1]
	set arg(interval_stride) 1
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
proc diffusion_coefficient::set_status {msg} {
    variable status_text $msg
    update
    puts "$msg"
}


# Number of dimensions
proc diffusion_coefficient::nd {} {
    variable arg
    set alongx $arg(alongx)
    set alongy $arg(alongy)
    set alongz $arg(alongz)
    set ND [expr $alongx+$alongy+$alongz]
    return $ND
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
# molecule, returns two lists with lag times and MSD (floats).
# Drift removal: see
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1303338/
proc diffusion_coefficient::compute_avg_msd {} {
    variable arg

    set from $arg(from)
    set to $arg(to)
    set step $arg(step)
    set interval_from $arg(interval_from)
    set interval_to $arg(interval_to)
    set interval_stride $arg(interval_stride)

    set alongx $arg(alongx); set alongy $arg(alongy); set alongz $arg(alongz)
    
    variable xt;   variable yt;   variable zt

    check_selection

    set_status "Initializing"
    set as [atomselect top $arg(selection)]
    set T [molinfo top get numframes]
    set N [$as num]


    if {$to=="last"}		{ set to [expr $N-1] }
    if {$interval_to=="last"}	{ set interval_to [expr $N-1] }

    # make three monster arrays x/y/z arranged for easy indexing
    # lindex $xt 4   returns the vector of all X's at time 4
    set xt {};    set yt {};     set zt {}
    for {set t 0} {$t<$T} {incr t} {
	$as frame $t
	set xvec [$as get x]
	set yvec [$as get y]
	set zvec [$as get z]
	if {$arg(remove_drift)==1} {
	    set xvec [veccenter $xvec]
	    set yvec [veccenter $yvec]
	    set zvec [veccenter $zvec]
	}
	lappend xt $xvec
	lappend yt $yvec
	lappend zt $zvec
    }
    $as delete

    # Form windows of varying sizes; ws=tau
    set tau_list {}
    set msd_list {}
    for {set ws $from} {$ws<=$to} {incr ws $step} {
	set msdavg 0
	set ns 0
	# and slide them
	for {set t0 $interval_from} \
	    {$t0<[expr $interval_to-$ws]} \
	    {incr t0 $interval_stride} {
		set t1 [expr $t0+$ws]
		set msd [msd_between  $t0 $t1]
		set msdavg [expr $msdavg+$msd]
		incr ns
	    }

	# convert frames into time units
	set tau [expr $ws*$arg(dt)]
	lappend tau_list $tau
	lappend msd_list [expr 1.*$msdavg/$ns]

	set_status [format "Computing: %2.0f%% done" \
		    [expr 100.*($ws-$from+1)/($to-$from+1)] ]
    }

    # return
    set_status "Ready"
    return [list $tau_list $msd_list]
}



# Gets tau, MSD(tau) and returns MSD(tau)/2/D/tau
proc diffusion_coefficient::msd_to_d {tau_list msd_list} {
    variable arg

    set dt $arg(dt)
    foreach tau $tau_list msd $msd_list {
	lappend d_list [expr $msd/2.0/$tau/[nd]]
    }
    return $d_list
}

