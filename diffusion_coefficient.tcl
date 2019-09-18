
# Core functions for VMD Diffusion Coefficient Plugin. VMD Diffusion
# Coefficient tool. Computes one, two or three-dimensional MSD-based
# diffusion coefficients of a chosen molecular species. See the
# LICENSE file for copyright information.

# Toni Giorgino, ISIB-Consiglio Nazionale delle Ricerche, 2012
# https://github.com/giorginolab/vmd_diffusion_coefficient/

package provide diffusion_coefficient 1.2

namespace eval ::diffusion_coefficient:: {
    # Variables matching command line options. 
    variable arg
    variable arg_defaults {
	selection "water and noh"
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
	fitD            -
    }
    array set arg $arg_defaults

    # List of args in "preferred" order
    variable arg_list {selection dt alongx alongy alongz remove_drift 
	from to step  interval_from interval_to interval_stride  }

    # Status text, bound by the GUI, otherwise unused
    variable status_text

    # If set the computations stops
    variable abort_flag 0
}


# User-accessible proc
proc diffusion_coefficient { args } {
    return [eval ::diffusion_coefficient::diffusion_coefficient $args] }


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
-fitD range   Compute D by a linear fit of MSD over the specified
              range. Returns a list of {D D_err S S_err} where
              D is the MSD slope divided by 2D; S is the MSD intercept; 
              and _err are the respective standard errors.

See https://github.com/giorginolab/vmd_diffusion_coefficient/ for
definitions. If you don't understand what MSD(tau) is, don't use this
tool.

Toni Giorgino, National Research Council of Italy.

Options (with defaults):"
    foreach k $arg_list {
	puts "   -$k [list $arg($k)]"
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


proc ::diffusion_coefficient::get_cli_equivalent {} {
    variable arg
    variable arg_list
    
    foreach k $arg_list {
	append cli "-$k [list $arg($k)] "
    }
    return $cli
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


    if { ($arg(msd)!="-") + ($arg(d)!="-") + ($arg(fitD)!="-") != 1 } {
	error "Exactly one of -msd, -d, or -fitD must be given"
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
    } elseif { $arg(fitD)=="range" } {
	lassign [compute_avg_msd] tlist msdlist
	set dfit [msd_fit $tlist $msdlist]
	return $dfit
    } else {
	error "Unknown invocation type."
    }

}


# Save to file
proc ::diffusion_coefficient::save_to_file {fn} {
    lassign [compute_avg_msd] tlist msdlist
    set dlist [msd_to_d $tlist $msdlist]
    set dfit [msd_fit $tlist $msdlist]

    set fp [open $fn w]
    puts $fp "# Created with the Diffusion Coefficient Tool for VMD."
    puts $fp "# Options:   [get_cli_equivalent]"
    puts $fp [format "# Fit value for D (first two columns): %.4g \xB1 %.4g \xC5\xB2/ns (intercept %.4g \xB1 %.4g \xC5\xB2)" {*}$dfit]
    puts $fp "# Tau MSD D"
    foreach tau $tlist msd $msdlist d $dlist {
	puts $fp [format "%.4f %.4f %.4f" $tau $msd $d]
    }
    close $fp
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


# Requires atomselect, axis ("x", "y" or "z"), start and end frames.
# Returns the squared distance moved along that axis
proc diffusion_coefficient::delta2_between {as axis t0 t1} {
    variable arg
    $as frame $t0
    set v0 [$as get $axis]
    $as frame $t1
    set v1 [$as get $axis]
    if {$arg(remove_drift)==1} {
	set v0 [veccenter $v0]
	set v1 [veccenter $v1]
    }
    set dv [vecsub $v1 $v0]
    set dv2 [vecdot $dv $dv]
    return $dv2
}


# Uses the atomselect
proc diffusion_coefficient::msd_between {as t0 t1} {
    variable arg
    set alongx $arg(alongx) 
    set alongy $arg(alongy)
    set alongz $arg(alongz)

    set N [$as num]

    set dx2 0
    set dy2 0
    set dz2 0
    
    if {$alongx==1} {
	set dx2 [delta2_between $as x $t0 $t1]
    }
    if {$alongy==1} {
	set dy2 [delta2_between $as y $t0 $t1]
    }
    if {$alongz==1} {
	set dz2 [delta2_between $as z $t0 $t1]
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
    variable abort_flag

    set from $arg(from)
    set to $arg(to)
    set step $arg(step)
    set interval_from $arg(interval_from)
    set interval_to $arg(interval_to)
    set interval_stride $arg(interval_stride)

    set alongx $arg(alongx); set alongy $arg(alongy); set alongz $arg(alongz)

    check_selection

    set_status "Initializing"
    
    set as [atomselect top $arg(selection)]
    
    set T [molinfo top get numframes]
    if {$to=="last"}		{ set to [expr $T-1] }
    if {$interval_to=="last"}	{ set interval_to [expr $T-1] }

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
		set msd [msd_between $as $t0 $t1]
		set msdavg [expr $msdavg+$msd]
		incr ns
	}

	# convert frames into time units
	set tau [expr $ws*$arg(dt)]
	lappend tau_list $tau
	lappend msd_list [expr 1.*$msdavg/$ns]

	set_status [format "Computing: %2.0f%% done" \
			[expr 100.*($ws-$from+1)/($to-$from+1)] ]

	if {$abort_flag == 1} {
	    set abort_flag 0
	    set_status "Aborted."
	    $as delete
	    error Aborted.
	}
    }

    $as delete

    # return
    set_status "Ready"
    return [list $tau_list $msd_list]
}



# Gets tau, MSD(tau) and returns MSD(tau)/2/D/tau
proc diffusion_coefficient::msd_to_d {tau_list msd_list} {
    variable arg

    set dt $arg(dt);		# not needed- we get real times in
    foreach tau $tau_list msd $msd_list {
	lappend d_list [expr $msd/2.0/$tau/[nd]]
    }
    return $d_list
}


# Gets tau, MSD(tau) and returns MSD(tau)/2/D/tau, intercept, and
# standard errors computed by linear fits.
proc diffusion_coefficient::msd_fit {tau_list msd_list} {
    variable arg

    set dt $arg(dt);		# not needed- we get real times in

    set fit [linear-model $tau_list $msd_list 1]
    set fit_slope [lindex $fit 1]
    set fit_int   [lindex $fit 0]
    set fit_slope_se [lindex $fit 7]
    set fit_int_se [lindex $fit 5]

    set out [list \
     	     [expr $fit_slope/2.0/[nd] ] \
     	     [expr $fit_slope_se/2.0/[nd] ] \
     	     $fit_int \
     	     $fit_int_se]

    set_status [format "Fit value for D: %.4g \xB1 %.4g \xC5\xB2/ns (intercept %.4g \xB1 %.4g \xC5\xB2)" {*}$out]
    
    return $out
}






## From TCLLIB, https://github.com/tcltk/tcllib/blob/master/modules/math/statistics.tcl
#
# linear-model
#    Determine the coefficients for a linear regression between
#    two series of data (the model: Y = A + B*X)
#
# Arguments:
#    xdata        Series of independent (X) data
#    ydata        Series of dependent (Y) data
#    intercept    Whether to use an intercept or not (optional)
#
# Result:
#    List of the following items:
#    - (Estimate of) Intercept A
#    - (Estimate of) Slope B
#    - Standard deviation of Y relative to fit
#    - Correlation coefficient R2
#    - Number of degrees of freedom df
#    - Standard error of the intercept A
#    - Significance level of A
#    - Standard error of the slope B
#    - Significance level of B
#
#
proc diffusion_coefficient::linear-model { xdata ydata {intercept 1} } {
   variable TOOFEWDATA

   if { [llength $xdata] < 3 } {
      return -code error -errorcode ARG "$TOOFEWDATA: not enough independent data"
   }
   if { [llength $ydata] < 3 } {
      return -code error -errorcode ARG "$TOOFEWDATA: not enough dependent data"
   }
   if { [llength $xdata] != [llength $ydata] } {
      return -code error -errorcode ARG "$TOOFEWDATA: number of dependent data differs from number of independent data"
   }

   set sumx  0.0
   set sumy  0.0
   set sumx2 0.0
   set sumy2 0.0
   set sumxy 0.0
   set df    0
   foreach x $xdata y $ydata {
      if { $x != "" && $y != "" } {
         set sumx  [expr {$sumx+$x}]
         set sumy  [expr {$sumy+$y}]
         set sumx2 [expr {$sumx2+$x*$x}]
         set sumy2 [expr {$sumy2+$y*$y}]
         set sumxy [expr {$sumxy+$x*$y}]
         incr df
      }
   }

   if { $df <= 2 } {
      return -code error -errorcode ARG "$TOOFEWDATA: too few valid data"
   }
   if { $sumx2 == 0.0 } {
      return -code error -errorcode ARG "$TOOFEWDATA: independent values are all the same"
   }

   #
   # Calculate the intermediate quantities
   #
   set sx  [expr {$sumx2-$sumx*$sumx/$df}]
   set sy  [expr {$sumy2-$sumy*$sumy/$df}]
   set sxy [expr {$sumxy-$sumx*$sumy/$df}]

   #
   # Calculate the coefficients
   #
   if { $intercept } {
      set B [expr {$sxy/$sx}]
      set A [expr {($sumy-$B*$sumx)/$df}]
   } else {
      set B [expr {$sumxy/$sumx2}]
      set A 0.0
   }

   #
   # Calculate the error estimates
   #
   set stdevY 0.0
   set varY   0.0

   if { $intercept } {
      set ve [expr {$sy-$B*$sxy}]
      if { $ve >= 0.0 } {
         set varY [expr {$ve/($df-2)}]
      }
   } else {
      set ve [expr {$sumy2-$B*$sumxy}]
      if { $ve >= 0.0 } {
         set varY [expr {$ve/($df-1)}]
      }
   }
   set seY [expr {sqrt($varY)}]

   if { $intercept } {
      set R2    [expr {$sxy*$sxy/($sx*$sy)}]
      set seA   [expr {$seY*sqrt(1.0/$df+$sumx*$sumx/($sx*$df*$df))}]
      set seB   [expr {sqrt($varY/$sx)}]
      set tA    {}
      set tB    {}
      if { $seA != 0.0 } {
         set tA    [expr {$A/$seA*sqrt($df-2)}]
      }
      if { $seB != 0.0 } {
         set tB    [expr {$B/$seB*sqrt($df-2)}]
      }
   } else {
      set R2    [expr {$sumxy*$sumxy/($sumx2*$sumy2)}]
      set seA   {}
      set tA    {}
      set tB    {}
      set seB   [expr {sqrt($varY/$sumx2)}]
      if { $seB != 0.0 } {
         set tB    [expr {$B/$seB*sqrt($df-1)}]
      }
   }

   #
   # Return the list of parameters
   #
   return [list $A $B $seY $R2 $df $seA $tA $seB $tB]
}
