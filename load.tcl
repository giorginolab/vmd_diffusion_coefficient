# Startup helper.


set here [file dirname [file normalize [info script]]]

# source [file join $diffusion_coefficient_gui::SCRIPTDIR diffusion_coefficient_gui_ui.tcl]

puts "diffusion_coefficient) Adding directory $here to auto_path and registering the menu"

lappend auto_path $here
package require diffusion_coefficient_gui
::diffusion_coefficient_gui::register_menu

