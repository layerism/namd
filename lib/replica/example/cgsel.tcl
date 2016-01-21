set paramfile "vdw.inp"
set psffile 	"his_glu_sol_ion.psf"
set pdbfile 	"his_glu_sol_ion.pdb"

#mol delete 	top
mol new 		$pdbfile
mol addfile $psffile


set all [atomselect top all]
set s 	[atomselect top "water"]
set p 	[atomselect top "protein"]



set namelist 	[$all get name]
set indlist		""
set ind				0
foreach iname $namelist {
	if {![regexp "^H" $iname]} {
		set indlist "$indlist$ind "	
	}
	incr ind
}
set cg [atomselect top "index $indlist"]
#set cg [atomselect top all]
#set cg [atomselect top "name CA N O C CB CG OH2 SOD CLA"]
#------------------build partition pdf-------------------------
$all 	set beta 0
$p 		set beta 0
$cg 	set beta 1

$all writepdb cg.pdb


#------------------read aa psf and param file-------------------------



set psf 		[open $psffile r]
while {[gets $psf line]!=-1} {
	set len [llength $line]	
	if { $len==9 && [regexp "^\[A-Z]" [lindex $line 3]] } {
		set ind    	[expr [lindex $line 0]-1]
		set charge  [lindex $line 6]
		set type 		[lindex $line 5]
		set aatype($ind) 			$type
		set aacharge($ind)		$charge
	}
}
close $psf



set param 	[open $paramfile r]
while {[gets $param line]!=-1} {
	#puts $line
	if {[regexp "^\[A-Z]" [lindex $line 0]] && [regexp "\[0-9]" [lindex $line 1]] && \
			[regexp "\[0-9]" [lindex $line 2]]} {		
		set type				[lindex $line 0]		
		set vdwEps			[lindex $line 2]	
		set vdwSgm			[lindex $line 3]	
		set aavdwEps($type)	$vdwEps
		set aavdwSgm($type)	[expr $vdwSgm*1.0]
		#puts $line
	}
	if {[regexp "\[0-9]" [lindex $line 4]] && [regexp "\[0-9]" [lindex $line 5]] && \
			[regexp "\[0-9]" [lindex $line 6]]} {		
		set type				[lindex $line 0]		
		set vdwEps			[lindex $line 5]	
		set vdwSgm			[lindex $line 6]	
		#set aavdwEps($type)	$vdwEps
		#set aavdwSgm($type)	$vdwSgm
		#puts $line
	}
}
close $param


#------------------build cg param-------------------------



set aalist 	[$all get index]
set cglist 	[$cg get index]

if {0} {
	set cgparam 	[open "zzzz.inp" w]
	foreach i $aalist {
		set itype				$aatype($i)
		set icharge 		$aacharge($i)
		set ivdwEps			$aavdwEps($itype)
		set ivdwSgm			$aavdwSgm($itype)	
		puts $cgparam [format "%-5d %-5s %+2.5f   %+2.5f   %+2.5f" [expr $i+1] $itype $icharge $ivdwEps $ivdwSgm]
	}
	close $cgparam
}


foreach cgi $cglist {

	set ipatch				"$cgi "
	set cgiatom 			[atomselect top "index $cgi"]		
	set cgibondlist 	[lindex [$cgiatom getbonds] 0]

	foreach nbi $cgibondlist {			
		if {[lsearch $cglist $nbi]==-1} {
			set ipatch 		"$ipatch$nbi "
		}
	}

	set picharge			0.0
	foreach pi $ipatch {
		set type				$aatype($pi)
		set charge			$aacharge($pi)
		set	picharge		[expr $picharge+$charge]
	}

	set pivdwA				0.0
	set pivdwB				0.0
	#puts $ipatch
	for {set i 0} {$i<[llength $ipatch]} {incr i} {
		for {set j $i} {$j<[llength $ipatch]} {incr j} { 		
			set pi					[lindex $ipatch $i]
			set pj					[lindex $ipatch $j]
			set itype				$aatype($pi)
			set jtype				$aatype($pj)
			set ivdwEps			$aavdwEps($itype)
			set ivdwSgm			$aavdwSgm($itype)
			set jvdwEps			$aavdwEps($jtype)
			set jvdwSgm			$aavdwSgm($jtype)
			set Aij					[expr sqrt(abs($ivdwEps)*abs($jvdwEps))*pow($ivdwSgm+$jvdwSgm,12)]
			set Bij					[expr sqrt(abs($ivdwEps)*abs($jvdwEps))*pow($ivdwSgm+$jvdwSgm,6)]			
			set pivdwA			[expr $pivdwA+$Aij]
			set pivdwB			[expr $pivdwB+$Bij]
			#puts "$itype $aavdwEps(SOD)"
		}
	}
	
	set pivdwEps			[expr -$pivdwB*$pivdwB/$pivdwA]
	set pivdwSgm			[expr pow($pivdwA/$pivdwB,1/6.0)/2.0]


	set itype							$aatype($cgi)
	set cgcharge($cgi)		$picharge
	set cgvdwEps($itype)	$pivdwEps
	set cgvdwSgm($itype)	$pivdwSgm		

	mol delete $cgiatom
	
}

set cgparam 	[open "cgparam.inp" w]
foreach i $aalist {
	set itype				$aatype($i)
	
	if {[lsearch $cglist $i]!=-1} {
		set icharge 		$cgcharge($i)
		set ivdwEps			$cgvdwEps($itype)
		set ivdwSgm			$cgvdwSgm($itype)
		puts $cgparam [format "%-5d %+1.5f   %+1.5f   %+1.5f" [expr $i+1] $icharge $ivdwEps $ivdwSgm]
	} else {
		puts $cgparam [format "%-5d %+1.5f   %+1.5f   %+1.5f" [expr $i+1] 0.0 0.0 0.0]
	}

}
close $cgparam


mol delete top





