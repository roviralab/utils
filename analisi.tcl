mol new 6jg1_g4g3g.top type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
#Loop over the folders
for {set i 1} {$i <= 15} {incr i} {
   mol addfile din_$i.nc type netcdf waitfor all
}
set totalframes [molinfo top get numframes]

#distances
set outfile [open dist.dat a+]
puts $outfile "#f \t h2oh_nuco2 \t c1_nuco1 \t o4_hab \t h3oh_hid207 \t h4oh_lys206 \t h4oh_asp95o1 \t h6oh_asp95o2 \t h3oh_4gb_o5 \t arg158_o2oh \t oab_h6oh \t tyr253_h6oh \t arg291_o6oh \t dih_conf \t ch2oh \t comw434_com4gb \t comw286_com4gb"
set h2oh  [[atomselect top "resname 0GB and name H2O"] get index]
set nuco2 [[atomselect top "resid 285 and name OD2"]   get index]
set c1    [[atomselect top "resname 0GB and name C1"]  get index]
set nuco1 [[atomselect top "resid 285 and name OD1"]   get index]
set o4    [[atomselect top "resname 4GB and name O4"]   get index]
set hab   [[atomselect top "resid 491 and name HE2"]   get index]
set hid207   [[atomselect top "resid 207 and name NE2"]   get index]
set h3oh  [[atomselect top "resname 0GB and name H3O"] get index]
set lys206   [[atomselect top "resid 206 and name NZ"]   get index]
set h4oh  [[atomselect top "resname 0GB and name H4O"] get index]
set asp95o1   [[atomselect top "resid 95 and name OD1"]   get index]
set asp95o2   [[atomselect top "resid 95 and name OD2"]   get index]
set h6oh  [[atomselect top "resname 0GB and name H6O"] get index]
set o5  [[atomselect top "resname 0GB and name O5"] get index]
set h3oh_4gb  [[atomselect top "resname 4GB and name H3O"] get index]
set arg158 [[atomselect top "resid 158 and name HH12"] get index]
set o2oh  [[atomselect top "resname 0GB and name O2"] get index]
set oab   [[atomselect top "resid 491 and name OE1"]   get index]
set tyr253   [[atomselect top "resid 253 and name OH"]   get index]
set arg291   [[atomselect top "resid 291 and name HH12"]   get index]
set o6oh4gb  [[atomselect top "resname 4GB and name O6"] get index]
set h6oh4gb  [[atomselect top "resname 4GB and name H6O"] get index]
set c2       [[atomselect top "resname 0GB and name C2"] get index]
set c5       [[atomselect top "resname 0GB and name C5"] get index]
set c4	     [[atomselect top "resname 0GB and name C4"] get index]
set c6	     [[atomselect top "resname 0GB and name C6"] get index]
set o6	     [[atomselect top "resname 0GB and name O6"] get index]

for {set f 1} {$f < $totalframes} {incr f} {
      set comw434 [measure center [atomselect top "sidechain and resid 434" frame $f] weight mass]
      set comw286 [measure center [atomselect top "sidechain and resid 286" frame $f] weight mass]
      set com4gb [measure center [atomselect top "resname 4GB" frame $f] weight mass]
      set h2oh_nuco2    [measure bond "$h2oh $nuco2"  frame $f]
      set c1_nuco1      [measure bond "$c1 $nuco1"  frame $f]
      set o4_hab        [measure bond "$o4 $hab"  frame $f]
      set h3oh_hid207        [measure bond "$h3oh $hid207"  frame $f]
      set h4oh_lys206        [measure bond "$h4oh $lys206"  frame $f]
      set h4oh_asp95o1        [measure bond "$h4oh $asp95o1"  frame $f]
      set h6oh_asp95o2        [measure bond "$h6oh $asp95o2"  frame $f]
      set h3oh_4gb_o5        [measure bond "$h3oh_4gb $o5"  frame $f]
      set arg158_o2oh        [measure bond "$arg158 $o2oh"  frame $f]
      set oab_h6oh4gb        [measure bond "$oab $h6oh4gb"  frame $f]
      set tyr253_h6oh4gb        [measure bond "$tyr253 $h6oh4gb"  frame $f]
      set arg291_o6oh4gb        [measure bond "$arg291 $o6oh4gb"  frame $f]
      set dih_conf		[measure dihed "$c2 $c1 $o5 $c5" frame $f]
      set ch2oh          [measure dihed "$c4 $c5 $c6 $o6" frame $f]
      set comw434_com4gb [veclength [vecsub $comw434 $com4gb]]
      set comw286_com4gb [veclength [vecsub $comw286 $com4gb]]
      puts $outfile "$f \t $h2oh_nuco2 \t $c1_nuco1 \t $o4_hab \t $h3oh_hid207 \t $h4oh_lys206 \t $h4oh_asp95o1 \t $h6oh_asp95o2 \t $h3oh_4gb_o5 \t $arg158_o2oh \t $oab_h6oh4gb \t $tyr253_h6oh4gb \t $arg291_o6oh4gb \t $dih_conf \t $ch2oh \t $comw434_com4gb \t $comw286_com4gb"
}

close $outfile

#rmsd
set outfile [open rmsd.dat a+]

set wholeSystem [atomselect top "all"]

set protBB    [atomselect top "(protein and backbone)"]
set refprotBB [atomselect top "(protein and backbone)" frame 0]

set substrate [atomselect top "resname 0GB 4GB 3GB ROH"]
set refsubstrate [atomselect top "resname 0GB 4GB 3GB ROH" frame 0]

set subs_0gb [atomselect top "resname 0GB"]
set refsubs_0gb [atomselect top "resname 0GB" frame 0]

set subs_4gb [atomselect top "resname 4GB"]
set refsubs_4gb [atomselect top "resname 4GB" frame 0]

set subs_3gb [atomselect top "resname 3GB ROH"]
set refsubs_3gb [atomselect top "resname 3GB ROH" frame 0]

set active_site [atomselect top "resid 58 99 148 162 210 211 254 257 289 290 291 295 320 434 438 495"]
set refactive_site [atomselect top "resid 58 99 148 162 210 211 254 257 289 290 291 295 320 434 438 495" frame 0]

set trp_pos1 [atomselect top "resid 438 290"]
set reftrp_pos1 [atomselect top "resid 438 290" frame 0]

for {set f 1} {$f < $totalframes} {incr f} {
   $protBB frame $f
   set align_matrix [measure fit $protBB $refprotBB]
   $wholeSystem frame $f
   $wholeSystem move $align_matrix

   $substrate frame $f
   $subs_0gb frame $f
   $subs_4gb frame $f
   $subs_3gb frame $f
   $active_site frame $f
   $trp_pos1 frame $f

   set rmsdprotBB [measure rmsd $refprotBB $protBB]
   set rmsdsubstrate [measure rmsd $refsubstrate $substrate]
   set rmsdsubs_0gb [measure rmsd $refsubs_0gb $subs_0gb]
   set rmsdsubs_4gb [measure rmsd $refsubs_4gb $subs_4gb]
   set rmsdsubs_3gb [measure rmsd $refsubs_3gb $subs_3gb]
   set rmsdactive_site [measure rmsd $refactive_site $active_site]
   set rmsdtrp_pos1 [measure rmsd $reftrp_pos1 $trp_pos1]

   puts $outfile "$f \t $rmsdprotBB \t $rmsdsubstrate \t $rmsdsubs_0gb \t $rmsdsubs_4gb \t $rmsdsubs_3gb \t $rmsdactive_site \t $rmsdtrp_pos1"
}

close $outfile

set last_f [expr ([molinfo top get numframes] - 1)]
animate write dcd substrate.dcd beg 1 end $last_f waitfor all sel $substrate top
animate write pdb substrate.pdb beg 1 end 1 waitfor all sel $substrate top
#
quit
