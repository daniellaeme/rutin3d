# load complex
load rutin_mpges1_complex.pdb, complex
hide everything

# separate selections
select protein, polymer.protein
select gsh, resn GSH
select rutin, resn UNL

# show protein as translucent grey cartoon to give context without visual noise
show cartoon, protein
color gray90, protein
set cartoon_transparency, 0.4, protein

# show Glutathione cofactor as sticks with lightorange carbons (matches PLIP)
show sticks, gsh
color lightorange, gsh and name C*
set stick_radius, 0.22, gsh

# show the docked Rutin ligand with high-contrast, vibrant forest-green carbons
show sticks, rutin
color forest, rutin and name C*
set stick_radius, 0.32, rutin

# select interacting residues identified in PLIP (Glu66, Leu69, Arg70, Arg73, Asn74, Glu77, Arg126)
select plip_residues, resi 66+69+70+73+74+77+126 and protein
show sticks, plip_residues
color marine, plip_residues and name C*
set stick_radius, 0.22, plip_residues

# label the C-alpha of interacting residues
label plip_residues and name CA, "%s-%s" % (resn, resi)
set label_size, 10
set label_color, black
set label_font_id, 7

# hide all standard non-polar hydrogens to eliminate spatial clutter
hide sticks, element H

# draw thin yellow dashed lines for non-covalent polar contacts
# Arg126 to Rutin Hydrogen Bonds
distance hb_arg126_a, (resi 126 and name NH1), (rutin and id 1342), 3.2, mode=2
distance hb_arg126_b, (resi 126 and name NH2), (rutin and id 1342), 3.0, mode=2

# Arg73 to Rutin Hydrogen Bond

distance hb_arg73, (resi 73 and name NH1), (rutin and id 1312), 3.7, mode=2

# Arg70 to Rutin Hydrogen Bond

distance hb_arg70, (resi 70 and name NH1), (rutin and id 1344), 3.2, mode=2

# Glu66 to Rutin Hydrogen Bond

distance hb_glu66, (resi 66 and name OE1), (rutin and id 1344), 3.7, mode=2

color yellow, hb_*
set dash_width, 2.0
set dash_gap, 0.2
hide labels, hb_*

# Configure professional lighting, background and view rendering
bg_color white
set ray_opaque_background, on
set depth_cue, on
set spec_count, 1
set shininess, 10

# Center camera on the interfacial active site binding cleft
zoom gsh or rutin, 8
center gsh or rutin

print "PyMOL Visual Setup Complete. Ready for ray 1200, 900!"

# Render the scene with professional lighting
# set ray_trace_mode, 1
# ray 1200, 900
# png ../results/rutin_docking_final.png