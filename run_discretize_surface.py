vdw_alteration = 1.8
vdw_minimum = 1.2
vdw_maximum = 4.2
x_coord = 15.243419647216797
y_coord = 8.344330787658691
z_coord = 4.119806289672852
radius = 8.0
basename = '3w32-benzene-PARTIAL-DATA-'
count = 0

cmd.pseudoatom('reference',pos="[%.2f,%.2f,%.2f]"%(x_coord,y_coord,z_coord))
cmd.alter(selection='resname ALP',expression="vdw = q - %.2f"%vdw_alteration)
cmd.color(color='magenta',selection='all')
cmd.set(name='surface_quality',value=0) # For fingerprint generation which only references the geometric center, consider making this value higher
cmd.set(name='surface_type',value=1)
objects = []
objects = cmd.get_object_list('all')
for obj in objects:
        count += 1
        print(f"Operating on: {obj}")
        cmd.hide('all')
        cmd.create(name='hydro_coloring',selection="(bychain resname ALP) within %.2f of reference"%radius)
        cmd.create(name='vdw_coloring',selection="(bychain resname ALP) within %.2f of reference"%radius)
        cmd.color(color='white',selection='hydro_coloring and elem O')
        cmd.color(color='black',selection='hydro_coloring and elem C')
        cmd.spectrum(expression='vdw',palette='rainbow',selection='vdw_coloring',minimum="%.2f"%vdw_minimum,maximum="%.2f"%vdw_maximum)
        cmd.show(representation='surface',selection='hydro_coloring')
        cmd.show(representation='surface',selection='vdw_coloring')
        basename + str(count) + ".wrl"
        cmd.save(filename="%s"%output_filename,selection='enabled')
