vdw_alteration = 1.8
vdw_minimum = 1.2
vdw_maximum = 4.2
x_coord = -19.91500473022461
y_coord = 5.902616500854492
z_coord = 9.151772499084473
radius = 5.0
output_filename = 'test.wrl'

cmd.pseudoatom('reference',pos="[%.2f,%.2f,%.2f]"%(x_coord,y_coord,z_coord))
cmd.alter(selection='resname ALP',expression="vdw = q - %.2f"%vdw_alteration)
cmd.color(color='magenta',selection='all')
cmd.set(name='surface_quality',value=0) # For fingerprint generation which only references the geometric center, consider making this value higher
cmd.set(name='surface_type',value=1)
objects = []
objects = cmd.get_object_list('all')
for obj in objects:
        print(f"Operating on: {obj}")
        cmd.hide('all')
        cmd.create(name='hydro_coloring',selection="(bychain resname ALP) within %.2f of reference"%radius)
        cmd.create(name='vdw_coloring',selection="(bychain resname ALP) within %.2f of reference"%radius)
        cmd.color(color='white',selection='hydro_coloring and elem O')
        cmd.color(color='black',selection='hydro_coloring and elem C')
        cmd.spectrum(expression='vdw',palette='rainbow',selection='vdw_coloring',minimum="%.2f"%vdw_minimum,maximum="%.2f"%vdw_maximum)
        cmd.show(representation='surface',selection='hydro_coloring')
        cmd.show(representation='surface',selection='vdw_coloring')
        cmd.save(filename="%s"%output_filename,selection='enabled')
