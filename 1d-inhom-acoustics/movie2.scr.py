t0 = 110
t1 = 300
dt = 10
pz = '0.01'

ofl = open('movie2', 'w')

ofl.write('''#### gnu.wave ####

set grid

se xlabel "x"
se ylabel "y"

se xr [0:2*pi]
se yr [0:2*pi]

p "Initial_Grid_Data.dat" u 2:3:7 w ima notitle
pause 0.1

#set terminal pngcairo size 1500,1500 enhanced font 'Verdana,9'\n\n''')

for i in range(t0, t1 + 1, dt):
    #ofl.write('''se output "mag_field_%i.png"\n'''%(i))
    ofl.write('''p "./fort.%i" u 2:3:7 w ima notitle\npause %s\n''' % (i, pz))

ofl.close()

# ffmpeg -framerate 20 -pattern_type glob -i 'mag_field_*.png' -c:v libx264 -pix_fmt yuv420p mag_field_movie.mp4
