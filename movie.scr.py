t0 = 101
t1 = 384
dt = 1
pz = '0.05'

ofl = open('movie', 'w')

ofl.write('''#### gnu.wave ####

se xr[-30:30]
se yr[-30:30]
se zr[0:50]

#set terminal pngcairo size 1500,1500 enhanced font 'Verdana,9'\n\n''')

for i in range(t0, t1 + 1, dt):
    #ofl.write('''se output "Tracer_p_%i.png"\n'''%(i))
    ofl.write('''sp "fort.%i" u 3:4:5 w l\npause %s\n''' % (i, pz))

ofl.close()

# ffmpeg -framerate 20 -pattern_type glob -i 'Tracer_p_*.png' -c:v libx264 -pix_fmt yuv420p Tracer_p_movie.mp4
