from pylab import *

# parameters
lshaft = 5.  # [C]
lwing  = 3.5 # [C]
gwing  = 15*pi/180  # [deg]
r = 0.5      # [C]
c = 1	     # [-]

# load wing section
h105 = loadtxt('h105.dat',skiprows=1)



# generating line y,z,th,c
gen = [ (0,0,0,c) ]
[ gen.append( ( r+r*cos(th) , -lshaft+r*sin(th) , th+pi , c ) ) for th in linspace(pi,3*pi/2+gwing,10) ] 
gen.append( ( r+lwing*cos(gwing),-r-lshaft+lwing*sin(gwing),pi/2,0.5) )

a = asarray(gen)

nsect = h105.shape[0]
ngen = a.shape[0]


with open('workfile.ply', 'w') as f:

  header = '''ply
format ascii 1.0           
comment made by Gianluca Meneghello
comment this file is a Test Foil
element vertex %i
property float x           
property float y           
property float z           
element face %i             
property list uchar int vertex_index 
end_header\n'''

  #f.write(header%(nsect*ngen,(nsect-1)*(ngen-1)))
  f.write(header%(nsect*ngen,(nsect-1)*(ngen-1)*2))


  # write vertices
  for yc,zc,th,c in a:
    for x,y in h105:
      f.write("%f %f %f\n"%(c*x,yc+c*y*cos(th),zc+c*y*sin(th)))

  # write faces

  for j in xrange(ngen-1):
    ns = j*nsect
    for i in xrange(nsect-1):
      #f.write('4 %i %i %i %i\n'%(ns+i,ns+i+1,ns+nsect+i,ns+nsect+i+1))
      f.write('3 %i %i %i\n'%(ns+i,ns+i+1,ns+nsect+i+1))
      f.write('3 %i %i %i\n'%(ns+i+1,ns+nsect+i+1,ns+nsect+i))

  f.write('\n')



