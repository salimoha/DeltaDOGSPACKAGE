from pylab import *
import fileinput,re

# parameters
lshaft = 5.  # [C]
lwing  = 3.5 # [C]
gwing  = 15*pi/180  # [deg]
r = 0.5      # [C]
c = 1	     # [-]

filename = 'workfile.ply'

# load wing section
h105 = loadtxt('h105.dat',skiprows=1)[:-1,:]   # remove the last point to avoid redundant point



# generating line y,z,th,c
gen = [ (0,0,0,c) ]
[ gen.append( ( r+r*cos(th) , -lshaft+r*sin(th) , th+pi , c ) ) for th in linspace(pi,3*pi/2+gwing,10) ] 
gen.append( ( r+lwing*cos(gwing),-r-lshaft+lwing*sin(gwing),pi/2,0.5) )


a = asarray(gen)
plot(a[:,0],a[:,1])
show()

nsect = h105.shape[0]
ngen = a.shape[0]

print "nsect,ngen = %i %i"%(nsect,ngen)

with open(filename, 'w') as f:

  header = '''ply
format ascii 1.0
comment made by Gianluca Meneghello
comment this file is a Test Foil
element vertex
property float x
property float y
property float z
element face
property list uchar int vertex_index
end_header\n'''

  #f.write(header%(nsect*ngen,(nsect-1)*(ngen-1)))
  f.write(header)


  # write vertices
  nvertices = 0
  for yc,zc,th,c in a:
    for x,y in h105:
      nvertices = nvertices+1
      f.write("%f %f %f\n"%(c*x-c*0.25,yc+c*y*cos(th),zc+c*y*sin(th)))

  # write faces
  nfaces = 0

  # root lid
  ns = 0
  for i in xrange(nsect-1):
    ii,jj,kk = ns+i,ns+i+1,ns+nsect-1-i
    if (ii != jj and jj != kk and ii != kk):
      nfaces = nfaces+1
      f.write('3 %i %i %i\n'%(ii,jj,kk))

  # foil surface
  for j in xrange(ngen-1):
    ns = j*nsect
    for i in xrange(nsect-1):
      nfaces = nfaces+2
      #f.write('4 %i %i %i %i\n'%(ns+i,ns+i+1,ns+nsect+i,ns+nsect+i+1))
      f.write( '3 %i %i %i\n'%( ns+(i+1)%(nsect-1)       , ns+(i)%(nsect-1) , ns+nsect+(i+1)%(nsect-1) ) )
      f.write( '3 %i %i %i\n'%( ns+nsect+(i+1)%(nsect-1) , ns+(i)%(nsect-1) , ns+nsect+(i)%(nsect-1)   ) )

  # tip lid
  ns = (nsect)*(ngen-1)
  for i in xrange(nsect-1):
    ii,jj,kk = ns+i+1,ns+i,ns+nsect-1-i
    if (ii != jj and jj != kk and ii != kk):
      nfaces = nfaces+1
      f.write('3 %i %i %i\n'%(ii,jj,kk))

  f.write('\n')


# update number of vertices and number of faces 
for line in fileinput.input(filename, inplace=True):
    print re.sub("element vertex","element vertex %i"%nvertices,line),
for line in fileinput.input(filename, inplace=True):
    print re.sub("element face","element face %i"%nfaces,line),

