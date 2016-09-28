## bezier curves for chord and angle of attack distribution

from pylab import *
import fileinput,re
from scipy import interpolate
from bezier import bezier_curve

# parameters
#lshaft = 5.  # [C]
#lwing  = 3.5 # [C]
#gwing  = 15*pi/180  # [deg]
#r = 0.5      # [C]
#c = 1	     # [-]

filename = 'workfile.ply'  # output ply file

h105 = loadtxt('h105.dat',skiprows=1)[:-1,:]   # load wing section, remove the last point to avoid redundant point

rr = array( [ 5    , .5   , 3    , .5   ] )
tt = array( [ pi/4 , pi/4 , pi/4 , pi/3 ] )
smax =  sum(rr*tt)   # maximum of spanwise coordinate

# quantities specified at end of sections
#ss = zeros(len(rr)+1)   # spanwise coordinate at chord locations
#for i in xrange(len(rr)):
#  ss[i+1] = ss[i] + rr[i]*tt[i]
#
#ss = ss/smax
#sscc = ssaa = ss.copy()
#cc = array( [ 1 , 1. , .5 , .5 , .1 ] )
#aa = array( [ 0 , 0  , 0  , 0  , 0  ] )*pi/180
#ww = array( [ 1 , .5 , .5 , .5 , 1  ] )  # weigths for interpolation0

# quantities specified as beizer curves
bezier_cc = array([[0,1],[1.,2],[1,.1]])  # control points for chord distribution
sscc,cc = bezier_curve(bezier_cc)
bezier_aa = array([[0,0],[.5,0],[1,0]])	# control points for angle of attack distribution
ssaa,aa = bezier_curve(bezier_aa)

cfun = interpolate.interp1d(sscc,cc)   # interpolated chord funcion
afun = interpolate.interp1d(ssaa,aa)   # interpolated chord funcion

# figure
figure('spec')
xx = linspace(0,1,1000)
cf = cfun(xx)
subplot(211), title('angle of attack')
af = afun(xx)
plot(xx,af,'b',label='$\alpha$')
plot(ssaa,aa,'.b')
plot(bezier_aa[:,0],bezier_aa[:,1],'ob')
autoscale(tight=False)
subplot(212), title('chord')
plot(xx,cf,'r',label='chord')
plot(sscc,cc,'r.')
plot(bezier_cc[:,0],bezier_cc[:,1],'or')
show()

# generatice
gen = [ ]   # list containing generating quantities

s0 = 0 
yc,zc = -rr[0],0
ths = 0; 
r,dth = rr[0],tt[0]
the = ths+dth
for th in linspace(ths,the,int(r*dth/0.05),endpoint=False):
  s = s0 + r*(th-ths)/smax
  gen.append( ( -0.25*cfun(s) , yc+r*cos(th) , zc+r*sin(th) , th+pi , cfun(s) , -afun(s) ) ) 

for i in xrange(1,len(tt)):
  s0 = s0 + r*(the-ths)/smax
  yc,zc = yc+(r-rr[i])*cos(the),zc+(r-rr[i])*sin(the)
  ths = the; 
  r,dth = rr[i],tt[i]
  the = ths+dth
  for th in linspace( ths , the , ceil(r*dth/0.05) , endpoint=(i == (len(rr)-1)) ):
    s = s0 + r*(th-ths)/smax
    gen.append( ( -0.25*cfun(s) , yc+r*cos(th) , zc+r*sin(th) , th+pi , cfun(s) , -afun(s) ) ) 
# end of generating line

# generate the ply file
gen = asarray(gen)

nsect = h105.shape[0]
ngen = gen.shape[0]

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
  for xc,yc,zc,th,c,a in gen:
    Rth= array([[1,0,0],[0,cos(th),-sin(th)],[0,sin(th),cos(th)]])
    Ra = array([[cos(a),-sin(a),0],[sin(a),cos(a),0],[0,0,1]])
    R = Rth.dot(Ra)
    for x,y in h105:
      x = xc+x
      #x , y , z = c*x , c*y*cos(th) , c*y*sin(th)  # rotation for the spanwise axis
      #x , y , z = ( xc+x ) , -( yc+y ) , -( zc+z )     # translation
      x,y,z = R.dot(array([c*x,c*y,0]))
      x , y , z = ( x ) , -( yc+y ) , -( zc+z )     # translation
      nvertices = nvertices+1
      f.write("%f %f %f\n"%(x,y,z))

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
    for i in xrange(nsect):
      nfaces = nfaces+2
      #f.write('4 %i %i %i %i\n'%(ns+i,ns+i+1,ns+nsect+i,ns+nsect+i+1))
      f.write( '3 %i %i %i\n'%( ns+(i+1)%(nsect)       , ns+(i)%(nsect) , ns+nsect+(i+1)%(nsect) ) )
      f.write( '3 %i %i %i\n'%( ns+nsect+(i+1)%(nsect) , ns+(i)%(nsect) , ns+nsect+(i)%(nsect)   ) )

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

