# REMARKS:
#
# in def generatrix(...), the line 
# if zc+r*sin(th) >= 0:  
# is used to avoid the foil piercing the surface. This has to be removed before running NFA
#
# flipud is used before interpolate.interp1d to have the incrasing ordering of the interpolated data. This is for backward compatibility

from pylab import *
from subprocess import Popen,STDOUT,PIPE
from re import sub
import fileinput
from scipy import interpolate
from bezier import bezier_curve
from smallestenclosingcircle import make_circle 

############################################
## FUNCTION TO DEFINE THE FOIL GENERATRIX ##
############################################
def generatrix(rr,tt,bezier_cc,bezier_aa,Ssurf,ds=0.05,figflg=False):
  '''This function defines the generatrix for the foil. The generatrix is
  composed of multiple circular sections specified with the rr and tt lists.
  The chord and angle of attach distribution along the curvilinear coordinate
  [0:1] produced by the generatrix are obtained by interpolating bezier curves
  defined by points contained in bezier_cc and bezier_aa respectively. 
  Inputs:
   - rr : the radii of curvatures for each section
   - tt : the arc (in radians) for each section
   - bezier_cc : a set of bezier curves control points for the chord
     distribution in the form [ [ x1,y1 ] , [ x2,y2] , .... , [xn,yn] ]
     defining a bezier curve for 0 < x < 1
   - bezier_aa : a set of bezier curves control points for the angle of attack
     distribution in the form [ [ x1,y1 ] , [ x2,y2] , .... , [xn,yn] ]
     defining a bezier curve for 0 < x < 1
   - figure=False : set to true if you want a plot of the generatrix
  '''
  rr = asarray(rr)
  tt = asarray(tt)
  bezier_cc = asarray(bezier_cc)
  bezier_aa = asarray(bezier_aa)

  smax =  sum(rr*tt)   								# maximum of spanwise coordinate
  s0 = rr[0]*tt[0]/smax ; bezier_cc[:,0] = bezier_cc[:,0]*(1-s0)+s0             # rescale bezier_cc so that the first section is constnat chord
  sscc,cc = bezier_curve(bezier_cc); sscc = flipud(sscc); cc = flipud(cc); cfun = interpolate.interp1d(sscc,cc,bounds_error=False,fill_value=1.)   	# interpolated chord function (first section is unitary chord)
  ssaa,aa = bezier_curve(bezier_aa); ssaa = flipud(ssaa); aa = flipud(aa); afun = interpolate.interp1d(ssaa,aa,bounds_error=False,fill_value=0.)   	# interpolated angle of attack function
  
  gen = [ ]   									# list containing generating quantities (one element per section)
  yzc = [ ]

  for i in xrange(0,len(tt)):
    if i == 0:
      s0 = 0 
      yc,zc = -rr[0],0
      ths = 0; 
    else:
      s0 = s0 + r*(the-ths)/smax
      yc,zc = yc+(r-rr[i])*cos(the),zc+(r-rr[i])*sin(the)
      ths = the; 

    r,dth = rr[i],tt[i]
    the = ths+dth
    nn = ceil(r*dth/ds)

    yzc.append( ( array([yc+r*cos(ths),yc,yc+r*cos(the) ]) , array([zc+r*sin(ths),zc,zc+r*sin(the) ]) ) )

    for th in linspace( ths , the , nn , endpoint=(i == (len(rr)-1)) ):
      s = min( s0 + r*(th-ths)/smax , 1. )
      gen.append( ( 0 , yc+r*cos(th) , zc+r*sin(th) , th+pi , cfun(s) , afun(s) , s ) ) 

  gen = asarray(gen)

  # computing the foil surface with current configuration and rescaling the chord to have the correct surface
  cc       = (gen[1:,4]+gen[:-1,4])/2   # average chord at each section
  ds       = diff(gen[:,6])*smax	# length of the section
  Stmp     = cc.dot(ds);		# the current surface of the foil
  gen[:,4] = gen[:,4]*Ssurf/Stmp	# the rescaled chord


  if figflg:
    ion()
    figure('spec',figsize=(10,10))
    xx = linspace(0,1,1000)

    ax = subplot2grid((3,2), (0,0)); 
    #ax.set_title('yz view'); 
    ax.axis('equal'); 
    ax.grid('on')
    x,y,r = make_circle(-gen[:,1:3]); ttt = linspace(0,2*pi,100); plot(x+r*cos(ttt),y+r*sin(ttt),'g:'); 
    ax.set_xlim(x-r,x+r)
    ax.set_ylim(y-r,y+r)
    [ plot(-x,-y,'b') for x,y in yzc]
    ax.plot(-gen[:,1],-gen[:,2],'r',lw=3); 
    ax.set_xlabel('y')
    ax.set_ylabel('z')
 
    ax = subplot2grid((3,2), (0,1)); 
    #ax.set_title('xz view'); 
    ax.plot(-gen[:,0],-gen[:,2],'r',lw=3); axis('equal'); grid('on')
    ax.yaxis.tick_right()
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.yaxis.set_label_position("right")
   
    ax = subplot2grid((3,2), (1,0)); 
    #ax.set_title('chord'); 
    ax.plot(xx,cfun(xx),'r',lw=3); 
    ax.plot([bezier_cc[-1,0],1],[bezier_cc[-1,1],0],'r',lw=3)
    ax.plot(bezier_cc[:,0],bezier_cc[:,1],'ob'); 
    ax.grid('on'); 
    ax.axis([-.05,1.1,0,2.])
    ax.plot(bezier_cc[:,0],bezier_cc[:,1])
    ax.set_xlabel('s/smax')
    ax.set_ylabel('chord')
    [ ax.annotate("%s"%i,xy=(bezier_cc[i][0],bezier_cc[i][1]),textcoords = "offset points",xytext = (+7,+7), ha = "left", va = "bottom") for i in xrange(len(bezier_cc))]

    ax = subplot2grid((3,2), (1,1)); 
    #ax.set_title('alpha'); 
    ax.plot(xx,afun(xx)*180/pi,'r',lw=3); 
    ax.plot(bezier_aa[:,0],bezier_aa[:,1]*180/pi,'ob'); 
    ax.grid('on'); 
    ax.axis([-.05,1.05,-5,5])
    ax.plot(bezier_aa[:,0],bezier_aa[:,1]*180/pi)
    ax.yaxis.tick_right()
    ax.set_xlabel('s/smax')
    ax.set_ylabel('alpha')
    [ ax.annotate("%s"%i,xy=(bezier_aa[i][0],bezier_aa[i][1]*180/pi),textcoords = "offset points",xytext = (0,+10), ha = "center", va = "bottom") for i in xrange(len(bezier_aa))]
    ax.yaxis.set_label_position("right")

   
    ax = subplot2grid((3,2), (2,0),colspan=2); ax.fill_between(smax*gen[:,6],gen[:,4]*0.25,-gen[:,4]*0.75); ax.axis('equal'); 
    stations = [sum(rr[:i]*tt[:i]) for i in xrange(len(rr))]
    stations.append(smax)
    ax.set_xticks(stations); ax.grid('on')
    text(.1,-.27,"Area = %.3f [m^2]"%Ssurf)
    ax.set_xlabel('s')
    ax.set_ylabel('x')
    draw(); 

  return gen
 
########################################
## FUNCTION TO WRITE THE PLY GEOMETRY ##
########################################
def writePLY(gen,afileName):
  """ 
  Writes a PLY file for the foil described in the gen array and with
  airfoil afileName (in XFOIL format, e.g. 'h105.dat'). 
  The output is saved in foil.ply
  """

  afile = loadtxt(afileName,skiprows=1)[:-1,:]   	# load wing section, remove the last point to avoid redundancy
 
  ngen , nsect = gen.shape[0] , afile.shape[0]		# count how many section and how many points per section

  ## START WRITING THE FILE ##
  with open('foil.ply', 'w') as f:

    ## HEADER OF THE PLY FILE ##
    header = '''ply
format ascii 1.0
comment made by Gianluca Meneghello (gianluca.meneghello@gmail.com)
comment this file is a Test Foil
element vertex
property float x
property float y
property float z
element face
property list uchar int vertex_index
end_header\n'''

    f.write(header)

    ## WRITE VERTICES ##
    nvertices = 0
    for xc,yc,zc,th,c,a,s in gen:
      a = -a                                                         # the angle of attack has to be changed because the geometry will be mirrored
      Rth= array([[1,0,0],[0,cos(th),-sin(th)],[0,sin(th),cos(th)]]) # rotation by the generatrix angle
      Ra = array([[cos(a),-sin(a),0],[sin(a),cos(a),0],[0,0,1]])     # rotation by the angle of attach (around the geneatrix line)
      R = Rth.dot(Ra)                                                # composite rotation
      for x,y in afile:
        x = x-0.25                                                   # center the profile on quarter chord before rotation and scaling
        x,y,z = R.dot(array([c*x,c*y,0]))                            # scale and rotate
        x , y , z = ( xc+x ) , -( yc+y ) , -( zc+z )                 # translation along the generatrix (minus signs if for mirroring)
        f.write("%f %f %f\n"%(x,y,z))                                # write the vertex to file
        nvertices = nvertices+1

    ## WRITE FACES (i.e. the connection between vertices) ##
    nfaces = 0

    # write the root lid (foil has to be closed)
    ns = 0 							     # index of the first point of the lid
    for i in xrange(nsect-1):
      ii,jj,kk = ns+i,ns+i+1,ns+nsect-1-i
      if (ii != jj and jj != kk and ii != kk):
        nfaces = nfaces+1
        f.write('3 %i %i %i\n'%(ii,jj,kk))

    # write the foil surface: we write two triangles at a time, corresponding to a square in the geometry 
    for j in xrange(ngen-1):
      ns = j*nsect
      for i in xrange(nsect):
        nfaces = nfaces+2
        f.write( '3 %i %i %i\n'%( ns+(i+1)%(nsect)       , ns+(i)%(nsect) , ns+nsect+(i+1)%(nsect) ) )
        f.write( '3 %i %i %i\n'%( ns+nsect+(i+1)%(nsect) , ns+(i)%(nsect) , ns+nsect+(i)%(nsect)   ) )

    # write the tip lid (foil has to be closed)
    ns = (nsect)*(ngen-1)    				 	     # index of the first point of the lid
    for i in xrange(nsect-1):
      ii,jj,kk = ns+i+1,ns+i,ns+nsect-1-i
      if (ii != jj and jj != kk and ii != kk):
        nfaces = nfaces+1
        f.write('3 %i %i %i\n'%(ii,jj,kk))

    f.write('\n')


  ## UPDATE NUMBER OF VERTICES AND NUMBER OF FACES IN THE HEADER ##
  for line in fileinput.input("foil.ply", inplace=True):
      print sub("element vertex","element vertex %i"%nvertices,line),
  for line in fileinput.input("foil.ply", inplace=True):
      print sub("element face","element face %i"%nfaces,line),

  return 0


########################################
## FUNCTION TO WRITE THE AVL GEOMETRY ##
########################################
def writeAVL(gen,afileName,Nchord=12,Nspan=101):
  """ 
  Writes an AVL input file for the foil described in the gen array and with
  airfoil afileName (in XFOIL format, e.g. 'h105.dat'). Nchord and
  Nspansise the number of horseshoe vortices in the chordwise and spanwise
  direction respectively (see AVL manual)
  """
  
  ## HEADER FOR THE AVL FILE ##
  header="""\
Test foil
#Mach
 0.0
#IYsym   IZsym   Zsym
 0       -1       0.0
#Sref    Cref    Bref
 1        1       1
#Xref    Yref    Zref
0.0      0.0     0.0
#CDp (already included in xfoil data (see CDCL entry), zero here)
0.0\n""" 

  ## HEADER FOR EACH SURFACE IN THE GEOMETRY ##
  ## (Nchord and Nspan have to be specified) ##
  surface_header = """\
SURFACE
foil
%i       1.0       %i      1.0  #Nchord  Cspace    Nspan   Sspace
COMPONENT 
1
SCALE
1.0  1.0  1.0
TRANSLATE
0.0  0.0  0.0
"""

  ## WRITING THE AVL FILE ##
  with open("foil.avl", "w") as text_file:   # write avl input geometry
    text_file.write(header)
    text_file.write(surface_header%(Nchord,Nspan))
    for xc,yc,zc,th,c,a,s in gen:
      text_file.write("SECTION\n")
      text_file.write("#Xle     Yle       Zle      Chord    Ainc      [ Nspan Sspace ]\n")
      text_file.write("%f %f %f %f %f\n"%(xc-0.25*c,-yc,-zc,c,a*180/pi)) #Xle    Yle    Zle     Chord   Ainc 
      text_file.write("CDCL\n0 0.00594117901910205 0.25 0.00543696634753934 0.5 0.0061644624025087\n")   # data obtained by fitting xfoil results in h105_n4.polar for cl [0.0:0.5]
      text_file.write("AFILE\n%s\n\n"%afileName)
    text_file.write("\n")
  text_file.close()

  return 0

###########################################################
## FUNCTION TO WRITE AVL FILE GIVEN A LIST OF PARAMETERS ##
###########################################################
def preAVL((th0,th1,th2,l0,l1,l2,S,xc,dyc,ctip,aroot,xa,ya,atip),Nchord=12,Nspan=101):
  ''' This function is to be used within an optimization loop. 
  To use this functions, set parameters from the input into the
  rr,tt,bezier_cc,bezier_aa arrays.  Define the rr,tt,bezier_cc,bezier_aa
  quantities by themselves otherwise'''

  afileName = 'h105.dat'                                                 # wing section file (XFOIL format)

  # optimum with these parameters, upwind working conditions.
  # Angles are Alpha =  2.95938 and Beta  =  0.70687 (for CL and CY as in upwind working condition)
  #p1 = array([0.1000,1.0505,0.2838,0.5000,0.2600,1.2312,0.3875]);    th0,th1,th2,l0,l1,l2,S = p1
  #p2 = array([ xc,dyc,ctip , aroot,xa,ya,atip ])

  tt        = array( [ th0,th1,th2 ] )                                                   # radii of curvatures
  ll        = array( [ l0 ,l1 ,l2 ] )                                                   # arch lenght of the foil section
  Ssurf     = S                                                                       # foil planform area
  rr        = ll/tt                                                                   # angle of the foil section
  bezier_cc = array([ [ 0 ,  1     ],[ xc ,  1   ],[ 1. , ctip+dyc ],[ 1. , ctip ] ]) # control points for chord distribution
  bezier_aa = array([ [ 0 ,  aroot ],[ xa ,   ya ],[ 1. , atip   ] ])                 # control points for angle of attack distribution

  gen = generatrix(rr,tt,bezier_cc,bezier_aa,Ssurf,ds=0.1 ,figflg=False)               # generate the generatrix line

  writeAVL(gen,afileName,Nchord,Nspan)                                      # write the AVL input file

  return gen

###########################
## PARSING OF AVL OUTPUT ##
###########################
def postAVL(lines):
  '''This function is to be used within an optimization loop. 
  It parses AVL output as contained in the lines input. 
  It looks for aerodynamic coefficients and returns the efficiency of the foil'''

  check = False
  # parse avl output
  for line in lines: 
    #print line,
    if "CXtot" in line: CXtot = float(line.split()[2])
    if "CYtot" in line: CYtot = float(line.split()[2])
    if "CZtot" in line: CZtot = float(line.split()[2])
    if "CLtot" in line: CLtot = float(line.split()[2])
    if "CDtot" in line: CDtot = float(line.split()[2])
    if "CDvis" in line: CDvis,CDind = float(line.split()[2]) , float(line.split()[5])
    if "CLff"  in line: 
      CLff,CDff   = float(line.split()[2]) , float(line.split()[5])
      check = True
    if "CYff"  in line: CYff,e      = float(line.split()[2]) , float(line.split()[5])

  if check:
    if (CLff+CYff) != 0:
      effteo = (CLff+CYff)
      eff = (CDvis+CDff)/(CLff + CYff)
      print '%5f %5f %5f %5f %5f %5f # AVL CLtot, CYtot, CDVis, CDff, e, eff'%(CLtot,CYtot,CDvis,CDff,e,eff)
    else:
      eff = 1.
      print '# AVL --- Lift is zero, I cannot compute inverse of efficiency'
   
  else:
    eff = 1.
    print '# AVL failed'



  return eff

#################################################
## TEST AVL CALLING SEQUENCE (BASED ON preAVL) ##
#################################################
def testAVL(params,afileName):
  gen = preAVL(params);
  writePLY(gen,afileName)
  p   = Popen("./avl foil.avl < avl_script | grep -e 'CXtot\|CYtot\|CZtot\|CLtot\|CDtot\|CDvis\|CLff\|CYff'", shell=True, stdout=PIPE, stderr=STDOUT)
  p.wait();  
  if p.returncode == 0:
    eff = postAVL( p.stdout.readlines() )
    return eff
  else:
    print '''
    An error as occurred during the avl computation. Return code is %i. 
    I do not have all the output. Try to run
        ./avl foil.avl < avl_script
    to see where it failed\n'''





##########
## MAIN ##
##########
if __name__ == "__main__":

  ########################################
  ## definition of construction section ##  
  ########################################
  afileName = 'h105.dat'                                          # wing section file (XFOIL format)

  tt = array( [0.1000,1.0505,0.2838] )
  ll = array( [0.5000,0.2600,1.2312] )
  rr = ll/tt

  bezier_cc = array([[0,1],[0,1.],[1.,1.],[1,1.]])               # control points for chord distribution
  bezier_aa = array([[0,0],[.5,0],[1,0]]) # control points for angle of attack distribution

  gen = generatrix(rr,tt,bezier_cc,bezier_aa,Ssurf=0.3875)                     # producing the generatrix

  ##########################################
  ## print geometry description to screen ##
  ##########################################
  #print '''FOIL: profile is %s'''%afileName
  #print 'radii [chord lengths]:',''.join('%5f  '%k for k in rr)
  #print 'arc   [radians]      :',''.join('%5f  '%k for k in tt)
  #print 'chord control points :',''.join('(%3f,%3f) '%tuple(k) for k in bezier_cc)
  #print 'alpha control points :',''.join('(%3f,%3f) '%tuple(k) for k in bezier_aa)

  ####################
  ## write PLY file ## 
  ####################
  print 'writing PLY file...',
  writePLY(gen,afileName)			# write ply file
  print 'done.'

  ####################
  ## write AVL file ## 
  ####################
  print 'writing AVL file...',
  writeAVL(gen,afileName,Nchord = 12 , Nspan=101) # write avl file
  print 'done.'

  ###############################################
  ## try interface with avl, if avl is present ##
  ###############################################
#  print 'running AVL example...'
#
#  # optimization for constant chord
#  #    th0    , th1    , th2    , l0     , l1     , l2     , Ssurf
#  #    0.1000 , 1.0505 , 0.2838 , 0.5000 , 0.2600 , 1.2312 , 0.3875 	optimum
#  #    0.1    , pi/4   , 0.1    , 0.5    , 0.1    , 0.5    , 0.2    	lower bound
#  #    pi/2   , pi/2   , pi/2   , 2      , 0.5    , 2      , 0.5    	upper bound
#  # sum(th) < pi/2				
#  # sum(ll) < 2
#  #
#  # optimization starting from the optimium for constant chord:
#  #    xc1 , yc2-yc3 , yc3    , a0     , xa1 , a1     , a2
#  #    0   , 1.6     , 0.8512 , 0.0052 , 0.2 , 0.0340 , -0.0541   optimum
#  # lower bound
#  # upper bound
#  # no additional constraints
#
#  params = array([0 , 0 , 1 , 0 , 0.5 , 0 , 0])

