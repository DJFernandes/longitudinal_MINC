import pyminc.volumes.factory as pm
import numpy as np
import os
import time
from scipy import interpolate
import sys
import argparse
import bisect

def flcslist(s):
    try:
        sspl = map(float, s.split(','))
        return sspl
    except:
        raise argparse.ArgumentTypeError("MUST BE A COMMA SEPERATED LIST OF NUMBERS")

def strcslist(s):
    try:
        sspl = s.split(',')
        return sspl
    except:
        raise argparse.ArgumentTypeError("MUST BE A COMMA SEPERATED LIST OF FILES")

parser = argparse.ArgumentParser(description='Spline Interpolate xfm files')

parser.add_argument('timepoints', 
                    help="corresponding timepoints of the xfmfiles (comma seperated list). The last timepoint is the reference timepoint.", 
                    type=flcslist, nargs=1)

parser.add_argument('xfmfiles', 
                    help="xfmfiles to be splined (comma seperated list)", 
                    type=strcslist, nargs=1)

parser.add_argument('interp_time', 
                    help="Timepoints for which to find interpolated xfm (comma seperated list)", 
                    type=flcslist, nargs=1)

parser.add_argument('outdir', 
                    help="Output Directory", 
                    type=str, nargs=1)

args = parser.parse_args()

tknots=args.timepoints[0]
xfms=args.xfmfiles[0]
tinterp=args.interp_time[0]
outdir=args.outdir[0]

if len(xfms)!=(len(tknots)-1):
   print "number of xfmfiles must be ONE LESS than number of timepoints"
   print "THe Last Timepoint is the reference timepoint"
   sys.exit()


comm="%"+time.strftime("%d/%m/%Y")+" "+time.strftime("%H:%M:%S")+">>> "+'python xfmspline.py "'+",".join(map(str,tknots))+'" "'+",".join(xfms)+'" "'+",".join(map(str,tinterp))+'" '+outdir

#tknots=[65,36,29]
#xfms=["diff_transforms/p65_to_p36/p65_to_p36_CM.xfm","diff_transforms/p36_to_p29/p36_to_p29_CM.xfm"]
#tinterp=[29,33]
#outdir="splinetest"

if not os.path.exists(outdir):
    os.makedirs(outdir)

outdir=outdir+"/"

#for leading zeroes
def zpad(val, n):
    bits = val.split('.')
    return "%s.%s" % (bits[0].zfill(n), bits[1])

# Function that applies xfm-notation 3x4 matrix (Rotation | Translation) (arr) to 3-component vector a
def afftrans(a,arr): 
   return np.dot(arr[:,0:3],(a+arr[:,3]))

#Read the affine component of xfm
def read_affine_xfm(xfmfile):
   with open(xfmfile) as f:
      content = f.readlines()
   startline=content.index("Linear_Transform =\n")+1 ; stopline=startline+3
   content[startline-1]="EXTRACTED_Linear_Transform =\n"
   content[stopline-1]=content[stopline-1].replace(';','')
   xfmarr=np.genfromtxt(content[startline:stopline])
   return xfmarr


#this function reads xfmfiles are outputs a vector field encoding the transformation
def xfm2vectorField(xfmfile):
   with open(xfmfile) as f:
     content = f.readlines()   

   #Extract All of Grid Transforms
   search='Displacement_Volume = '
   linestxt=filter(lambda x:search in x, content)
   
   NUMtrans=len(linestxt)   #Number of grid transforms
   
   gridfilesbase=[linetxt.replace(search,"").replace(';\n',"") for linetxt in linestxt]
   xfmgrids=gridfilesbase
   if len(os.path.dirname(xfmfile))!=0:
      xfmgrids=[os.path.dirname(xfmfile)+'/'+gridfilebase for gridfilebase in gridfilesbase]
   
   xfmarrs = [None] * len(linestxt)
   for i in range(NUMtrans):
      startline=content.index("Linear_Transform =\n")+1 ; stopline=startline+3
      content[startline-1]="EXTRACTED_Linear_Transform =\n"
      content[stopline-1]=content[stopline-1].replace(';','')
      xfmarr=np.genfromtxt(content[startline:stopline])
      xfmarrs[i]=xfmarr
 
   for i in range(NUMtrans):
      print "concatenating transform "+str(i+1)+" of "+str(NUMtrans)
      if i==0:
        v=np.apply_along_axis(afftrans, 0, pm.volumeFromFile(xfmgrids[i]).data ,xfmarrs[i])
        continue
      v=v+np.apply_along_axis(afftrans, 0, pm.volumeFromFile(xfmgrids[i]).data ,xfmarrs[i])

   return v

#this function reads xfmfile. Returns first nlin grid
def xfm2vfirstnlingrid(xfmfile):
   with open(xfmfile) as f:
     content = f.readlines() 
   search='Displacement_Volume = '
   linetxt=filter(lambda x:search in x, content)[0]
   linenum=content.index(linetxt)
   gridfilebase=linetxt.replace(search,"").replace(';\n',"")
   if len(os.path.dirname(xfmfile))==0:
    xfmgrid=gridfilebase
   else:
    xfmgrid=os.path.dirname(xfmfile)+'/'+gridfilebase
   return xfmgrid


#replace last number of occurances of string
def rreplace(s, old, new):
   li = s.rsplit(old, 1)
   return new.join(li)

def xfmwrite(xfmfile,likefile,outaffinearr,outnlingrid,comment):
   outgridfile=rreplace(os.path.basename(xfmfile),'.xfm','_grid_0.mnc')
   if len(os.path.dirname(xfmfile))!=0:
      outgridfilefull=os.path.dirname(xfmfile)+"/"+outgridfile
   else:
      outgridfilefull=outgridfile

   outgrid = pm.volumeLikeFile(likefile, outgridfilefull)
   outgrid.data=outnlingrid
   outgrid.writeFile()
   outgrid.closeVolume()

   outaffinearrlist=map(str,outaffinearr.tolist())

   newcontent=['MNI Transform File']
   newcontent.append("%"+comment)
   newcontent.append("")
   newcontent.append('Transform_Type = Grid_Transform;')
   newcontent.append('Displacement_Volume = '+outgridfile+';')
   newcontent.append('Transform_Type = Linear;')
   newcontent.append('Linear_Transform =')
   newcontent.append(outaffinearrlist[0].replace(",","").replace("[","").replace("]",""))
   newcontent.append(outaffinearrlist[1].replace(",","").replace("[","").replace("]",""))
   newcontent.append(outaffinearrlist[2].replace(",","").replace("[","").replace("]","")+";")

   with open(xfmfile,'w') as f:
       for item in newcontent:
            f.write(item+"\n")

# Read all xfms
numKNOTS=len(xfms)

for i in range(numKNOTS):
   print "reading "+xfms[i]
   v=xfm2vectorField(xfms[i])
   
   #Read first volume and initialize empty matrix (flattened and one for each vector component)
   if i==0:
      voldims=(np.prod(v.shape[1:4]),numKNOTS+1)
      x_series=np.zeros(voldims)
      y_series=np.zeros(voldims)
      z_series=np.zeros(voldims)

   #broadcast vector field into empty matrices
   x_series[:,i]=np.ravel(v[0,:,:,:])
   y_series[:,i]=np.ravel(v[1,:,:,:])
   z_series[:,i]=np.ravel(v[2,:,:,:])

# The final timepoint is the reference timepoint
# Thus, its vector transformation is IDENTITY (AKA 0 everywhere)

print x_series.shape

def myInterpFunc(yknots):
   return interpolate.interp1d(tknots, yknots,kind='cubic')(tinterp)

print "interpolating x components"
interp_x_series=np.apply_along_axis(myInterpFunc, 1, x_series)
print "interpolating y components"
interp_y_series=np.apply_along_axis(myInterpFunc, 1, y_series)
print "interpolating z components"
interp_z_series=np.apply_along_axis(myInterpFunc, 1, z_series)

likefile=xfm2vfirstnlingrid(xfms[-1])

voldims=tuple(list(pm.volumeFromFile(likefile).data.shape)[1:4])


def get_reference(x):
    stk=sorted(tknots)
    ntrv = bisect.bisect_right(stk,x)
    if ntrv>(len(tknots)-1):
       ntrv=len(tknots)-1
    return stk[ntrv]


for i in range(len(tinterp)):
   vectransgrid=np.empty(tuple([3]+list(voldims)))
   vectransgrid[0,:,:,:]=np.reshape(interp_x_series[:,i],voldims)
   vectransgrid[1,:,:,:]=np.reshape(interp_y_series[:,i],voldims)
   vectransgrid[2,:,:,:]=np.reshape(interp_z_series[:,i],voldims)
   
   outxfmname=outdir+'p'+zpad(str(float(get_reference(tinterp[i]))),2).replace(".","d")+"_to_"+'p'+zpad(str(float(tinterp[i])),2).replace(".","d")+".xfm"

   xfmwrite(outxfmname,likefile,np.asarray([[1,0,0,0],[0,1,0,0],[0,0,1,0]]),vectransgrid,comm)


