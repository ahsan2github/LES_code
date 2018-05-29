#!/usr/bin/env python3
import numpy as np
import os.path
import os, string, ast
import struct
import subprocess as sp

varDic=[]
vars=[  'Von','Sc','Co','nnn','pi','model','averaging','fgr',\
        'tfr','mom_nodes','scl_nodes','FFT_FLAG','nsteps','hprocs','dt',\
        'Nx','Ny','Nz','aNx','l_r','z_i','l_z','sponge','z_d',\
        'rlx_time','Ugal','Vgal','Ug','Vg',"u_star",\
        'press_cor','press_step','u_avg','S_advec','scalarcount',\
        'f_c','M_advec','scalarFlags','nhrs','startUTC',\
        'surfaceFlags','surfaceFluxes','scalarScales',\
        'dsdtHomogeneous','scalarAdvecFlags','inversion',\
        'theta_0','Ri_flag',"c_count",'p_count','cs_count',\
        'nframe','sframe','framestep','n_tslevels','timeseriesfreq',\
        'inituA','initsA','ruler','npart','start_release',\
        'partstep','freq','nr','c_flag','Cd',\
        'pressureScale','densityAir','Cp_air','densityWater','soilLevels','zt',\
        'latentHeatWater','heatCapWater','waterGasConst','moistureCriteria',\
        'temperatureCriteria','tempFluxCriteria','maxFluxIterations',\
        'maxTempIterations','convFactor','endConstSEB','updateFreqSEB',\
        'integrateSoilDiffFreq','albedoFlag','radiationFlag','stepsPerRadVal',\
        'SB_constant','solarIrradiance','latitude','longitude','day',\
        'emissivity','evolflag','amr_crit','min','max','amr_cnt', 'filtLevelsLambda']


def readInput(file_in, vars = vars):
    with open(file_in, mode='r') as f:
        rbound=3
        dicvar={}
        lines = f.readlines()
        #print lines
        for x in reversed(lines):
                str1=x.replace(',','')        # str1 is a string, remove all commas(,)
                lst=str1.split()   # this is a list NOT A STRING
                for item in reversed(vars):
                    vkey=item
                    tf=lst[0:rbound].count(vkey)
                    if tf>0 and vkey=='Von':
                            dicvar.update({'kappa':float(str1.split()[0].replace('d','0'))})
                            vars=[item for item in vars if item!=vkey]
                    elif tf>0 and vkey=='min':
                            dicvar.update({'minGrid':float(str1.split()[0].replace('d','0'))})
                            vars=[item for item in vars if item!=vkey]
                    elif tf>0 and vkey=='max':
                            dicvar.update({'maxGrid':float(str1.split()[0].replace('d','0'))})
                            #var.remove(vkey)
                    elif tf>0:
                          #  print vkey
                            temp=lst[0].replace('d','0')
                            number2store=ast.literal_eval(temp)
                            #print vkey,':',number2store
                            dicvar.update({vkey:number2store})
                            vars=[item for item in vars if item!=vkey]
                    else:
                            pass
                    lines = [item for item in lines if item!= x]
        dicvar.update({'dx' : dicvar['z_i'] * dicvar['pi'] * 2.0 / dicvar['Nx']})
        dicvar.update({'dy' : dicvar['z_i'] * dicvar['pi'] * 2.0 * dicvar['l_r'] / dicvar['Ny']})
        dicvar.update({'dz' : dicvar['l_z']  / (dicvar['Nz']-1)})

    f.close()
    return dicvar

inparam = readInput('../input/LESinputs.txt')
dt = inparam['dt']
Nx = inparam['Nx']
Ny = inparam['Ny']
Nz = inparam['Nz']
u_star = inparam['u_star']
Ugal = inparam['Ugal']
Vgal = inparam['Vgal']
aNx = inparam['aNx']
z_i = inparam['dx']
l_z = inparam['dy']
dz = inparam['dz']
resDir = os.getcwd()
baseDir= ".."
# read break file to get tend
lines = np.loadtxt("../checkpoints/break", comments="#", delimiter="\n", unpack=False)
# get only the last line
tend = int(lines[-1]); 
ttt = int(lines[0]);

# get all frame list
uframeDir = baseDir + "/output/u_frame/"
if uframeDir[-1] != '/':
    uframeDir = uframeDir + '/'
os.chdir(uframeDir)
#  get a list of available frames
uflist = sp.check_output('ls')
uflist = uflist.split()

vframeDir = baseDir + "/output/v_frame/"
if vframeDir[-1] != '/':
    vframeDir = vframeDir + '/'
os.chdir(vframeDir)
#  get a list of available frames
vflist = sp.check_output('ls')
vflist = vflist.split()

wframeDir = baseDir + "/output/w_frame/"
if wframeDir[-1] != '/':
    wframeDir = wframeDir + '/'
os.chdir(wframeDir)
#  get a list of available frames
wflist = sp.check_output('ls')
wflist = wflist.split() 

lambda2FrameDir = baseDir + "/output/lambda2_frame/"
if lambda2FrameDir[-1] != '/':
    lambda2FrameDir = lambda2FrameDir + '/'
os.chdir(lambda2FrameDir)
#  get a list of available frames
lambda2flist = sp.check_output('ls')
lambda2flist = lambda2flist.split() 
print("List of Lambda2 Frames")
print(lambda2flist)
os.chdir(resDir)
nocf = (ttt-inparam['sframe']) / inparam['framestep']
if(nocf < 0) :
    print("No frames outputted yet !!")
else :    
    for i  in np.arange(len(uflist)):
        fout = open(baseDir +'/output/' + 'data_tstep_' + str(inparam['sframe']+i*inparam['framestep'])+'.xmf', 'w')
        fout.write(r"<?xml version=\"1.0\" ?>"); fout.write();
        fout.write(r"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>"); fout.write();
        fout.write(r"<Xdmf Version=\"2.0\">"); fout.write();
        fout.write(r"  <Domain>"); fout.write();
        fout.write(r"    <Grid Name=\"mesh\" GridType=\"Uniform\">")
        fout.write(r"      <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\" " + \
                            str(Nx) + str(r" ") +str(Ny) + str(r" ") + str(Nz) + "\"" +"\>")
        fout.write(r"      <Geometry GeometryType=\"Origin_DxDyDz\">")
        fout.write(r"        <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" \
                                Format=\"Xml\">")
        fout.write(r"      " + str(0.0) + str(0.0) + str(0.0))
        fout.write(r"        </DataItem>")
        fout.write(r"        <DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"Xml\">")
        fout.write(r"        " + str(dx) + str(r" ") + str(dy)+ str(r" ") + str(dz))
        fout.write(r"        </DataItem>")
        fout.write(r"      </Geometry>")
    # u-frame
        fout.write(r"      <Attribute Name = \"u\" AttributeType = \"Scalar\"  Center=\"Node\">")
        fout.write(r"        <DataItem Dimensions =\" " + str(Nx) + r" " + str(Ny) + r"  " + str(Nz) + \
                             r"\" "+ r"NumberType = \"Float\" " + r"Precision = \"8\"" +\
                            r"Format =\"Binary\"" + ">")
        fout.write(r"        " + uframeDir + uflist[i])
        fout.write(r"        </DataItem>")
        fout.write(r"      </Attribute>")    
    # v-frame
        fout.write(r"      <Attribute Name = \"v\" AttributeType = \"Scalar\"  Center=\"Node\">")
        fout.write(r"        <DataItem Dimensions =\" " + str(Nx) + r" " + str(Ny) + r"  " + str(Nz) + \
                             r"\" "+ r"NumberType = \"Float\" " + r"Precision = \"8\"" +\
                            r"Format =\"Binary\"" + ">")
        fout.write(r"        " + vframeDir + vflist[i])
        fout.write(r"        </DataItem>")
        fout.write(r"      </Attribute>")    
    # w-frame
        fout.write(r"      <Attribute Name = \"w\" AttributeType = \"Scalar\"  Center=\"Node\">")
        fout.write(r"        <DataItem Dimensions =\" " + str(Nx) + r" " + str(Ny) + r"  " + str(Nz) + \
                             r"\" "+ r"NumberType = \"Float\" " + r"Precision = \"8\"" +\
                            r"Format =\"Binary\"" + ">")
        fout.write(r"        " + wframeDir + wflist[i])
        fout.write(r"        </DataItem>")
        fout.write(r"      </Attribute>") 
  
    # lambda2-frame
        fout.write(r"      <Attribute Name = \"lambda2\" AttributeType = \"Scalar\"  Center=\"Node\">")
        fout.write(r"        <DataItem Dimensions =\" " + str(Nx) + r" " + str(Ny) + r"  " + str(Nz) + \
                             r"\" "+ r"NumberType = \"Float\" " + r"Precision = \"8\"" +\
                            r"Format =\"Binary\"" + ">")
        fout.write(r"        " + lambda2FrameDir + lambda2flist [i])
        fout.write(r"        </DataItem>")
        fout.write(r"      </Attribute>") 
        for j in np.arange(inparam['filtLevelsLambda']):        
  # lambda2-frame-filtered
            fout.write(r"      <Attribute Name = \"lambda2Filtered\" AttributeType = \"Scalar\"  Center=\"Node\">")
            fout.write(r"        <DataItem Dimensions =\" " + str(Nx) + r" " + str(Ny) + r"  " + str(Nz) + \
                             r"\" "+ r"NumberType = \"Float\" " + r"Precision = \"8\"" +\
                            r"Format =\"Binary\"" + ">")
            fout.write(r"        " + lambda2FrameDir + str(lambda2flist [i+len(uflist)+j]))
            fout.write(r"        </DataItem>")
            fout.write(r"      </Attribute>")  
        
        fout.write(r"    </Grid>")
        fout.write(r"  </Domain>")
        fout.write(r"</Xdmf>")
        fout.close()
