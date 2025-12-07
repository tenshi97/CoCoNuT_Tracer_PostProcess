import os
import fnmatch
import h5py

tbounce = 0.
t_alt = 0.0
t_sav = 0.0
clear = 0

def animatio(pattern, incremental_time=0.):
    str0='aaaaaaaaaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbbbbbbbbccccccccccccccccccccccccddddddddddddddddddddddddddeeeeeeeeeeeeeeeeeeeeeeeee'
    str1='abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwzabcdefghijklmnopqrstuvwxyzabcdfghijklmnopqrstuvwxyz'

    t_alt = 0.0
    t_sav = 0.0
    clear = 0

    files = []
    searchdir_str = '/media/tm314/Seagate Backup Plus Drive/s3.5_envel/'
    searchdir = os.listdir(searchdir_str)
    # pattern = 's3.5_envel_new.o*'

    for entry in searchdir:
        if fnmatch.fnmatch(entry, pattern):
            files.append(entry)

    nfiles = len(files)

    f = open('files_3d.txt', 'w')

    for finx in range(nfiles):
        print(files[finx])
        for y_string in ['Yin', 'Yang']:
            xdfm_process(files[finx], 0.0, t_alt, t_sav, clear, searchdir_str, pattern, y_string)
        print('File processed: ', files[finx])

    f.close()

def xdfm_process(data, tbounce, t_alt, t_sav, clear, searchdir, pattern, yinyang=''):
    data = searchdir+data
    file_handle = h5py.File(data)
    nslices = len(file_handle.keys())
    group_name = list(file_handle.keys())[-1]
    if group_name == 'info':
        nslices = nslices-1

    iter = -1
    file_finished = (nslices == 0)

    while not file_finished:
        iter += 1
        print('Data ', data, iter)

        file_finished = (iter >= nslices-2)
        slice = iter

        print(nslices, slice, iter)
        group_name = list(file_handle.keys())[slice] + '/' + yinyang
        group = list(file_handle.keys())[slice]
        print('GROUPS = ',group_name, group)

        dataset = file_handle['/'+group+'/time']
        time = dataset[()] - tbounce

        dataset = file_handle['/'+group+'/nstep']
        nstep = dataset[()]

        dataset = file_handle['/'+group+'/xzn']
        xzn = dataset[:]
        qx = dataset.size

        dataset = file_handle['/'+group+'/yzn']
        yzn = dataset[:]
        qy = dataset.size

        dataset = file_handle['/'+group+'/zzn']
        zzn = dataset[:]
        qz = dataset.size

        qn = 21

        if time-t_sav < 1.2e3 and clear == 0:
            print(50,' rm ', data)
            clear = 1
        else:
            clear = 0
            t_sav = time

        if time > t_alt + 0.0000:
            t_alt = time
            print(group_name)
            g = open(pattern[:-3]+'_'+yinyang+'_{:d}.xdmf'.format(nstep), 'w')
            write = lambda x: g.write(x + '\n')
            write('<?xml version="1.0" ?>')
            write("<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd'>")
            write('<Xdmf Version="2.0">')
            write('  <Domain>')

            write('    <Grid GridType="Uniform" Name="hydro-grid">')
            write('      <Time Value="{:14.12f}"/>'.format(time))
            write('      <Topology NumberOfElements="{:3d} {:3d} {:3d} " TopologyType="3DRectMesh"/>'.format(qz, qy, qx))
            write('      <Geometry GeometryType="VXVYVZ">')
            write('	<DataItem Dimensions="{:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qx) + data + ':' + group + '/xzn</DataItem>')
            write('	<DataItem Dimensions="{:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qy) + data + ':' + group + '/yzn</DataItem>')
            write('	<DataItem Dimensions="{:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz) + data + ':' + group + '/zzn</DataItem>')
            write('      </Geometry>')


            write('      <Attribute AttributeType="Scalar" Center="Node" Name="den">')
            write('	<DataItem Dimensions="{:3d} {:3d} {:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz, qy, qx) + data + ':' + group_name + '/den</DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="sto">')
            write('	<DataItem Dimensions="{:3d} {:3d} {:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz, qy, qx) + data + ':' + group_name + '/sto</DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="pre">')
            write('	<DataItem Dimensions="{:3d} {:3d} {:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz, qy, qx) + data + ':' + group_name + '/pre</DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="qen">')
            write('	<DataItem Dimensions="{:3d} {:3d} {:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz, qy, qx) + data + ':' + group_name + '/qen</DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="vex">')
            write('	<DataItem Dimensions="{:3d} {:3d} {:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz, qy, qx) + data + ':' + group_name + '/vex</DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="vey">')
            write('	<DataItem Dimensions="{:3d} {:3d} {:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz, qy, qx) + data + ':' + group_name + '/vey</DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="vez">')
            write('	<DataItem Dimensions="{:3d} {:3d} {:3d}" Format="HDF" NumberType="Float" Precision="8">'.format(qz, qy, qx) + data + ':' + group_name + '/vez</DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="Y_e">')
            write('	<DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            {:2d} 0 0 0'.format(qn-1))
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_n">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            0 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')


            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_p">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            1 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_d">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            2 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_t">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            3 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_3He">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            4 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_alpha">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            2 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')


            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_C">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            5 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Nitrogen">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            6 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')


            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_O">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            7 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Ne">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            8 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Mg">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            9 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Si">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            10 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_S">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            11 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Ar">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            12 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Ca">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            13 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Ti">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            14 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Cr">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            15 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Fe52">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            16 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Fe54">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            17 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Ni56">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            18 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Fe56">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            19 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Fe60">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            20 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('      <Attribute AttributeType="Scalar" Center="Node" Name="X_Ni62">')
            write('	    <DataItem ItemType ="HyperSlab" Dimensions = "{:3d} {:3d} {:3d}" Type = "HyperSlab">'.format(qz, qy, qx))
            write('          <DataItem Dimensions = "3 4" Format = "XML">')
            write('            21 0 0 0')
            write('            1 1 1 1')
            write('            1 {:3d} {:3d} {:3d}'.format(qz, qy, qx))
            write('          </DataItem>')
            write('          <DataItem Name = "xnu" Dimensions = "{:2d} {:3d} {:3d} {:3d}" Format = "HDF5">'.format(qn, qz, qy, qx) + data + ':' + group_name + '/xnu</DataItem>')
            write('        </DataItem>')
            write('      </Attribute>')

            write('    </Grid>')

            write('  </Domain>')
            write('</Xdmf>')

        g.close()
