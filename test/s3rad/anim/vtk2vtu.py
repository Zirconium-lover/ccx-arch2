#!/usr/bin/env python3
"""Convert CalculiX DE1 legacy ASCII VTK frames to compressed XML VTU + a PVD
collection.  Geometry is the deformed configuration written by the solver;
cells are the live C3D4 elements, so terminal deletion shows as material
disappearing.  Cell arrays kept: DE1_D (damage) and MATERIAL_ID."""
import sys,os,re,zlib,struct,base64

BLOCK=32768
def enc(raw):
    blocks=[raw[i:i+BLOCK] for i in range(0,len(raw),BLOCK)] or [b'']
    comp=[zlib.compress(b,6) for b in blocks]
    hdr=struct.pack('<3I',len(blocks),BLOCK,len(blocks[-1]))+struct.pack('<%dI'%len(blocks),*[len(c) for c in comp])
    return hdr+b''.join(comp)

def read_vtk(path):
    txt=open(path,'r',errors='replace').read()
    tm=re.search(r'time=([0-9.eE+-]+)',txt)
    inc=re.search(r'inc=(\d+)',txt)
    time=float(tm.group(1)) if tm else 0.0
    incn=int(inc.group(1)) if inc else 0
    tok=txt.split()
    def idx(word,after=0):
        for i in range(after,len(tok)):
            if tok[i]==word: return i
        return -1
    i=idx('POINTS'); npts=int(tok[i+1])
    p=i+3
    pts=[float(x) for x in tok[p:p+3*npts]]
    i=idx('CELLS',p); ncell=int(tok[i+1]); nsize=int(tok[i+2])
    c=i+3
    flat=[int(x) for x in tok[c:c+nsize]]
    conn=[];off=[];k=0;acc=0
    while k<len(flat):
        n=flat[k]; conn.extend(flat[k+1:k+1+n]); acc+=n; off.append(acc); k+=1+n
    i=idx('CELL_TYPES',c+nsize); types=[int(x) for x in tok[i+2:i+2+ncell]]
    arrays={}
    j=idx('CELL_DATA',i)
    while True:
        s=idx('SCALARS',j+1)
        if s<0: break
        name=tok[s+1]; dtype=tok[s+2]
        lt=idx('LOOKUP_TABLE',s)          # SCALARS name type [n] LOOKUP_TABLE default
        d=lt+2
        vals=tok[d:d+ncell]
        arrays[name]=(dtype,vals)
        j=d+ncell
    return time,incn,npts,pts,ncell,conn,off,types,arrays

def write_vtu(path,npts,pts,ncell,conn,off,types,arrays,keep):
    payload=b''; offsets=[]; 
    def add(raw):
        nonlocal payload
        offsets.append(len(payload)); payload+=enc(raw)
    add(struct.pack('<%df'%len(pts),*pts))                      # 0 points
    add(struct.pack('<%di'%len(conn),*conn))                    # 1 connectivity
    add(struct.pack('<%di'%len(off),*off))                      # 2 offsets
    add(struct.pack('<%dB'%len(types),*types))                  # 3 types
    meta=[]
    for name in keep:
        if name not in arrays: continue
        dtype,vals=arrays[name]
        if dtype in ('int','long'):
            add(struct.pack('<%di'%ncell,*[int(v) for v in vals])); meta.append((name,'Int32'))
        else:
            add(struct.pack('<%df'%ncell,*[float(v) for v in vals])); meta.append((name,'Float32'))
    with open(path,'wb') as f:
        f.write(b'<?xml version="1.0"?>\n')
        f.write(b'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">\n')
        f.write(('<UnstructuredGrid>\n<Piece NumberOfPoints="%d" NumberOfCells="%d">\n'%(npts,ncell)).encode())
        f.write(('<Points><DataArray type="Float32" NumberOfComponents="3" format="appended" offset="%d"/></Points>\n'%offsets[0]).encode())
        f.write(('<Cells><DataArray type="Int32" Name="connectivity" format="appended" offset="%d"/>'%offsets[1]).encode())
        f.write(('<DataArray type="Int32" Name="offsets" format="appended" offset="%d"/>'%offsets[2]).encode())
        f.write(('<DataArray type="UInt8" Name="types" format="appended" offset="%d"/></Cells>\n'%offsets[3]).encode())
        f.write(b'<CellData>\n')
        for n,(name,t) in enumerate(meta):
            f.write(('<DataArray type="%s" Name="%s" format="appended" offset="%d"/>\n'%(t,name,offsets[4+n])).encode())
        f.write(b'</CellData>\n</Piece>\n</UnstructuredGrid>\n<AppendedData encoding="raw">\n_')
        f.write(payload)
        f.write(b'\n</AppendedData>\n</VTKFile>\n')

if __name__=='__main__':
    src,dst,stride=sys.argv[1],sys.argv[2],int(sys.argv[3]) if len(sys.argv)>3 else 1
    os.makedirs(dst,exist_ok=True)
    files=sorted(f for f in os.listdir(src) if re.match(r'.*\.de1\.\d+\.vtk$',f))
    files=files[::stride]
    ent=[]
    for n,fn in enumerate(files):
        t,inc,npts,pts,ncell,conn,off,types,arr=read_vtk(os.path.join(src,fn))
        out='s3rad_%04d.vtu'%n
        write_vtu(os.path.join(dst,out),npts,pts,ncell,conn,off,types,arr,['DE1_D','MATERIAL_ID'])
        ent.append((t,out,inc,ncell))
        if n%25==0: print("  %s -> %s  theta=%.6f cells=%d"%(fn,out,t,ncell),flush=True)
    with open(os.path.join(dst,'s3rad.pvd'),'w') as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n')
        for t,o,inc,nc in ent: f.write('<DataSet timestep="%.9f" group="" part="0" file="%s"/>\n'%(t,o))
        f.write('</Collection>\n</VTKFile>\n')
    print("frames=%d  theta %.6f -> %.6f  cells %d -> %d"%(len(ent),ent[0][0],ent[-1][0],ent[0][3],ent[-1][3]))
