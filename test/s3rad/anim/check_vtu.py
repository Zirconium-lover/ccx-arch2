#!/usr/bin/env python3
"""Read back a compressed appended-raw VTU and verify it against the source VTK."""
import sys,re,zlib,struct
def parse_vtu(path):
    raw=open(path,'rb').read()
    m=re.search(rb'<AppendedData encoding="raw">\s*_',raw)
    base=m.end()
    hdr=raw[:m.start()].decode('utf8',errors='replace')
    npts=int(re.search(r'NumberOfPoints="(\d+)"',hdr).group(1))
    ncell=int(re.search(r'NumberOfCells="(\d+)"',hdr).group(1))
    arrs={}
    for dm in re.finditer(r'<DataArray type="(\w+)"(?: NumberOfComponents="(\d+)")? Name="([^"]*)" format="appended" offset="(\d+)"',hdr):
        arrs[dm.group(3)]=(dm.group(1),int(dm.group(4)))
    pm=re.search(r'<Points><DataArray type="(\w+)" NumberOfComponents="3" format="appended" offset="(\d+)"',hdr)
    arrs['__points__']=(pm.group(1),int(pm.group(2)))
    def read(off):
        p=base+off
        nb,bs,last=struct.unpack_from('<3I',raw,p); p+=12
        cs=struct.unpack_from('<%dI'%nb,raw,p); p+=4*nb
        out=b''
        for c in cs:
            out+=zlib.decompress(raw[p:p+c]); p+=c
        return out
    fmt={'Float32':('f',4),'Int32':('i',4),'UInt8':('B',1)}
    out={}
    for name,(t,off) in arrs.items():
        b=read(off); ch,sz=fmt[t]
        out[name]=struct.unpack('<%d%s'%(len(b)//sz,ch),b)
    return npts,ncell,out

def parse_vtk(path):
    tok=open(path,errors='replace').read().split()
    def idx(w,a=0):
        for i in range(a,len(tok)):
            if tok[i]==w: return i
        return -1
    i=idx('POINTS'); npts=int(tok[i+1]); pts=[float(x) for x in tok[i+3:i+3+3*npts]]
    j=idx('CELLS',i); ncell=int(tok[j+1])
    k=idx('CELL_DATA',j); s=idx('SCALARS',k)
    lt=idx('LOOKUP_TABLE',s); vals=[float(x) for x in tok[lt+2:lt+2+ncell]]
    return npts,ncell,pts,vals,tok[s+1]

vtu,vtk=sys.argv[1],sys.argv[2]
npts,ncell,arrs=parse_vtu(vtu)
n2,c2,pts,vals,first=parse_vtk(vtk)
ok=True
def chk(label,a,b,tol=0.0):
    global ok
    good = (abs(a-b)<=tol) if isinstance(a,float) else (a==b)
    print("  %-28s vtu=%-16s vtk=%-16s %s"%(label,a,b,"ok" if good else "MISMATCH"))
    ok = ok and good
chk("points",npts,n2); chk("cells",ncell,c2)
P=arrs['__points__']
chk("point coords count",len(P),len(pts))
mx=max(abs(P[i]-pts[i]) for i in range(len(pts)))
print("  %-28s max|dx| = %.3e (float32 round-off)"%("coordinate round-trip",mx))
ok = ok and mx < 1e-4
D=arrs.get(first)
if D:
    md=max(abs(D[i]-vals[i]) for i in range(ncell))
    print("  %-28s max|d| = %.3e"%("%s round-trip"%first,md))
    ok = ok and md < 1e-6
conn=arrs['connectivity']; off=arrs['offsets']; ty=arrs['types']
chk("connectivity length",len(conn),4*ncell)
chk("offsets last",off[-1],4*ncell)
chk("all cell types = 10 (tetra)",all(t==10 for t in ty),True)
chk("node ids 0-based in range",max(conn)<npts and min(conn)>=0,True)
print("VERDICT:","PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
