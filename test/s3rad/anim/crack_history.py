#!/usr/bin/env python3
"""One VTU holding every element the run destroyed, with the grip displacement
at which it died.  Threshold on theta_deleted in ParaView to scrub the whole
crack history from a single small file."""
import re,zlib,struct,sys
BLOCK=32768
def enc(raw):
    blocks=[raw[i:i+BLOCK] for i in range(0,len(raw),BLOCK)] or [b'']
    comp=[zlib.compress(b,6) for b in blocks]
    return struct.pack('<3I',len(blocks),BLOCK,len(blocks[-1]))+struct.pack('<%dI'%len(comp),*[len(c) for c in comp])+b''.join(comp)

deck,dmg,out=sys.argv[1],sys.argv[2],sys.argv[3]
txt=open(deck,errors='replace').read()
nodes={}
for m in re.finditer(r'\*Node\b[^\n]*\n(.*?)(?=\n\*)',txt,re.S|re.I):
    for ln in m.group(1).strip().split('\n'):
        p=[x.strip() for x in ln.split(',') if x.strip()]
        if len(p)>=4:
            try: nodes[int(p[0])]=(float(p[1]),float(p[2]),float(p[3]))
            except: pass
elems={}
for m in re.finditer(r'\*Element,\s*Type=([A-Za-z0-9]+)[^\n]*\n(.*?)(?=\n\*)',txt,re.S|re.I):
    et=m.group(1).upper()
    for ln in m.group(2).strip().split('\n'):
        p=[x.strip() for x in ln.split(',') if x.strip()]
        if len(p)<2: continue
        try: elems[int(p[0])]=(et,[int(x) for x in p[1:]])
        except: pass
recs=[]
for l in open(dmg):
    if l.startswith('#'): continue
    f=l.split(); recs.append((int(f[0]),int(f[2]),float(f[3]),int(f[5])))
recs=[r for r in recs if r[0] in elems and elems[r[0]][0]=='C3D4']
print("destroyed C3D4 elements: %d"%len(recs))
used=sorted({n for e,_,_,_ in recs for n in elems[e][1]})
remap={n:i for i,n in enumerate(used)}
pts=[]
for n in used: pts.extend(nodes[n])
conn=[];off=[];acc=0;th=[];inc=[];mat=[]
for e,i,t,mt in recs:
    for n in elems[e][1]: conn.append(remap[n])
    acc+=4; off.append(acc); th.append(t); inc.append(i); mat.append(mt)
nc=len(recs)
payload=b'';offs=[]
def add(raw):
    global payload
    offs.append(len(payload)); payload+=enc(raw)
add(struct.pack('<%df'%len(pts),*pts))
add(struct.pack('<%di'%len(conn),*conn))
add(struct.pack('<%di'%len(off),*off))
add(struct.pack('<%dB'%nc,*([10]*nc)))
add(struct.pack('<%df'%nc,*th))
add(struct.pack('<%di'%nc,*inc))
add(struct.pack('<%di'%nc,*mat))
with open(out,'wb') as f:
    f.write(b'<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">\n')
    f.write(('<UnstructuredGrid>\n<Piece NumberOfPoints="%d" NumberOfCells="%d">\n'%(len(used),nc)).encode())
    f.write(('<Points><DataArray type="Float32" NumberOfComponents="3" format="appended" offset="%d"/></Points>\n'%offs[0]).encode())
    f.write(('<Cells><DataArray type="Int32" Name="connectivity" format="appended" offset="%d"/><DataArray type="Int32" Name="offsets" format="appended" offset="%d"/><DataArray type="UInt8" Name="types" format="appended" offset="%d"/></Cells>\n'%(offs[1],offs[2],offs[3])).encode())
    f.write(('<CellData Scalars="theta_deleted">\n<DataArray type="Float32" Name="theta_deleted" format="appended" offset="%d"/>\n<DataArray type="Int32" Name="increment" format="appended" offset="%d"/>\n<DataArray type="Int32" Name="material_id" format="appended" offset="%d"/>\n</CellData>\n'%(offs[4],offs[5],offs[6])).encode())
    f.write(b'</Piece>\n</UnstructuredGrid>\n<AppendedData encoding="raw">\n_')
    f.write(payload); f.write(b'\n</AppendedData>\n</VTKFile>\n')
print("wrote %s  points=%d cells=%d  theta %.6f -> %.6f"%(out,len(used),nc,min(th),max(th)))
