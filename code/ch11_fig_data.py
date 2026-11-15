import json, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import ch11_icosian as C
PH=(1+5**0.5)/2
f=lambda z: z.a+z.b*PH            # 観察者側の表示：小数にするのはここだけ
V=[[f(c)/2 for c in u] for u in C.U]
decs=[]
for d in C.D:
    x=min(d); decs.append([C.T[x][p] for p in C.pw])
# 底の点：x a x̄（a は g の虚部）
a=(C.ZERO,)+tuple(C.U[C.g][1:])
base=[]
for d in decs:
    x=C.U[d[0]]; b=C.mul(C.mul(x,a),C.conj(x)); base.append([f(c) for c in b[1:]])
near=[[i,j] for i in range(12) for j in C.nb[i] if i<j]
orth=[C.ortho[i][0] for i in range(12)]
json.dump({"v":[[round(x,5) for x in p] for p in V],"dec":decs,"base":[[round(x,5) for x in p] for p in base],"near":near,"orth":orth},open("ch11_fig.json","w"))
print(len(decs),len(near),orth)
