"""第18章：単位を二つ並べる ── 13の番地に a＋bφ を置くと26次元
13 の番地 i＝0〜12 に、それぞれ一つの数 z_i ＝ a_i ＋ b_i φ（a, b は整数）を置く。
整数で並べると (a_0, b_0, a_1, b_1, …, a_12, b_12) の 26 個。すべて整数の行列と整数の比較だけで判定する。
基準（走らせる前に固定。結果を見てから変えない）
  18-1 次元：13 個の a＋bφ の並びと、整数 26 個の並びは一対一。1 と φ を各番地に置いた 26 本の整数の列は、
       26×26 の行列として行列式 ±1（26 本で、整数の格子 Z^26 をちょうど張る）。
  18-2 φ倍の行列：φ·(a＋bφ) ＝ b ＋ (a＋b)φ なので、各番地で (a,b) → (b, a＋b)。これを 26×26 の整数行列 M に書き、
       26 本の基底と検査用の整数の列で、M を掛けた結果が番地ごとに φ を掛けた結果と一致する。
  18-3 可逆性：det M ＝ ±1。
  18-4 φ の関係：M² ＝ M ＋ I が整数行列として成り立つ。
  18-5 幾何の包含：正二十面体の頂点12個（(0,±1,±φ) の巡回、共通倍率 1）と、600胞体の頂点120個
       （単位イコシアン、成分を2倍して a＋bφ の整数にする：共通倍率 2）を、番地 0〜2 と 0〜3 に置き、
       整数 26 成分の列として書ける。像はどれも異なる（12個・120個）。
       計量は Z^26 の整数の内積とは別に、番地ごとの積の和 B(u,v) ＝ Σ z_i·w_i（a＋bφ の整数）として入れる。
       この計量で、正二十面体は12頂点の長さ²が1種類・最短の組が30、600胞体は120頂点の長さ²が1種類・最短の組が720。
  負の対照（記録。5項目の合否には含めない）
       各番地で (a,b)→(b,a) と入れ替えるだけの行列 S では S² ＝ I で、S² ＝ S ＋ I は成り立たない。
       正二十面体の像を Z^26 の整数の内積で測ると、最短の組は 30 にならない（計量を入れ替えると形は保たれない）。
"""
import itertools, collections, io, contextlib
from zphi import Z
P=Z(0,1); N=13; D=2*N
def to_int(zs):            # 13 個の a＋bφ → 26 個の整数
    out=[]
    for z in zs: out+= [z.a, z.b]
    return out
def to_phi(v): return [Z(v[2*i],v[2*i+1]) for i in range(N)]
def det(M):
    M=[r[:] for r in M]; n=len(M); prev=1; sgn=1
    for k in range(n-1):
        if M[k][k]==0:
            sw=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if sw is None: return 0
            M[k],M[sw]=M[sw],M[k]; sgn=-sgn
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sgn*M[n-1][n-1]
def matvec(M,v): return [sum(M[i][j]*v[j] for j in range(D)) for i in range(D)]
def matmul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(D)) for j in range(D)] for i in range(D)]
I=[[1 if i==j else 0 for j in range(D)] for i in range(D)]
# 18-1
cols=[]
for i in range(N):
    for unit in (Z(1),P):
        zs=[Z(0)]*N; zs[i]=unit; cols.append(to_int(zs))
B26=[[cols[j][i] for j in range(D)] for i in range(D)]
d1=det(B26)
tests=[[(7*i*i+3*j-5)%11-5 for j in range(D)] for i in range(5)]   # 検査用の整数の列（決まった式で作る）
round_ok= all(to_int(to_phi(v))==v for v in tests)
r1= abs(d1)==1 and round_ok and len(cols)==26
print("18-1 1 と φ を各番地に置いた26本の行列式",d1,"／整数26個 ⇄ a＋bφ 13個の往復",round_ok,"→","OK" if r1 else "NG")
# 18-2
M=[[0]*D for _ in range(D)]
for i in range(N):
    a,b=2*i,2*i+1
    M[a][b]=1            # 新しい a ＝ b
    M[b][a]=1; M[b][b]=1 # 新しい b ＝ a＋b
def phi_times(v): return to_int([P*z for z in to_phi(v)])
chk= all(matvec(M,v)==phi_times(v) for v in cols+tests)
r2= chk and all(isinstance(x,int) for r in M for x in r)
print("18-2 M を掛けた結果と、番地ごとに φ を掛けた結果の一致（基底26本＋検査5本）",chk,"→","OK" if r2 else "NG")
# 18-3
dM=det(M); r3= abs(dM)==1
print("18-3 det M ＝",dM,"→","OK" if r3 else "NG")
# 18-4
M2=matmul(M,M); MI=[[M[i][j]+I[i][j] for j in range(D)] for i in range(D)]
r4= M2==MI
print("18-4 M² ＝ M ＋ I",r4,"→","OK" if r4 else "NG")
# 18-5
def Bform(u,v):
    s=Z()
    for x,y in zip(to_phi(u),to_phi(v)): s=s+x*y
    return s
ico=set()
for s1,s2 in itertools.product((1,-1),repeat=2):
    for c in range(3):
        v=[Z(0)]*3; v[(c+1)%3]=Z(s1); v[(c+2)%3]=P*s2; ico.add(tuple(v))
ICO=[to_int(list(v)+[Z(0)]*(N-3)) for v in sorted(ico,key=repr)]
with contextlib.redirect_stdout(io.StringIO()):
    import ch11_icosian as C
CELL=[to_int(list(u)+[Z(0)]*(N-4)) for u in C.U]           # C.U は成分を2倍した単位イコシアン
def shape(V,form):
    norms=collections.Counter(repr(form(v,v)) for v in V)
    d2=collections.Counter()
    for u,v in itertools.combinations(V,2):
        w=[x-y for x,y in zip(u,v)]; d2[form(w,w)]+=1
    mn=min(d2); return norms,d2[mn],mn
nI,eI,mI=shape(ICO,Bform); nC,eC,mC=shape(CELL,Bform)
ints_ok= all(isinstance(x,int) for V in (ICO,CELL) for v in V for x in v) and all(len(v)==26 for V in (ICO,CELL) for v in V)
distinct= len(set(map(tuple,ICO)))==12 and len(set(map(tuple,CELL)))==120
r5= ints_ok and distinct and len(nI)==1 and eI==30 and len(nC)==1 and eC==720
print("18-5 整数26成分で書ける",ints_ok,"／像はすべて異なる",distinct)
print("     正二十面体（倍率1）：長さ²",dict(nI),"／最短の組",eI,"（長さ²",repr(mI),"）")
print("     600胞体（倍率2）：長さ²",dict(nC),"／最短の組",eC,"（長さ²",repr(mC),"） → ","OK" if r5 else "NG")
# 負の対照（記録）
S=[[0]*D for _ in range(D)]
for i in range(N): S[2*i][2*i+1]=1; S[2*i+1][2*i]=1
S2=matmul(S,S); SI=[[S[i][j]+I[i][j] for j in range(D)] for i in range(D)]
dot=lambda u,v: sum(x*y for x,y in zip(u,v))
nZ,eZ,mZ=shape(ICO,dot)
print("負の対照（記録） 入れ替え S：S²＝I",S2==I,"／S²＝S＋I",S2==SI,"／正二十面体を Z^26 の内積で測ると 長さ²",dict(nZ),"最短の組",eZ)
print("第18章の検定 総合（18-1〜18-5）","OK" if (r1 and r2 and r3 and r4 and r5) else "NG")
