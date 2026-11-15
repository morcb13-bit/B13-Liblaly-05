"""第20章：52次元 ── 13の番地に Z[ζ₁₀] の数を置く
13 の番地 i＝0〜12 に、それぞれ一つの数 z ＝ a＋bζ＋cζ²＋dζ³（a〜d は整数、ζ は36°の回転、ζ⁴ ＝ −1＋ζ−ζ²＋ζ³）を置く。
整数で並べると 13×4 ＝ 52 個。ζ を掛けることは、番地ごとの 4×4 の整数の区画 B（1→ζ→ζ²→ζ³→ζ⁴）で書ける。
基準（走らせる前に固定。結果を見てから変えない）
  20-1 次元：13×4 ＝ 52。各番地に 1, ζ, ζ², ζ³ を置いた 52 本の整数の列の行列式が ±1。
  20-2 回転行列：ζ を掛けることが 52×52 の整数の行列 R で書け、番地ごとに ζ を掛けた結果と一致する。
  20-3 10回転：R^10 ＝ I（整数の行列として）。
  20-4 固定方向：固有値 1 を持たない（det(R − I) ≠ 0、整数の行列式）。
  20-5 不変平面：どの番地の 4次元の区画も、36° の平面と 108° の平面の二つに分かれる。
       A₊ ＝ B² − φB ＋ I（36°：2cos36° ＝ φ）、A₋ ＝ B² − (1−φ)B ＋ I（108°：2cos108° ＝ 1−φ）を a＋bφ の行列として作り、
       ・A₊·A₋ ＝ A₋·A₊ ＝ 0
       ・A₊、A₋ の 3×3 の小行列式がすべて 0、2×2 の小行列式に 0 でないものがある（どちらも階数 2）
       ・A₊ − A₋ ＝ (1−2φ)B で、det B ≠ 0 かつ 1−2φ ≠ 0（二つの平面は 0 だけを共有）
       を確かめる。A₊ が 0 にする 2次元が 36° の平面、A₋ が 0 にする 2次元が 108° の平面。
       ［知識］2cos36° ＝ φ、2cos108° ＝ 1−φ。
  負の対照（記録。5項目の合否には含めない）
       12 回で一周する区画（x⁴ − x² ＋ 1、30°）で同じ行列を作ると、R^10 ≠ I で、B＋B⁻¹ の二乗は B＋B⁻¹＋I にならない。
"""
import itertools
from zphi import Z
P=Z(0,1); N=13; D=52
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
def mm(A,B):
    n=len(A); return [[sum(A[i][t]*B[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
def eye(n): return [[1 if i==j else 0 for j in range(n)] for i in range(n)]
def companion(c):          # 列ベクトル (a,b,c,d) に x を掛ける区画。x⁴ ＝ c0 + c1 x + c2 x² + c3 x³
    B=[[0]*4 for _ in range(4)]
    for j in range(3): B[j+1][j]=1
    for i in range(4): B[i][3]=c[i]
    return B
B=companion([-1,1,-1,1])                 # ζ⁴ ＝ −1 + ζ − ζ² + ζ³
def block_diag(B):
    R=[[0]*D for _ in range(D)]
    for k in range(N):
        for i in range(4):
            for j in range(4): R[4*k+i][4*k+j]=B[i][j]
    return R
# 20-1
cols=[]
for k in range(N):
    for p in range(4):
        v=[0]*D; v[4*k+p]=1; cols.append(v)
d1=det([[cols[j][i] for j in range(D)] for i in range(D)])
r1= D==13*4 and abs(d1)==1
print("20-1 13×4 ＝",D,"／1,ζ,ζ²,ζ³ を各番地に置いた52本の行列式",d1,"→","OK" if r1 else "NG")
# 20-2
R=block_diag(B)
def zeta_times(v):         # 番地ごとに ζ を掛ける（多項式として掛けて ζ⁴ を戻す）
    out=[]
    for k in range(N):
        a,b,c,d=v[4*k:4*k+4]
        # ζ·(a + bζ + cζ² + dζ³) ＝ aζ + bζ² + cζ³ + dζ⁴
        out+= [-d, a+d, b-d, c+d]
    return out
tests=[[(5*i*i+7*j+3)%9-4 for j in range(D)] for i in range(5)]
mv=lambda M,v:[sum(M[i][j]*v[j] for j in range(D)) for i in range(D)]
r2= all(mv(R,v)==zeta_times(v) for v in cols+tests) and all(isinstance(x,int) for r in R for x in r)
print("20-2 R を掛けた結果と、番地ごとに ζ を掛けた結果の一致（基底52本＋検査5本）",r2,"→","OK" if r2 else "NG")
# 20-3（区画の累乗から R の累乗を作り、52×52 でも確かめる）
pw=[eye(4)]
for k in range(10): pw.append(mm(pw[-1],B))
R10=block_diag(pw[10])
first=[k for k in range(1,11) if pw[k]==eye(4)]
r3= R10==eye(D)
print("20-3 R^10 ＝ I",r3,"／R^k ＝ I となる最初の k",first[0] if first else None,"／R^5 ＝ −I",pw[5]==[[-x for x in r] for r in eye(4)],"→","OK" if r3 else "NG")
# 20-4
RI=[[R[i][j]-(1 if i==j else 0) for j in range(D)] for i in range(D)]
d4=det(RI)
r4= d4!=0
print("20-4 det(R − I) ＝",d4,"（区画では",det([[B[i][j]-(1 if i==j else 0) for j in range(4)] for i in range(4)]),"）→","OK" if r4 else "NG")
# 20-5
def zmat(M): return [[Z(x) for x in r] for r in M]
def zmm(A,B): return [[sum((A[i][t]*B[t][j] for t in range(4)),Z()) for j in range(4)] for i in range(4)]
def zadd(*Ms): return [[sum((M[i][j] for M in Ms),Z()) for j in range(4)] for i in range(4)]
def zscale(M,s): return [[x*s for x in r] for r in M]
def zdet(M):
    n=len(M)
    if n==1: return M[0][0]
    s=Z()
    for j in range(n):
        if M[0][j]==Z(0): continue
        sub=[r[:j]+r[j+1:] for r in M[1:]]
        t=M[0][j]*zdet(sub); s= s+t if j%2==0 else s-t
    return s
def minors(M,k):
    return [zdet([[M[r][c] for c in cs] for r in rs]) for rs in itertools.combinations(range(4),k) for cs in itertools.combinations(range(4),k)]
Bz=zmat(B); B2=zmm(Bz,Bz); Iz=zmat(eye(4))
Ap=zadd(B2,zscale(Bz,Z(0)-P),Iz)              # B² − φB + I
Am=zadd(B2,zscale(Bz,P-Z(1)),Iz)              # B² − (1−φ)B + I
Zero=[[Z()]*4 for _ in range(4)]
prod_ok= zmm(Ap,Am)==Zero and zmm(Am,Ap)==Zero
rank2=lambda M: all(x==Z(0) for x in minors(M,3)) and any(x!=Z(0) for x in minors(M,2))
diff_ok= zadd(Ap,zscale(Am,Z(-1)))==zscale(Bz,Z(1)-P*2) and zdet(Bz)!=Z(0) and (Z(1)-P*2)!=Z(0)
# 補助：T ＝ B ＋ B⁻¹（B⁻¹ ＝ B⁹）が T² ＝ T ＋ I（第18章の φ 倍と同じ関係）
T=[[pw[1][i][j]+pw[9][i][j] for j in range(4)] for i in range(4)]
T2=mm(T,T); Tphi= T2==[[T[i][j]+(1 if i==j else 0) for j in range(4)] for i in range(4)]
r5= prod_ok and rank2(Ap) and rank2(Am) and diff_ok
print("20-5 A₊·A₋ ＝ A₋·A₊ ＝ 0",prod_ok,"／A₊ の階数 2",rank2(Ap),"／A₋ の階数 2",rank2(Am),"／二つの平面は 0 だけを共有",diff_ok,"→","OK" if r5 else "NG")
print("   補助：T ＝ B ＋ B⁻¹ は整数の行列で T² ＝ T ＋ I",Tphi,"／T ＝",T)
# 負の対照（記録）
C12=companion([-1,0,1,0])                # x⁴ ＝ −1 + x²（x⁴ − x² + 1 ＝ 0）
q=[eye(4)]
for k in range(12): q.append(mm(q[-1],C12))
T12=[[q[1][i][j]+q[11][i][j] for j in range(4)] for i in range(4)]
print("負の対照（記録） 30° の区画：R^10 ＝ I",q[10]==eye(4),"／R^12 ＝ I",q[12]==eye(4),"／(B＋B⁻¹)² ＝ B＋B⁻¹＋I",mm(T12,T12)==[[T12[i][j]+(1 if i==j else 0) for j in range(4)] for i in range(4)])
print("第20章の検定 総合（20-1〜20-5）","OK" if (r1 and r2 and r3 and r4 and r5) else "NG")
