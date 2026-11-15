"""第15章：偶数の次元と奇数の次元 ── 裏返せるか、真ん中に頂点があるか
n 次元の立方体（頂点は ±1 の列）を n＝1〜13 で数える。すべて整数。
改訂の記録：第4部の範囲を11次元から13次元に改めたため、同じ基準で n＝13 まで延ばした（旧版は ch15_parity_v1.py／_v1_out.txt）。
基準（走らせる前に決めた）
  P1 立方体の外側（n 次元の面そのものを除くすべての面）の数を、頂点から符号を交互に付けて足すと、
     偶数の n で 0、奇数の n で 2。
  P2 立方体を自分に重ねる回転（±1 を置いた置換の行列で行列式 ＋1）のうち、動かない向きを持たないものは、
     偶数の n にだけあり、奇数の n には一つもない。すべての頂点を向かいへ写す −1 倍は、偶数の n でだけ回転になる。
     n≤7 は行列を全部作って行列式で判定し、n≤13 は置換の輪（巡回）ごとの符号の積で数える。
     ［数え方］動かない向きがないのは、置換のどの輪でも符号の積が −1 のとき。n≤7 で二つの数え方が一致すること。
  P3 長い対角線に垂直で中心を通る面は、偶数の n では C(n, n/2) 個の頂点を通り、奇数の n ではどの頂点も通らない。
     奇数の n では、その面は辺の中点を C(n,(n−1)/2)·(n+1)/2 個通る。
  負の対照 奇数の n にも、動かない向きを持たない「行列式 −1 のもの」（鏡に映す写し）はある（−1 倍がその一つ）。
"""
import itertools
from math import comb
N=13
# P1 面を全部数える（自由に動く座標の組と、止める座標の符号）
def faces(n,k): return sum(1 for free in itertools.combinations(range(n),k) for _ in itertools.product((1,-1),repeat=n-k))
F={n:[faces(n,k) for k in range(n+1)] for n in range(1,N+1)}
bd={n:sum((-1)**k*F[n][k] for k in range(n)) for n in range(1,N+1)}
p1= all(bd[n]==(0 if n%2==0 else 2) for n in bd)
print("P1 外側の面を交互に足した数",bd,"→","OK" if p1 else "NG")
for n in (6,7,13): print(f"   {n}次元の立方体 面の数",F[n])
# P2 直接（n≤7）
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
def direct(n):
    rot=0; rotfree=0; reffree=0
    for perm in itertools.permutations(range(n)):
        for sg in itertools.product((1,-1),repeat=n):
            M=[[0]*n for _ in range(n)]
            for i in range(n): M[i][perm[i]]=sg[i]
            d=det(M)
            A=[[M[i][j]-(1 if i==j else 0) for j in range(n)] for i in range(n)]
            free=det(A)!=0
            if d==1:
                rot+=1; rotfree+=free
            elif free: reffree+=1
    return rot,rotfree,reffree
# P2 置換の輪で数える（n≤13）
def cycles(perm):
    seen=[False]*len(perm); out=[]
    for i in range(len(perm)):
        if not seen[i]:
            L=0; j=i
            while not seen[j]: seen[j]=True; j=perm[j]; L+=1
            out.append(L)
    return out
def bycycle(n):
    # 符号の選び方 2^n のうち、どの輪でも符号の積が −1 になるもの：輪ごとに 2^(L−1) 通り。
    # 行列式＝置換の符号 × 符号の積 ＝ Π(−1)^(L−1) × Π(−1) ＝ (−1)^n
    rotfree=0; reffree=0
    from collections import Counter
    # 置換を輪の長さの組ごとに数える（全部の置換をたどるのは n≤9、それより上は輪の型ごとの個数を整数で出す）
    def parts(n,mx):
        if n==0: yield []; return
        for k in range(min(n,mx),0,-1):
            for p in parts(n-k,k): yield [k]+p
    from math import factorial
    for p in parts(n,n):
        c=Counter(p); cnt=factorial(n)
        for L,m in c.items(): cnt//= (L**m)*factorial(m)
        ways=1
        for L in p: ways*=2**(L-1)
        detv=1
        for L in p: detv*=(-1)**(L-1)*(-1)
        if detv==1: rotfree+=cnt*ways
        else: reffree+=cnt*ways
    return rotfree,reffree
D={n:direct(n) for n in range(1,8)}
C={n:bycycle(n) for n in range(1,N+1)}
agree= all(D[n][1]==C[n][0] and D[n][2]==C[n][1] for n in range(1,8))
rf={n:C[n][0] for n in C}
neg1={n:(1 if n%2==0 else -1) for n in range(1,N+1)}      # −1 倍の行列式 (−1)^n
p2= agree and all((rf[n]>0)==(n%2==0) for n in rf) and all((neg1[n]==1)==(n%2==0) for n in neg1)
print("P2 動かない向きのない回転",rf,"／n≤7 で行列式の数え方と一致",agree,"／−1 倍の行列式",neg1,"→","OK" if p2 else "NG")
# P3
mid={}; edgemid={}
for n in range(1,N+1):
    V=list(itertools.product((1,-1),repeat=n))
    mid[n]=sum(1 for v in V if sum(v)==0)
    # 中心の面（座標の和 0）が通る辺の中点：辺の両端の和が +1 と −1
    edgemid[n]=sum(1 for v in V for i in range(n) if v[i]==1 and sum(v)==1)
p3= all(mid[n]==(comb(n,n//2) if n%2==0 else 0) for n in mid) and all((edgemid[n]==(comb(n,(n-1)//2)*(n+1)//2)) for n in edgemid if n%2==1) \
    and all(edgemid[n]==0 for n in edgemid if n%2==0)
print("P3 中心の面が通る頂点",mid,"／辺の中点",edgemid,"→","OK" if p3 else "NG")
neg= all(C[n][1]>0 for n in range(1,N+1) if n%2==1)
print("負の対照 奇数の n の、動かない向きのない鏡の写し（行列式 −1）",{n:C[n][1] for n in C if n%2==1},"→","OK" if neg else "NG")
print("第15章の検定 総合","OK" if (p1 and p2 and p3 and neg) else "NG")
