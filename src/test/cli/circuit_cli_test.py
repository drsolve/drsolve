"""Independent expanded-polynomial comparisons for circuit specialization."""
import os,re,random,subprocess,sys,tempfile
from pathlib import Path
binary=str(Path(sys.argv[1] if len(sys.argv)>1 else './drsolve').resolve())

def add(a,b,p):
 c=dict(a)
 for k,v in b.items():c[k]=(c.get(k,0)+v)%p
 return {k:v for k,v in c.items() if v}
def mul(a,b,p):
 c={}
 for (i,j),u in a.items():
  for (k,l),v in b.items():c[i+k,j+l]=(c.get((i+k,j+l),0)+u*v)%p
 return {k:v for k,v in c.items() if v}
def pretty(a):
 return '+'.join(str(v)+''.join('*'+name+(f'^{e}' if e!=1 else '') for name,e in zip('xy',k) if e) for k,v in sorted(a.items(),reverse=True)) or '0'
def monic(s,p):
 s=re.sub(r'\s+','',s);coeff={}
 for term in s.replace('-','+-').split('+'):
  if not term:continue
  m=re.fullmatch(r'([+-]?\d+)?\*?([xy])?(?:\^(\d+))?',term)
  assert m,(s,term)
  v=int(m[1]) if m[1] else 1;d=int(m[3]) if m[3] else 1 if m[2] else 0
  coeff[d]=(coeff.get(d,0)+v)%p
 coeff={k:v for k,v in coeff.items() if v}
 if coeff:
  inv=pow(coeff[max(coeff)],-1,p);coeff={k:v*inv%p for k,v in coeff.items()}
 return coeff

def pairs(log):
 return {tuple(map(int,m)) for m in re.findall(r'x = (\d+), y = (\d+) \[PASS\]',log)}

rng=random.Random(77)
with tempfile.TemporaryDirectory() as d:
 d=Path(d);path=d/'input.drc';out=d/'output';comparisons=0;extension=False
 def run(args):return subprocess.run([binary,*args],capture_output=True,text=True,timeout=120)
 for p,degree in [(17,3),(101,2),(65537,3),(2130706433,3),(18446744069414584321,2)]:
  for eliminate in ['x','y']:
   lines=['DRCIRCUIT 1',f'prime {p}',f'eliminate {eliminate}'];polys=[{(1,0):1},{(0,1):1}]
   def emit(line,poly):lines.append(line);polys.append(poly);return len(polys)-1
   def power(a,e):
    poly={(0,0):1}
    for _ in range(e):poly=mul(poly,polys[a],p)
    return emit(f'pow {a} {e}',poly)
   x=power(0,degree);y=power(1,degree)
   xy=emit('mul 0 1',mul(polys[0],polys[1],p))
   shared=emit(f'add {x} {y}',add(polys[x],polys[y],p))
   zero=emit(f'sub {xy} {xy}',{})
   const=emit('const -123456789123456789123456789',{(0,0):(-123456789123456789123456789)%p})
   one=power(zero,0)
   # Repeated references, shared subexpressions and zero powers exercise liveness.
   outputs=[]
   for k in range(2):
    terms=[(1,shared),(k+1,xy),(rng.randrange(1,p),0),(rng.randrange(1,p),1),(1,one),(1,zero),(1,const),(2,const)]
    poly={(0,0):k+1}
    for coeff,ref in terms:poly=add(poly,{e:coeff*v%p for e,v in polys[ref].items()},p)
    outputs.append(emit('lin '+str(k+1)+''.join(f' {a} {b}' for a,b in terms),poly))
   lines.append('outputs '+' '.join(map(str,outputs)));path.write_text('\n'.join(lines)+'\n')
   proc=run(['--circuit',str(path),'--threads','2','--resultant-only','-o',str(out)])
   assert proc.returncode==0,proc.stderr+proc.stdout
   circuit=monic(out.read_text().split('Resultant:\n')[1],p)
   extension|=bool(re.search(r'field degree ([2-9]|\d\d)',proc.stdout))
   if p==17 and eliminate=='x':
    # Exercise buffer growth, CRLF, blank lines and EOF without a newline.
    padded=lines.copy();padded[1]+=' '*10000+'# long line'
    path.write_bytes(('\r\n\r\n'+'\r\n'.join(padded)).encode())
    portable=run(['--circuit',str(path),'--threads','2','--resultant-only','-o',str(out)])
    assert portable.returncode==0,portable.stderr+portable.stdout
    assert circuit==monic(out.read_text().split('Resultant:\n')[1],p)
   equations=','.join(pretty(polys[i]) for i in outputs)
   ref=run([equations,eliminate,str(p),'--threads','1','--resultant-only','--no-bezout-block','-o',str(out)])
   assert ref.returncode==0,ref.stderr
   assert circuit==monic(out.read_text().split('Resultant:\n')[1],p),(p,eliminate,circuit,out.read_text())
   proc=run(['--circuit',str(path),'--threads','2','-o',str(out)])
   assert proc.returncode==0,proc.stderr+proc.stdout
   if p<=101:
    expected={(x,y) for x in range(p) for y in range(p) if all(sum(v*pow(x,i,p)*pow(y,j,p) for (i,j),v in polys[k].items())%p==0 for k in outputs)}
    assert pairs(proc.stdout)==expected,(proc.stdout,expected)
   else:
    ref=run([equations,str(p),'--threads','1','-v','2','-o',str(out)])
    assert ref.returncode==0 and pairs(proc.stdout)==pairs(ref.stdout),(proc.stdout,ref.stdout,ref.stderr)
   comparisons+=1
 # Several adaptive rounds must reuse read-only inverse-Hankel samples.
 # Compare against a fully expanded independent resultant, including retries
 # from a user-specified tiny budget.
 lines=['DRCIRCUIT 1','prime 2130706433','eliminate x','const 1']
 next_id=3;xp=[2,0];yp=[2,1]
 for e in range(2,17):
  xp.append(next_id);lines.append(f'pow 0 {e}');next_id+=1
  yp.append(next_id);lines.append(f'pow 1 {e}');next_id+=1
 terms=[]
 for i in range(17):
  for j in range(17-i):
   if i and j:
    ref_id=next_id;lines.append(f'mul {xp[i]} {yp[j]}');next_id+=1
   else:ref_id=xp[i] if i else yp[j]
   terms.append(((i,j),ref_id))
 outputs=[];equations=[]
 for k in range(2):
  poly={e:rng.randrange(1,2130706433) for e,_ in terms}
  lines.append('lin 0'+''.join(f' {poly[e]} {ref_id}' for e,ref_id in terms))
  outputs.append(next_id);next_id+=1;equations.append(pretty(poly))
 lines.append('outputs '+' '.join(map(str,outputs)))
 path.write_text('\n'.join(lines)+'\n')
 ref=run([','.join(equations),'x','2130706433',
          '--threads','1','--resultant-only','--no-bezout-block','-o',str(out)])
 assert ref.returncode==0,ref.stderr
 expected=monic(out.read_text().split('Resultant:\n')[1],2130706433)
 for flags in ([],['--bezout-initial-points','1']):
  proc=subprocess.run([binary,'--circuit',str(path),'--threads','2',
                       '--resultant-only','-o',str(out),*flags],
                      capture_output=True,text=True,timeout=120,
                      env=dict(os.environ,DRSOLVE_PROFILE_BEZOUT='1'))
  assert proc.returncode==0,proc.stdout+proc.stderr
  assert expected==monic(out.read_text().split('Resultant:\n')[1],2130706433)
  if flags:
   assert re.search(r'Circuit attempt [2-9]:.*reused=[1-9]',proc.stderr),proc.stderr
  else:
   assert 'Circuit attempt 1: points=86 status=accepted' in proc.stderr,proc.stderr

 # Shared power gates must reduce the startup estimate without relying on
 # Poseidon names or round metadata. An unused high-degree gate has no effect.
 for a,b in ((3,3),(2,5)):
  for eliminate in ('x','y'):
   path.write_text(f'''DRCIRCUIT 1
prime 65537
eliminate {eliminate}
add 0 1
pow 2 {a}
add 3 0
pow 4 {b}
lin 1 1 5 1 0
lin 2 1 5 1 1
pow 0 1000
outputs 6 7
''')
   proc=run(['--circuit',str(path),'--threads','2','--resultant-only','-o',str(out)])
   assert proc.returncode==0,proc.stdout+proc.stderr
   assert f'resultant_degree_estimate={a*b} ' in proc.stdout,proc.stdout
   got=monic(out.read_text().split('Resultant:\n')[1],65537)
   common=f'((x+y)^{a}+x)^{b}'
   ref=run([f'{common}+x+1,{common}+y+2',eliminate,'65537',
            '--threads','1','--resultant-only','--no-bezout-block','-o',str(out)])
   assert ref.returncode==0 and got==monic(out.read_text().split('Resultant:\n')[1],65537),(proc.stdout,ref.stdout,ref.stderr)

 # A cubic identity reduces to a quadratic: degree certification must fold
 # the exact cancellation, rather than guess the smaller bound from samples.
 path.write_text("""DRCIRCUIT 1
prime 101
add 0 1
pow 2 3
pow 0 3
pow 1 3
pow 0 2
pow 1 2
mul 6 1
mul 0 7
lin 0 1 3 -1 4 -1 5 -3 8 -3 9 1 6 1 7
lin 1 1 10 1 0
outputs 10 11
""")
 proc=run(['--circuit',str(path),'--threads','2','--resultant-only','-o',str(out)])
 assert proc.returncode==0,proc.stderr
 got=monic(out.read_text().split('Resultant:\n')[1],101)
 ref=run(['x^2+y^2,x^2+y^2+x+1','y','101','--resultant-only','--no-bezout-block','-o',str(out)])
 assert ref.returncode==0 and got==monic(out.read_text().split('Resultant:\n')[1],101)
 assert extension,'extension verification was not exercised'
 for text in [
  'DRCIRCUIT 1\nprime 15\noutputs 0 1',
  'DRCIRCUIT 1\nprime 101\nadd 0 2\noutputs 2 2',
  'DRCIRCUIT 1\nprime 101\npow 1 2\noutputs 1 2',
  'DRCIRCUIT 1\nprime 101\nmul 0 1\noutputs 2 2',
  'DRCIRCUIT 1\nprime 101\npow 1 18446744073709551615\noutputs 2 2',
  'DRCIRCUIT 1\nprime 101\nadd 0 1\noutputs 2 2', # identically zero resultant: must not claim no roots
  'DRCIRCUIT 1\nprime 101\nadd 0 1\noutputs 2 2\nconst 1',
 ]:
  path.write_text(text);proc=run(['--circuit',str(path),'--threads','2'])
  assert proc.returncode!=0 and 'No solutions' not in proc.stdout,(text,proc.stdout,proc.stderr)
 print(f'Circuit: {comparisons} expanded resultant/solution comparisons; extension checks and rejected inputs PASS')
