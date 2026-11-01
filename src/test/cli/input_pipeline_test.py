"""Input pipeline regressions checked against exhaustive solutions over F_101."""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary=str(Path(sys.argv[1] if len(sys.argv)>1 else './drsolve').resolve())
expected={(x,y) for x in range(101) for y in range(101)
          if (x*x+y*y-13)%101==0 and (x-y+1)%101==0}
with tempfile.TemporaryDirectory() as tmp:
    tmp=Path(tmp);count=0
    def run(text, x='x', y='y', cutoff=None):
        global count
        src=tmp/'input.dr';src.write_bytes(text.encode())
        env=dict(os.environ)
        if cutoff is not None:env['DRSOLVE_BEZOUT_HGCD_CUTOFF']=str(cutoff)
        proc=subprocess.run([binary,'-f',str(src),'-o',str(tmp/'out.dr'),
                             '--threads','1','-v','2'],capture_output=True,text=True,check=True,env=env)
        found=set()
        for line in proc.stdout.splitlines():
            if '[PASS]' in line:
                values=dict(re.findall(r'(\w+)\s*=\s*(\d+)',line))
                if x in values and y in values:found.add((int(values[x]),int(values[y])))
        assert found==expected,(text[:100],found,proc.stdout,proc.stderr)
        count+=1
        return proc
    for nl in ['\n','\r\n','\r']:
        for x,y in [('x','y'),('x_1','_y'),('x','x_1')]:
            for expr in [f'{x}^2+{y}^2-13',f'({x}+{y})^2-2*{x}*{y}-13']:
                run(nl.join(['# comment','101',expr,f'{x}-{y}+1',' # options: --method 1','']),x,y)
    for cutoff in [0,32,64,1024,'invalid']:
        proc=run('101\nx^2+y^2-13\nx-y+1\n',cutoff=cutoff)
        if cutoff=='invalid':assert 'Invalid DRSOLVE_BEZOUT_HGCD_CUTOFF' in proc.stderr
    big=10**180+123456789
    assert big%101
    run(f'101\n{big}*x^2+{big}*y^2-{13*big}, x-y+1\n')
    # Signs, repeated variables/monomials, cancelling coefficients and spaces.
    run('101\n  x*x + y^2 + 101*x - 13 , (x-y)+1\n')
    # More than one input-buffer chunk, CRLF split over a chunk boundary,
    # and options after polynomial data. No trailing newline required.
    run(' '*65535+'\r\n101\r\n'+'0+'*40000+'x^2+y^2-13\r\nx-y+1')
    run('101\n'+'0+'*40000+'x^2+y^2-13\nx-y+1\n'+' '*65540+'# options: --method 1\n')
    src=tmp/'input.dr';src.write_text('101\n'+'0+'*40000+'x+y\nx-y\n# options: --method 99\n')
    proc=subprocess.run([binary,'-f',str(src),'-o',str(tmp/'out.dr')],capture_output=True,text=True)
    assert proc.returncode!=0 and 'invalid --method' in proc.stderr,proc
    print(f'Input pipeline: {count} exhaustive solution comparisons; late option validation PASS')
