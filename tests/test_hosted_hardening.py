"""Hosted-mode isolation and error-sanitization gates. Offline."""
import os, sys
os.environ["OLIGOFORGE_HOSTED"] = "1"
os.environ["OLIGOFORGE_ALLOW_SERVER_STORAGE"] = "0"
os.environ["OLIGOFORGE_ALLOW_SHARED_CONDITIONS"] = "0"
sys.path.insert(0, ".")

from fastapi.testclient import TestClient
import app
from oligoforge import thermo as T

fails=[]
def check(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + (("  ["+str(detail)+"]") if detail and not cond else ""))
    if not cond: fails.append(name)

c=TestClient(app.app)
h=c.get('/healthz')
check('health hosted flags', h.json().get('hosted_mode') and not h.json().get('server_storage_enabled'))
check('security headers present', h.headers.get('x-content-type-options')=='nosniff' and h.headers.get('x-frame-options')=='DENY')
check('shared conditions disabled', c.post('/api/conditions',json={'mv_conc':80}).status_code==403)
check('server projects disabled', c.post('/api/project/save',json={'name':'x','assays':[]}).status_code==403)
check('local BLAST path disabled', c.post('/api/blast',json={'seq':'ACGTACGTACGT','mode':'local','db_path':'/tmp/private'}).status_code==403)

# Full server diagnostics stay in logs; response must not reveal exception text.
orig=app.RPT.build
try:
    app.RPT.build=lambda *a,**k: (_ for _ in ()).throw(RuntimeError('/tmp/private/APIKEY=SECRET'))
    r=c.post('/api/report',json={'panel':[]}).json()
finally:
    app.RPT.build=orig
check('hosted exception sanitized', r.get('error')=='report failed' and 'SECRET' not in str(r), r)

# Validation errors omit the rejected input value.
secret='SECRET_'+'X'*200
r=c.post('/api/qc',json={'seq':secret,'role':{'bad':'shape'}})
body=r.json()
check('validation response does not echo payload', r.status_code==422 and secret not in str(body), body)

# Request-local orthopanel conditions must not mutate shared process state.
before=(dict(T.COND),T.ANNEAL_C)
payload={'candidates':[{'name':'a','seq':'ACGATCAGTTGCATCAGGTA'}], 'mv_conc':123, 'dv_conc':4, 'anneal_c':55}
r=c.post('/api/orthogonal-panel',json=payload).json()
after=(dict(T.COND),T.ANNEAL_C)
check('orthopanel reports request conditions', r.get('conditions',{}).get('mv_conc')==123)
check('orthopanel does not mutate global conditions', before==after,(before,after))

# Strict order path rejects opaque/malformed modification syntax.
r=c.post('/api/order_csv',json={'oligos':[{'name':'x','kind':'probe_lna','seq':'ACGT+NACGT'}],'gblocks':[]}).json()
check('malformed modified order rejected', bool(r.get('error')),r)

print(('FAIL: '+', '.join(fails)) if fails else 'ALL HOSTED HARDENING ASSERTS PASS')
raise SystemExit(1 if fails else 0)
