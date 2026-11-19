from pathlib import Path
import json
import numpy as np
from worker import Simulation,main

root=Path.cwd();cfg=json.loads((root/'config.json').read_text());v=json.loads((root/'validation.json').read_text())
assert v['passed'] and v['initial']==cfg['initial'] and v['selected_divisor']==cfg['dt_divisor']
if not (root/'restart_gate.json').exists():
 s=Simulation(cfg,root/'gate')
 for _ in range(3):s.advance()
 s.save();s.advance();expected=s.h.copy()
 restored=Simulation(cfg,root/'gate');restored.advance();error=float(np.max(abs(expected-restored.h)))
 assert error<1e-10
 (root/'restart_gate.json').write_text(json.dumps(dict(max_field_error=error,passed=True)))
raise SystemExit(main())
