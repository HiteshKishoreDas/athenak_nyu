"""Independent exact arithmetic demonstrations; not native runtime tests."""
from fractions import Fraction as F
import json

stage_flux = [F(1), F(2), F(4)]
true_rk3 = sum(w*f for w,f in zip([F(1,6),F(1,6),F(2,3)],stage_flux))
base_saved = sum(stage_flux)/3
assert true_rk3 == F(19,6) and base_saved == F(7,3)
key_a = 104729*7919 + 0*104729
key_b = 0*7919 + 7919*104729
assert key_a == key_b
# The first face interval covers [0,.8), the second [.8,1): exact face bias.
requested = [F(4,5),F(4,5)]
realized = [F(4,5),F(1,5)]
assert sum(requested) > 1 and sum(realized) == 1
result = dict(status='source_arithmetic_demonstrated_not_native_runtime',
    rk3=dict(stage_flux=list(map(str,stage_flux)),correct=str(true_rk3),
             baseline_stage_average=str(base_saved)),
    rng_key_alias=dict(tag_cycle_pairs=[[104729,0],[0,7919]],key=key_a),
    overoutflow=dict(requested=list(map(str,requested)),
                    cumulative_uniform_draw_realized=list(map(str,realized))))
print(json.dumps(result,indent=2))
