"""Meaningful offline checks of fixture and legacy output-reader failure handling."""
from pathlib import Path
import struct
import tempfile

from common import edit_input, read_thermo
from qualification import rh_states


def rejected(path):
    try:
        read_thermo(path, allow_identical_replay=True)
    except ValueError:
        return
    raise AssertionError('malformed file accepted')


def main():
    # New sections remain safely editable: no repeated header/overriding old value.
    text = edit_input("<job>\nbasename = fixture\n", {
        "output4/file_type": "rst", "output4/dt": -1})
    text = edit_input(text, {"output4/dt": 1})
    assert text.count("<output4>") == 1 and "dt = 1" in text and "dt = -1" not in text

    text='<mesh>\nnx1 = 64\n<time>\ntlim = .25\n'
    changed=edit_input(text,{'mesh/nx1':128,'mesh/nx2':4,'output1/dt':1})
    assert 'nx1 = 128' in changed and 'nx2 = 4' in changed and '<output1>\ndt = 1' in changed
    for mach in (1.01,1.2,2,5,100):
        for gamma in (1.4,5/3):
            rh_states(mach,gamma)
            rh_states(mach,gamma,speed=-.5,reflected=True)
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/'sample.thp'
        header=struct.pack('=32s4i',b'ATHK_PRTCL_THERMO_HISTORY',2,8,1,3)+b'rho'
        block=struct.pack('=16s5id',b'ATHKTHPBLK',2,1,4,5,7,.125)
        values=struct.pack('=4i5d',7,123,1,0,.125,.2,.3,.4,1)
        path.write_bytes(header+block+values)
        rows,times=read_thermo(path)
        assert rows[7,123]['fields']=={'rho':1} and times=={7:.125}
        path.write_bytes(header+block+values+block+values)
        try:read_thermo(path)
        except ValueError:pass
        else:raise AssertionError('strict reader accepted duplicate identity')
        rows,times=read_thermo(path,allow_identical_replay=True)
        assert len(rows)==1 and rows.identical_replay_rows==1
        changed_values=struct.pack('=4i5d',7,123,1,0,.125,.21,.3,.4,1)
        path.write_bytes(header+block+values+block+changed_values);rejected(path)
        path.write_bytes(header+block+values[:-1]);rejected(path)
        path.write_bytes(header+block[:-1]);rejected(path)
        path.write_bytes(b'bad');rejected(path)
    print('offline MC fixture/RH/reader adversarial checks passed')


if __name__=='__main__':main()
