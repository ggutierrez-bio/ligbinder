import os
from typing import List
import subprocess


class AmberMDEngine:
    def __init__(
        self,
        path,
        crd_file: str = "initial.rst7",
        top_file: str = "top.prmtop",
        trj_file: str = "traj.nc",
        rst_file: str = "final.rst7",
        log_file: str = "md.log",
        ref_file: str = "ref.crd",
        inp_file: str = "md.in",
        steps = 250000,
        tstep = 4.0,
        use_gpu=True
    ) -> None:

        self.crd_file = os.path.join(path, crd_file)
        self.top_file = os.path.join(path, top_file)
        self.trj_file = os.path.join(path, trj_file)
        self.rst_file = os.path.join(path, rst_file)
        self.log_file = os.path.join(path, log_file)
        self.ref_file = os.path.join(path, ref_file)
        self.inp_file = os.path.join(path, inp_file)
        
        self.steps = steps
        self.tstep = tstep

        self.use_gpu = use_gpu
        self.binary = "pmemd.cuda" if self.use_gpu else "sander"

    def write_input(self):
        interval=self.steps//10
        lines = [
            "#  Constant Volume",
            "&cntrl",
            "ntx=5, irest=0, iwrap=1,",
            f"ntxo=2, ntpr={interval}, ntwx={interval}, ntwv=0, ntwe=0, ioutfm=1,",
            f"nstlim={self.steps}, dt={self.tstep/1000},",
            "ntc=2, ntf=2,",
            "ntb=1, cut=9.0,",
            "ntt=3, gamma_ln=4.0, ig=-1,",
            "temp0=300, ",
            "&end",
        ]
        with open(self.inp_file,'w') as file:
            file.writelines(lines)

    def run(self):
        self.write_input()
        return subprocess.run(self._get_command(), check=True)
    
    def _get_command(self) -> List[str]:
        return [
            self.binary,
            "-O",
            f"-i {self.inp_file}",
            f"-o {self.log_file}",
            f"-p {self.top_file}",
            f"-c {self.crd_file}",
            f"-r {self.rst_file}",
            f"-x {self.trj_file}"
        ]