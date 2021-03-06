import os
from typing import List, Optional
import subprocess
import logging

from ligbinder.settings import SETTINGS


logger = logging.getLogger(__name__)


class AmberMDEngine:
    def __init__(
        self,
        path,
        crd_file: Optional[str] = None,
        top_file: Optional[str] = None,
        trj_file: Optional[str] = None,
        rst_file: Optional[str] = None,
        log_file: Optional[str] = None,
        ref_file: Optional[str] = None,
        inp_file: Optional[str] = None,
        steps=250000,
        tstep=4.0,
        use_gpu=True,
        use_hmr=True,
        apply_restraints=True,
        restraint_force=1.0
    ) -> None:

        self.crd_file = os.path.join(
            path, crd_file if crd_file is not None else SETTINGS["md"]["crd_file"]
        )
        self.trj_file = os.path.join(
            path, trj_file if trj_file is not None else SETTINGS["md"]["trj_file"]
        )
        self.top_file = os.path.join(
            path, top_file if top_file is not None else SETTINGS["md"]["top_file"]
        )
        self.rst_file = os.path.join(
            path, rst_file if rst_file is not None else SETTINGS["md"]["rst_file"]
        )
        self.log_file = os.path.join(
            path, log_file if log_file is not None else SETTINGS["md"]["log_file"]
        )
        self.ref_file = os.path.join(
            path, ref_file if ref_file is not None else SETTINGS["md"]["ref_file"]
        )
        self.inp_file = os.path.join(
            path, inp_file if inp_file is not None else SETTINGS["md"]["inp_file"]
        )

        self.steps = steps
        self.tstep = tstep

        self.use_gpu = use_gpu
        self.binary = "pmemd.cuda" if self.use_gpu else "sander"
        self.apply_restraints = apply_restraints
        self.restraint_force = restraint_force
        self.use_hmr = use_hmr

    def write_input(self):
        interval = self.steps // 10
        restraints = f'restraint_wt={self.restraint_force}, restraintmask=\'{SETTINGS["system"]["restraint_mask"]}\','
        restraints = restraints if self.apply_restraints else ""
        lines = [
            "#  Constant Volume",
            "&cntrl",
            "ntx=1, irest=0, iwrap=1,",
            f"ntxo=2, ntpr={interval}, ntwx={interval}, ntwv=0, ntwe=0, ioutfm=1,",
            f"nstlim={self.steps}, dt={self.tstep/1000},",
            "ntc=2, ntf=2,",
            "ntb=1, cut=9.0,",
            "ntt=3, gamma_ln=4.0, ig=-1,",
            "temp0=300,",
            restraints,
            "&end",
            ""
        ]
        lines = [line for line in lines if line is not None]
        msg = "\n".join(lines)
        with open(self.inp_file, "w") as file:
            file.write(msg)

    def run(self):
        self.write_input()
        cmd = self._get_command()
        logger.info('Running md engine')
        logger.debug(f'{" ".join(cmd)}')
        return subprocess.run(self._get_command(), check=True)

    def _get_command(self) -> List[str]:
        command = [
            self.binary,
            "-O",
            "-i", f"{self.inp_file}",
            "-o", f"{self.log_file}",
            "-p", f"{self.top_file}",
            "-c", f"{self.crd_file}",
            "-x", f"{self.trj_file}",
            "-r", f"{self.rst_file}"
        ]
        command += ["-ref", f"{self.ref_file}"] if self.apply_restraints else []
        return command
