"""An FHI-aims jobflow runner."""

from __future__ import annotations

import json
import logging
import os
import shlex
import subprocess
from os.path import expandvars
from typing import TYPE_CHECKING, Any

from ase.calculators.aims import Aims, AimsProfile
from ase.calculators.socketio import SocketIOCalculator
from custodian import Custodian
from custodian.aims.handlers import RoundChargeErrorHandler, UnconvergedErrorHandler
from custodian.aims.jobs import AimsJob
from custodian.aims.validators import AimsOutputValidator
from emmet.core.types.enums import ValueEnum
from monty.json import MontyDecoder
from pymatgen.io.ase import AseAtomsAdaptor

from atomate2 import SETTINGS

if TYPE_CHECKING:
    from collections.abc import Sequence

    from custodian.custodian import ErrorHandler, Validator
    from pymatgen.core import Molecule, Structure

    from atomate2.aims.schemas.task import AimsTaskDoc

DEFAULT_HANDLERS = (UnconvergedErrorHandler(), RoundChargeErrorHandler())
_DEFAULT_VALIDATORS = (AimsOutputValidator(),)

logger = logging.getLogger(__name__)


class JobType(ValueEnum):
    """
    Type of FHI-aims job.

    - ``DIRECT``: Run FHI-aims without using custodian.
    - ``NORMAL``: Normal custodian :obj:`.AimsJob`.
    """

    DIRECT = "direct"
    NORMAL = "normal"


def run_aims(
    aims_cmd: str = None,
    job_type: JobType | str = JobType.DIRECT,
    max_errors: int = SETTINGS.AIMS_CUSTODIAN_MAX_ERRORS,
    scratch_dir: str = SETTINGS.CUSTODIAN_SCRATCH_DIR,
    handlers: Sequence[ErrorHandler] = DEFAULT_HANDLERS,
    validators: Sequence[Validator] = _DEFAULT_VALIDATORS,
    aims_job_kwargs: dict[str, Any] = None,
    custodian_kwargs: dict[str, Any] = None,
) -> None:
    """
    Run FHI-aims.

    Supports running FHI-aims with or without custodian (see :obj:`JobType`).

    Parameters
    ----------
    aims_cmd : str
        The command used to run FHI-aims (defaults to SETTINGS.AIMS_CMD).
    job_type : str or .JobType
        The job type. Defaults to DIRECT, i.e. without custodian.
    max_errors : int
        The maximum number of errors allowed by custodian.
    scratch_dir : str
        The scratch directory used by custodian.
    handlers : list of .ErrorHandler
        The error handlers used by custodian.
    validators : list of .Validator
        The validators handlers used by custodian.
    aims_job_kwargs : dict
        Keyword arguments that are passed to :obj:`.AimsJob`.
    custodian_kwargs : dict
        Keyword arguments that are passed to :obj:`.Custodian`.
    """
    if aims_cmd is None:
        aims_cmd = SETTINGS.AIMS_CMD

    aims_cmd = expandvars(aims_cmd)

    if job_type == JobType.DIRECT:
        logger.info(f"Running command: {aims_cmd}")
        return_code = subprocess.call(["/bin/bash", "-c", aims_cmd], env=os.environ)
        logger.info(f"{aims_cmd} finished running with return code: {return_code}")
        return

    aims_job_kwargs = aims_job_kwargs or {}
    custodian_kwargs = custodian_kwargs or {}

    # AimsJob redirects the output itself: drop the "> aims.out" of AIMS_CMD
    split_aims_cmd = shlex.split(aims_cmd)
    if ">" in split_aims_cmd:
        idx = split_aims_cmd.index(">")
        aims_job_kwargs.setdefault("output_file", split_aims_cmd[idx + 1])
        split_aims_cmd = split_aims_cmd[:idx]

    if job_type == JobType.NORMAL:
        jobs = [AimsJob(split_aims_cmd, **aims_job_kwargs)]
    else:
        raise ValueError(f"Unsupported {job_type=}")

    custodian_manager = Custodian(
        handlers,
        jobs,
        validators=validators,
        max_errors=max_errors,
        scratch_dir=scratch_dir,
        **custodian_kwargs,
    )

    logger.info("Running FHI-aims using custodian.")
    custodian_manager.run()


def should_stop_children(
    task_document: AimsTaskDoc,
    handle_unsuccessful: bool | str = True,
) -> bool:
    """
    Decide whether child jobs should continue.

    Parameters
    ----------
    task_document : .TaskDoc
        An FHI-aims task document.
    handle_unsuccessful : bool or str
        This is a three-way toggle on what to do if your job looks OK, but is actually
        not converged (either electronic or ionic):

        - `True`: Mark job as completed, but stop children.
        - `False`: Do nothing, continue with workflow as normal.
        - `"error"`: Throw an error.

    Returns
    -------
    bool
        Whether to stop child jobs.
    """
    if task_document.state == "successful":
        return False

    if isinstance(handle_unsuccessful, bool):
        return handle_unsuccessful

    if handle_unsuccessful == "error":
        raise RuntimeError("Job was not successful (not converged)!")

    raise RuntimeError(f"Unknown option for handle_unsuccessful: {handle_unsuccessful}")


def run_aims_socket(
    structures_to_calculate: list[Structure | Molecule], aims_cmd: str = None
) -> None:
    """Use the ASE interface to run FHI-aims from the socket.

    Parameters
    ----------
    structures_to_calculate: list[Structure or Molecule]
        The list of structures to run scf calculations on
    aims_cmd: str
        The aims command to use (defaults to SETTINGS.AIMS_CMD).
    """
    if aims_cmd is None:
        aims_cmd = SETTINGS.AIMS_CMD

    ase_adaptor = AseAtomsAdaptor()
    atoms_to_calculate = [
        ase_adaptor.get_atoms(structure) for structure in structures_to_calculate
    ]

    with open("parameters.json") as param_file:
        parameters = json.load(param_file, cls=MontyDecoder)

    del_keys = []
    for key, val in parameters.items():
        if val is None:
            del_keys.append(key)

    for key in del_keys:
        parameters.pop(key)

    if aims_cmd:
        command = aims_cmd
    elif "command" in parameters:
        command = parameters["command"]
    else:
        command = SETTINGS.AIMS_CMD

    parameters["profile"] = AimsProfile(command=command)

    calculator = Aims(**parameters)
    port = parameters["use_pimd_wrapper"][1]
    atoms = atoms_to_calculate[0].copy()

    with SocketIOCalculator(calc=calculator, port=port) as calc:
        for atoms_calc in atoms_to_calculate:
            # Delete prior calculation results
            calc.results.clear()

            # Reset atoms information to the new cell
            atoms.info = atoms_calc.info
            atoms.cell = atoms_calc.cell
            atoms.positions = atoms_calc.positions

            calc.calculate(atoms, system_changes=["positions", "cell"])

        calc.close()
