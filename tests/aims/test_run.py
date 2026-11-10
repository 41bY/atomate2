import shutil
from pathlib import Path

from pyfhiaims.control.control import AimsControl

from atomate2.aims.run import JobType, run_aims
from atomate2.aims.schemas.task import AimsTaskDoc


def test_run_aims_custodian(test_dir):
    ref_dir = test_dir / "aims" / "static-si"
    for file in (ref_dir / "inputs").iterdir():
        shutil.copy(file, ".")

    # zcat stands in for FHI-aims, printing the output of a converged run
    run_aims(f"zcat {ref_dir}/outputs/aims.out.gz > aims.out", job_type=JobType.NORMAL)

    assert Path("control.in.orig").is_file()
    task_doc = AimsTaskDoc.from_directory(Path.cwd())
    assert task_doc.custodian[0]["corrections"] == []
    # set for ScfMonitorHandler, one of the default handlers
    parameters = AimsControl.from_file("control.in").parameters
    assert parameters["sc_iter_limit"] == "1000"
    assert parameters["elsi_restart"] == "read_and_write 100"
