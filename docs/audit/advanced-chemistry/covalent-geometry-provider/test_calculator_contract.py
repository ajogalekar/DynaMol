"""No-compute coverage for native worker versus supervisor classification."""
import pytest
from calculator import is_quantum_worker

@pytest.mark.parametrize('command',[
 '/runtime/python /project/backend/qm_worker.py request.json output',
 '/runtime/python3.12 -u /cache/qm_worker.py request.json output',
 '/runtime/python -X dev /project/covalent_reference_worker.py input output',
 '/runtime/python /project/qualify_calculator.py',
 '/runtime/python /project/geometric_adapter.py --worker request.json output',
])
def test_actual_quantum_workers_are_counted(command):
    assert is_quantum_worker(command)

@pytest.mark.parametrize('command',[
 '/runtime/python /project/test_geometric_adapter.py',
 '/runtime/python /project/test_geometric_adapter.py /project/geometric_adapter.py --worker',
 '/runtime/python /project/geometric_adapter.py --preflight request.json',
 '/runtime/python -m pytest /project/test_geometric_adapter.py',
 '/runtime/python -c "print(\'qm_worker.py\')"',
 '/bin/bash -lc "/runtime/python /project/qm_worker.py input output"',
 '/runtime/python -',
])
def test_idle_supervisors_and_test_scripts_are_not_quantum_workers(command):
    assert not is_quantum_worker(command)

@pytest.mark.parametrize('extra_worker,expected_rejection',[(False,False),(True,True)])
def test_execution_guard_ignores_supervisor_but_counts_true_worker(tmp_path,monkeypatch,extra_worker,expected_rejection):
    import calculator
    rows=['101 /runtime/python /project/qm_worker.py input output',
          '102 /runtime/python /project/qm_worker.py another output',
          '103 /runtime/python /project/test_geometric_adapter.py']
    if extra_worker:rows.append('104 /runtime/python /project/geometric_adapter.py --worker input output')
    monkeypatch.setattr(calculator.subprocess,'check_output',lambda *args,**kwargs:'\n'.join(rows))
    (tmp_path/'build/advanced-chemistry').mkdir(parents=True)
    if expected_rejection:
        with pytest.raises(ValueError,match='No free quantum execution slot'):
            with calculator.execution_guard(tmp_path):pytest.fail('Fourth recognized worker admitted')
    else:
        with calculator.execution_guard(tmp_path) as record:
            assert len(record['active_other_quantum_workers'])==2
            assert record['native_threads']=={'tblite_openmp_max_threads':1,'tblite_openblas_threads':1}
