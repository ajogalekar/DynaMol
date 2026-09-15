"""Write a fresh reproducible synthetic integration audit; no native imports."""
import ast
import copy
import difflib
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
ORIGINAL_SHA='ddb67a1f8e9be6ed05c3c95e9864e20b79f0233604e27926e30e5ffc09295bc1'
HELPER_SHA='b3a0543b5b2394373e7a787ad01ecd77865cc507f85bec8ec99efe6efed26528'


def main():
    out=HERE/(sys.argv[1] if len(sys.argv)>1 else 'mock-audit')
    out.mkdir(exist_ok=False)
    original=(ROOT/'backend/qm_worker.py').read_bytes()
    candidate=(HERE/'qm_worker.py').read_bytes()
    helper=(HERE/'qm_evaluation_journal.py').read_bytes()
    assert hashlib.sha256(original).hexdigest()==ORIGINAL_SHA
    assert hashlib.sha256(helper).hexdigest()==HELPER_SHA
    before,after=ast.parse(original),ast.parse(candidate)
    reduced=copy.deepcopy(after)
    additions=[node for node in reduced.body if isinstance(node,ast.FunctionDef) and node.name=='_create_evaluation_journal']
    assert len(additions)==1
    reduced.body.remove(additions[0])
    allowed={
        'record_evaluation = _create_evaluation_journal(data, elements, initial_bohr, mean_field, out, report, lib.param.BOHR)',
        'evaluation_record = record_evaluation(env)',
        "if evaluation_record is not None:\n    step['evaluation_record'] = evaluation_record",
    }
    removed=[]

    class RemoveHooks(ast.NodeTransformer):
        def visit_Assign(self,node):
            if ast.unparse(node) in allowed:
                removed.append(ast.unparse(node));return None
            return self.generic_visit(node)

        def visit_If(self,node):
            if ast.unparse(node) in allowed:
                removed.append(ast.unparse(node));return None
            return self.generic_visit(node)

    reduced=RemoveHooks().visit(reduced)
    assert len(removed)==3 and set(removed)==allowed
    assert ast.dump(before,include_attributes=False)==ast.dump(reduced,include_attributes=False)
    (out/'source.diff').write_text(''.join(difflib.unified_diff(original.decode().splitlines(True),
        candidate.decode().splitlines(True),fromfile='backend/qm_worker.py',tofile='candidate/qm_worker.py')))
    ast_review={'original_worker_sha256':ORIGINAL_SHA,'new_top_level_functions':['_create_evaluation_journal'],
        'added_original_worker_statements':removed,'original_ast_identical_after_removing_only_observer_and_hooks':True,
        'numerical_method_and_acceptance_code_changes':False}
    (out/'ast-review.json').write_text(json.dumps(ast_review,indent=2)+'\n')
    result=subprocess.run([sys.executable,str(HERE/'test_integration.py'),'-v'],capture_output=True,text=True,timeout=30)
    (out/'tests.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0 and 'Ran 20 tests' in result.stderr
    spec=importlib.util.spec_from_file_location('integration_fixture',HERE/'test_integration.py')
    fixture=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixture)
    case=fixture.IntegrationTests();case.setUp()
    try:
        callback=fixture.existing_callback(case.record,case.data,case.report,lambda *a,**kw:None)
        callback(case.env)
        destination=out/'synthetic-evaluation';destination.mkdir()
        shutil.copytree(case.out/'evaluation-journal',destination/'evaluation-journal')
        (destination/'input.json').write_text(json.dumps(case.input))
        (destination/'report.json').write_text(json.dumps(case.report,indent=2)+'\n')
        (destination/'README.txt').write_text('Synthetic mock data only. No energy or gradient was computed. Not a molecular model.\n')
    finally:
        case.doCleanups()
    report={'schema_version':1,'status':'mock_integration_passed','test_count':20,'test_exit_code':result.returncode,
        'candidate_worker_sha256':hashlib.sha256(candidate).hexdigest(),'original_worker_sha256':ORIGINAL_SHA,
        'journal_helper_sha256':HELPER_SHA,'journal_helper_byte_identical_to_reviewed_35_test_version':True,
        'test_source_sha256':hashlib.sha256((HERE/'test_integration.py').read_bytes()).hexdigest(),
        'production_worker_changed':False,'native_qm_calculations':0,'native_qualification':'not_performed',
        'active_integration':False,'accepted':False,'physical_acceptance':False,
        'scope':'Completed callback observer only; candidate ready for independent review before later small native qualification'}
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    files=sorted(p for p in out.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files += [HERE/name for name in ['qm_worker.py','qm_evaluation_journal.py','test_integration.py','run_mock_audit.py','README.md']]
    manifest={'schema_version':1,'artifacts':[{'file':str(p.relative_to(HERE)),
        'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in files]}
    (out/'artifact-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
