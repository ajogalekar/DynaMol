"""Build native RESP with its four large work matrices on the heap.

The numerical source and original dimension limits stay unchanged. This avoids
the macOS dyld failure caused by approximately 2 GB of static COMMON storage.
Use a separately extracted, pinned AmberTools source tree; validation is separate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess


MODULE = '''module resp_heap_workspace
  implicit none
  double precision, allocatable :: apot(:,:), awt(:,:)
  double precision, allocatable :: a(:,:), awork(:,:)
contains
  subroutine initialize_resp_heap(n)
    integer, intent(in) :: n
    integer :: allocation_status
    allocate(apot(n,n), awt(n,n), a(n,n), awork(n,n), &
             source=0.0d0, stat=allocation_status)
    if (allocation_status /= 0) then
      write(*,*) 'RESP: failed to allocate the charge-fitting workspace'
      stop 1
    end if
  end subroutine initialize_resp_heap
end module resp_heap_workspace
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def transform(source):
    """Change storage declarations only, with per-routine exact scoping."""
    lines = source.splitlines(keepends=True)
    statements = []
    for line in lines:
        if len(line) > 5 and line[0] not in 'cC*!#' and line[5] not in ' 0':
            if not statements:
                raise ValueError('Continuation without preceding statement')
            statements[-1].append(line)
        else:
            statements.append([line])
    routines, current = [], []
    for statement in statements:
        first = statement[0]
        if re.match(r'^\s*(program|subroutine)\s+', first, re.I):
            if current:
                routines.append(current)
            current = []
        current.append(statement)
    if current:
        routines.append(current)
    output, edits = [], []
    pattern = re.compile(r'\b(apot|awt|a|awork)\s*\(maxq,maxq\)\s*,\s*', re.I)
    for routine in routines:
        imported, rebuilt = set(), []
        for statement in routine:
            text = ''.join(statement)
            if re.match(r'^\s*common\s*/(?:ESPCOM|CALCUL|worker)/', statement[0], re.I):
                logical = statement[0][6:].strip()
                logical += ' '.join(line[6:].strip() for line in statement[1:])
                names = {match.group(1).lower() for match in pattern.finditer(logical)}
                if names:
                    logical = pattern.sub('', logical)
                    imported.update(names)
                    text = '      ' + logical + '\n'
                    edits.append({'routine': routine[0][0].strip(), 'matrices': sorted(names)})
            rebuilt.append(text)
        if imported:
            uses = sorted(imported)
            if re.match(r'^\s*program\s+resp\s*$', routine[0][0], re.I):
                uses.append('initialize_resp_heap')
            rebuilt.insert(1, '      use resp_heap_workspace, only: ' + ', '.join(uses) + '\n')
        output.extend(rebuilt)
    result = ''.join(output)
    if result.count('      call filein\n') != 1:
        raise ValueError('Unexpected RESP main entry; source version must be reviewed')
    # Preserve static COMMON's initial zero values before any numerical work.
    result = result.replace('      do jn = 1,maxq\n', '      call initialize_resp_heap(maxq)\n      do jn = 1,maxq\n', 1)
    if not edits or any(name not in {entry for row in edits for entry in row['matrices']}
                        for name in ('apot', 'awt', 'a', 'awork')):
        raise ValueError('Expected four RESP work matrices were not found')
    return result, edits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True, help='AmberTools/src/etc')
    parser.add_argument('--compiler', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, compiler, output = args.source.resolve(), args.compiler.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = {name: source / name for name in ('resp.F', 'lapack.F', 'limits.h')}
    rewritten, edits = transform(paths['resp.F'].read_text())
    (output / 'resp.F').write_text(rewritten)
    (output / 'resp_heap_workspace.f90').write_text(MODULE)
    for name in ('lapack.F', 'limits.h'):
        (output / name).write_bytes(paths[name].read_bytes())
    command = [str(compiler), '-O0', '-cpp', '-std=legacy', '-fallow-argument-mismatch',
               '-ffixed-line-length-none', 'resp_heap_workspace.f90', 'resp.F', 'lapack.F', '-o', 'resp']
    completed = subprocess.run(command, cwd=output, capture_output=True, text=True)
    (output / 'build.log').write_text(completed.stdout + completed.stderr)
    report = {'source': {name: {'path': str(path), 'sha256': sha(path)} for name, path in paths.items()},
              'builder_sha256': sha(__file__), 'compiler': str(compiler),
              'compiler_version': subprocess.run([str(compiler), '--version'], capture_output=True, text=True).stdout,
              'command': command, 'returncode': completed.returncode, 'storage_edits': edits,
              'scope': 'Only four maxq×maxq COMMON matrices become zero-initialized allocatable arrays; fitting equations, native LAPACK, numeric constants and limits remain unchanged.',
              'validated': False}
    if completed.returncode == 0:
        report['binary_sha256'] = sha(output / 'resp')
    (output / 'build.json').write_text(json.dumps(report, indent=2) + '\n')
    if completed.returncode:
        raise RuntimeError('RESP build failed; inspect build.log')
    print(output / 'resp')


if __name__ == '__main__':
    main()
