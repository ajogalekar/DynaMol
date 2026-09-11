// Run: node docs/audit/check_plots.cjs. Analytic plot checks, no browser required.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const crypto = require('node:crypto');
const { execFileSync } = require('node:child_process');
const root = path.resolve(__dirname, '../..');
const temp = fs.mkdtempSync(path.join(os.tmpdir(), 'dynamol-plot-review-'));
const bundle = path.join(temp, 'plot.cjs');
execFileSync(path.join(root, 'frontend/node_modules/.bin/esbuild'), [path.join(root, 'frontend/src/components/PlotPanel.tsx'), '--bundle', '--platform=node', '--format=cjs', '--outfile=' + bundle], { stdio: 'pipe' });
const { createPlotData } = require(bundle);
const result = [];
const measurement = (values, times, kind = 'distance') => ({ id: 'fixture', kind, atoms: [], label: 'Analytic fixture', unit: kind === 'dihedral' ? '°' : 'Å', values, times_ps: times, color: '#fff' });
function check(name, action) {
  try { result.push({ name, status: 'pass', evidence: action() }); }
  catch (error) { result.push({ name, status: 'fail', error: String(error) }); }
}
check('Nonuniform timestamps determine x positions and nearest-time seek', () => {
  const plot = createPlotData(measurement([1, 2, 3], [0, 1, 10]));
  assert.ok(Math.abs((plot.x(1) - plot.x(0)) / (plot.x(2) - plot.x(0)) - .1) < 1e-10);
  assert.equal(plot.nearestIndex(.5), 1);
  assert.equal(plot.nearestIndex(.7), 2);
  return { middleFrameFraction: .1, seekAt50Percent: 1, seekAt70Percent: 2 };
});
check('Torsion mean and segments respect the branch cut', () => {
  const plot = createPlotData(measurement([179, -179], [0, 1], 'dihedral'));
  assert.ok(Math.abs(Math.abs(plot.mean) - 180) < 1e-10);
  assert.equal(plot.segments.length, 2);
  return { circularMeanDegrees: plot.mean, segments: plot.segments.length };
});
check('A zero circular resultant has no reported mean', () => {
  const plot = createPlotData(measurement([0, 180], [0, 1], 'dihedral'));
  assert.equal(plot.mean, null);
  return { circularMean: null };
});
check('Nonmonotonic times explicitly fall back to frame index', () => {
  const plot = createPlotData(measurement([1, 2, 3], [0, 2, 1]));
  assert.equal(plot.timestampsUsable, false);
  assert.ok(plot.axisCaption.startsWith('Frame index'));
  assert.equal(plot.nearestIndex(.5), 1);
  return { caption: plot.axisCaption };
});
check('Unknown physical time uses a frame label', () => {
  const plot = createPlotData(measurement([1, 2, 3], [0, 5, 10]), 'frame');
  assert.equal(plot.axisCaption, 'Original frame index (time not supplied)');
  return { caption: plot.axisCaption };
});
check('Single-frame plot remains finite and seekable', () => {
  const plot = createPlotData(measurement([1], [2]));
  assert.ok(Number.isFinite(plot.x(0)));
  assert.equal(plot.nearestIndex(.9), 0);
  assert.equal(plot.ticks.length, 1);
  return { ticks: plot.ticks.length };
});
fs.rmSync(temp, { recursive: true, force: true });
const report = {
  recorded_at: new Date().toISOString(),
  scope: 'Analytic chart calculations only; these tests do not assess MD convergence.',
  node_version: process.version,
  source_sha256: crypto.createHash('sha256').update(fs.readFileSync(path.join(root, 'frontend/src/components/PlotPanel.tsx'))).digest('hex'),
  checks: result,
  passed: result.filter(item => item.status === 'pass').length,
  failed: result.filter(item => item.status === 'fail').length,
};
fs.writeFileSync(path.join(__dirname, 'plot-checks.json'), JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify(report, null, 2));
process.exitCode = report.failed ? 1 : 0;
