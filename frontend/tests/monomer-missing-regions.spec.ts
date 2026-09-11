import { test } from './testWorkspace';
import { expect, type Locator, type Page } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

const fixtures = fileURLToPath(new URL('../../docs/audit/preparation-fixtures/', import.meta.url));

async function openStudio(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  return page.getByRole('dialog', { name: 'Set molecules in motion.' });
}

async function upload(page: Page, studio: Locator, buffer: Buffer, label: string) {
  const received = page.waitForResponse(
    (response) =>
      response.url().endsWith('/api/structures/upload') && response.request().method() === 'POST',
  );
  await studio
    .getByLabel('Upload simulation structure', { exact: true })
    .setInputFiles({ name: `${label}-${Date.now()}.pdb`, mimeType: 'chemical/x-pdb', buffer });
  const response = await received;
  expect(response.ok(), await response.text()).toBe(true);
  const dataset = await response.json();
  await expect(studio.locator('.source-summary strong')).toHaveText(dataset.name);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  return dataset;
}

async function twoChains() {
  const text = await readFile(`${fixtures}/six_residues_intact.pdb`, 'utf8');
  const atoms = text.split('\n').filter((line) => line.startsWith('ATOM  '));
  const chainB = atoms.map((line, index) => {
    const serial = String(index + atoms.length + 2).padStart(5);
    const x = (Number(line.slice(30, 38)) + 35).toFixed(3).padStart(8);
    return (
      line.slice(0, 6) + serial + line.slice(11, 21) + 'B' + line.slice(22, 30) + x + line.slice(38)
    );
  });
  return Buffer.from(
    [
      'SEQRES   1 A    6  MET GLN ILE PHE VAL LYS',
      'SEQRES   1 B    6  MET GLN ILE PHE VAL LYS',
      ...atoms,
      'TER',
      ...chainB,
      'TER',
      'END',
      '',
    ].join('\n'),
  );
}

function gate() {
  let release!: () => void;
  const promise = new Promise<void>((resolve) => {
    release = resolve;
  });
  return { release, promise };
}

test('missing-atom and known-loop choices are visible before Prep while detailed options stay closed', async ({
  page,
}, info) => {
  const studio = await openStudio(page);
  await upload(
    page,
    studio,
    await readFile(`${fixtures}/six_residues_known_gap.pdb`),
    'visible-loop',
  );
  const repairs = studio.getByRole('region', { name: 'Missing atoms & residues' });
  const addAtoms = repairs.getByRole('checkbox', { name: 'Add missing heavy atoms', exact: true });
  const build = repairs.getByRole('checkbox', {
    name: 'Build supported missing loops / residues',
    exact: true,
  });
  await expect(addAtoms).toBeVisible();
  await expect(addAtoms).toBeChecked();
  await expect(build).toBeEnabled();
  await expect(build).not.toBeChecked();
  await expect(
    studio.getByRole('button', { name: 'Preparation options', exact: true }),
  ).toHaveAttribute('aria-expanded', 'false');
  await expect(studio.locator('.prep-options')).toHaveCount(0);
  await expect(repairs.locator('.missing-residues')).toContainText('ILE');
  await expect(repairs).toContainText('can build a starting model');
  await expect(repairs).toContainText('6 residues per internal gap and 12 residues total');
  const beforePrep = await studio.evaluate(
    (element) =>
      !!(
        element
          .querySelector('.missing-structure-setup')!
          .compareDocumentPosition(element.querySelector('.prep-action-row')!) &
        Node.DOCUMENT_POSITION_FOLLOWING
      ),
  );
  expect(beforePrep).toBe(true);
  await build.check();
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  await studio.getByRole('button', { name: 'Preparation options', exact: true }).click();
  await expect(
    studio.getByRole('checkbox', { name: 'Add missing heavy atoms', exact: true }),
  ).toHaveCount(1);
  await expect(
    studio.getByRole('checkbox', { name: 'Build supported missing loops / residues', exact: true }),
  ).toHaveCount(1);
  await repairs.scrollIntoViewIfNeeded();
  await page.screenshot({ path: info.outputPath('missing-regions-before-prep.png') });
});

test('one-click suggested chain and alternate chain extraction are reversible and preserve the original structure', async ({
  page,
  request,
}, info) => {
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  const studio = await openStudio(page);
  const source = await upload(page, studio, await twoChains(), 'two-chain-source');
  const choicesResponse = await request.get(`/api/datasets/${source.id}/monomers`);
  expect(choicesResponse.ok(), await choicesResponse.text()).toBe(true);
  const choices = await choicesResponse.json();
  expect(choices.chains).toHaveLength(2);
  const chooser = studio.getByRole('region', { name: 'Choose one protein chain' });
  await expect(chooser.getByRole('combobox', { name: 'Protein chain for monomer' })).toHaveValue(
    String(choices.recommended_chain_index),
  );
  await expect(chooser).toContainText(
    'chain count alone does not establish the biological assembly',
  );
  const tracking = studio.getByRole('region', { name: 'Track during simulation', exact: true });
  await tracking.getByRole('button', { name: /^Track during simulation/ }).click();
  for (const chain of choices.chains) {
    const atoms = source.atoms
      .filter(
        (atom: { name: string; chain: string }) =>
          atom.name === 'CA' && atom.chain === chain.chain_id,
      )
      .slice(0, 2);
    expect(atoms).toHaveLength(2);
    await tracking.getByRole('button', { name: 'Add measurement', exact: true }).click();
    const editor = tracking.locator('.tracking-editor');
    await editor
      .locator('.measure-types')
      .getByRole('button', { name: /Distance$/ })
      .click();
    for (const atom of atoms) {
      await editor
        .getByRole('textbox', { name: 'Find an atom', exact: true })
        .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
      await editor
        .locator('.atom-search-results')
        .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
        .click();
    }
    await expect(editor).toHaveCount(0);
  }
  await expect
    .poll(
      async () => (await (await request.get('/api/workspace')).json()).state.measurements.length,
    )
    .toBe(2);
  const measurementsBefore = (await (await request.get('/api/workspace')).json()).state
    .measurements;
  // Seed only the history response to exercise dismissal without running preparation.
  const previousPreparation = {
    id: 'monomer-history-fixture',
    name: 'Completed preparation fixture',
    engine: 'preparation',
    status: 'completed',
    stage: 'Complete',
    progress: 100,
    completed_steps: 5,
    total_steps: 5,
    elapsed_seconds: 1,
    logs: [],
    config: { dataset_id: source.id },
    dataset_id: source.id,
    created_at: new Date().toISOString(),
  };
  await page.route('**/api/jobs', async (route) => {
    if (route.request().method() !== 'GET') return route.continue();
    const response = await route.fetch();
    await route.fulfill({ response, json: [previousPreparation, ...(await response.json())] });
  });
  await expect(studio.locator('.structure-job-monitor')).toContainText(
    'Protein preparation complete',
  );
  let nativePreparations = 0;
  page.on('request', (request) => {
    if (
      request.method() === 'POST' &&
      (request.url().endsWith('/api/preparations') || request.url().endsWith('/api/jobs'))
    )
      nativePreparations++;
  });
  let firstOutput: any;
  for (const chainIndex of [
    choices.recommended_chain_index,
    choices.chains.find(
      (chain: { index: number }) => chain.index !== choices.recommended_chain_index,
    ).index,
  ]) {
    if (firstOutput)
      await chooser
        .getByRole('combobox', { name: 'Protein chain for monomer' })
        .selectOption(String(chainIndex));
    const extracted = page.waitForResponse(
      (response) =>
        response.url().endsWith(`/api/datasets/${source.id}/monomer`) &&
        response.request().method() === 'POST',
    );
    await chooser.getByRole('button', { name: 'Use one monomer', exact: true }).click();
    const response = await extracted;
    expect(response.ok(), await response.text()).toBe(true);
    expect(response.request().postDataJSON()).toEqual({
      chain_index: chainIndex,
      keep_associated_molecules: true,
    });
    const result = await response.json();
    expect(result.monomer_selection.chain_index).toBe(chainIndex);
    expect(result.parent_dataset_id).toBe(source.id);
    expect(result.n_atoms).toBe(source.n_atoms / 2);
    expect(result.preparation).toBeUndefined();
    expect(result.solvation).toBeUndefined();
    await expect(page.locator('.structure-card h2')).toHaveText(result.name);
    await expect(studio.locator('.monomer-result')).toContainText('selected');
    await expect(chooser).toHaveCount(0);
    await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
    await expect(studio.locator('.structure-job-monitor')).toHaveCount(0);
    await expect(page.locator('.measurement-tab')).toHaveCount(firstOutput ? 0 : 1);
    await expect(
      page
        .getByRole('status')
        .filter({ hasText: 'Measurements outside the retained monomer were removed' }),
    ).toBeVisible();
    await expect
      .poll(async () => (await (await request.get('/api/workspace')).json()).state.dataset_id)
      .toBe(result.id);
    const stateAfter = (await (await request.get('/api/workspace')).json()).state;
    if (!firstOutput) {
      const retainedBefore = measurementsBefore.find((measurement: { atoms: number[] }) =>
        measurement.atoms.every((index) =>
          result.monomer_selection.retained_atom_indices.includes(index),
        ),
      );
      expect(retainedBefore).toBeTruthy();
      expect(stateAfter.measurements).toHaveLength(1);
      expect(stateAfter.measurements[0].id).toBe(retainedBefore.id);
      expect(stateAfter.measurements[0].atoms).toEqual(
        retainedBefore.atoms.map((index: number) =>
          result.monomer_selection.retained_atom_indices.indexOf(index),
        ),
      );
      expect(stateAfter.measurements[0].trackDuringRun).toBe(true);
    } else expect(stateAfter.measurements).toHaveLength(0);
    await expect(
      studio.locator('.job-card').filter({ hasText: previousPreparation.name }),
    ).toHaveCount(1);
    const original = await (await request.get(`/api/datasets/${source.id}`)).json();
    expect(original.n_atoms).toBe(source.n_atoms);
    expect(original.atoms).toEqual(source.atoms);
    if (!firstOutput) {
      firstOutput = result;
      await studio.getByRole('button', { name: 'Restore full structure', exact: true }).click();
      await expect(page.locator('.structure-card h2')).toHaveText(source.name);
      await expect(chooser).toBeVisible();
    }
  }
  await studio.locator('.monomer-result').scrollIntoViewIfNeeded();
  await page.screenshot({ path: info.outputPath('one-chain-selected.png') });
  expect(nativePreparations).toBe(0);
  expect(errors).toEqual([]);
});

test('pending monomer extraction locks Start and cannot overwrite a new setup', async ({
  page,
}) => {
  const studio = await openStudio(page);
  const source = await upload(page, studio, await twoChains(), 'stale-monomer');
  const chooser = studio.getByRole('region', { name: 'Choose one protein chain' });
  await expect(chooser).toBeVisible();
  await expect(studio.getByRole('button', { name: 'Start simulation', exact: true })).toBeEnabled();
  const held = gate();
  let responseReady!: () => void;
  const ready = new Promise<void>((resolve) => {
    responseReady = resolve;
  });
  await page.route(`**/api/datasets/${source.id}/monomer`, async (route) => {
    const result = await route.fetch();
    responseReady();
    await held.promise;
    await route.fulfill({ response: result });
  });
  await chooser.getByRole('button', { name: 'Use one monomer', exact: true }).click();
  await ready;
  await expect(
    chooser.getByRole('button', { name: 'Selecting monomer…', exact: true }),
  ).toBeDisabled();
  await expect(
    studio.getByRole('button', { name: 'Start simulation', exact: true }),
  ).toBeDisabled();
  await expect(studio.getByRole('button', { name: 'GROMACS', exact: true })).toBeDisabled();
  await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const finished = page.waitForResponse((response) =>
    response.url().endsWith(`/api/datasets/${source.id}/monomer`),
  );
  held.release();
  await finished;
  await expect(chooser).toBeVisible();
  await expect(page.locator('.structure-card h2')).toHaveText(source.name);
  await expect(studio.locator('.monomer-result')).toHaveCount(0);
  await expect(studio.getByRole('button', { name: 'Start simulation', exact: true })).toBeEnabled();
});
