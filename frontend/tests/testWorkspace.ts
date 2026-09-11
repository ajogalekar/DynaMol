import { test as base } from '@playwright/test';

/** Older UI regressions explicitly start from the bundled demonstration now
 * that normal application startup restores a user's last scene. Never enable
 * this fixture against a real user's data root. */
export const test = base.extend<{ isolatedWorkspace: void }>({
  isolatedWorkspace: [
    async ({ request, baseURL }, use) => {
      if (
        process.env.DYNAMOL_TEST_ISOLATED !== '1' ||
        !baseURL ||
        ['8765', '5173', '4173'].includes(new URL(baseURL).port)
      ) {
        throw new Error(
          'These regressions reset the workspace. Start an isolated DYNAMOL_DATA_DIR server on a dedicated test port, then set DYNAMOL_TEST_ISOLATED=1 and DYNAMOL_BASE_URL. User-facing ports 8765/5173/4173 are protected.',
        );
      }
      const response = await request.get('/api/datasets/demo');
      if (!response.ok())
        throw new Error('Copy the bundled demo into the isolated test data root first.');
      const demo = await response.json();
      const ca = demo.atoms.filter(
        (atom: { name: string; category: string }) =>
          atom.name === 'CA' && atom.category === 'protein',
      );
      const atoms = [
        ca[Math.floor(ca.length * 0.27)].index,
        ca[Math.floor(ca.length * 0.67)].index,
      ];
      const measured = await request.post(`/api/datasets/${demo.id}/measurements`, {
        data: { kind: 'distance', atoms },
      });
      if (!measured.ok()) throw new Error(`Demo measurement failed: ${await measured.text()}`);
      const measurement = {
        ...(await measured.json()),
        id: 'audit-default-distance',
        color: '#74dfc5',
        visible: false,
        label: `${demo.atoms[atoms[0]].residue}${demo.atoms[atoms[0]].resid} ↔ ${demo.atoms[atoms[1]].residue}${demo.atoms[atoms[1]].resid}`,
      };
      const saved = await request.post('/api/workspace', {
        data: {
          state: {
            version: 1,
            dataset_id: demo.id,
            measurements: [measurement],
            active_measurement: measurement.id,
          },
        },
      });
      if (!saved.ok())
        throw new Error(`Could not initialize isolated test workspace: ${await saved.text()}`);
      await use();
    },
    { auto: true },
  ],
});
