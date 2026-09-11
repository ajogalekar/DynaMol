import { useEffect, useRef, useState } from 'react';
import { api } from './api';
import type { Job, RunMeasurementSnapshot } from './types';

interface FollowedRun {
  id: string;
  source: string;
  automatic: boolean;
}
const key = 'dynamol.followed-run.v1';
function readFollowed(): FollowedRun | null {
  try {
    const saved = JSON.parse(localStorage.getItem(key) ?? 'null');
    return saved && typeof saved.id === 'string' && typeof saved.source === 'string'
      ? { ...saved, automatic: saved.automatic === true }
      : null;
  } catch {
    return null;
  }
}

/** Run data stays separate from the starting structure's saved-frame analysis. */
export function useSimulationTracking(options: {
  jobs: Job[];
  datasetId?: string;
  ready: boolean;
  onComplete: (job: Job, snapshot: RunMeasurementSnapshot) => Promise<boolean>;
  onNotice: (message: string) => void;
}) {
  const [followed, setFollowed] = useState(readFollowed);
  const [snapshot, setSnapshot] = useState<RunMeasurementSnapshot | null>(null);
  const [pollError, setPollError] = useState('');
  const [result, setResult] = useState<Job | null>(null);
  const latest = useRef(options);
  latest.current = options;
  const currentFollow = useRef(followed);
  currentFollow.current = followed;
  function update(next: FollowedRun | null) {
    currentFollow.current = next;
    setFollowed(next);
    try {
      if (next) localStorage.setItem(key, JSON.stringify(next));
      else localStorage.removeItem(key);
    } catch {
      /* Tracking works even when local persistence is unavailable. */
    }
  }
  const job = options.jobs.find((j) => j.id === followed?.id);
  useEffect(() => {
    if (!followed) return;
    const id = followed.id;
    let stopped = false;
    let timer = 0;
    const controller = new AbortController();
    async function poll() {
      try {
        const next = await api.jobMeasurements(id, controller.signal);
        if (stopped || currentFollow.current?.id !== id) return;
        setSnapshot(next);
        setPollError('');
        const state = latest.current;
        const job = state.jobs.find((j) => j.id === id);
        if (job && state.ready && job.status === 'completed' && job.dataset_id) {
          const shouldOpen =
            currentFollow.current.automatic && state.datasetId === currentFollow.current.source;
          const opened = shouldOpen ? await state.onComplete(job, next) : false;
          if (stopped || currentFollow.current?.id !== id) return;
          if (!opened) {
            setResult(job);
            state.onNotice(`${job.name} is ready to explore. Your current scene was kept.`);
          }
          update(null);
          return;
        }
        if (job && ['failed', 'cancelled', 'interrupted'].includes(job.status)) {
          // Retain the last saved values; resume can follow the same job again.
          return;
        }
      } catch (error) {
        if (stopped) return;
        setPollError(`Live measurements are temporarily unavailable: ${(error as Error).message}`);
      }
      if (!stopped) timer = window.setTimeout(() => void poll(), 1500);
    }
    void poll();
    return () => {
      stopped = true;
      controller.abort();
      window.clearTimeout(timer);
    };
  }, [followed?.id, job?.status]);
  return {
    job,
    snapshot,
    pollError,
    result,
    showing: !!followed?.automatic && options.datasetId === followed.source,
    follow(job: Job) {
      if (!['openmm', 'gromacs'].includes(job.engine)) return;
      setSnapshot(null);
      setPollError('');
      setResult(null);
      update({ id: job.id, source: String(job.config.dataset_id), automatic: true });
    },
    keepCurrentScene() {
      if (currentFollow.current) update({ ...currentFollow.current, automatic: false });
    },
    dismissResult() {
      setResult(null);
    },
  };
}
