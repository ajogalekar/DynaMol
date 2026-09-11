import { useEffect, useRef, useState } from 'react';
import {
  ArrowUpRight,
  Check,
  CircleAlert,
  Clock3,
  LoaderCircle,
  Square,
  Terminal,
  X,
} from 'lucide-react';
import type { Job } from '../types';

interface StructureJobMonitorProps {
  complexPreparation?: boolean;
  job?: Job;
  submitting?: 'preparation' | 'solvation' | null;
  loadingResult?: boolean;
  loadError?: string;
  resultInView?: boolean;
  onCancel?: () => void;
  onOpenResult?: () => void;
  onViewLog?: () => void;
  onDismiss?: () => void;
}

const terminalStatuses = new Set(['completed', 'failed', 'cancelled', 'interrupted']);

function elapsedLabel(seconds: number) {
  const whole = Math.max(0, Math.floor(seconds));
  if (whole < 60) return `${whole}s`;
  if (whole < 3600) return `${Math.floor(whole / 60)}m ${whole % 60}s`;
  return `${Math.floor(whole / 3600)}h ${Math.floor((whole % 3600) / 60)}m`;
}

export default function StructureJobMonitor({
  complexPreparation = false,
  job: previousJob,
  submitting,
  loadingResult = false,
  loadError,
  resultInView = false,
  onCancel,
  onOpenResult,
  onViewLog,
  onDismiss,
}: StructureJobMonitorProps) {
  const job = submitting ? undefined : previousJob;
  const [now, setNow] = useState(Date.now);
  const active = !!job && !terminalStatuses.has(job.status);
  const ticking = active || !!submitting;
  const clock = useRef<{ key: string; startedAt: number } | null>(null);
  const clockKey = job?.id ?? (submitting ? `starting-${submitting}` : 'idle');
  if (clock.current?.key !== clockKey) {
    const createdAt = job ? Date.parse(job.created_at) : NaN;
    const currentTime = Date.now();
    clock.current = {
      key: clockKey,
      startedAt: Number.isFinite(createdAt)
        ? Math.min(createdAt, currentTime)
        : currentTime - Math.max(0, job?.elapsed_seconds ?? 0) * 1000,
    };
  }
  useEffect(() => {
    if (!ticking) return;
    setNow(Date.now());
    const timer = window.setInterval(() => setNow(Date.now()), 1000);
    return () => window.clearInterval(timer);
  }, [clockKey, ticking]);

  if (!job && !submitting && !loadError) return null;

  const operation = submitting ?? (job?.engine === 'solvation' ? 'solvation' : 'preparation');
  const water = operation === 'solvation';
  const complex = job ? job.config.complex === true : complexPreparation;
  const subject = complex ? 'Complex' : 'Protein';
  const complete = job?.status === 'completed';
  const failed = job?.status === 'failed';
  const stopped = job?.status === 'cancelled' || job?.status === 'interrupted';
  const attention = failed || !!loadError;
  const spinning = ticking || loadingResult;
  const title = submitting
    ? water
      ? 'Starting explicit-water setup…'
      : `Starting ${complex ? 'complex' : 'protein'} preparation…`
    : complete
      ? water
        ? 'Explicit water ready'
        : `${subject} preparation complete`
      : failed
        ? water
          ? 'Explicit-water setup failed'
          : `${subject} preparation failed`
        : stopped
          ? `${water ? 'Explicit-water setup' : `${subject} preparation`} ${job?.status}`
          : active
            ? water
              ? 'Building explicit water'
              : `Preparing your ${complex ? 'complex' : 'protein'}`
            : 'Preparation needs attention';
  const stage = submitting
    ? 'Starting a background worker. You can keep exploring the viewer.'
    : loadingResult
      ? 'Opening the prepared structure in the viewer…'
      : complete
        ? resultInView
          ? water
            ? 'The explicit-water structure is in the viewer.'
            : 'Prepared structure is in the viewer.'
          : 'Open the prepared structure in the viewer.'
        : (job?.stage ?? 'The request could not finish.');
  const totalStages = Math.max(0, Math.floor(job?.total_steps ?? 0));
  const completedStages = Math.min(totalStages, Math.max(0, Math.floor(job?.completed_steps ?? 0)));
  const progress = totalStages > 0 ? Math.round((completedStages / totalStages) * 100) : 0;
  const elapsed = ticking
    ? Math.max(job?.elapsed_seconds ?? 0, (now - clock.current.startedAt) / 1000)
    : (job?.elapsed_seconds ?? 0);
  const error = loadError || (failed || job?.status === 'interrupted' ? job?.error : undefined);
  const Icon = spinning
    ? LoaderCircle
    : attention
      ? CircleAlert
      : complete
        ? Check
        : stopped
          ? Square
          : CircleAlert;

  return (
    <section
      className={`structure-job-monitor${complete && !attention ? ' is-complete' : ''}${attention ? ' needs-attention' : ''}`}
      aria-label={water ? 'Explicit-water setup monitor' : `${subject} preparation monitor`}
    >
      <div className="structure-job-monitor-heading">
        <span className="structure-job-monitor-icon" aria-hidden="true">
          <Icon size={17} className={spinning ? 'structure-job-monitor-spinner' : undefined} />
        </span>
        <div
          className="structure-job-monitor-status"
          role="status"
          aria-live="polite"
          aria-atomic="true"
        >
          <strong>{title}</strong>
          <span>{stage}</span>
        </div>
        {!ticking && onDismiss && (
          <button
            type="button"
            className="structure-job-monitor-dismiss"
            onClick={onDismiss}
            aria-label={water ? 'Dismiss explicit-water status' : 'Dismiss preparation status'}
          >
            <X size={14} />
          </button>
        )}
      </div>

      {job && (active || complete) && totalStages > 0 && (
        <div className="structure-job-monitor-progress">
          <div
            className="structure-job-monitor-track"
            role="progressbar"
            aria-label={water ? 'Explicit-water setup progress' : `${subject} preparation progress`}
            aria-valuemin={0}
            aria-valuemax={100}
            aria-valuenow={progress}
            aria-valuetext={`${completedStages} of ${totalStages} stages complete`}
          >
            <span style={{ width: `${progress}%` }} />
          </div>
          <div className="structure-job-monitor-progress-label">
            <span>
              Stage progress · {completedStages}/{totalStages} stages
            </span>
            <b>{progress}%</b>
          </div>
        </div>
      )}

      {error && (
        <details className="structure-job-monitor-error" role="alert">
          <summary>
            <span>{error}</span>
            <b>Details</b>
          </summary>
          <p>{error}</p>
        </details>
      )}

      <div className="structure-job-monitor-footer">
        {job || submitting ? (
          <span
            className="structure-job-monitor-elapsed"
            aria-label={`Elapsed time: ${elapsedLabel(elapsed)}`}
          >
            <Clock3 size={11} aria-hidden="true" /> {elapsedLabel(elapsed)} elapsed
          </span>
        ) : (
          <span />
        )}
        <div className="structure-job-monitor-actions">
          {active && onCancel && (
            <button
              type="button"
              onClick={onCancel}
              aria-label={
                water
                  ? 'Cancel explicit-water setup'
                  : `Cancel ${complex ? 'complex' : 'protein'} preparation`
              }
            >
              <Square size={10} aria-hidden="true" /> Cancel
            </button>
          )}
          {job && onViewLog && (
            <button
              type="button"
              onClick={onViewLog}
              aria-label={
                water
                  ? 'View explicit-water setup log'
                  : `View ${complex ? 'complex' : 'protein'} preparation log`
              }
            >
              <Terminal size={12} aria-hidden="true" /> View log
            </button>
          )}
          {complete && !resultInView && onOpenResult && (
            <button
              type="button"
              className="structure-job-monitor-open"
              onClick={onOpenResult}
              disabled={loadingResult}
            >
              <ArrowUpRight size={12} aria-hidden="true" /> Open result
            </button>
          )}
        </div>
      </div>
    </section>
  );
}
