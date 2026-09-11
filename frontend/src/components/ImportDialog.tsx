import { useRef, useState } from 'react';
import { ArrowUpFromLine, FileBox, FileVideo, FolderOpen, LoaderCircle, X } from 'lucide-react';
import { api } from '../api';
import type { Dataset } from '../types';
export default function ImportDialog({
  onClose,
  onLoaded,
}: {
  onClose: () => void;
  onLoaded: (dataset: Dataset) => void;
}) {
  const [topology, setTopology] = useState<File | null>(null),
    [trajectory, setTrajectory] = useState<File | null>(null),
    [stride, setStride] = useState(1),
    [interval, setInterval] = useState(''),
    [busy, setBusy] = useState(false),
    [error, setError] = useState('');
  const topInput = useRef<HTMLInputElement>(null),
    trajInput = useRef<HTMLInputElement>(null);
  function assign(files: File[]) {
    for (const f of files) {
      if (/\.(pdb|cif|mmcif|gro|h5)$/i.test(f.name)) setTopology(f);
      else setTrajectory(f);
    }
  }
  async function submit(e: React.FormEvent) {
    e.preventDefault();
    if (!topology) return;
    setBusy(true);
    setError('');
    const form = new FormData();
    form.append('topology', topology);
    if (trajectory) form.append('trajectory', trajectory);
    form.append('stride', String(stride));
    if (interval) form.append('frame_interval_ps', interval);
    try {
      onLoaded(await api.upload(form));
    } catch (e) {
      setError((e as Error).message);
      setBusy(false);
    }
  }
  return (
    <div
      className="modal-backdrop"
      onMouseDown={(e) => {
        if (e.target === e.currentTarget && !busy) onClose();
      }}
    >
      <section
        className="modal import-modal"
        role="dialog"
        aria-modal="true"
        aria-labelledby="import-title"
      >
        <div className="modal-header">
          <div className="eyebrow">YOUR NEXT DISCOVERY</div>
          <button
            className="icon-button"
            aria-label="Close import"
            disabled={busy}
            onClick={onClose}
          >
            <X size={18} />
          </button>
        </div>
        <h2 id="import-title">Bring your molecules.</h2>
        <p className="modal-subtitle">Open a structure on its own, or pair it with a trajectory.</p>
        <form onSubmit={submit}>
          <div
            className="dropzone"
            onDragOver={(e) => e.preventDefault()}
            onDrop={(e) => {
              e.preventDefault();
              assign(Array.from(e.dataTransfer.files));
            }}
          >
            <div className="drop-icon">
              <ArrowUpFromLine size={25} />
            </div>
            <strong>Drop your molecular files here</strong>
            <span>or choose them below. Files stay on this computer.</span>
          </div>
          <div className="file-pair">
            <button
              type="button"
              className={`file-slot ${topology ? 'chosen' : ''}`}
              onClick={() => topInput.current?.click()}
            >
              <FileBox size={23} />
              <b>{topology?.name ?? 'Choose structure'}</b>
              <small>PDB, mmCIF, GRO, HDF5</small>
              <FolderOpen size={14} />
            </button>
            <button
              type="button"
              className={`file-slot ${trajectory ? 'chosen' : ''}`}
              onClick={() => trajInput.current?.click()}
            >
              <FileVideo size={23} />
              <b>{trajectory?.name ?? 'Add trajectory'}</b>
              <small>XTC, DCD, TRR, NetCDF, H5, PDB</small>
              <FolderOpen size={14} />
            </button>
          </div>
          <input
            hidden
            ref={topInput}
            type="file"
            accept=".pdb,.cif,.mmcif,.gro,.h5"
            onChange={(e) => setTopology(e.target.files?.[0] ?? null)}
          />
          <input
            hidden
            ref={trajInput}
            type="file"
            accept=".xtc,.dcd,.trr,.nc,.netcdf,.h5,.pdb,.xyz,.mdcrd,.crd"
            onChange={(e) => setTrajectory(e.target.files?.[0] ?? null)}
          />
          <div className="form-row">
            <label>
              Read every{' '}
              <div className="input-unit">
                <input
                  type="number"
                  min="1"
                  max="1000"
                  value={stride}
                  onChange={(e) => setStride(Number(e.target.value))}
                  required
                />
                <span>frames</span>
              </div>
            </label>
            <label>
              Time between original frames{' '}
              <div className="input-unit">
                <input
                  type="number"
                  min="0.000001"
                  step="any"
                  placeholder="From file"
                  value={interval}
                  onChange={(e) => setInterval(e.target.value)}
                />
                <span>ps</span>
              </div>
            </label>
          </div>
          <p className="form-note">
            Topology and trajectory must have exactly the same atom order. For files without
            readable timing, supply a frame interval to get a physical time axis. Increase the
            stride for large trajectories.
          </p>
          {error && (
            <div className="error-box" role="alert">
              {error}
            </div>
          )}
          <div className="modal-actions">
            <button type="button" className="text-button" disabled={busy} onClick={onClose}>
              Cancel
            </button>
            <button className="primary-button" disabled={!topology || busy}>
              {busy ? (
                <>
                  <LoaderCircle size={16} className="spin" /> Reading molecules…
                </>
              ) : (
                <>
                  <FolderOpen size={16} /> Open in DynaMol
                </>
              )}
            </button>
          </div>
        </form>
      </section>
    </div>
  );
}
