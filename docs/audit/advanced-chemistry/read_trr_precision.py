"""Read the first modern TRR frame without downcasting its stored precision.

Layout follows pinned GROMACS2025.4 trrio.cpp and gmxfio_xdr.cpp. This audit
reader rejects obsolete payload blocks and requires coordinates, box and forces.
"""
from pathlib import Path
import struct
import numpy as np


def read_first_trr(path):
    data = Path(path).read_bytes()
    offset = 0

    def unpack(fmt):
        nonlocal offset
        size = struct.calcsize(fmt)
        if offset + size > len(data):
            raise ValueError('Truncated TRR header')
        value = struct.unpack_from(fmt, data, offset)
        offset += size
        return value

    magic, declared, length = unpack('>3i')
    if magic != 1993 or not 0 < length <= declared <= 256:
        raise ValueError('Invalid TRR magic/string header')
    version = data[offset:offset+length].rstrip(b'\0')
    offset += (length + 3) // 4 * 4
    if version != b'GMX_trn_file':
        raise ValueError('Unrecognized TRR version')
    fields = unpack('>13i')
    ir, energy, box_size, virial, pressure, topology, symmetry, x_size, v_size, f_size, n, step, nre = fields
    if any([ir, energy, topology, symmetry]) or n <= 0 or nre != 0:
        raise ValueError('Unsupported legacy TRR blocks or atom inventory')
    precision = box_size // 9
    if precision not in (4, 8) or box_size != 9 * precision:
        raise ValueError('Invalid stored TRR real precision')
    vector_size = n * 3 * precision
    if x_size != vector_size or f_size != vector_size or v_size not in (0, vector_size):
        raise ValueError('TRR coordinate/force payload length mismatch')
    if virial not in (0, box_size) or pressure not in (0, box_size):
        raise ValueError('TRR matrix length mismatch')
    time, lambda_value = unpack('>2d' if precision == 8 else '>2f')

    def array(byte_count, shape):
        nonlocal offset
        if offset + byte_count > len(data):
            raise ValueError('Truncated TRR array')
        result = np.frombuffer(data, dtype='>f8' if precision == 8 else '>f4',
                               count=byte_count // precision, offset=offset).astype(np.float64).reshape(shape)
        offset += byte_count
        if not np.isfinite(result).all():
            raise ValueError('Nonfinite TRR array')
        return result

    box = array(box_size, (3, 3))
    if virial:
        array(virial, (3, 3))
    if pressure:
        array(pressure, (3, 3))
    coordinates = array(x_size, (n, 3))
    if v_size:
        array(v_size, (n, 3))
    forces = array(f_size, (n, 3))
    if not np.isfinite([time, lambda_value]).all():
        raise ValueError('Nonfinite frame time/lambda')
    return {'stored_precision_bytes': precision, 'step': step, 'time_ps': time,
            'coords_nm': coordinates, 'box_nm': box, 'forces_kj_mol_nm': forces,
            'first_frame_bytes': offset, 'file_bytes': len(data)}
