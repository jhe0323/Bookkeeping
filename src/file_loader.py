# file_loader.py
# -*- coding: utf-8 -*-
"""
Utilities for loading LUH2 states/transitions and PFT maps, and normalizing their
latitude/longitude conventions for downstream simulation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import netCDF4 as nc
import numpy as np


@dataclass
class GridMeta:
    lat: np.ndarray
    lon: np.ndarray
    lat_desc: bool
    lon_0_360: bool
    lat_name: str
    lon_name: str
    internal_lat: np.ndarray
    internal_lon: np.ndarray

class _MemoryDim:
    def __init__(self, size: int):
        self.size = int(size)

class _MemoryDataset:
    """
    A light in-memory dataset wrapper.
    Provides .variables and .dimensions similar to netCDF4.Dataset.
    """
    def __init__(self, variables: dict, dimensions: dict):
        self.variables = variables
        self.dimensions = {k: _MemoryDim(v) for k, v in dimensions.items()}
        

def _guess_lat_lon_names(ds: nc.Dataset) -> Tuple[str, str]:
    candidates_lat = ("lat", "latitude", "y")
    candidates_lon = ("lon", "longitude", "x")
    lat_name = next((n for n in candidates_lat if n in ds.variables), None)
    lon_name = next((n for n in candidates_lon if n in ds.variables), None)
    if lat_name is None or lon_name is None:
        raise KeyError("Cannot find lat/lon coordinate variables.")
    return lat_name, lon_name

def _as_1d_coord(v, coord_type: str = "lat") -> np.ndarray:
    arr = np.asarray(v[:], dtype=float)
    if arr.ndim == 1:
        arr_1d = arr
    elif arr.ndim == 2:
        if coord_type == "lat":
            col0 = arr[:, 0]
            if np.allclose(arr, col0[:, None], equal_nan=True):
                arr_1d = col0
            else:
                row0 = arr[0, :]
                arr_1d = row0 if np.allclose(arr, row0[None, :], equal_nan=True) else col0
        else:
            row0 = arr[0, :]
            if np.allclose(arr, row0[None, :], equal_nan=True):
                arr_1d = row0
            else:
                col0 = arr[:, 0]
                arr_1d = col0 if np.allclose(arr, col0[:, None], equal_nan=True) else row0
    else:
        raise ValueError(f"Unsupported coordinate ndim={arr.ndim}")

    arr_1d = np.asarray(arr_1d, dtype=float).ravel()
    return arr_1d[np.isfinite(arr_1d)]

def _is_descending(a: np.ndarray) -> bool:
    return bool(np.all(np.diff(a) < 0)) if a.size >= 2 else False

def _is_lon_0_360(lon: np.ndarray) -> bool:
    lon = lon[np.isfinite(lon)]
    return bool(lon.min() >= 0 and lon.max() > 180) if lon.size > 0 else False

def _wrap_lon_180(lon: np.ndarray) -> np.ndarray:
    return ((lon + 180.0) % 360.0) - 180.0

def _optimize_indices_to_slice(idx_array):
    """
    Convert a regular NumPy index array to a native slice.
    This avoids expensive fancy indexing in netCDF4 reads.
    """
    if not isinstance(idx_array, np.ndarray):
        return idx_array

    if idx_array.size == 0:
        return slice(0, 0)

    if idx_array.size == 1:
        i = int(idx_array[0])
        return slice(i, i + 1)

    step = int(idx_array[1]) - int(idx_array[0])

    if np.all(np.diff(idx_array) == step):
        if step > 0:
            return slice(int(idx_array[0]), int(idx_array[-1]) + 1, step)
        else:
            stop = int(idx_array[-1]) - 1
            if stop < 0:
                stop = None
            return slice(int(idx_array[0]), stop, step)

    return idx_array


class FileLoader:
    def load_luh2_dataset(
        self,
        path: str,
        lat_slice: slice = None,
        lon_slice: slice = None,
    ) -> Tuple[_MemoryDataset, GridMeta]:
        
        ds = nc.Dataset(path, "r")
        try:
            lat_name, lon_name = _guess_lat_lon_names(ds)
            lat_raw_full = _as_1d_coord(ds.variables[lat_name], coord_type="lat")
            lon_raw_full = _as_1d_coord(ds.variables[lon_name], coord_type="lon")

            is_desc = _is_descending(lat_raw_full)
            lon_0_360 = _is_lon_0_360(lon_raw_full)
            nlat_full = lat_raw_full.size
            nlon_full = lon_raw_full.size

            if lat_slice is None: lat_slice = slice(0, nlat_full)
            if lon_slice is None: lon_slice = slice(0, nlon_full)

            # 1. 纬度切片映射
            if is_desc:
                i0_asc = lat_slice.start if lat_slice.start is not None else 0
                i1_asc = lat_slice.stop if lat_slice.stop is not None else nlat_full
                lat_indices = np.arange(nlat_full)[slice(nlat_full - i1_asc, nlat_full - i0_asc)]
            else:
                lat_indices = np.arange(nlat_full)[lat_slice]

            lat_raw = lat_raw_full[lat_indices]
            lat_meta = lat_raw[::-1] if is_desc else lat_raw

            # 2. 经度切片映射
            lon_res = abs(float(np.median(np.diff(lon_raw_full))))
            nlon_internal = int(round(360.0 / lon_res))
            lon_internal_full = np.linspace(-180.0 + 0.5 * lon_res, 180.0 - 0.5 * lon_res, nlon_internal)
            target_lon_180 = lon_internal_full[lon_slice]
            lon_raw_180 = _wrap_lon_180(lon_raw_full)

            lon_indices = []
            for x in target_lon_180:
                jj = int(np.nanargmin(np.abs(lon_raw_180 - x)))
                diff = abs(lon_raw_180[jj] - x)
                if diff > 0.51 * lon_res:
                    raise ValueError(
                        f"Cannot map internal lon {x} to source lon. "
                        f"nearest={lon_raw_full[jj]}, wrapped={lon_raw_180[jj]}, diff={diff}"
                    )
                lon_indices.append(jj)
            
            lon_indices = np.asarray(lon_indices, dtype=np.int64)
            lon_raw = lon_raw_full[lon_indices]

            lat_res = abs(float(np.median(np.diff(lat_raw_full))))
            internal_lat_full = np.linspace(-90.0 + 0.5 * lat_res, 90.0 - 0.5 * lat_res, nlat_full)
            meta = GridMeta(
                lat=lat_meta, lon=lon_raw, lat_desc=is_desc, lon_0_360=lon_0_360,
                lat_name=lat_name, lon_name=lon_name,
                internal_lat=internal_lat_full[lat_slice], internal_lon=target_lon_180,
            )

            variables = {}
            dimensions = {}

            for dim_name, dim in ds.dimensions.items():
                if dim_name == lat_name: dimensions[dim_name] = len(lat_indices)
                elif dim_name == lon_name: dimensions[dim_name] = len(lon_indices)
                else: dimensions[dim_name] = len(dim)

            # 生成优化后的切片
            lat_slice_opt = _optimize_indices_to_slice(lat_indices)
            lon_slice_opt = _optimize_indices_to_slice(lon_indices)

            for name, var in ds.variables.items():
                if name == lat_name:
                    variables[name] = np.asarray(lat_raw, dtype=np.float32)
                    continue
                if name == lon_name:
                    variables[name] = np.asarray(lon_raw, dtype=np.float32)
                    continue

                index = []
                for d in var.dimensions:
                    if d == lat_name:
                        index.append(lat_slice_opt)
                    elif d == lon_name:
                        index.append(lon_slice_opt)
                    else:
                        index.append(slice(None))

                try:
                    arr = var[tuple(index)]
                except Exception as e:
                    raise RuntimeError(f"Failed reading variable {name}") from e

                # 降精度并去除掩码
                if hasattr(arr, "filled"):
                    arr = arr.filled(0.0)
                arr = np.asarray(arr, dtype=np.float32)

                variables[name] = arr

            out = _MemoryDataset(variables=variables, dimensions=dimensions)
            
            print(
                f"[LOAD] {path}: "
                f"lat={len(lat_indices)}, lon={len(lon_indices)}, "
                f"lat_slice_opt={lat_slice_opt}, lon_slice_opt={lon_slice_opt}, "
                f"lon_raw_first_last=({float(lon_raw[0])}, {float(lon_raw[-1])}), "
                f"lon_internal_first_last=({float(target_lon_180[0])}, {float(target_lon_180[-1])}), "
                f"vars={len(variables)}",
                flush=True,
            )

            return out, meta

        finally:
            ds.close()

    def load_pft_map(
        self,
        pft_path: str,
        pft_var: Optional[str] = None,
        *,
        target_lat_asc: np.ndarray,
        target_lon_180: np.ndarray,
        target_lon_res: float,
        expected_n_pft: Optional[int] = None,
        lat_slice: slice = None,
        lon_slice: slice = None,
    ) -> np.ndarray:
        """Load a dominant or fractional PFT map and return PFT fractions.

        The returned array always has shape ``(lat, lon, pft)``. Two input
        layouts are supported:

        1. dominant map: ``(lat, lon)`` with 1-based PFT codes;
        2. fractional map: any 3-D ordering of ``pft``, ``lat`` and ``lon``.

        Fractional values are cleaned and normalized independently in each
        grid cell. A coarser PFT map may be copied to a finer model grid by
        nearest cell centre. A finer PFT map is not silently aggregated to a
        coarser model grid because that operation should be area weighted.

        ``lat_slice`` and ``lon_slice`` are retained for call compatibility;
        alignment is performed directly from the supplied target coordinates.
        """
        del lat_slice, lon_slice  # target coordinates already describe the local slice

        target_lat_asc = np.asarray(target_lat_asc, dtype=float).ravel()
        target_lon_180 = _wrap_lon_180(
            np.asarray(target_lon_180, dtype=float).ravel()
        )
        if target_lat_asc.size == 0 or target_lon_180.size == 0:
            raise ValueError("Target PFT grid is empty.")

        pft_dim_aliases = {
            "pft", "pfts", "npft", "veget", "veg", "vegtype",
            "vegetation", "class", "classes", "type", "types",
        }
        time_aliases = {"time", "time_counter", "year"}
        lat_aliases = {"lat", "latitude", "y"}
        lon_aliases = {"lon", "longitude", "x"}

        def nearest_indices(source, target, *, periodic=False):
            source = np.asarray(source, dtype=float).ravel()
            target = np.asarray(target, dtype=float).ravel()
            if source.size == 0:
                raise ValueError("Empty source coordinate.")

            out = np.empty(target.size, dtype=np.int64)
            for k, value in enumerate(target):
                if periodic:
                    diff = np.abs(((source - value + 180.0) % 360.0) - 180.0)
                else:
                    diff = np.abs(source - value)
                out[k] = int(np.nanargmin(diff))
            return out

        def coord_resolution(values):
            values = np.asarray(values, dtype=float).ravel()
            if values.size < 2:
                return np.nan
            unique = np.unique(np.round(values, 10))
            if unique.size < 2:
                return np.nan
            return abs(float(np.median(np.diff(np.sort(unique)))))

        def choose_variable(ds):
            if pft_var is not None:
                if pft_var not in ds.variables:
                    raise KeyError(
                        f"PFT variable '{pft_var}' not found in {pft_path}. "
                        f"Available variables: {list(ds.variables)}"
                    )
                return pft_var

            preferred = [
                "pft_fraction", "pft_frac", "pft_fractions",
                "vegetation_fraction", "vegetfrac", "maxvegetfrac",
                "pft", "PFT", "dominant_pft", "pft_map",
                "pft_dominant", "dominantPFT",
            ]
            for name in preferred:
                if name in ds.variables:
                    return name

            candidates = []
            for name, candidate in ds.variables.items():
                dims = {d.lower() for d in candidate.dimensions}
                non_time_ndim = sum(d.lower() not in time_aliases for d in candidate.dimensions)
                has_lat = bool(dims & lat_aliases)
                has_lon = bool(dims & lon_aliases)
                if has_lat and has_lon and non_time_ndim in (2, 3):
                    candidates.append(name)

            if len(candidates) == 1:
                return candidates[0]
            raise KeyError(
                "Cannot infer the PFT data variable. Pass pft_var explicitly. "
                f"Candidate variables: {candidates}; all variables: {list(ds.variables)}"
            )

        ds = nc.Dataset(pft_path, "r")
        try:
            chosen_var = choose_variable(ds)
            var = ds.variables[chosen_var]
            original_dims = list(var.dimensions)

            lat_name, lon_name = _guess_lat_lon_names(ds)
            src_lat = _as_1d_coord(ds.variables[lat_name], coord_type="lat")
            src_lon_raw = _as_1d_coord(ds.variables[lon_name], coord_type="lon")
            src_lon_180 = _wrap_lon_180(src_lon_raw)

            src_lat_res = coord_resolution(src_lat)
            src_lon_res = coord_resolution(src_lon_180)
            target_lat_res = coord_resolution(target_lat_asc)
            target_lon_res_eff = (
                float(target_lon_res)
                if np.isfinite(target_lon_res) and target_lon_res > 0
                else coord_resolution(target_lon_180)
            )

            # Do not disguise a required area-weighted aggregation as nearest-neighbour.
            if (
                np.isfinite(src_lat_res)
                and np.isfinite(target_lat_res)
                and src_lat_res < target_lat_res * (1.0 - 1e-6)
            ) or (
                np.isfinite(src_lon_res)
                and np.isfinite(target_lon_res_eff)
                and src_lon_res < target_lon_res_eff * (1.0 - 1e-6)
            ):
                raise ValueError(
                    "The PFT map is finer than the model grid. Pre-aggregate PFT "
                    "fractions with area weighting before loading. "
                    f"source_res=({src_lat_res}, {src_lon_res}), "
                    f"target_res=({target_lat_res}, {target_lon_res_eff})"
                )

            lat_indices = nearest_indices(src_lat, target_lat_asc, periodic=False)
            lon_indices = nearest_indices(src_lon_180, target_lon_180, periodic=True)

            lat_indexer = _optimize_indices_to_slice(lat_indices)
            lon_indexer = _optimize_indices_to_slice(lon_indices)

            # Coordinate variable names and data dimension names are not always identical.
            dims_lower = [d.lower() for d in original_dims]
            lat_axis_original = next(
                (ax for ax, d in enumerate(dims_lower) if d in lat_aliases),
                None,
            )
            lon_axis_original = next(
                (ax for ax, d in enumerate(dims_lower) if d in lon_aliases),
                None,
            )
            if lat_axis_original is None or lon_axis_original is None:
                raise ValueError(
                    f"PFT variable '{chosen_var}' must contain lat/lon dimensions. "
                    f"dimensions={original_dims}"
                )

            index = []
            remaining_dims = []
            for d in original_dims:
                dl = d.lower()
                if dl in time_aliases:
                    index.append(0)
                elif dl in lat_aliases:
                    index.append(lat_indexer)
                    remaining_dims.append(d)
                elif dl in lon_aliases:
                    index.append(lon_indexer)
                    remaining_dims.append(d)
                else:
                    index.append(slice(None))
                    remaining_dims.append(d)

            arr = var[tuple(index)]
            if hasattr(arr, "filled"):
                arr = arr.filled(0.0)
            arr = np.asarray(arr)

        finally:
            ds.close()

        dims_lower = [d.lower() for d in remaining_dims]
        lat_axis = next((ax for ax, d in enumerate(dims_lower) if d in lat_aliases), None)
        lon_axis = next((ax for ax, d in enumerate(dims_lower) if d in lon_aliases), None)
        if lat_axis is None or lon_axis is None:
            raise ValueError(
                f"Cannot infer lat/lon axes after reading PFT map. "
                f"shape={arr.shape}, dimensions={remaining_dims}"
            )

        nlat_t = target_lat_asc.size
        nlon_t = target_lon_180.size

        if arr.ndim == 2:
            # Legacy dominant map -> one-hot fraction cube.
            if expected_n_pft is None:
                raise ValueError(
                    "expected_n_pft is required when loading a dominant PFT map."
                )
            arr = np.moveaxis(arr, (lat_axis, lon_axis), (0, 1))
            if arr.shape != (nlat_t, nlon_t):
                raise ValueError(
                    f"Dominant PFT shape mismatch: arr={arr.shape}, "
                    f"target=({nlat_t}, {nlon_t})"
                )

            dominant = np.asarray(arr, dtype=np.float64)
            rounded = np.rint(dominant)
            integer_like = np.isfinite(dominant) & (np.abs(dominant - rounded) <= 1e-5)
            codes = rounded.astype(np.int64, copy=False)
            valid = integer_like & (codes >= 1) & (codes <= int(expected_n_pft))

            fractions = np.zeros(
                (nlat_t, nlon_t, int(expected_n_pft)),
                dtype=np.float32,
            )
            rows, cols = np.nonzero(valid)
            fractions[rows, cols, codes[rows, cols] - 1] = 1.0
            mode = "dominant->fraction"

        elif arr.ndim == 3:
            pft_axis = next(
                (ax for ax, d in enumerate(dims_lower) if d in pft_dim_aliases),
                None,
            )
            if pft_axis is None and expected_n_pft is not None:
                candidates = [
                    ax for ax, size in enumerate(arr.shape)
                    if ax not in (lat_axis, lon_axis) and size == int(expected_n_pft)
                ]
                if len(candidates) == 1:
                    pft_axis = candidates[0]

            if pft_axis is None:
                remaining_axes = [ax for ax in range(3) if ax not in (lat_axis, lon_axis)]
                if len(remaining_axes) == 1:
                    pft_axis = remaining_axes[0]

            if pft_axis is None:
                raise ValueError(
                    f"Cannot infer PFT axis. shape={arr.shape}, dimensions={remaining_dims}"
                )

            fractions = np.moveaxis(
                arr,
                (lat_axis, lon_axis, pft_axis),
                (0, 1, 2),
            ).astype(np.float32, copy=False)

            if fractions.shape[:2] != (nlat_t, nlon_t):
                raise ValueError(
                    f"Fractional PFT shape mismatch: arr={fractions.shape}, "
                    f"target=({nlat_t}, {nlon_t}, npft)"
                )
            if expected_n_pft is not None and fractions.shape[2] != int(expected_n_pft):
                raise ValueError(
                    f"PFT map contains {fractions.shape[2]} classes, but the current "
                    f"parameter configuration contains {expected_n_pft}."
                )

            fractions = np.where(np.isfinite(fractions), fractions, 0.0)
            fractions[fractions < 0.0] = 0.0
            sums = fractions.sum(axis=2, keepdims=True, dtype=np.float64)
            fractions = np.divide(
                fractions,
                sums,
                out=np.zeros_like(fractions),
                where=sums > 0.0,
            )
            mode = "fraction"

        else:
            raise ValueError(
                f"Unsupported PFT variable rank after removing time: ndim={arr.ndim}, "
                f"shape={arr.shape}, dimensions={remaining_dims}. "
                "Expected a 2-D dominant map or a 3-D fractional map."
            )

        zero_cells = int(np.count_nonzero(fractions.sum(axis=2) <= 0.0))
        print(
            f"[PFT] {pft_path}: var={chosen_var}, mode={mode}, "
            f"shape={fractions.shape}, zero_fraction_cells={zero_cells}",
            flush=True,
        )
        return fractions
