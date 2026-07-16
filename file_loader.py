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
        lat_slice: slice = None,
        lon_slice: slice = None,
    ) -> np.ndarray:
        
        ds = nc.Dataset(pft_path, "r")
        try:
            if pft_var is None:
                common = ["maxvegetfrac", "pft", "PFT", "dominant_pft", "pft_map", "pft_dominant", "dominantPFT"]
                pft_var = next((n for n in common if n in ds.variables), None)

            var = ds.variables[pft_var]
            dim_names = list(var.dimensions)

            lat_name = next((n for n in ("lat", "latitude", "y") if n in ds.variables), None)
            lon_name = next((n for n in ("lon", "longitude", "x") if n in ds.variables), None)

            lat = np.asarray(ds.variables[lat_name][:], dtype=float).ravel()
            lon = np.asarray(ds.variables[lon_name][:], dtype=float).ravel()
            src_lat_desc = bool(lat.size >= 2 and lat[0] > lat[-1])

            if lat_slice is None: lat_slice = slice(0, lat.size)
            if lon_slice is None: lon_slice = slice(0, lon.size)

            if src_lat_desc:
                i0_asc = lat_slice.start if lat_slice.start is not None else 0
                i1_asc = lat_slice.stop if lat_slice.stop is not None else lat.size
                pft_lat_indices = np.arange(lat.size)[slice(lat.size - i1_asc, lat.size - i0_asc)]
            else:
                pft_lat_indices = np.arange(lat.size)[lat_slice]

            lon_res = abs(float(np.median(np.diff(lon))))
            nlon_internal = int(round(360.0 / lon_res))
            lon_internal_full = np.linspace(-180.0 + 0.5 * lon_res, 180.0 - 0.5 * lon_res, nlon_internal)
            target_lon = lon_internal_full[lon_slice]
            lon_180 = _wrap_lon_180(lon)

            pft_lon_indices = []
            for x in target_lon:
                jj = int(np.nanargmin(np.abs(lon_180 - x)))
                diff = abs(lon_180[jj] - x)
                if diff > 0.51 * lon_res:
                    raise ValueError(
                        f"Cannot map PFT internal lon {x} to source lon. "
                        f"nearest={lon[jj]}, wrapped={lon_180[jj]}, diff={diff}"
                    )
                pft_lon_indices.append(jj)
            pft_lon_indices = np.asarray(pft_lon_indices, dtype=np.int64)

            # 生成优化后的切片
            lat_slice_opt = _optimize_indices_to_slice(pft_lat_indices)
            lon_slice_opt = _optimize_indices_to_slice(pft_lon_indices)

            index = []
            for d in dim_names:
                dl = d.lower()
                if dl in ("time", "time_counter"):
                    index.append(0)
                elif d == lat_name or dl in ("lat", "latitude", "y"):
                    index.append(lat_slice_opt)
                elif d == lon_name or dl in ("lon", "longitude", "x"):
                    index.append(lon_slice_opt)
                else:
                    index.append(slice(None))

            arr = var[tuple(index)]
            
            if hasattr(arr, "filled"):
                arr = arr.filled(0.0)
            arr = np.asarray(arr, dtype=np.float32)
            
            dim_names = [d for d in dim_names if d.lower() not in ("time", "time_counter")]

        finally:
            ds.close()

        if arr.ndim == 3:
            dims_lower = [d.lower() for d in dim_names]
            pft_axis, lat_axis, lon_axis = None, None, None
            for ax, d in enumerate(dims_lower):
                if d in ("veget", "pft", "class", "classes"): pft_axis = ax
                elif d in ("lat", "latitude", "y"): lat_axis = ax
                elif d in ("lon", "longitude", "x"): lon_axis = ax

            if pft_axis is None:
                candidate_axes = [ax for ax, n in enumerate(arr.shape) if n == 15]
                if len(candidate_axes) == 1: pft_axis = candidate_axes[0]

            if lat_axis is None or lon_axis is None or pft_axis is None:
                raise ValueError(
                    f"Cannot infer axes for PFT map. shape={arr.shape}, dims={dim_names}"
                )

            arr = np.moveaxis(arr, (lat_axis, lon_axis, pft_axis), (0, 1, 2))
            if src_lat_desc: arr = np.flip(arr, axis=0)

        nlat_t, nlon_t = target_lat_asc.size, target_lon_180.size
        
        if arr.ndim == 3:
            if arr.shape[0] != nlat_t or arr.shape[1] != nlon_t:
                raise ValueError(
                    f"PFT local shape mismatch: arr={arr.shape}, "
                    f"target=({nlat_t}, {nlon_t}, npft)"
                )
                
            arr = np.where(np.isfinite(arr), arr, 0.0)
            arr[arr < 0.0] = 0.0
            s = arr.sum(axis=2, keepdims=True)
            ok = s > 0.0
            arr = np.divide(arr, s, out=np.zeros_like(arr), where=ok)

        return arr