"""Shared helpers for XClone model pipelines (RDR, BAF, Combine)."""

import os
from datetime import datetime, timezone
from typing import Optional, Tuple

import xclone


def configure_warnings(ignore_warnings: bool) -> None:
    """Silence warnings when requested."""
    if ignore_warnings:
        import warnings

        warnings.filterwarnings("ignore")


def load_config(module_name: str, config_file):
    """Resolve the configuration object and default to provided module settings."""
    from .._config import XCloneConfig

    if config_file is None:
        print(
            "Model configuration file not specified.\n"
            f"Default settings in XClone-{module_name} will be used."
        )
        return XCloneConfig(module=module_name)
    return config_file


def resolve_output_dirs(out_dir: Optional[str]) -> Tuple[str, str, str]:
    """Return base, data, and plot output directories, creating them if needed."""
    base_out_dir = out_dir or os.path.join(os.getcwd(), "XCLONE_OUT")
    out_data_dir = os.path.join(base_out_dir, "data")
    out_plot_dir = os.path.join(base_out_dir, "plot")

    xclone.al.dir_make(out_data_dir)
    xclone.al.dir_make(out_plot_dir)
    return base_out_dir, out_data_dir, out_plot_dir


def write_adata_safe(adata, filepath: str, label: str) -> None:
    """Write AnnData safely, surfacing a concise warning if it fails."""
    try:
        adata.write(filepath)
    except Exception as exc:  # pylint: disable=broad-except
        print(f"[XClone Warning] {exc}")
    else:
        print(f"[XClone hint] {label} saved in {os.path.dirname(filepath)}.")


def log_duration(logger, start_time: datetime, label: str) -> None:
    """Log elapsed time in seconds with a module label."""
    time_passed = datetime.now(timezone.utc) - start_time
    logger.info("%s finished (%d seconds)", label, time_passed.total_seconds())
