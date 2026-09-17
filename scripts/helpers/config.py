import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict

import yaml


@dataclass
class Style:
    dpi: int
    axis_label_fontsize: int
    axis_label_fontweight: str
    tick_fontsize: int
    tick_fontweight: str
    kpam_color: str
    lit_color: str
    pyruvylation_color: str
    acetylation_color: str
    sgnh_domain_color: str
    gray_color: str
    both_color: str
    no_ecod_color: str
    ssrbh_color: str
    sslbh_color: str
    af3_set_colors: Dict[str, str] = field(default_factory=dict)


class Config:
    _CONFIG_PATH = Path(__file__).resolve().parents[2] / "config" / "config.yml"

    def __init__(self):
        with open(self._CONFIG_PATH) as f:
            data = yaml.safe_load(f)
        paths = data["paths"]
        # THESIS_INPUT_DIR / THESIS_OUTPUT_DIR override config.yml, e.g. to write a trial
        # run somewhere else without editing the file.
        self.input_dir = Path(os.environ.get("THESIS_INPUT_DIR") or paths["input_dir"]).expanduser()
        self.output_dir = Path(os.environ.get("THESIS_OUTPUT_DIR") or paths["output_dir"]).expanduser()
        # gwas_path is optional: the GWAS data ships inside the input folder.
        self.gwas_path = Path(
            os.environ.get("THESIS_GWAS_PATH") or paths.get("gwas_path") or self.input_dir / "data-gwas"
        ).expanduser()
        if not self.input_dir.is_dir():
            raise FileNotFoundError(
                f"input_dir does not exist: {self.input_dir}\n"
                f"Set paths.input_dir in {self._CONFIG_PATH} to the unpacked Figshare input."
            )
        self.output_dir.mkdir(parents=True, exist_ok=True)

        s = data.get("style", {})
        self.style = Style(
            dpi=s.get("dpi", 300),
            axis_label_fontsize=s.get("axis_label_fontsize", 12),
            axis_label_fontweight=s.get("axis_label_fontweight", "bold"),
            tick_fontsize=s.get("tick_fontsize", 8),
            tick_fontweight=s.get("tick_fontweight", "bold"),
            kpam_color=s.get("kpam_color", "#1f77b4"),
            lit_color=s.get("lit_color", "#d62728"),
            pyruvylation_color=s.get("pyruvylation_color", "#3283bf"),
            acetylation_color=s.get("acetylation_color", "#d45515"),
            sgnh_domain_color=s.get("sgnh_domain_color", "#c9a227"),
            gray_color=s.get("gray_color", "#bfbfbf"),
            both_color=s.get("both_color", "#7A3FD4"),
            no_ecod_color=s.get("no_ecod_color", "#ffffff"),
            ssrbh_color=s.get("ssrbh_color", "#1f77b4"),
            sslbh_color=s.get("sslbh_color", "#9467bd"),
            af3_set_colors=data.get("af3_set_colors", {}),
        )
