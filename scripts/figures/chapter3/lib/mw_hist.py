from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_mw_hist(input_csv: Path, output_path: Path, style=None) -> None:
    dpi = style.dpi if style else 200
    label_fs = style.axis_label_fontsize if style else 12
    label_fw = style.axis_label_fontweight if style else "bold"
    tick_fs = style.tick_fontsize if style else 8
    tick_fw = style.tick_fontweight if style else "bold"

    df = pd.read_csv(input_csv)

    series = pd.to_numeric(df["mw_struct_kda"], errors="coerce").dropna()

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.hist(series, bins=40, edgecolor="black")
    ax.set_xlabel("mw_struct_kda", fontsize=label_fs, fontweight=label_fw)
    ax.set_ylabel("count", fontsize=label_fs, fontweight=label_fw)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(tick_fs)
        label.set_fontweight(tick_fw)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
