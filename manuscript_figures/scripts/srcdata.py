"""Collect per-panel Source Data blocks for one figure (read by build_source_data.py)."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "build" / "source_data"


class SourceData:
    def __init__(self, fig):
        self.dir = ROOT / fig
        self.dir.mkdir(parents=True, exist_ok=True)
        for old in self.dir.glob("*.csv"):
            old.unlink()
        self.blocks = []

    def add(self, panel, desc, df):
        """panel e.g. 'Fig. 6d'; desc = one-line description incl. units; df = tidy DataFrame."""
        name = f"{len(self.blocks):02d}_{panel.replace('Fig. ', '').replace(' ', '_').replace('/', '-')}.csv"
        df.to_csv(self.dir / name, index=False)
        self.blocks.append({"panel": panel, "desc": desc, "csv": name})

    def save(self):
        (self.dir / "manifest.json").write_text(json.dumps(self.blocks, indent=1))
