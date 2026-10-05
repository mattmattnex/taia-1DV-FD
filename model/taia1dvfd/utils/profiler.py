import time
import torch
from collections import defaultdict


class StepProfiler:
    """Lightweight per-section timer with running statistics."""

    def __init__(self):
        self._t0 = None
        self._last = None
        self.times = defaultdict(float)        # this step's section times (ms)
        self.totals = defaultdict(float)       # cumulative across all steps
        self.counts = defaultdict(int)
        self.peaks = defaultdict(float)        # slowest occurrence per section

    def start(self):
        self._t0 = time.perf_counter()
        self._last = self._t0
        self.times.clear()

    def mark(self, name):
        """Record elapsed time since the previous mark under `name`."""
        now = time.perf_counter()
        dt_ms = (now - self._last) * 1000.0
        self._last = now
        self.times[name] = dt_ms
        self.totals[name] += dt_ms
        self.counts[name] += 1
        self.peaks[name] = max(self.peaks[name], dt_ms)

    @property
    def total_ms(self):
        return sum(self.times.values())

    def avg(self, name):
        n = self.counts[name]
        return self.totals[name] / n if n else 0.0

    def format_step(self, kt, nt):
        parts = " | ".join(f"{name}={dt:6.2f}" for name, dt in self.times.items())
        return (f"[{kt:>6}/{nt}]  total={self.total_ms:7.2f}ms   {parts}")

    def format_summary(self):
        lines = ["", f"{'section':<8}{'avg(ms)':>10}{'peak(ms)':>10}{'calls':>8}"]
        lines.append("-" * 36)
        for name in self.totals:
            lines.append(f"{name:<8}{self.avg(name):>10.2f}"
                         f"{self.peaks[name]:>10.2f}{self.counts[name]:>8}")
        return "\n".join(lines)