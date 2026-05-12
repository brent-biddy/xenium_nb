from contextlib import contextmanager
import time

# Module-level list accumulates (label, elapsed) tuples across all timer calls
# within a notebook execution. Persists for the lifetime of the Jupyter kernel.
_timings = []


@contextmanager
def timer(label):
    """Context manager that times a block and prints elapsed time inline."""
    start = time.perf_counter()
    yield
    elapsed = time.perf_counter() - start
    _timings.append((label, elapsed))
    minutes, seconds = divmod(elapsed, 60)
    if minutes > 0:
        print(f"[{label}] {int(minutes)}m {seconds:.1f}s")
    else:
        print(f"[{label}] {seconds:.2f}s")


def write_timings_tsv(path):
    """Write recorded timings to a TSV file with columns: step, seconds.

    Appends a final 'Total' row summing all recorded elapsed times.
    """
    with open(path, "w") as f:
        f.write("step\tseconds\n")
        total = 0.0
        for label, elapsed in _timings:
            f.write(f"{label}\t{elapsed:.4f}\n")
            total += elapsed
        f.write(f"Total\t{total:.4f}\n")


def timing_summary(path=None):
    """Print a formatted table of all recorded timings and a total.

    If path is provided, also write the timings to a TSV file.
    """
    if not _timings:
        print("No timings recorded.")
        return
    col = max(max(len(label) for label, _ in _timings), 4)
    print(f"\n{'Step':<{col}}  {'Time':>10}")
    print("-" * (col + 13))
    total = 0
    for label, elapsed in _timings:
        minutes, seconds = divmod(elapsed, 60)
        time_str = f"{int(minutes)}m {seconds:.1f}s" if minutes > 0 else f"{seconds:.2f}s"
        print(f"{label:<{col}}  {time_str:>10}")
        total += elapsed
    print("-" * (col + 13))
    minutes, seconds = divmod(total, 60)
    total_str = f"{int(minutes)}m {seconds:.1f}s" if minutes > 0 else f"{seconds:.2f}s"
    print(f"{'Total':<{col}}  {total_str:>10}")
    if path is not None:
        write_timings_tsv(path)
