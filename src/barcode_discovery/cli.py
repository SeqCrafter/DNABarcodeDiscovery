from __future__ import annotations

import argparse

from . import align, config, evaluate, scan


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        prog="barcode-discovery",
        description="Pi sliding-window + sequence-cluster-ratio barcode screening.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("align", help="step 0 (optional): MSA with pymafft (Rust MAFFT)")
    m.add_argument("fasta", help="unaligned multi-FASTA")
    m.add_argument("-o", "--out", default="aligned.fasta")
    m.add_argument("--strategy", default="auto",
                   choices=["auto", "fftns2", "fftnsi", "ginsi", "linsi", "einsi"])
    m.add_argument("--maxiterate", type=int, default=0)
    m.add_argument("-t", "--threads", type=int, default=0, help="0 = all cores")
    m.add_argument("--adjust-direction", action="store_true",
                   help="detect and reverse-complement sequences on the opposite strand")

    s = sub.add_parser("scan", help="step 1: sliding-window Pi -> window_pi.tsv + html plot")
    s.add_argument("fasta", help="MAFFT-aligned multi-FASTA")
    s.add_argument("-o", "--out", default="results")
    s.add_argument("-w", "--window", type=int, default=200, help="window size")
    s.add_argument("-s", "--step", type=int, default=10, help="step size")
    s.add_argument("--pi-threshold", type=float, default=0.05,
                   help="threshold line drawn on the plot")
    s.add_argument("--gap-mode", choices=["excluded", "considered"], default="excluded",
                   help="DnaSP 'Gaps in Sliding Window': excluded = window/step counted "
                        "in valid sites (default), considered = in alignment columns")

    i = sub.add_parser("init-config", help="write an empty step-2 config template")
    i.add_argument("-o", "--out", default="config.toml")
    i.add_argument("--force", action="store_true", help="overwrite existing file")

    e = sub.add_parser("evaluate", help="step 2: cluster ratio from a config file")
    e.add_argument("config", help="TOML config (see init-config)")

    a = p.parse_args(argv)
    try:
        if a.cmd == "align":
            align.run(a.fasta, a.out, a.strategy, a.maxiterate, a.threads,
                      a.adjust_direction)
        elif a.cmd == "scan":
            scan.run(a.fasta, a.out, a.window, a.step, a.pi_threshold,
                     a.gap_mode == "excluded")
        elif a.cmd == "init-config":
            config.write_template(a.out, a.force)
            print(f"Template written to {a.out}")
        else:
            cfg = config.load(a.config)
            cands = evaluate.run(cfg)
            print(f"{'rank':<5}{'name':<14}{'aln_region':<22}{'Pi':<10}{'ratio%':<9}pass")
            for c in cands:
                reg = ";".join(f"{x}-{y}" for x, y in c.aln_ranges)
                print(f"{c.rank:<5}{c.name:<14}{reg:<22}{c.pi:<10.5f}"
                      f"{c.cluster_ratio:<9.2f}{'PASS' if c.passed else 'FAIL'}")
            print(f"\nResults written to {cfg.out}/barcode_candidates.tsv")
    except (ValueError, FileExistsError, FileNotFoundError) as exc:
        p.exit(1, f"error: {exc}\n")


if __name__ == "__main__":
    main()
