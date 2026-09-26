"""premise_bench: the PREMISE comparative benchmark and parameter ablation.

Run `python3 -m premise_bench --help` inside `nix develop .#benchmark` for the commands.

Layout, by stage:
    config, splits, utils, metrics    configuration, the split registry, shared helpers, metrics
    runner, pipeline                  child processes and the end-to-end driver (`run`)
    methods/                          one module per compared tool: build, classify, load
    data/                             reference cleaning, decoys, synthetic and real samples
    evaluate/                         truth and scoring
    ablation/                         the one-at-a-time parameter sweep
"""
