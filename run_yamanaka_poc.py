"""Command-line entry point for the Stage 2.11 Yamanaka POC."""
from dynamics.yamanaka_poc import run


if __name__ == "__main__":
    import json
    print(json.dumps(run(), indent=2, default=str))
