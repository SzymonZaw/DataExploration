import json
from dynamics.yamanaka_trajectory_poc import run

if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))
