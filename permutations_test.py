from tabulate import tabulate
import math

def count_paths(N, M):
    """
    Count number of walks of length N on {0,…,M}, starting at 0,
    with steps +1 or -1, that never go below 0 or above M.
    
    Returns a list counts[0..M] of path‐counts ending at each position.
    """
    # dp_prev[k] = # ways to be at k after the previous step
    dp_prev = [0] * (M+1)
    dp_prev[0] = 1  # at step 0, we're at position 0 in exactly 1 way

    for _ in range(N):
        dp_curr = [0] * (M+1)
        for k in range(M+1):
            # step up from k–1 → k
            if k > 0:
                dp_curr[k] += dp_prev[k-1]
            # step down from k+1 → k
            if k < M:
                dp_curr[k] += dp_prev[k+1]
        dp_prev = dp_curr

    return dp_prev


def show_all(N_max, M):
    """
    For each move count 0..N_max, compute counts[0..M] and print a table where
    for each step k you have two adjacent columns: the comma-formatted count and its log2.
    """
    # Build interleaved headers: ["Moves", "Step 0", "log2(0)", "Step 1", "log2(1)", …]
    headers = ["Moves"]
    for k in range(M+1):
        headers += [f"Step {k}", f"log2(Step {k})"]

    table = []
    for n in range(N_max + 1):
        counts = count_paths(n, M)
        row = [n]
        for c in counts:
            # formatted count
            fmt = f"{c:,}"
            # log2 or “-inf”
            lg = f"{math.log2(c):.2f}" if c > 0 else "-inf"
            row += [fmt, lg]
        table.append(row)

    print(tabulate(table, headers, tablefmt="github"))

if __name__ == "__main__":
    N_max = 60   # up to 10 moves
    M =60        # steps 0..3
    show_all(N_max, M)

# # Example usage:
# if __name__ == "__main__":
#     N = 50    # total moves
#     M = 30    # steps range 0..3
#     counts = count_paths(N, M)
#     for target_step, num_paths in enumerate(counts):
#         print(f"→ step {target_step}: {num_paths:,} paths")