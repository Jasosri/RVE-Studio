from rvestudio.validation import run_all14_benchmark
b=run_all14_benchmark()
print(b['summary'].to_string(index=False))
print(f"Paper agreement: {b['within_tolerance']}/{b['total']} within +/-{b['tolerance_pct']}%; max: {b['max_abs_paper_pct']:.4f}%")
